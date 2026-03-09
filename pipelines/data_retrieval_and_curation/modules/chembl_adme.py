import aiohttp
import asyncio
import nest_asyncio
import os
import json

from chembl_webresource_client.new_client import new_client

from multiprocessing import Pool, cpu_count, get_context
from multiprocessing.pool import ThreadPool
import pandas as pd
from tqdm import tqdm

from pipeline.task_registry import register_task

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

@register_task(
    "retrieve_chembl_assays",
    category="ADME",
    description="Retrieve assay metadata from CHEMBL API"
)
def retrieve_chembl_assays(config, data=None):

    import requests
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from tqdm import tqdm
    import pandas as pd

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/assay.json"
    LIMIT = config.get("limit", 1000)

    params = {
        "assay_type__in": "A,P",
        "limit": LIMIT
    }

    r = requests.get(BASE_URL, params=params)
    meta = r.json()

    total = meta["page_meta"]["total_count"]
    pages = (total // LIMIT) + 1

    def fetch_page(page):

        p = {
            "assay_type__in": "A,P",
            "limit": LIMIT,
            "offset": page * LIMIT
        }

        r = requests.get(BASE_URL, params=p)
        return r.json()["assays"]

    results = []

    with ThreadPoolExecutor(max_workers=16) as executor:

        futures = [executor.submit(fetch_page, p) for p in range(pages)]

        for f in tqdm(as_completed(futures), total=pages, desc="Downloading assays"):
            results.extend(f.result())

    assay_df = pd.DataFrame(results)

    logger.info(f"Retrieved {len(assay_df)} ADME assays")

    return {"assays": assay_df}

@register_task(
    "filter_assays_by_keywords",
    category="ADME",
    description="Filter CHEMBL assays using keyword matching"
)
def filter_assays_by_keywords(config, assays=None):

    if isinstance(assays, dict):
        assays = assays["assays"]

    assay_keywords = config.get("assay_keywords", {})

    assays["description"] = assays["description"].fillna("").str.lower()

    filtered = []

    for category, keywords in assay_keywords.items():

        mask = assays["description"].apply(
            lambda x: any(k in x for k in keywords)
        )

        tmp = assays[mask].copy()
        tmp["category"] = category

        filtered.append(tmp)

    filtered_df = pd.concat(filtered)

    assay_ids = {
        cat: filtered_df.loc[
            filtered_df["category"] == cat,
            "assay_chembl_id"
        ].unique().tolist()
        for cat in assay_keywords.keys()
    }

    return {
        "assay_ids": assay_ids,
        "filtered_assays": filtered_df
    }

@register_task(
    "retrieve_assay_activities",
    category="ADME",
    description="Retrieve activity records for selected assays"
)
def retrieve_assay_activities(config, assay_ids=None):

    if isinstance(assay_ids, dict) and "assay_ids" in assay_ids:
        assay_ids = assay_ids["assay_ids"]

    BATCH_SIZE = config.get("batch_size", 50)
    MAX_CONCURRENT_REQUESTS = config.get("max_concurrent_requests", 10)

    RELEVANT_TYPES = {
        "papp", "apparent permeability", "permeability",
        "solubility", "logp", "logd", "efflux ratio"
    }

    def build_url(batch):
        base = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        assay_filter = "&".join([f"assay_chembl_id__in={aid}" for aid in batch])
        fields = (
            "activity_id,assay_chembl_id,molecule_chembl_id,"
            "standard_type,standard_units,relation,standard_value,target_chembl_id"
        )
        url = f"{base}?{assay_filter}&only={fields}&limit=1000"
        return url

    async def fetch_batch(session, batch):

        url = build_url(batch)

        for attempt in range(3):
            try:
                async with session.get(url) as response:

                    if response.status == 200:
                        data = await response.json()
                        acts = data.get("activities", [])

                        return [
                            a for a in acts
                            if a.get("standard_type", "").lower().strip() in RELEVANT_TYPES
                        ]

                    await asyncio.sleep(2 ** attempt)

            except Exception:
                await asyncio.sleep(2 ** attempt)

        return []

    async def fetch_all_batches(batches):

        connector = aiohttp.TCPConnector(limit_per_host=MAX_CONCURRENT_REQUESTS)

        async with aiohttp.ClientSession(connector=connector) as session:

            tasks = [fetch_batch(session, b) for b in batches]

            results = []

            for f in tqdm(
                asyncio.as_completed(tasks),
                total=len(tasks),
                desc="Fetching activities"
            ):
                results.extend(await f)

            return results

    async def run_retrieval():

        all_results = {}

        for cat, ids in assay_ids.items():

            logger.info(f"Fetching {cat} activities")

            assay_batches = [
                ids[i:i+BATCH_SIZE]
                for i in range(0, len(ids), BATCH_SIZE)
            ]

            acts = await fetch_all_batches(assay_batches)

            logger.info(f"Retrieved {len(acts)} {cat} activity records")

            all_results[cat] = acts

        return all_results

    results = asyncio.run(run_retrieval())

    return {"activities": results}

nest_asyncio.apply()  # nested asyncio loops work in notebooks or pipelines

@register_task(
    "enrich_activity_metadata",
    category="ADME",
    description="Fetch molecule, assay and target metadata in parallel using aiohttp"
)
def enrich_activity_metadata(config, activities=None):
    """
    Enrich per-assay activity records with molecule SMILES, assay description/type, and target organism.

    Saves per-assay JSONL files to OUTPUT_DIR so progress can be resumed.
    """
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    MAX_CONCURRENT = config.get("max_concurrent", 20)
    BATCH_SIZE = config.get("batch_size", 100)
    OUTPUT_DIR = config.get("output", {}).get("directory", "outputs/adme")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    RELEVANT_FIELDS = [
        "activity_id",
        "assay_chembl_id",
        "molecule_chembl_id",
        "standard_type",
        "standard_units",
        "relation",
        "standard_value",
        "target_chembl_id"
    ]

    if activities is None:
        raise ValueError("No activities provided to enrich_activity_metadata()")

    # unwrap if wrapped by previous task
    if isinstance(activities, dict) and "activities" in activities:
        activities = activities["activities"]

    async def fetch_json(session, url, max_retries=3, base_delay=1.0):
        for attempt in range(max_retries + 1):
            try:
                async with session.get(url) as response:
                    if response.status == 200:
                        return await response.json()
                    elif response.status in {429, 500, 502, 503, 504}:
                        await asyncio.sleep(base_delay * (2 ** attempt))
                    else:
                        return None
            except Exception:
                await asyncio.sleep(base_delay * (2 ** attempt))
        return None

    async def enrich_activity_api(session, act):
        """Enrich a single activity record safely with molecule, assay, and target info."""
        record = {field: act.get(field) for field in RELEVANT_FIELDS}
        record.update({
            "smiles": None,
            "assay_description": None,
            "assay_type": None,
            "target_organism": None,
            "error": None
        })

        try:
            if record.get("molecule_chembl_id"):
                mol_url = f"{BASE_URL}/molecule/{record['molecule_chembl_id']}.json"
                mol_json = await fetch_json(session, mol_url)
                record["smiles"] = (mol_json or {}).get("molecule_structures", {}).get("canonical_smiles")

            if record.get("assay_chembl_id"):
                assay_url = f"{BASE_URL}/assay/{record['assay_chembl_id']}.json"
                assay_json = await fetch_json(session, assay_url)
                if assay_json:
                    record["assay_description"] = assay_json.get("description")
                    record["assay_type"] = assay_json.get("assay_type")

            if record.get("target_chembl_id"):
                target_url = f"{BASE_URL}/target/{record['target_chembl_id']}.json"
                target_json = await fetch_json(session, target_url)
                if target_json:
                    record["target_organism"] = target_json.get("organism")

        except Exception as e:
            record["error"] = str(e)

        return record

    def load_existing_records(assay_name):
        """Load previously saved JSONL records for a given assay."""
        enriched_cat = []
        outfile = os.path.join(OUTPUT_DIR, f"{assay_name}.json")
        if os.path.exists(outfile):
            with open(outfile, "r") as f:
                for line in f:
                    try:
                        enriched_cat.append(json.loads(line))
                    except Exception:
                        continue
        return enriched_cat

    async def enrich_category(assay_name, acts):
        """Enrich a single assay's activities asynchronously, in batches."""
        
        acts = [a for a in acts if isinstance(a, dict)]

        enriched_cat = load_existing_records(assay_name)
        processed_ids = {r["activity_id"] for r in enriched_cat}

        remaining_acts = [a for a in acts if a.get("activity_id") not in processed_ids]
        total_acts = len(remaining_acts)
        pbar = tqdm(total=total_acts, desc=assay_name, unit="activity")
        outfile = os.path.join(OUTPUT_DIR, f"{assay_name}.jsonl")

        connector = aiohttp.TCPConnector(limit_per_host=MAX_CONCURRENT)
        async with aiohttp.ClientSession(connector=connector) as session:
            with open(outfile, "a") as f_out:
                for i in range(0, total_acts, BATCH_SIZE):
                    batch = remaining_acts[i:i + BATCH_SIZE]
                    tasks = [enrich_activity_api(session, act) for act in batch]
                    for fut in asyncio.as_completed(tasks):
                        enriched = await fut
                        enriched_cat.append(enriched)
                        f_out.write(json.dumps(enriched) + "\n")
                        f_out.flush()
                        pbar.update(1)

        pbar.close()
        return enriched_cat

    async def run_all():
        results = {}
        for assay_name, acts in activities.items():
            results[assay_name] = await enrich_category(assay_name, acts)
        return results

    enriched_data = asyncio.run(run_all())
    
    return enriched_data