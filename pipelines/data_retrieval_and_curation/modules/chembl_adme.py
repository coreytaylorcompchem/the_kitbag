import aiohttp
import asyncio
import nest_asyncio
import os
import json
import requests

from concurrent.futures import ThreadPoolExecutor, as_completed

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
    targets = config.get("targets", {})

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

    def build_url(batch):
        assay_filter = f"assay_chembl_id__in={','.join(batch)}"
        fields = "activity_id,assay_chembl_id,molecule_chembl_id,standard_type,standard_units,relation,standard_value,target_chembl_id"
        url = f"{BASE_URL}?{assay_filter}&only={fields}&limit=1000"
        return url

    async def fetch_batch(session, batch):
        url = build_url(batch)
        results = []

        while url:
            async with session.get(url) as response:
                data = await response.json()
                results.extend(data.get("activities", []))

                next_url = data["page_meta"].get("next")
                if next_url and next_url.startswith("/"):
                    next_url = "https://www.ebi.ac.uk" + next_url

                url = next_url

        return results

    async def fetch_target_page(session, target_id, offset):

        params = {
            "target_chembl_id": target_id,
            "assay_type": "B",
            "limit": 1000,
            "offset": offset
        }

        async with session.get(BASE_URL, params=params) as response:
            data = await response.json()
            return data.get("activities", [])

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

    async def fetch_target_activities(target_name, target_id):

        connector = aiohttp.TCPConnector(limit_per_host=MAX_CONCURRENT_REQUESTS)

        async with aiohttp.ClientSession(connector=connector) as session:

            # First request to determine total pages
            params = {
                "target_chembl_id": target_id,
                "assay_type": "B",
                "limit": 1000
            }

            async with session.get(BASE_URL, params=params) as response:
                meta = await response.json()

            total = meta["page_meta"]["total_count"]
            pages = (total // 1000) + 1

            logger.info(f"Retrieved {total} activity records")

            offsets = [i * 1000 for i in range(pages)]

            tasks = [
                fetch_target_page(session, target_id, offset)
                for offset in offsets
            ]

            results = []

            for f in tqdm(
                asyncio.as_completed(tasks),
                total=len(tasks),
                desc=f"Downloading {target_name}"
            ):
                results.extend(await f)

            return results

    async def run_retrieval():

        all_results = {}

        # ---- ADME assays via keyword matching ----
        for cat, ids in assay_ids.items():

            logger.info(f"Fetching {cat} activities")

            assay_batches = [
                ids[i:i+BATCH_SIZE]
                for i in range(0, len(ids), BATCH_SIZE)
            ]

            acts = await fetch_all_batches(assay_batches)

            logger.info(f"Retrieved {len(acts)} {cat} activity records")

            all_results[cat] = acts

        # ---- Target-based retrieval (e.g. hERG) ----
        for target_name, target_id in targets.items():

            logger.info(f"Fetching {target_name} activities via target {target_id}")

            acts = await fetch_target_activities(target_name, target_id)

            logger.info(f"Retrieved {len(acts)} {target_name} activity records")

            all_results[target_name] = acts

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

    Uses in-memory caching so each molecule/assay/target is only fetched once.
    """

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    MAX_CONCURRENT = config.get("max_concurrent", 30)
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

    if isinstance(activities, dict) and "activities" in activities:
        activities = activities["activities"]

    # CACHES 
    mol_cache = {}
    assay_cache = {}
    target_cache = {}

    # Robust async fetch 
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

        record = {field: act.get(field) for field in RELEVANT_FIELDS}

        record.update({
            "smiles": None,
            "assay_description": None,
            "assay_type": None,
            "target_organism": None,
            "error": None
        })

        try:

            mol_id = record.get("molecule_chembl_id")
            assay_id = record.get("assay_chembl_id")
            target_id = record.get("target_chembl_id")

            # MOLECULE 
            if mol_id:
                if mol_id not in mol_cache:
                    mol_url = f"{BASE_URL}/molecule/{mol_id}.json"
                    mol_json = await fetch_json(session, mol_url)
                    mol_cache[mol_id] = (
                        (mol_json or {})
                        .get("molecule_structures", {})
                        .get("canonical_smiles")
                    )
                record["smiles"] = mol_cache[mol_id]

            # ASSAY 
            if assay_id:
                if assay_id not in assay_cache:
                    assay_url = f"{BASE_URL}/assay/{assay_id}.json"
                    assay_json = await fetch_json(session, assay_url)

                    assay_cache[assay_id] = {
                        "description": (assay_json or {}).get("description"),
                        "assay_type": (assay_json or {}).get("assay_type"),
                    }

                record["assay_description"] = assay_cache[assay_id]["description"]
                record["assay_type"] = assay_cache[assay_id]["assay_type"]

            # TARGET 
            if target_id:
                if target_id not in target_cache:
                    target_url = f"{BASE_URL}/target/{target_id}.json"
                    target_json = await fetch_json(session, target_url)

                    target_cache[target_id] = (target_json or {}).get("organism")

                record["target_organism"] = target_cache[target_id]

        except Exception as e:
            record["error"] = str(e)

        return record

    def load_existing_records(assay_name):

        enriched_cat = []
        outfile = os.path.join(OUTPUT_DIR, f"{assay_name}.jsonl")

        if os.path.exists(outfile):
            with open(outfile, "r") as f:
                for line in f:
                    try:
                        enriched_cat.append(json.loads(line))
                    except Exception:
                        continue

        return enriched_cat

    async def enrich_category(assay_name, acts):

        acts = [a for a in acts if isinstance(a, dict)]

        enriched_cat = load_existing_records(assay_name)
        processed_ids = {r["activity_id"] for r in enriched_cat}

        remaining_acts = [
            a for a in acts if a.get("activity_id") not in processed_ids
        ]

        total_acts = len(remaining_acts)

        pbar = tqdm(total=total_acts, desc=assay_name, unit="activity")

        outfile = os.path.join(OUTPUT_DIR, f"{assay_name}.jsonl")
        connector = aiohttp.TCPConnector(limit_per_host=MAX_CONCURRENT)

        async with aiohttp.ClientSession(connector=connector) as session:
            with open(outfile, "a") as f_out:
                for i in range(0, total_acts, BATCH_SIZE):

                    batch = remaining_acts[i:i + BATCH_SIZE]

                    tasks = [
                        enrich_activity_api(session, act)
                        for act in batch
                    ]

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
            logger.info(f"Enriching {assay_name}")
            results[assay_name] = await enrich_category(assay_name, acts)
        return results

    logger.info("Supplementing ADME activities with assay metadata")

    enriched_data = asyncio.run(run_all())

    return enriched_data