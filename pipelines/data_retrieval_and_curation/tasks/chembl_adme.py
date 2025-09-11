from pipeline.task_registry import register_task
from chembl_webresource_client.new_client import new_client
from multiprocessing import Pool, cpu_count, get_context
from multiprocessing.pool import ThreadPool
import pandas as pd
from tqdm import tqdm

# Readout name variants mapping
ADME_READOUT_ALIASES = {
    "LogD": ["LogD", "Log D", "logD", "log D"],
    "mhCLint": ["CLint", "Clint", "Microsomal clearance", "Microsomal CLint", "Intrinsic clearance"],
    "Solubility": ["Solubility", "Aqueous solubility", "Water solubility"],
    "Caco-2 Permeability": [
        "Caco-2 Permeability", "Caco-2", "Apparent permeability", "Papp", "Caco2 Permeability"
    ],
    "hPPB": ["Plasma protein binding", "PPB", "Protein binding", "Plasma Binding"]
}

def fetch_single_alias_task(args):
    alias, readout_key, base_filters = args
    filters = base_filters.copy()
    filters["standard_type"] = alias

    print(f"Fetching for alias: '{alias}' (readout: {readout_key}) with filters: {filters}")
    query = new_client.activity.filter(**filters)

    rows = []
    for item in query:
        rows.append(item)

    if rows:
        df = pd.DataFrame(rows)
        df["readout"] = readout_key
        return df
    else:
        print(f"No records for alias '{alias}'")
        return pd.DataFrame()

@register_task("retrieve_chembl_adme_data")
def retrieve_chembl_adme_data(config, data=None):
    readouts = config.get("readout", ["LogD"])
    relation = config.get("relation")
    assay_type = config.get("assay_type")
    readout_filters = config.get("filters", {})

    if isinstance(readouts, str):
        readouts = [readouts]

    results = []

    for readout_key in tqdm(readouts, desc="Querying ADME Readouts (types)"):
        aliases = ADME_READOUT_ALIASES.get(readout_key, [readout_key])
        base_filters = {
            "assay_type": assay_type,
        }

        if relation:
            base_filters["standard_relation"] = relation

        # Add any user-specified filters for this readout
        base_filters.update(readout_filters.get(readout_key, {}))

        # Prepare alias fetch tasks
        tasks = [(alias, readout_key, base_filters) for alias in aliases]

        # Use multiprocessing to fetch aliases in parallel
        num_workers = min(len(tasks), max(1, cpu_count() - 2))
        print(f"Parallel fetching with {num_workers} workers for '{readout_key}' aliases...")

        # Spawn context to allow nested multiprocessing
        with ThreadPool(num_workers) as pool:    alias_dfs = list(tqdm(pool.imap(fetch_single_alias_task, tasks), total=len(tasks), desc=f"Fetching {readout_key}"))
        readout_df = pd.concat([df for df in alias_dfs if not df.empty], ignore_index=True)

        if readout_df.empty:
            print(f"❌ No data found for readout: {readout_key}")
        else:
            print(f"Retrieved {len(readout_df)} rows for readout: {readout_key}")
            results.append(readout_df)

    if results:
        return {"df": pd.concat(results, ignore_index=True)}
    else:
        return {"df": pd.DataFrame()}
    
# from chembl_webresource_client.new_client import new_client
# import pandas as pd

# # Aliases to check
# aliases = ["LogD", "Log D", "logD", "log D", "LogD7.4", "LogD7_4", "pLogD", "logD (pH 7.4)"]

# results_summary = {}

# for alias in aliases:
#     print(f"\nChecking alias: '{alias}'")
#     query = new_client.activity.filter(standard_type=alias)
    
#     # Fetch the first N items (you can adjust the cap here)
#     try:
#         sample_data = list(query[:100])
#         count = len(sample_data)
#         results_summary[alias] = count
#         print(f" → Found {count} records for standard_type = '{alias}'")

#         if count > 0:
#             df = pd.DataFrame(sample_data)
#             print("Sample data:")
#             print(df[["molecule_chembl_id", "standard_type", "standard_value", "assay_type", "standard_units"]].head())
    
#     except Exception as e:
#         print(f"Error fetching alias '{alias}': {e}")
#         results_summary[alias] = 0

# # Print final summary
# print("\nSummary of found records:")
# for alias, count in results_summary.items():
#     print(f" - {alias!r}: {count} records")
