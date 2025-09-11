from pipeline.task_registry import register_task
from chembl_webresource_client.new_client import new_client
import pandas as pd
from tqdm import tqdm

# Readout name variants mapping
ADME_READOUT_ALIASES = {
    "LogD": ["LogD7.4", "LogD", "Log D", "logD", "log D"],
    "CLint": ["CLint", "Clint", "Microsomal clearance", "Microsomal CLint", "Intrinsic clearance"],
    "Solubility": ["Solubility", "Aqueous solubility", "Water solubility"],
    "Caco-2 Permeability": [
        "Caco-2 Permeability", "Caco-2", "Apparent permeability", "Papp", "Caco2 Permeability"
    ],
    "PPB": ["Plasma protein binding", "PPB", "Protein binding", "Plasma Binding"]
}

@register_task("retrieve_chembl_adme_data")
def retrieve_chembl_adme_data(config, data=None):
    readouts = config.get("readout", ["LogD"])
    relation = config.get("relation")
    assay_type = config.get("assay_type", "A")
    readout_filters = config.get("filters", {}) 

    if isinstance(readouts, str):
        readouts = [readouts]

    activity = new_client.activity
    all_dfs = []

    for key in tqdm(readouts, desc="Querying ADME Readouts (types)"):
        aliases = ADME_READOUT_ALIASES.get(key, [key])
        combined_df = pd.DataFrame()

        for alias in aliases:
            # Base filters
            filters = {
                "standard_type": alias,
                "assay_type": assay_type
            }

            if relation:
                filters["standard_relation"] = relation

            # 🔍 Merge in user-defined YAML filters for this readout (e.g. CLint → Homo sapiens)
            filters.update(readout_filters.get(key, {}))

            print(f"🔍 Fetching for alias: '{alias}' (readout: {key}) with filters: {filters}")
            query = activity.filter(**filters)

            rows = []
            for item in tqdm(query, desc=f" ↳ Records for '{alias}'", leave=False):
                rows.append(item)

            if rows:
                df = pd.DataFrame(rows)
                df["readout"] = key
                combined_df = pd.concat([combined_df, df], ignore_index=True)

        if combined_df.empty:
            print(f"❌ No data found for readout: {key}")
        else:
            print(f"✅ Retrieved {len(combined_df)} rows for readout: {key}")
            all_dfs.append(combined_df)

    if all_dfs:
        final = pd.concat(all_dfs, ignore_index=True)
        return {"df": final}
    else:
        return {"df": pd.DataFrame()}
