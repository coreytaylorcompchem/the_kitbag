import os
import pandas as pd
from pipeline.task_registry import get_task
from pipeline.workflow_registry import register_workflow
from tqdm import tqdm

@register_workflow("chembl_adme_data", description="Retrieve ADME data from ChEMBL")
def run_chembl_adme_workflow(config):
    output_cfg = config.get("output", {})
    output_dir = output_cfg.get("directory", "outputs/adme_data")
    os.makedirs(output_dir, exist_ok=True)

    data = None
    for step in config.get("workflow", []):
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found.")

        data = task_func(config, data)

    if isinstance(data, dict) and "df" in data:
        df = data["df"]
    else:
        df = data

    if df is not None and not df.empty:
        file_path = os.path.join(output_dir, output_cfg.get("filename", "combined_adme_data.csv"))
        df.to_csv(file_path, index=False)
        print(f"✅ ADME data saved to {file_path}")
    else:
        print("⚠️ No ADME data retrieved.")

    return df