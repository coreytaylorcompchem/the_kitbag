from pipeline.workflow_registry import register_workflow
from pipeline.task_registry import get_task

from pathlib import Path

@register_workflow(
    'chembl_bioactivity_single_target',
    description="Retrieve and clean ChEMBL bioactivities and compounds."
)
def run_chembl_bioactivity_workflow(config):
    steps = config.get('workflow', [])
    data = None

    for step in steps:
        print(f"\n>>> Running step: {step}")
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found in registry.")
        data = task_func(config, data)

        # Logging shapes
        if isinstance(data, dict) and "df" in data and hasattr(data["df"], "shape"):
            print(f"[{step}] → Output shape: {data['df'].shape}")
        elif hasattr(data, "shape"):
            print(f"[{step}] → Output shape: {data.shape}")
        else:
            print(f"[{step}] → Output data type: {type(data)}")

    # Save final output
    output_cfg = config.get("output", {})
    out_dir = output_cfg.get("directory", "outputs/single_target")
    out_file = output_cfg.get("filename", "output.csv")
    overwrite = output_cfg.get("overwrite", False)

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_path = Path(out_dir) / out_file

    # Determine DataFrame to save
    if isinstance(data, dict) and "df" in data:
        df_to_save = data["df"]
    elif hasattr(data, "to_csv"):
        df_to_save = data
    else:
        print("❌ No DataFrame to save.")
        return data

    if out_path.exists() and not overwrite:
        print(f"❌ File already exists at {out_path} and overwrite is False. Skipping save.")
    else:
        df_to_save.to_csv(out_path, index=False)
        print(f"✅ Saved output to {out_path}")

    return data

