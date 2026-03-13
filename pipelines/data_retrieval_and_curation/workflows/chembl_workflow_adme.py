import os
import pandas as pd
from copy import deepcopy
from pipeline.task_registry import get_task
from workflows import register_workflow

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)


@register_workflow(
    "chembl_adme_data",
    description="Retrieve and construct multi-task ADME dataset from ChEMBL."
)
def run_chembl_mtl_workflow(config):
    local_config = deepcopy(config)
    data = None

    for step in local_config.get("workflow", []):
        task_func = get_task(step)

        if not task_func:
            raise ValueError(f"Task '{step}' not found in registry.")

        logger.info(f"Running task: {step}")

        data = task_func(local_config, data)

    if isinstance(data, dict) and "df" in data:
        df = data["df"]
    elif isinstance(data, pd.DataFrame):
        df = data
    else:
        raise ValueError("Final workflow output did not return a dataframe.")

    output_cfg = local_config.get("output", {})
    out_dir = output_cfg.get("directory", "outputs/adme")
    filename = output_cfg.get("filename", "chembl_mtl_dataset.csv")

    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, filename)
    df.to_csv(output_path, index=False)

    logger.info(f"Saved dataset to: {output_path}")

    return {"df": df, "path": output_path}