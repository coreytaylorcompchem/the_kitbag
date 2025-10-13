import os
import pandas as pd

from workflows import register_workflow
from pipeline.task_registry import get_task
from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def run_single_target(config):
    uniprot_id = config.get("uniprot_id")
    base_output_dir = config.get("output", {}).get("directory", "outputs/multi_target_pdb")
    target_output_dir = os.path.join(base_output_dir, uniprot_id)
    os.makedirs(target_output_dir, exist_ok=True)
    config["output"]["directory"] = target_output_dir

    steps = config.get("workflow", [])
    data = None
    for step in steps:
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found.")
        data = task_func(config, data)

    pdb_paths = sorted([
        os.path.join(target_output_dir, f)
        for f in os.listdir(target_output_dir)
        if f.endswith("_std.pdb")
    ])

    aligned_paths = data if isinstance(data, list) else []

    df_summary = generate_target_summary(uniprot_id, pdb_paths, aligned_paths, output_dir=target_output_dir)

    return {
        "uniprot_id": uniprot_id,
        "df": df_summary
    }

def generate_target_summary(uniprot_id, pdb_paths, aligned_paths, output_dir):
    summary = {
        "uniprot_id": uniprot_id,
        "num_structures": len(pdb_paths),
        "num_aligned": len(aligned_paths),
        "aligned_files": ";".join(os.path.basename(p) for p in aligned_paths),
        "output_dir": output_dir
    }
    return pd.DataFrame([summary])

@register_workflow("pdb_multi_target", description="Process PDB structures for multiple UniProt targets.")
def run_pdb_multi_target(config):

    # --- NEW: Optional InterPro pre-analysis stage ---
    if config.get("interpro_id"):
        logger.info("🔍 Running InterPro analysis to determine UniProt targets...")
        analysis_task = get_task("analyse_pdbs_by_interpro_accession")
        if not analysis_task:
            raise ValueError("Task 'analyse_pdbs_by_interpro_accession' not found.")
        data = analysis_task(config)
        if not data or "uniprots" not in data:
            raise ValueError("No UniProt IDs retrieved from InterPro analysis.")
        config["uniprot_ids"] = data["uniprots"]
        logger.info(f"✅ Retrieved {len(data['uniprots'])} UniProt targets from InterPro analysis.")
    # -------------------------------------------------
    
    if not config.get("uniprot_ids"):
        raise ValueError("❌ No UniProt IDs provided in config under 'uniprot_ids'")

    runner = ParallelWorkflowRunner(
        workflow_func=run_single_target,
        config=config,
        input_key="uniprot_ids",
        output_key="uniprot_id",
        filename_pattern="{uniprot_id}_aligned.csv",
        output_dir=config.get("output", {}).get("directory", "outputs/multi_target_pdb")
    )

    return runner.run()
