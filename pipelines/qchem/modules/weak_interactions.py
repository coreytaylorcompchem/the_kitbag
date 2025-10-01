import os
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from backends.utils.wfn_export import write_wfn_file

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    'nci',
    description="Performs NCI calculation, outputs cube files.",
    modifies_geometry=False,
    category='Wave function analysis'
)
def run(backend, xyz_file_or_wfn, step_config, global_config=None):
    global_output_cfg = global_config.get("output", {}) if global_config else {}
    output_dir = global_output_cfg.get("directory", "outputs/nci")
    overwrite = global_output_cfg.get("overwrite", False)
    os.makedirs(output_dir, exist_ok=True)

    if isinstance(xyz_file_or_wfn, str):
        base_name = os.path.splitext(os.path.basename(xyz_file_or_wfn))[0]
    else:
        base_name = step_config.get("name", "wavefunction")

    wfn_path = os.path.join(output_dir, f"{base_name}.molden")
    log_path = os.path.join(output_dir, f"{base_name}_nci.log")

    if os.path.exists(log_path) and not overwrite:
        logger.info(f"Skipping NCI analysis - {log_path} already exists.")
        return log_path

    # Handle input wavefunction
    if hasattr(xyz_file_or_wfn, 'molecule'):
        logger.debug("Received wavefunction object directly, writing to file.")
        write_wfn_file(xyz_file_or_wfn, wfn_path, os.path.basename(wfn_path))
    elif isinstance(xyz_file_or_wfn, str) and xyz_file_or_wfn.endswith(('.wfn', '.molden')):
        logger.debug(f"Using provided wavefunction file: {xyz_file_or_wfn}")
        wfn_path = xyz_file_or_wfn
    else:
        raise ValueError("NCI analysis task requires a .wfn or .molden file or wavefunction object.")

    # Prepare Multiwfn script as a string
    script = "\n".join([
        "20",   # Visual study of weak interaction
        "1",    # NCI analysis
        "3",    # High quality grid
        "1",    # Save the scatter graph to file
        "3",    # Output cube files to func1.cub and func2.cub in current folder
        "0",    # Back one menu
        "0",    # Back one menu
        "q", # Quit
    ]) + "\n"

    logger.info(f"Running Multiwfn on: {wfn_path}")
    backend.run_multiwfn(
        input_file=wfn_path,
        output_dir=output_dir,
        ncpu=step_config.get("ncpu", 1),
        log_path=log_path,
        script_input=script
    )

    return log_path
