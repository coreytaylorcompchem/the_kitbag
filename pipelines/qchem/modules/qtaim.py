import os
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from backends.utils.wfn_export import write_wfn_file

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    'qtaim',
    description="Performs QTAIM analysis on .molden formatted file.",
    modifies_geometry=False,
    category='Wave function analysis'
)
def run(backend, xyz_file_or_wfn, step_config, global_config=None):
    global_output_cfg = global_config.get("output", {}) if global_config else {}
    output_dir = global_output_cfg.get("directory", "outputs/qtaim")
    overwrite = global_output_cfg.get("overwrite", False)
    os.makedirs(output_dir, exist_ok=True)

    if isinstance(xyz_file_or_wfn, str):
        base_name = os.path.splitext(os.path.basename(xyz_file_or_wfn))[0]
    else:
        base_name = step_config.get("name", "wavefunction")

    wfn_path = os.path.join(output_dir, f"{base_name}.molden")
    log_path = os.path.join(output_dir, f"{base_name}_qtaim.log")

    if os.path.exists(log_path) and not overwrite:
        logger.info(f"Skipping QTAIM - {log_path} already exists.")
        return log_path

    # Handle input wavefunction
    if hasattr(xyz_file_or_wfn, 'molecule'):
        logger.debug("Received wavefunction object directly, writing to file.")
        write_wfn_file(xyz_file_or_wfn, wfn_path, os.path.basename(wfn_path))
    elif isinstance(xyz_file_or_wfn, str) and xyz_file_or_wfn.endswith(('.wfn', '.molden')):
        logger.debug(f"Using provided wavefunction file: {xyz_file_or_wfn}")
        wfn_path = xyz_file_or_wfn
    else:
        raise ValueError("QTAIM task requires a .wfn or .molden file or wavefunction object.")

    # Prepare Multiwfn script as a string
    script = "\n".join([
        "2",    # Topology analysis
        "2",    # Search CPs from nuclear positions
        "3",    # Search CPs from midpoint of atomic pairs
        "4",    # Search CPs from triangle center of three atoms
        "5",    # Search CPs from pyramid center of four atoms
        "8",    # Generate paths connecting (3,-3) and (3,-1) CPs
        "-4",   # Modify/export CPs
        "6",    # Export CPs as CPs.pdb
        "0",    # Return
        "-5",   # Modify/export paths
        "6",    # Export paths as paths.pdb
        "0",    # Return
        "7",    # Show real space values at CPs
        "0",    # Print all to CPProp.txt
        "-10",  # Return to main menu
        "q"     # Quit
    ]) + "\n"

    logger.info(f"Running Multiwfn on: {wfn_path}")
    backend.run_multiwfn(
        input_file=wfn_path,
        output_dir=output_dir,
        ncpu=step_config.get("ncpu", 1),
        log_path=log_path,
        script_input=script
    )

    return {
        "wfn_file": wfn_path,
        "log": log_path
    }
