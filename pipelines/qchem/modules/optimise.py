import os
import io
import sys
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    'optimise',
    description="Performs geometry optimization.",
    modifies_geometry=True,
    category='Optimization'
)
def run(backend, xyz_file, step_config, global_config=None):
    # Load global output settings
    global_output_cfg = global_config.get("output", {}) if global_config else {}
    output_dir = global_output_cfg.get("directory", "results")
    overwrite = global_output_cfg.get("overwrite", False)
    os.makedirs(output_dir, exist_ok=True)

    # Get base name for fallback filenames
    base_name = os.path.splitext(os.path.basename(xyz_file))[0]
    output_cfg = step_config.get("output", {})

    log_filename = output_cfg.get("log", f"{base_name}_opt.log")
    energy_filename = output_cfg.get("energy", f"{base_name}_opt_energy.txt")
    geometry_filename = output_cfg.get("geometry", f"{base_name}_opt.xyz")
    wfn_filename = output_cfg.get("wfn", f"{base_name}.wfn")

    log_path = os.path.join(output_dir, log_filename)
    energy_path = os.path.join(output_dir, energy_filename)
    geometry_path = os.path.join(output_dir, geometry_filename)
    wfn_path = os.path.join(output_dir, wfn_filename)

    # Checkpoint: skip if log exists and overwrite is false
    if os.path.exists(log_path) and not overwrite:
        logger.warning(f"Skipping - {log_path} already exists.")
        result_dict = {"wfn_file": wfn_path if os.path.exists(wfn_path) else None,
                       "geometry": None,
                       "energy": None,
                       "wfn": None}
        return result_dict

    # Capture stdout from Psi4 or backend output
    stdout_capture = io.StringIO()
    sys_stdout_backup = sys.stdout
    sys.stdout = stdout_capture

    try:
        # Pass the full log_path to backend optimise call if supported
        # Adjust this call if backend.optimise() accepts log_path parameter
        # At present, Orca and Psi4 does.
        result = backend.optimise(xyz_file, step_config, log_path=log_path)
    finally:
        sys.stdout = sys_stdout_backup

    log_output = stdout_capture.getvalue()

    # Save the stdout capture (complementary to backend log file)
    with open(log_path, "a") as f:
        f.write("\n[Stdout Capture]\n")
        f.write(log_output)
    logger.info(f"Saved log output to: {log_path}")

    energy = None
    final_geom = None
    wfn = None

    if isinstance(result, dict):
        energy = result.get("energy")
        wfn = result.get("wfn")
        if wfn:
            try:
                wfn.to_file(wfn_path)
                logger.info(f"Saved wavefunction file to: {wfn_path}")
            except Exception as e:
                logger.warning(f"Failed to write .wfn file: {e}")

        if wfn:
            final_geom = wfn.molecule().save_string_xyz()
    else:
        energy = result

    if energy is not None:
        energy_value = energy["energy"] if isinstance(energy, dict) and "energy" in energy else energy
        with open(energy_path, "w") as f:
            f.write(f"Energy: {energy_value:.10f} Hartree\n")
        logger.info(f"Saved energy to: {energy_path}")

    if final_geom and isinstance(final_geom, str):
        with open(geometry_path, "w") as f:
            f.write(final_geom)
        logger.info(f"Saved geometry to: {geometry_path}")

    # Return dict for downstream tasks
    return {
        "energy": energy,
        "method": step_config.get("method"),
        "basis": step_config.get("basis"),
        "wfn_file": wfn_path if wfn else None,
        "wfn": wfn,
        "geometry": final_geom
    }
