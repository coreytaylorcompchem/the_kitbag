import os
import io
import sys
from pipeline.task_registry import register_task

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    'single_point',
    description="Performs a single-point energy calculation using the selected backend.",
    modifies_geometry=False,
    category='Energy'
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

    log_filename = output_cfg.get("log", f"{base_name}_spe.log")
    energy_filename = output_cfg.get("energy", f"{base_name}_energy.txt")
    geometry_filename = output_cfg.get("geometry", f"{base_name}_spe.xyz")

    log_path = os.path.join(output_dir, log_filename)
    energy_path = os.path.join(output_dir, energy_filename)
    geometry_path = os.path.join(output_dir, geometry_filename)

    # Checkpoint: skip if log exists and overwrite is false
    if os.path.exists(log_path) and not overwrite:
        logger.warning(f"Skipping – {log_path} already exists.")
        return log_path

    # Capture output printed to stdout
    stdout_capture = io.StringIO()
    sys_stdout_backup = sys.stdout
    sys.stdout = stdout_capture

    try:
        # Pass the full log_path to the backend
        result = backend.single_point(xyz_file, step_config, log_path=log_path)
    finally:
        sys.stdout = sys_stdout_backup

    log_output = stdout_capture.getvalue()

    # Save full stdout log
    with open(log_path, "a") as f:
        f.write("\n[Stdout Capture]\n")
        f.write(log_output)
    logger.info(f"Saved log output to: {log_path}")

    # Save energy and geometry if returned
    if isinstance(result, tuple):
        energy, final_geom = result
    else:
        energy = result
        final_geom = None

    if energy is not None:
        with open(energy_path, "w") as f:
            f.write(f"Energy: {energy['energy']:.10f} Hartree\n")
        logger.info(f"Saved energy to: {energy_path}")

    if final_geom and isinstance(final_geom, str):
        with open(geometry_path, "w") as f:
            f.write(final_geom)
        logger.info(f"Saved geometry to: {geometry_path}")

    return log_path
