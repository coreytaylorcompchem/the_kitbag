import os
import io
import sys
from pipeline.task_registry import register_task

from pipeline.logger import setup_logger
logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    'mesp',
    description="Calculates Molecular Electrostatic Potential (MESP) and outputs cube files for visualization.",
    modifies_geometry=False,
    category='Property'
)
def run(backend, xyz_or_wfn, step_config, global_config=None):
    global_output_cfg = global_config.get("output", {}) if global_config else {}
    output_dir = global_output_cfg.get("directory", "results")
    overwrite = global_output_cfg.get("overwrite", False)
    os.makedirs(output_dir, exist_ok=True)

    # Determine base name for files
    if hasattr(xyz_or_wfn, 'molecule'):
        base_name = step_config.get("base_name", "calc")
    else:
        base_name = os.path.splitext(os.path.basename(xyz_or_wfn))[0]

    output_cfg = step_config.get("output", {})

    log_filename = output_cfg.get("log", f"{base_name}_mesp.log")
    mesp_cube_filename = output_cfg.get("mesp_cube", f"{base_name}_mesp.cube")
    density_cube_filename = output_cfg.get("density_cube", f"{base_name}_density.cube")

    log_path = os.path.join(output_dir, log_filename)
    mesp_cube_path = os.path.join(output_dir, mesp_cube_filename)
    density_cube_path = os.path.join(output_dir, density_cube_filename)

    if os.path.exists(log_path) and not overwrite:
        logger.warning(f"Skipping – {log_path} already exists.")
        return log_path

    # Safeguard: fallback to last method/basis if available
    if global_config:
        method_fallback = global_config.get("last_method")
        basis_fallback = global_config.get("last_basis")
        if "method" not in step_config and method_fallback:
            step_config["method"] = method_fallback
        if "basis" not in step_config and basis_fallback:
            step_config["basis"] = basis_fallback

    stdout_capture = io.StringIO()
    sys_stdout_backup = sys.stdout
    sys.stdout = stdout_capture

    try:
        backend.mesp(
            xyz_or_wfn,
            step_config,
            mesp_cube_path=mesp_cube_path,
            density_cube_path=density_cube_path,
            log_path=log_path
        )
    finally:
        sys.stdout = sys_stdout_backup

    log_output = stdout_capture.getvalue()

    with open(log_path, "a") as f:
        f.write("\n[Stdout Capture]\n")
        f.write(log_output)
    logger.info(f"Saved log output to: {log_path}")

    return log_path
