import os
import io
import sys
import re
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def print_progress(iteration, width=30):
    """
    Very simple text progress bar (in-place).
    """
    filled = int(width * iteration / 100)
    bar = "[" + "#" * filled + "-" * (width - filled) + "]"
    sys.stdout.write(f"\r{bar} {iteration}%")
    sys.stdout.flush()


def ensure_job_directory(base_output_dir, base_name):
    """
    Makes a subdirectory inside the output directory for each optimisation.
    Returns the path to it.
    """
    job_dir = os.path.join(base_output_dir, base_name)
    os.makedirs(job_dir, exist_ok=True)
    return job_dir


@register_task(
    'optimise',
    description="Performs geometry optimization.",
    modifies_geometry=True,
    category='Optimization'
)
def run(backend, xyz_file, step_config, global_config=None):

    # Load global output settings
    global_output_cfg = global_config.get("output", {}) if global_config else {}
    base_output_dir = global_output_cfg.get("directory", "results")
    overwrite = global_output_cfg.get("overwrite", False)
    os.makedirs(base_output_dir, exist_ok=True)

    # Extract the optimisation name (without .xyz)
    base_name = os.path.splitext(os.path.basename(xyz_file))[0]

    # Make a unique directory for this calculation
    job_dir = ensure_job_directory(base_output_dir, base_name)

    # Filenames inside job_dir
    output_cfg = step_config.get("output", {})

    log_filename = output_cfg.get("log", f"{base_name}.log")
    energy_filename = output_cfg.get("energy", f"{base_name}_energy.txt")
    geometry_filename = output_cfg.get("geometry", f"{base_name}_opt.xyz")
    wfn_filename = output_cfg.get("wfn", f"{base_name}.wfn")

    log_path = os.path.join(job_dir, log_filename)
    energy_path = os.path.join(job_dir, energy_filename)
    geometry_path = os.path.join(job_dir, geometry_filename)
    wfn_path = os.path.join(job_dir, wfn_filename)

    # Skip if exists and not overwriting
    if os.path.exists(log_path) and not overwrite:
        logger.warning(f"Skipping - {log_path} already exists.")
        result_dict = {
            "wfn_file": wfn_path if os.path.exists(wfn_path) else None,
            "geometry": None,
            "energy": None,
            "wfn": None
        }
        return result_dict

    # Capture backend stdout (ORCA, Psi4, etc.)
    stdout_capture = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = stdout_capture

    try:
        # Run backend optimisation
        result = backend.optimise(
            xyz_file,
            step_config,
            log_path=log_path   # ORCA/Psi4 should print here
        )
    finally:
        sys.stdout = old_stdout

    # The captured output may contain iteration info
    captured_text = stdout_capture.getvalue()

    # --- PROGRESS BAR PARSING --------------------------------------------
    # Try to count iterations from backend output
    iterations = 0

    # ORCA emits "GEOMETRY OPTIMIZATION CYCLE N"
    orca_cycles = re.findall(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)", captured_text)
    if orca_cycles:
        iterations = max(map(int, orca_cycles))

    # Psi4 emits "@OptKing Iteration"
    psi4_cycles = re.findall(r"@OptKing\s+Iteration\s+(\d+)", captured_text)
    if psi4_cycles:
        iterations = max(iterations, max(map(int, psi4_cycles)))

    # If we didn't find cycles, assume a single step
    if iterations == 0:
        iterations = 1

    # Print final progress bar to 100%
    print_progress(100)
    print("")  # newline

    # Append captured stdout to the log file
    with open(log_path, "a") as f:
        f.write("\n[Stdout Capture]\n")
        f.write(captured_text)
    logger.info(f"Saved log output to: {log_path}")

    # Extract ORCA/Psi4 result from backend
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
        # Backend returned only energy
        energy = result

    # Save energy
    if energy is not None:
        val = energy["energy"] if isinstance(energy, dict) and "energy" in energy else energy
        with open(energy_path, "w") as f:
            f.write(f"Energy: {val:.10f} Hartree\n")
        logger.info(f"Saved energy: {energy_path}")

    # Save geometry
    if final_geom:
        with open(geometry_path, "w") as f:
            f.write(final_geom)
        logger.info(f"Saved geometry: {geometry_path}")

    # Return dict for downstream tasks
    return {
        "energy": energy,
        "method": step_config.get("method"),
        "basis": step_config.get("basis"),
        "wfn_file": wfn_path if wfn else None,
        "wfn": wfn,
        "geometry": final_geom
    }
