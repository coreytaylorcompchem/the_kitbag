import os
from task_registry import register_task, get_task_metadata

@register_task(
    'optimise',
    description="Performs geometry optimization using the selected backend.",
    modifies_geometry=True,
)
def run(backend, xyz_file, config):
    output_dir = config.get("output", {}).get("directory", "results")
    os.makedirs(output_dir, exist_ok=True)

    base_name = os.path.splitext(os.path.basename(xyz_file))[0]
    default_filename = f"{base_name}_optimised.xyz"
    custom_filenames = config.get("output", {}).get("filenames", {})
    filename = custom_filenames.get("optimise", default_filename)

    output_path = os.path.join(output_dir, filename)

    # Checkpoint: does the file already exist? If so, skip.
    overwrite = config.get("output", {}).get("overwrite", False)
    if os.path.exists(output_path) and not overwrite:
        print(f"[optimise] Skipping – {output_path} already exists.")
        return output_path

    # run the opt with the specified backend
    optimised_xyz = backend.optimise(xyz_file, config)

    with open(output_path, "w") as f:
        f.write(optimised_xyz)

    print(f"[optimise] Saved optimised geometry to: {output_path}")
    return output_path