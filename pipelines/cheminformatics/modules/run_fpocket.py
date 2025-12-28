import os
import subprocess
import pandas as pd
from pathlib import Path
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def run_fpocket_on_structure(pdb_file, input_dir, output_dir):
    """
    Run fpocket on a PDB structure file after creating a symbolic link to the output directory.
    The output pocket descriptors are captured and saved in a CSV file.
    """
    try:
        # Ensure input and output directories exist
        input_dir = Path(input_dir).resolve()  # Absolute path for input
        output_dir = Path(output_dir).resolve()  # Absolute path for output

        # Create the symbolic link in the output directory pointing to the pdb file
        pdb_name = pdb_file.stem  # PDB file name without extension
        linked_pdb = output_dir / f"{pdb_name}.pdb"

        # Remove any existing symlink with the same name
        if linked_pdb.exists():
            linked_pdb.unlink()

        # Create the symbolic link using absolute path
        os.symlink(pdb_file.resolve(), linked_pdb)  # Correct symlink using absolute path
        logger.info(f"Symbolically linked {pdb_file} to {linked_pdb}")

        # Run fpocket on the symbolic link and output descriptors to a file
        output_prefix_dir = output_dir / f"{pdb_name}_out"  # Subdirectory based on pdb name
        output_prefix_dir.mkdir(parents=True, exist_ok=True)

        # Command to run fpocket with -d flag to output pocket descriptors to stdout
        pocket_descriptor_file = output_prefix_dir / f"{pdb_name}_pocket_descriptors.csv"
        command = f"fpocket -f {linked_pdb} -o {output_prefix_dir} -d > {pocket_descriptor_file}"
        subprocess.run(command, shell=True, check=True)

        # Check if pocket descriptors were successfully generated
        if pocket_descriptor_file.exists():
            logger.info(f"Successfully ran fpocket on {pdb_file.name}. Pocket descriptors saved to {pocket_descriptor_file}")
        else:
            logger.error(f"Fpocket did not generate a valid pocket descriptor file for {pdb_file.name}")
            return None

        return pocket_descriptor_file

    except Exception as e:
        logger.error(f"Error running fpocket on {pdb_file.name}: {e}")
        return None


def run_fpocket_on_gpcr_structures(gpcr_pdb_dir, output_dir, max_pockets_per_structure):
    """
    Run fpocket on multiple GPCR structures (PDB files) from the input directory.
    It will create symbolic links and run fpocket, and then save the results to the output directory.
    """
    pdb_files = list(Path(gpcr_pdb_dir).glob("*.pdb"))  # Get all pdb files in the input directory

    # Run fpocket on each PDB file in the directory
    pocket_files = []
    for pdb_file in pdb_files:
        # Make sure output_dir is passed to the run_fpocket_on_structure function
        pocket_file = run_fpocket_on_structure(pdb_file, gpcr_pdb_dir, output_dir)  # Ensure the output_dir is passed here
        
        if pocket_file:
            pocket_files.append(pocket_file)
    
    return pocket_files


@register_task("run_fpocket", category="Pocket detection", description="Run fpocket on protein structures to detect binding pockets.")
def run_fpocket(config: dict, data: dict = None) -> dict:
    """
    The main task to run fpocket on protein structures.
    """
    params = config.get("run_fpocket", {})

    gpcr_pdb_dir = Path(params.get("input_file_gpcr_pdbs"))
    output_dir = Path(params.get("output_directory"))
    max_pockets_per_structure = params.get("max_pockets_per_structure", 100)

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Running fpocket on GPCR structures in {gpcr_pdb_dir}...")
    
    # Run fpocket on all GPCR PDBs
    pocket_files = run_fpocket_on_gpcr_structures(gpcr_pdb_dir, output_dir, max_pockets_per_structure)
    
    if not pocket_files:
        logger.error("No pocket files were generated.")
        return {}

    # Combine all pocket CSVs into one
    all_pockets_df = pd.concat([pd.read_csv(f) for f in pocket_files], ignore_index=True)
    
    # Save combined pocket CSV
    combined_pocket_file = output_dir / "pockets.csv"
    all_pockets_df.to_csv(combined_pocket_file, index=False)
    logger.info(f"Combined pocket data saved to {combined_pocket_file}")

    return {"pockets_file": str(combined_pocket_file)}
