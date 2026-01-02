import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

from pathlib import Path
from tqdm import tqdm

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def run_fpocket_on_structure(pdb_file, input_dir, output_dir):
    """
    Run fpocket on a PDB structure file.
    The output pocket descriptors are captured and saved in a CSV file.
    """
    try:
        # Annoying workaround: fpocket does not allow us to control the output dir 
        # So we need to symlink the pdb where we want to dump data (output dir from yaml).

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
        logger.debug(f"Symbolically linked {pdb_file} to {linked_pdb}")

        # Run fpocket on the symbolic link and output descriptors to a file
        output_prefix_dir = output_dir / f"{pdb_name}_out"  # Subdirectory -> pdb name
        output_prefix_dir.mkdir(parents=True, exist_ok=True)

        # Command to run fpocket with -d flag to output pocket descriptors to stdout
        pocket_descriptor_file = output_prefix_dir / f"{pdb_name}_pocket_descriptors.csv"
        command = f"fpocket -f {linked_pdb} -o {output_prefix_dir} -d > {pocket_descriptor_file}"
        subprocess.run(command, shell=True, check=True)

        # Check if pocket descriptors were successfully generated
        if pocket_descriptor_file.exists():
            logger.debug(f"Successfully ran fpocket on {pdb_file.name}.")
            logger.debug(f"Pocket descriptors saved to {pocket_descriptor_file}")
        else:
            logger.error(f"Fpocket did not generate a valid pocket descriptor file for {pdb_file.name}")
            return None

        return pocket_descriptor_file

    except Exception as e:
        logger.error(f"Error running fpocket on {pdb_file.name}: {e}")
        return None


def run_fpocket_on_gpcr_structures(
    gpcr_pdb_dir,
    output_dir,
    max_pockets_per_structure,
    n_jobs=1
):
    pdb_files = sorted(Path(gpcr_pdb_dir).glob("*.pdb"))

    pocket_files = []

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = {
            executor.submit(
                run_fpocket_on_structure,
                pdb_file,
                gpcr_pdb_dir,
                output_dir
            ): pdb_file
            for pdb_file in pdb_files
        }

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Running fpocket",
            unit="structure"
        ):
            pdb_file = futures[future]

            try:
                pocket_file = future.result()
                if pocket_file:
                    pocket_files.append(pocket_file)
            except Exception as e:
                logger.error(
                    f"Fpocket failed for {pdb_file.name}: {e}"
                )

    return pocket_files