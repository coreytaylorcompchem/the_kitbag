import os
import shutil
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from Bio.PDB import PDBParser, PDBIO, Select
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class ChainALigandSelect(Select):
    def __init__(self, remove_waters):
        self.remove_waters = remove_waters

    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        hetfield, _, _ = residue.id
        is_ligand = hetfield != ' ' and hetfield != 'W'
        is_water = residue.get_resname() == "HOH" or hetfield == 'W'
        return is_ligand or (residue.parent.id == "A" and not (self.remove_waters and is_water))

def process_structure(pdb_path, output_dir, remove_waters):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    try:
        structure_id = os.path.basename(pdb_path).split('.')[0].lower().replace("pdb", "")
        structure = parser.get_structure(structure_id, pdb_path)
        io.set_structure(structure)

        out_path = os.path.join(output_dir, f"{structure_id}_std.pdb")
        io.save(out_path, select=ChainALigandSelect(remove_waters))
        logger.debug(f"Standardised: {pdb_path} → {out_path}")
        return out_path
    except Exception as e:
        logger.warning(f"Failed to standardise {pdb_path}: {e}")
        return None

@register_task("standardise_pdbs", category="PDB", description="Standardise PDBs (Chain A + ligand, remove solvent).")
def standardise_pdbs(config, pdb_paths):
    output_dir = os.path.join(config.get("output", {}).get("directory", "outputs/pdb_standardised"))
    os.makedirs(output_dir, exist_ok=True)

    remove_waters = config.get("remove_waters", True)
    output_paths = []

    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(process_structure, path, output_dir, remove_waters): path
            for path in pdb_paths
        }
        for future in as_completed(futures):
            result = future.result()
            if result:
                output_paths.append(result)

    logger.info(f"{len(output_paths)} PDBs standardised successfully.")
    return output_paths