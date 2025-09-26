import os
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import Superimposer

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def get_chain_a_ca_atoms(structure):
    atoms = []
    try:
        for model in structure:
            for chain in model:
                if chain.id == "A":
                    for residue in chain:
                        if "CA" in residue:
                            atoms.append(residue["CA"])
                    break
            break  # only first model
    except Exception as e:
        logger.warning(f"Failed to extract CA atoms: {e}")
    return atoms

@register_task("align_structures", category="PDB", description="Align all PDBs to the first retrieved structure.")
def align_structures(config, pdb_paths):
    if not pdb_paths or len(pdb_paths) < 2:
        logger.warning("Not enough structures to align.")
        return pdb_paths

    output_dir = os.path.join(config.get("output", {}).get("directory", "outputs/aligned"))
    os.makedirs(output_dir, exist_ok=True)

    parser = PDBParser(QUIET=True)
    io = PDBIO()

    # Reference structure
    ref_path = pdb_paths[0]
    ref_structure = parser.get_structure("ref", ref_path)
    ref_atoms = get_chain_a_ca_atoms(ref_structure)

    aligned_paths = [ref_path]  # reference goes first

    for path in pdb_paths[1:]:
        try:
            moving_structure = parser.get_structure("mov", path)
            mov_atoms = get_chain_a_ca_atoms(moving_structure)

            if len(ref_atoms) != len(mov_atoms) or len(ref_atoms) < 3:
                logger.warning(f"Skipping alignment (CA mismatch): {path}")
                continue

            sup = Superimposer()
            sup.set_atoms(ref_atoms, mov_atoms)
            sup.apply(moving_structure.get_atoms()) 

            aligned_path = os.path.join(output_dir, os.path.basename(path).replace(".pdb", "_aligned.pdb"))
            io.set_structure(moving_structure)
            io.save(aligned_path)
            aligned_paths.append(aligned_path)

            logger.info(f"Aligned: {path} → {aligned_path}")

        except Exception as e:
            logger.warning(f"Alignment failed for {path}: {e}")

    return aligned_paths
