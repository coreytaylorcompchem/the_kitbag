import re
import subprocess
import shutil
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.neighbors import KDTree

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from modules.docking_tasks import get_protein_preparer

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


# ================================================================
# Utility: SDF → PDB
# ================================================================
def convert_sdf_to_pdb(sdf_file: Path, output_dir: Path) -> Path:
    """
    Convert an SDF to a PDB using OpenBabel.
    """
    pdb_file = output_dir / f"{sdf_file.stem}.pdb"

    result = subprocess.run(
        ["obabel", str(sdf_file), "-O", str(pdb_file)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        logger.error(f"Error converting SDF to PDB: {result.stderr}")
        raise RuntimeError(f"Failed to convert {sdf_file} to PDB.")

    logger.debug(f"Converted {sdf_file} → {pdb_file}")
    return pdb_file


# ================================================================
# Charge detection via SMILES generation
# ================================================================
def get_ligand_charge_from_smiles(pdb_path: Path) -> int:
    """
    Get ligand formal charge using OpenBabel protonation at pH 7.4.
    Extracts charges from bracketed SMILES atoms.
    """
    result = subprocess.run(
        ["obabel", str(pdb_path), "-osmi", "-p", "7.4"],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"Open Babel SMILES conversion failed:\n{result.stderr}")

    smiles = result.stdout.strip()
    if not smiles:
        raise RuntimeError("SMILES conversion returned empty output.")

    charges = 0
    for token in re.findall(r'\[.*?\]', smiles):
        for sign, magnitude in re.findall(r'([+-])(\d*)', token):
            mag = int(magnitude) if magnitude else 1
            charges += mag if sign == '+' else -mag

    return charges


# ================================================================
# RDKit ligand reconstruction (critical part)
# ================================================================
def rdkit_reconstruct_ligand(docked_sdf: Path, output_dir: Path) -> Path:
    """
    Simply load docked SDF and write to PDB.
    This preserves coordinates, atom order, and protonation.
    """
    output_dir.mkdir(exist_ok=True, parents=True)
    corrected_pdb = output_dir / f"{docked_sdf.stem}_corrected.pdb"

    mol = Chem.MolFromMolFile(str(docked_sdf), removeHs=False)
    if mol is None:
        raise RuntimeError(f"Cannot read docked SDF: {docked_sdf}")

    Chem.MolToPDBFile(mol, str(corrected_pdb))
    logger.debug(f"Copied docked SDF to PDB → {corrected_pdb}")

    return corrected_pdb


# ================================================================
# PDB → MOL2 via OpenBabel
# ================================================================
def pdb_to_mol2_with_obabel(pdb_file: Path, mol2_file: Path):
    """
    Convert PDB to GAFF-compatible MOL2 using OpenBabel.
    """
    result = subprocess.run(
        ["obabel", str(pdb_file), "-O", str(mol2_file), "-p", "7.4"],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"OpenBabel PDB→MOL2 failed:\n{result.stderr}")


# ================================================================
# Antechamber wrapper
# ================================================================
def parameterize_ligand(docked_sdf: Path, output_dir: Path):
    """
    Fully robust antechamber + parmchk2 ligand parameterization.
    """

    output_dir.mkdir(parents=True, exist_ok=True)

    ligand_name = docked_sdf.stem
    final_mol2 = output_dir / f"{ligand_name}.mol2"
    final_frcmod = output_dir / f"{ligand_name}.frcmod"

    # --- Step 1: RDKit SDF -> PDB (preserve coords / names) ---
    mol = Chem.MolFromMolFile(str(docked_sdf), removeHs=False)
    if mol is None:
        raise RuntimeError(f"Failed to read SDF: {docked_sdf}")

    rdkit_pdb = output_dir / f"{ligand_name}_rdkit.pdb"
    Chem.MolToPDBFile(mol, str(rdkit_pdb))

    # --- Step 2: Charge detection ---
    net_charge = get_ligand_charge_from_smiles(rdkit_pdb)

    # --- Step 3: Create tmpdir and stage input inside it ---
    tmpdir = output_dir / f"{ligand_name}_ante_tmp"
    tmpdir.mkdir(exist_ok=True)

    ante_input_pdb = tmpdir / "input.pdb"
    ante_output_mol2 = tmpdir / "ante.mol2"
    ante_output_frcmod = tmpdir / "ante.frcmod"

    shutil.copy(rdkit_pdb, ante_input_pdb)

    # --- Step 4: antechamber (NO PATHS — ONLY FILENAMES) ---
    cmd_ante = [
        "antechamber",
        "-i", "input.pdb",
        "-fi", "pdb",
        "-o", "ante.mol2",
        "-fo", "mol2",
        "-c", "bcc",
        "-nc", str(net_charge),
        "-at", "gaff2",
        "-j", "5",
        "-s", "2"
    ]

    result = subprocess.run(cmd_ante, cwd=tmpdir, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Antechamber failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    # --- Step 5: parmchk2 ---
    cmd_parmchk = [
        "parmchk2",
        "-i", "ante.mol2",
        "-f", "mol2",
        "-o", "ante.frcmod"
    ]
    result = subprocess.run(cmd_parmchk, cwd=tmpdir, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"parmchk2 failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    # --- Step 6: Copy the GOOD files out of tmpdir ---
    shutil.copy(ante_output_mol2, final_mol2)
    shutil.copy(ante_output_frcmod, final_frcmod)

    # --- Step 7: Cleanup ---
    shutil.rmtree(tmpdir, ignore_errors=True)

    return final_mol2, final_frcmod


# ================================================================
# tleap
# ================================================================
def prepare_topology_and_coordinates(receptor_pdb: Path, ligand_mol2: Path,
                                     ligand_frcmod: Path, output_dir: Path,
                                     force_field: str = 'leaprc.protein.ff14SB'):
    """
    Create Amber topology/coordinate files for receptor and ligand.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    receptor_top = output_dir / "receptor.prmtop"
    receptor_crd = output_dir / "receptor.inpcrd"
    ligand_top = output_dir / "ligand.prmtop"
    ligand_crd = output_dir / "ligand.inpcrd"

    tleap_input = f"""
source {force_field}
source leaprc.gaff2

receptor = loadPdb {receptor_pdb}
ligand = loadMol2 {ligand_mol2}

loadAmberParams {ligand_frcmod}

saveAmberParm receptor {receptor_top} {receptor_crd}
saveAmberParm ligand {ligand_top} {ligand_crd}

quit
"""

    script = output_dir / "tleap_input.in"
    with open(script, 'w') as f:
        f.write(tleap_input)

    result = subprocess.run(["tleap", "-f", str(script)], capture_output=True, text=True)

    logfile = output_dir / "tleap_debug.log"
    with open(logfile, "w") as f:
        f.write(result.stdout)
        f.write("\n\n--- STDERR ---\n\n")
        f.write(result.stderr)

    logger.debug("\n\nFULL TLEAP LOG SAVED TO:", logfile)

    if result.returncode != 0:
        logger.error(f"tleap failed:\n{result.stderr}")
        raise RuntimeError("tleap failed")

    logger.info("tleap topology generation successful.")
    return receptor_top, receptor_crd, ligand_top, ligand_crd


# ================================================================
# 3D-RISM
# ================================================================
def run_3d_rism_calculation(receptor_top: Path, receptor_crd: Path,
                            ligand_top: Path, ligand_crd: Path,
                            output_dir: Path, config):
    """
    Run AmberTools rism3d.
    """
    solvent = config['rism'].get('solvent', 'water')
    grid_spacing = config['rism'].get('grid_spacing', 0.5)
    probe_radius = config['rism'].get('probe_radius', 1.4)

    output_dir.mkdir(exist_ok=True)

    out_nc = output_dir / "rism_output.nc"

    cmd = [
        "rism3d",
        "-p", str(receptor_top),
        "-c", str(receptor_crd),
        "-l", str(ligand_top),
        "-x", str(ligand_crd),
        "-o", str(out_nc),
        "-solvent", solvent,
        "-grid", str(grid_spacing),
        "-radius", str(probe_radius),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"3D-RISM error:\n{result.stderr}")
        raise RuntimeError("3D-RISM failed")

    logger.info(f"3D-RISM complete → {out_nc}")
    return out_nc


# ================================================================
# MASTER TASK
# ================================================================
@register_task("run_3d_rism",
               category='Solvent modeling',
               description="Run 3D-RISM calculation using AmberTools.")
def run_3d_rism(backend, ligand, config, **kwargs):

    # ✓ Always use the prepared receptor, never the raw file
    preparer = get_protein_preparer(backend, config)
    if preparer.protonated_pdb is None:
        raise RuntimeError("Receptor was not prepared before RISM step. Run prepare_receptor_pdbqt first.")

    receptor_pdb = preparer.protonated_pdb  # <------ FIXED

    output_dir = Path(config['docking']['output_dir'])
    output_dir.mkdir(exist_ok=True)

    # GNINA docking output
    sdf_file = output_dir / f"{ligand['name']}_conf0_docked.sdf"
    if not sdf_file.exists():
        raise FileNotFoundError(f"Docked SDF missing: {sdf_file}")

    # --- Step 1: Convert docked SDF → PDB for tleap ---
    ligand_pdb = convert_sdf_to_pdb(sdf_file, output_dir)

    # --- Step 2: Ligand parameterization ---
    ligand_mol2, ligand_frcmod = parameterize_ligand(
        sdf_file,       # correct input: SDF file 
        output_dir,
    )

    # --- Step 3: tleap ---
    receptor_top, receptor_crd, ligand_top, ligand_crd = prepare_topology_and_coordinates(
        receptor_pdb, ligand_mol2, ligand_frcmod, output_dir
    )

    # --- Step 4: 3D-RISM ---
    rism_output = run_3d_rism_calculation(
        receptor_top, receptor_crd, ligand_top, ligand_crd, output_dir, config
    )

    return rism_output
