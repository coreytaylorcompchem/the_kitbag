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
    Parameterize ligand by converting docked SDF → PDB, then mol2 → antechamber, etc.
    """
    docked_sdf = docked_sdf.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(exist_ok=True, parents=True)

    corrected_pdb = rdkit_reconstruct_ligand(
        docked_sdf,
        output_dir
    )

    mol2_file = output_dir / f"{corrected_pdb.stem}.mol2"
    pdb_to_mol2_with_obabel(corrected_pdb, mol2_file)

    net_charge = get_ligand_charge_from_smiles(corrected_pdb)
    logger.debug(f"Detected ligand charge = {net_charge}")

    frcmod_file = output_dir / f"{corrected_pdb.stem}.frcmod"

    tmpdir = (output_dir / f"{corrected_pdb.stem}_ante_tmp")
    tmpdir.mkdir(exist_ok=True)

    cmd = [
        "antechamber",
        "-i", str(mol2_file),
        "-fi", "mol2",
        "-o", str(mol2_file),
        "-fo", "mol2",
        "-c", "bcc",
        "-nc", str(net_charge),
        "-at", "gaff2",
        "-j", "4",
        "-s", "2"
    ]

    result = subprocess.run(cmd, cwd=tmpdir, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Antechamber error: {result.stderr}")
        raise RuntimeError("Antechamber failed")

    result = subprocess.run(
        ["parmchk2", "-i", str(mol2_file), "-f", "mol2", "-o", str(frcmod_file)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        logger.error(f"parmchk2 error: {result.stderr}")
        raise RuntimeError("parmchk2 failed")

    logger.info(f"Generated MOL2 + FRCMOD: {mol2_file}, {frcmod_file}")
    shutil.rmtree(tmpdir, ignore_errors=True)

    return mol2_file, frcmod_file



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

    receptor_pdb = Path(config['protein']['pdb_path'])
    output_dir = Path(config['docking']['output_dir'])

    sdf_file = output_dir / f"{ligand['name']}_conf0_docked.sdf"
    original_smiles = ligand["smiles"]

    if not receptor_pdb.exists():
        raise FileNotFoundError(f"Receptor PDB missing: {receptor_pdb}")
    if not sdf_file.exists():
        raise FileNotFoundError(f"SDF missing: {sdf_file}")

    # Step 1: SDF → PDB (still needed for tleap)
    ligand_pdb = convert_sdf_to_pdb(sdf_file, output_dir)

    # Step 2: Parameterization (FIXED: pass SDF!)
    ligand_mol2, ligand_frcmod = parameterize_ligand(
        sdf_file,        # <---- FIXED
        output_dir,
    )

    # Step 3: tleap
    receptor_top, receptor_crd, ligand_top, ligand_crd = prepare_topology_and_coordinates(
        receptor_pdb, ligand_mol2, ligand_frcmod, output_dir
    )

    # Step 4: 3D-RISM
    rism_output = run_3d_rism_calculation(
        receptor_top, receptor_crd, ligand_top, ligand_crd, output_dir, config
    )

    return rism_output
