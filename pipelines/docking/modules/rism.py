import re
import subprocess
from pathlib import Path
import shutil
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def convert_sdf_to_pdb(sdf_file: Path, output_dir: Path) -> Path:
    """
    Convert the ligand from sdf to pdb format using Open Babel or any other conversion tool.
    """
    pdb_file = output_dir / f"{sdf_file.stem}.pdb"
    
    try:
        result = subprocess.run(
            ["obabel", str(sdf_file), "-O", str(pdb_file)],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            logger.error(f"Error converting SDF to PDB: {result.stderr}")
            raise RuntimeError(f"Failed to convert {sdf_file} to PDB.")
        
        logger.info(f"Converted {sdf_file} to {pdb_file}.")
        return pdb_file
    
    except Exception as e:
        logger.error(f"Error during SDF to PDB conversion: {e}")
        raise


def get_ligand_charge_from_smiles(pdb_path: Path) -> int:
    """
    Get formal charge by converting to SMILES and summing charge annotations.
    Uses Open Babel protonation at pH 7.4 for correct charge assignment.
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


def parameterize_ligand(ligand_pdb: Path, output_dir: Path):
    """
    Parameterise ligand using antechamber with explicit charge.
    Runs antechamber in a per-ligand working directory to avoid parallel collisions.
    """

    # Resolve to absolute paths
    ligand_pdb = ligand_pdb.resolve()
    output_dir = output_dir.resolve()

    net_charge = get_ligand_charge_from_smiles(ligand_pdb)
    logger.debug(f"Detected ligand charge = {net_charge} for {ligand_pdb.name}")

    mol2_file = (output_dir / f"{ligand_pdb.stem}.mol2").resolve()
    frcmod_file = (output_dir / f"{ligand_pdb.stem}.frcmod").resolve()

    # Create temp directory for antechamber
    tmpdir = (output_dir / f"{ligand_pdb.stem}_antechamber_tmp").resolve()
    tmpdir.mkdir(exist_ok=True, parents=True)

    # Run antechamber inside tmpdir
    cmd = [
        "antechamber",
        "-i", str(ligand_pdb),
        "-fi", "pdb",
        "-o", str(mol2_file),
        "-fo", "mol2",
        "-c", "bcc",
        "-nc", str(net_charge),
        "-at", "gaff2",
        "-j", "4",          # <--- IMPORTANT FIX
        "-s", "2"
    ]

    result = subprocess.run(cmd, cwd=tmpdir, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Antechamber error: {result.stderr}")
        raise RuntimeError(f"Failed to generate mol2: {ligand_pdb}")

    logger.info(f"Generated mol2 file: {mol2_file}")

    # Parmchk2 – safe to run in global dir
    result = subprocess.run(
        ["parmchk2", "-i", str(mol2_file), "-f", "mol2", "-o", str(frcmod_file)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        logger.error(f"parmchk2 error: {result.stderr}")
        raise RuntimeError(f"Failed to generate frcmod for {ligand_pdb}")

    logger.info(f"Generated frcmod: {frcmod_file}")

    # Cleanup temp folder
    try:
        shutil.rmtree(tmpdir)
    except Exception as e:
        logger.warning(f"Couldn't remove temp directory {tmpdir}: {e}")

    return mol2_file, frcmod_file


def prepare_topology_and_coordinates(receptor_pdb: Path, ligand_mol2: Path, ligand_frcmod: Path, output_dir: Path, force_field: str = 'ff14SB'):
    """
    Prepare the topology and coordinate files using tleap from AmberTools.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    receptor_top = output_dir / "receptor.prmtop"
    receptor_crd = output_dir / "receptor.inpcrd"
    ligand_top = output_dir / "ligand.prmtop"
    ligand_crd = output_dir / "ligand.inpcrd"

    tleap_input = f"""
source {force_field}
receptor = loadPdb {receptor_pdb}
ligand = loadMol2 {ligand_mol2}
loadAmberParams {ligand_frcmod}
saveAmberParm receptor {receptor_top} {receptor_crd}
saveAmberParm ligand {ligand_top} {ligand_crd}
quit
"""

    tleap_script = output_dir / "tleap_input.in"
    with open(tleap_script, 'w') as f:
        f.write(tleap_input)

    tleap_cmd = ["tleap", "-f", str(tleap_script)]
    result = subprocess.run(tleap_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"tleap failed:\n{result.stderr}")
        raise RuntimeError(f"Failed to prepare topology and coordinates for {receptor_pdb.name} and {ligand_mol2.name}")

    logger.info(f"Topologies and coordinates generated successfully.")
    return receptor_top, receptor_crd, ligand_top, ligand_crd


def run_3d_rism_calculation(receptor_top: Path, receptor_crd: Path, ligand_top: Path, ligand_crd: Path, output_dir: Path, config):
    """
    Runs 3D-RISM calculations using AmberTools with prepared topology and coordinate files.
    """
    solvent = config['rism'].get('solvent', 'water')
    grid_spacing = config['rism'].get('grid_spacing', 0.5)
    probe_radius = config['rism'].get('probe_radius', 1.4)

    output_dir.mkdir(parents=True, exist_ok=True)

    rism_output = output_dir / "rism_output.nc"

    rism_cmd = [
        "rism3d",
        "-p", str(receptor_top),
        "-c", str(receptor_crd),
        "-l", str(ligand_top),
        "-x", str(ligand_crd),
        "-o", str(rism_output),
        "-solvent", solvent,
        "-grid", str(grid_spacing),
        "-radius", str(probe_radius),
    ]

    result = subprocess.run(rism_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"3D-RISM calculation failed:\n{result.stderr}")
        raise RuntimeError(f"3D-RISM calculation failed for {receptor_top.name}")

    logger.info(f"3D-RISM calculation completed successfully. Output saved to: {rism_output}")
    return rism_output


@register_task("run_3d_rism", 
               category='Solvent modeling',
               description="Run 3D-RISM calculation using AmberTools.")
def run_3d_rism(backend, ligand, config, **kwargs):
    """
    Task to run 3D-RISM calculations and integrate high-probability waters into the receptor.
    """
    receptor_pdb = Path(config['protein']['pdb_path'])
    output_dir = Path(config['docking']['output_dir'])

    sdf_file = output_dir / f"{ligand['name']}_conf0_docked.sdf"

    if not receptor_pdb.exists():
        raise FileNotFoundError(f"Receptor PDB not found at: {receptor_pdb}")
    if not sdf_file.exists():
        raise FileNotFoundError(f"Ligand SDF not found at: {sdf_file}")

    ligand_pdb = convert_sdf_to_pdb(sdf_file, output_dir)

    ligand_mol2, ligand_frcmod = parameterize_ligand(ligand_pdb, output_dir)

    receptor_top, receptor_crd, ligand_top, ligand_crd = prepare_topology_and_coordinates(
        receptor_pdb, ligand_mol2, ligand_frcmod, output_dir,
        force_field='leaprc.protein.ff14SB'
    )
    
    rism_output = run_3d_rism_calculation(receptor_top, receptor_crd, ligand_top, ligand_crd, output_dir, config)
    
    return rism_output
