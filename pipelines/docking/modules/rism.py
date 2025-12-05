import subprocess
from pathlib import Path
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=True,  # Enable debug mode for better visibility
    simple_format=True
)

def convert_sdf_to_pdb(sdf_file: Path, output_dir: Path) -> Path:
    """
    Convert the ligand from sdf to pdb format using Open Babel or any other conversion tool.
    """
    pdb_file = output_dir / f"{sdf_file.stem}.pdb"
    
    try:
        # Assuming Open Babel is installed
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

def parameterize_ligand(ligand_pdb: Path, output_dir: Path):
    """
    Parameterize the ligand using antechamber and parmchk2 to generate mol2 and frcmod files.
    """
    mol2_file = output_dir / f"{ligand_pdb.stem}.mol2"
    frcmod_file = output_dir / f"{ligand_pdb.stem}.frcmod"

    # Step 1: Use antechamber to generate the mol2 file
    try:
        result = subprocess.run(
            ["antechamber", "-i", str(ligand_pdb), "-o", str(mol2_file), "-fi", "pdb", "-fo", "mol2", "-c", "bcc", "-s", "2"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            logger.error(f"Error running antechamber: {result.stderr}")
            raise RuntimeError(f"Failed to generate mol2 for ligand {ligand_pdb.name}")
        logger.info(f"Generated mol2 file: {mol2_file}")
    except Exception as e:
        logger.error(f"Error during ligand parameterization with antechamber: {e}")
        raise

    # Step 2: Use parmchk2 to generate the frcmod file for missing parameters
    try:
        result = subprocess.run(
            ["parmchk2", "-i", str(mol2_file), "-f", "mol2", "-o", str(frcmod_file)],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            logger.error(f"Error running parmchk2: {result.stderr}")
            raise RuntimeError(f"Failed to generate frcmod for ligand {ligand_pdb.name}")
        logger.info(f"Generated frcmod file: {frcmod_file}")
    except Exception as e:
        logger.error(f"Error during frcmod generation with parmchk2: {e}")
        raise

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

    # Correct tleap input to use mol2 and frcmod for ligand
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
    try:
        result = subprocess.run(tleap_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"tleap failed:\n{result.stderr}")
            raise RuntimeError(f"Failed to prepare topology and coordinates for {receptor_pdb.name} and {ligand_mol2.name}")
        logger.info(f"Topologies and coordinates generated successfully.")
    except Exception as e:
        logger.error(f"❌ Error during tleap preparation: {e}")
        raise

    return receptor_top, receptor_crd, ligand_top, ligand_crd

def run_3d_rism_calculation(receptor_top: Path, receptor_crd: Path, ligand_top: Path, ligand_crd: Path, output_dir: Path, config):
    """
    Runs 3D-RISM calculations using AmberTools with prepared topology and coordinate files.
    """
    # Retrieve the 3D-RISM-specific options from config
    solvent = config['rism'].get('solvent', 'water')
    grid_spacing = config['rism'].get('grid_spacing', 0.5)
    probe_radius = config['rism'].get('probe_radius', 1.4)

    # Make sure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    rism_output = output_dir / "rism_output.nc"

    rism_cmd = [
        "rism3d",  # Assuming this is the path to AmberTools' rism3d tool
        "-p", str(receptor_top),  # Path to prepared receptor topology file
        "-c", str(receptor_crd),  # Path to receptor coordinates file
        "-l", str(ligand_top),  # Path to ligand topology file
        "-x", str(ligand_crd),  # Path to ligand coordinates file
        "-o", str(rism_output),  # Output file for results
        "-solvent", solvent,  # Solvent type (from config)
        "-grid", str(grid_spacing),  # Grid spacing in Å (from config)
        "-radius", str(probe_radius),  # Probe radius for solvent calculation (from config)
    ]

    try:
        result = subprocess.run(rism_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"3D-RISM calculation failed:\n{result.stderr}")
            raise RuntimeError(f"3D-RISM calculation failed for {receptor_top.name}")
        logger.info(f"3D-RISM calculation completed successfully. Output saved to: {rism_output}")
        return rism_output
    except Exception as e:
        logger.error(f"❌ Error during 3D-RISM calculation: {e}")
        raise

@register_task("run_3d_rism", 
               category='Solvent modeling',
               description="Run 3D-RISM calculation using AmberTools.")
def run_3d_rism(backend, ligand, config, **kwargs):
    """
    Task to run 3D-RISM calculations and integrate high-probability waters into the receptor.
    """
    receptor_pdb = Path(config['protein']['pdb_path'])
    output_dir = Path(config['docking']['output_dir'])

    # Assuming ligand has already been docked and converted to PDB
    sdf_file = Path(config['docking']['output_dir']) / f"{ligand['name']}_conf0_docked.sdf"

    if not receptor_pdb.exists():
        raise FileNotFoundError(f"Receptor PDB not found at: {receptor_pdb}")
    if not sdf_file.exists():
        raise FileNotFoundError(f"Ligand SDF not found at: {sdf_file}")

    # Convert SDF to PDB for the ligand
    ligand_pdb = convert_sdf_to_pdb(sdf_file, output_dir)

    # Parameterize the ligand (generate mol2 and frcmod)
    ligand_mol2, ligand_frcmod = parameterize_ligand(ligand_pdb, output_dir)

    # Prepare topology and coordinate files for receptor and ligand
    receptor_top, receptor_crd, ligand_top, ligand_crd = prepare_topology_and_coordinates(
        receptor_pdb, ligand_mol2, ligand_frcmod, output_dir, force_field='leaprc.protein.ff14SB'  # Use the desired force field here
    )
    
    # Run the 3D-RISM calculation with the provided configuration options
    rism_output = run_3d_rism_calculation(receptor_top, receptor_crd, ligand_top, ligand_crd, output_dir, config)
    
    return rism_output
