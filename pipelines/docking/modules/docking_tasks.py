from pathlib import Path
from modules.vanilla_docking import LigandPreparer
from docking_task_registry import register_task

from rdkit import Chem
from tqdm import tqdm

@register_task("standardize_ligand", description="Standardize ligand from SMILES.", supported_backends=["gnina"])
def standardize_ligand(backend, ligand_info, config):
    lp = LigandPreparer(smiles=ligand_info['smiles'], name=ligand_info['name'])
    lp.standardize()
    backend.cache[ligand_info['name']] = lp

@register_task("generate_conformers", description="Generate RDKit conformers.")
def generate_conformers(backend, ligand_info, config):
    lp = backend.cache.get(ligand_info['name'])
    if not lp:
        raise ValueError("Ligand not found in cache. Did you run 'standardize_ligand'?")
    n_confs = config.get("n_conformers", 250)
    lp.generate_conformers(n_confs=n_confs)

@register_task("cluster_conformers", description="Cluster and select conformers.")
def cluster_conformers(backend, ligand_info, config):
    lp = backend.cache.get(ligand_info['name'])
    if not lp:
        raise ValueError("Ligand not found in cache.")

    docking_cfg = config.get("docking", {})
    final_n = docking_cfg.get("final_n_conformers", 5)
    rmsd_thresh = docking_cfg.get("rmsd_threshold", 0.75)
    min_gap = docking_cfg.get("min_energy_gap", 0.5)

    lp.cluster_and_select(
        final_n=final_n,
        rmsd_threshold=rmsd_thresh,
        min_energy_gap=min_gap
    )

@register_task("optimize_with_xtb", description="Optimize conformers using GFN1-xTB.")
def optimize_with_xtb(backend, ligand_info, config):
    lp = backend.cache.get(ligand_info['name'])
    if not lp:
        raise ValueError("Ligand not found in cache.")
    output_dir = Path(config['output_dir'])
    output_dir.mkdir(exist_ok=True, parents=True)
    lp.optimize_with_xtb(output_dir=output_dir)

@register_task("save_final_conformers", description="Save final conformers to SDF.")
def save_final_conformers(backend, ligand_info, config):
    lp = backend.cache.get(ligand_info['name'])
    if not lp:
        raise ValueError("Ligand not found in cache.")
    output_dir = Path(config['output_dir'])
    sdf_path = lp.save_final_conformers(output_dir)
    backend.cache[f"{ligand_info['name']}_sdf_path"] = sdf_path

@register_task("convert_to_pdbqt", description="Convert final conformers to PDBQT.")
def convert_to_pdbqt(backend, ligand_info, config):
    lp = backend.cache.get(ligand_info['name'])
    if not lp:
        raise ValueError("Ligand not found in cache.")
    
    output_dir = Path(config['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)

    docking_mode = config.get("docking", {}).get("docking_mode", "ensemble")

    # Use the mode from config, and let lp handle all file creation
    pdbqt_paths = lp.convert_to_pdbqt(output_dir=output_dir, mode=docking_mode)

    # Cache the generated PDBQT paths
    backend.cache[f"{ligand_info['name']}_pdbqt_path"] = pdbqt_paths

    print(f"[INFO] Converted {len(pdbqt_paths)} conformers to PDBQT for {ligand_info['name']}")


@register_task("dock", description="Run docking using Gnina backend.")
def dock(backend, ligand_info, config):
    output_dir = Path(config['output_dir'])
    ligand_name = ligand_info['name']
    
    pdbqt_paths = backend.cache.get(f"{ligand_name}_pdbqt_path")
    if pdbqt_paths is None:
        raise ValueError(f"PDBQT path not found for ligand '{ligand_name}'. Did you run 'convert_to_pdbqt'?")

    if not isinstance(pdbqt_paths, list):
        pdbqt_paths = [pdbqt_paths]  # Backward compatibility

    receptor_pdbqt = backend.cache.get("receptor_pdbqt")
    if receptor_pdbqt is None:
        raise ValueError("Receptor PDBQT path not found in backend cache.")

    docking_cfg = config.get("docking", {})
    docking_mode = docking_cfg.get("docking_mode", "ensemble")

    if 'center' not in docking_cfg or 'size' not in docking_cfg:
        raise ValueError("Docking 'center' and 'size' must be specified in config under 'docking'.")

    center = tuple(docking_cfg['center'])
    size = tuple(docking_cfg['size'])

    for i, pdbqt_path in enumerate(tqdm(pdbqt_paths, desc=f"[GNINA] Docking {ligand_name}", unit="conf")):
        output_path = output_dir / f"{ligand_name}_conf{i}_docked.sdf"
        backend.dock(
            receptor_path=receptor_pdbqt,
            ligand_path=pdbqt_path,
            output_path=output_path,
            center=center,
            size=size
        )
