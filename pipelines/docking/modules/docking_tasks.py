from pathlib import Path
from modules.vanilla_docking import LigandPreparer
from docking_task_registry import register_task

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
    final_n = config.get("final_n_conformers", 35)
    lp.cluster_and_select(final_n=final_n)

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
    docking_mode = config.get('docking_mode', 'ensemble')

    pdbqt_paths = lp.convert_to_pdbqt(output_dir=output_dir, mode=docking_mode)
    backend.cache[f"{ligand_info['name']}_pdbqt_paths"] = pdbqt_paths

@register_task("dock", description="Run docking using Gnina backend.")
def dock(backend, ligand_info, config):
    output_dir = Path(config['output_dir'])
    ligand_name = ligand_info['name']

    pdbqt_paths = backend.cache.get(f"{ligand_name}_pdbqt_paths")
    if not pdbqt_paths:
        raise ValueError(f"PDBQT paths not found for ligand '{ligand_name}'. Did you run 'convert_to_pdbqt'?")

    receptor_pdbqt = backend.cache.get("receptor_pdbqt")
    if receptor_pdbqt is None:
        raise ValueError("Receptor PDBQT path not found.")

    docking_cfg = config.get('docking', {})
    center = tuple(docking_cfg['center'])
    size = tuple(docking_cfg['size'])

    for pdbqt_path in pdbqt_paths:
        suffix = pdbqt_path.stem.split('_')[-1]  # conf0, conf1, etc.
        output_path = output_dir / f"{ligand_name}_{suffix}_docked.sdf"

        backend.dock(
            receptor_path=receptor_pdbqt,
            ligand_path=pdbqt_path,
            output_path=output_path,
            center=center,
            size=size
        )
