from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from pathlib import Path
import requests

import modeller

logger = setup_logger(__name__, debug_mode=True, simple_format=True)


def fetch_uniprot_fasta(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if not response.ok:
        raise ValueError(f"Could not fetch UniProt FASTA for {uniprot_id}")
    return ''.join(response.text.split('\n')[1:])


def group_missing_residues(missing_indices):
    if not missing_indices:
        return []

    groups = []
    start = prev = missing_indices[0]

    for idx in missing_indices[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            groups.append((start, prev))
            start = prev = idx
    groups.append((start, prev))
    return groups


from backends.modeller import ModellerBackend  # Adjust path as needed

@register_task("fix_residues", category="Protein modelling", description="Fix incomplete residues and extract Chain A")
def fix_residues(backend, config, **kwargs):
    protein_cfg = config.get("protein", {})
    pdb_path = protein_cfg["pdb_path"]
    uniprot_id = protein_cfg["uniprot_id"]
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    n_loop_models = config["backend"].get("n_loop_models", 2)
    loop_refinement = config["backend"].get("loop_refinement", {})

    backend = ModellerBackend(pdb_path, n_loop_models=n_loop_models, loop_refinement=loop_refinement, output_dir=output_dir)

    chain_id = "A"
    chain_selection = modeller.Selection(backend.model)
    for res in backend.model.residues:
        if res.chain == chain_id:
            chain_selection.add(res)

    if not chain_selection:
        raise ValueError(f"Chain {chain_id} not found in PDB")

    chain_a_pdb = output_dir / "chain_A_only.pdb"
    chain_selection.write(file=str(chain_a_pdb))
    logger.info(f"Chain A only PDB saved to: {chain_a_pdb}")

    # Reload chain A model into backend with updated output_dir
    backend = ModellerBackend(str(chain_a_pdb), n_loop_models=n_loop_models, loop_refinement=loop_refinement, output_dir=output_dir)

    backend.complete_missing_atoms(uniprot_id=uniprot_id)

    fixed_pdb = output_dir / "fixed.pdb"
    backend.write_pdb(str(fixed_pdb))
    logger.info(f"Fixed PDB saved to: {fixed_pdb}")

    backend.cache = {
        "chain_a_pdb": str(chain_a_pdb),
        "fixed_pdb": str(fixed_pdb),
    }

    return backend.cache


@register_task("fix_loops", category="Protein modelling", description="Model missing loops.")
def fix_loops(backend, config, **kwargs):
    protein_cfg = config.get("protein", {})
    uniprot_id = protein_cfg["uniprot_id"]
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    fixed_pdb = output_dir / "fixed.pdb"
    n_loop_models = config["backend"].get("n_loop_models", 2)

    backend = ModellerBackend(str(fixed_pdb), n_loop_models=n_loop_models, output_dir=output_dir)
    # For simplicity, loop ranges should be passed in config or calculated
    missing_residues = config["protein"].get("missing_residues", [])
    loop_ranges = group_missing_residues(missing_residues)

    backend.fix_loops(loop_ranges)

    loop_model_pdb = output_dir / "loop_model.pdb"
    backend.write_pdb(str(loop_model_pdb))
    logger.info(f"Loop modeled PDB saved to: {loop_model_pdb}")

    backend.cache = {"loop_model_pdb": str(loop_model_pdb)}

    return backend.cache


@register_task("refine_loops", category="Protein modelling", description="Refine loop structures.")
def refine_loops(backend, config, **kwargs):
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    loop_model_pdb = output_dir / "loop_model.pdb"
    n_refine_models = config["backend"].get("loop_refinement", {}).get("n_models", 1)

    backend = ModellerBackend(str(loop_model_pdb), loop_refinement={"n_loop_models": n_refine_models}, output_dir=output_dir)
    backend.refine_loops()

    refined_pdb = output_dir / "refined.pdb"
    backend.write_pdb(str(refined_pdb))
    logger.info(f"Refined PDB saved to: {refined_pdb}")

    backend.cache = {"refined_pdb": str(refined_pdb)}

    return backend.cache
