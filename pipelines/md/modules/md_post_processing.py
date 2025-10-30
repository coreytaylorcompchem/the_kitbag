import os
from modules.automated_analyses import HydrogenBondAnalysisTask
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


class MDPostProcessingWorkflow:
    """Handles post-processing analyses for MD workflows."""

    def __init__(self, config, context=None):
        self.config = config
        self.context = context or {}

    def post_processing(self):
        post_cfg = self.config.get("post_processing", {})
        solvent_cfg = post_cfg.get("solvent_analyses")

        if solvent_cfg is True:
            solvent_cfg = {"enabled": True}

        if not solvent_cfg.get("enabled", False):
            logger.info("Solvent analysis disabled — skipping post-processing.")
            return {}

        topology = solvent_cfg.get("topology")
        trajectory = solvent_cfg.get("trajectory")

        # Fallback: look for previous step outputs
        if not topology:
            prep_cfg = self.config.get("prepare_system", {})
            prep_out = prep_cfg.get("output_trajectory", "output")
            topology = os.path.join(prep_out, "topology.pdb")
        if not trajectory:
            prod_cfg = self.config.get("production", {})
            trajectory = prod_cfg.get("output_trajectory", "output/trajectory.dcd")
            output_dir = os.path.dirname(trajectory) or "output"

        # Check files exist
        if not os.path.exists(topology):
            raise FileNotFoundError(f"Topology file not found: {topology}")
        if not os.path.exists(trajectory):
            raise FileNotFoundError(f"Trajectory file not found: {trajectory}")

        logger.info(f"Loading topology: {topology}")
        logger.info(f"Loading trajectory: {trajectory}")

        task = HydrogenBondAnalysisTask(
            topology=topology,
            trajectory=trajectory,
            binding_site_resids=[],
            ligand_resname=solvent_cfg.get("ligand_resname", "UNK"),
            start=solvent_cfg.get("start", 0),
            stop=solvent_cfg.get("stop", -1),
            step=solvent_cfg.get("step", 1),
            water_resname=solvent_cfg.get("water_resname", "HOH"),
            output_dir=output_dir,
            top_n=solvent_cfg.get("top_n", 20),
        )

        results = task.run()
        return results