import os
from modules.automated_analyses import (
    HydrogenBondAnalysisTask,
    RMSDAnalysisTask,
    RMSFAnalysisTask,
)
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


class MDPostProcessingWorkflow:
    """Handles post-processing analyses for MD workflows."""

    def __init__(self, config, context=None):
        self.config = config
        self.context = context or {}

    def post_processing(self):
        post_cfg = self.config.get("post_processing", {})

        # --------------------------------------------------------------
        # General / shared parameters
        # --------------------------------------------------------------
        output_dir = post_cfg.get("output_dir", "output_postproc")
        ligand_resname = post_cfg.get("ligand_resname", "UNK")
        water_resname = post_cfg.get("water_resname", "HOH")
        start = post_cfg.get("start", 0)
        stop = post_cfg.get("stop", -1)
        step = post_cfg.get("step", 1)

        os.makedirs(output_dir, exist_ok=True)

        # Determine topology and trajectory
        topology = post_cfg.get("topology")
        trajectory = post_cfg.get("trajectory")

        # Fallback to earlier workflow outputs
        if not topology:
            prep_cfg = self.config.get("prepare_system", {})
            prep_out = prep_cfg.get("output_trajectory", "output")
            topology = os.path.join(prep_out, "topology.pdb")

        if not trajectory:
            prod_cfg = self.config.get("production", {})
            trajectory = prod_cfg.get("output_trajectory", os.path.join(output_dir, "trajectory.dcd"))

        # Sanity checks
        if not os.path.exists(topology):
            raise FileNotFoundError(f"Topology file not found: {topology}")
        if not os.path.exists(trajectory):
            raise FileNotFoundError(f"Trajectory file not found: {trajectory}")

        logger.info(f"Topology: {topology}")
        logger.info(f"Trajectory: {trajectory}")

        results = {}

        # ==============================================================
        # 1️⃣  Time-series analyses (RMSD, RMSF, interactions, etc.)
        # ==============================================================
        traj_analysis_cfg = post_cfg.get("trajectory_analysis", {})
        time_series = traj_analysis_cfg.get("time_series", [])

        # --- RMSD ---
        if "rmsd" in time_series:
            logger.info("Running RMSD analysis...")
            rmsd_outdir = os.path.join(output_dir, "rmsd")
            os.makedirs(rmsd_outdir, exist_ok=True)
            rmsd_task = RMSDAnalysisTask(
                topology=topology,
                trajectory=trajectory,
                start=start,
                stop=stop,
                step=step,
                output_dir=rmsd_outdir,
                ligand_resname=ligand_resname,
            )
            results["rmsd"] = rmsd_task.run()

        # --- RMSF ---
        if "rmsf" in time_series:
            logger.info("Running RMSF analysis...")
            rmsf_outdir = os.path.join(output_dir, "rmsf")
            os.makedirs(rmsf_outdir, exist_ok=True)
            rmsf_task = RMSFAnalysisTask(
                topology=topology,
                trajectory=trajectory,
                start=start,
                stop=stop,
                step=step,
                output_dir=rmsf_outdir,
            )
            results["rmsf"] = rmsf_task.run()

        # (future: interactions, etc.)
        # ==============================================================
        # 2️⃣  Solvent-mediated hydrogen bond analysis
        # ==============================================================
        solvent_cfg = post_cfg.get("solvent_analyses", {})
        if isinstance(solvent_cfg, bool):
            solvent_cfg = {"enabled": solvent_cfg}

        if solvent_cfg.get("enabled", False):
            logger.info("Running solvent hydrogen bond analysis...")

            solvent_outdir = os.path.join(output_dir, "solvent_hbonds")
            os.makedirs(solvent_outdir, exist_ok=True)

            task = HydrogenBondAnalysisTask(
                topology=topology,
                trajectory=trajectory,
                ligand_resname=ligand_resname,
                water_resname=water_resname,
                start=start,
                stop=stop,
                step=step,
                binding_site_cutoff=solvent_cfg.get("binding_site_cutoff", 5.0),
                top_n=solvent_cfg.get("top_n", 20),
                output_dir=solvent_outdir,
            )
            results["solvent_hbonds"] = task.run()
        else:
            logger.info("Solvent analysis disabled — skipping.")

        return results
