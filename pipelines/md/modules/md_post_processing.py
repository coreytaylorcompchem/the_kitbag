import os
import re
from natsort import natsorted  # pip install natsort
from modules.automated_analyses import (
    HydrogenBondAnalysisTask,
    RMSDAnalysisTask,
    RMSFAnalysisTask,
    InteractionFingerprintTask,
    ProteinLigandCommunityTask,        
    HydrationSiteEnergyTask,        
    TemporalMotifPersistenceTask,     
    NetworkEmbeddingAnalysisTask,    
    ProteinProteinNetworkEmbeddingTask,
)
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


class MDPostProcessingWorkflow:
    """Handles post-processing analyses for MD workflows."""

    def __init__(self, config, context=None):
        self.config = config
        self.context = context or {}

    def _collect_trajectories(self, traj_path):
        """Detect and return trajectory files in numerical order."""
        if os.path.isdir(traj_path):
            traj_dir = traj_path
        else:
            traj_dir = os.path.dirname(traj_path)
            if traj_dir == "":
                traj_dir = "."

        traj_base = os.path.splitext(os.path.basename(traj_path))[0]
        traj_ext = os.path.splitext(traj_path)[1]                   

        # Collect all trajs files
        pattern = re.compile(rf"{re.escape(traj_base)}[_-]?(\d+){re.escape(traj_ext)}$")
        all_files = [
            os.path.join(traj_dir, f)
            for f in os.listdir(traj_dir)
            if f.endswith(traj_ext) and pattern.match(f)
        ]

        if not all_files:
            # Fallback for single traj.dcd
            if os.path.exists(traj_path):
                return [traj_path]
            raise FileNotFoundError(f"No trajectory files found matching pattern {traj_base}*.{traj_ext}")

        # Numerical sort
        traj_files = natsorted(all_files)
        logger.info(f"Detected {len(traj_files)} trajectory files.")
        for f in traj_files:
            logger.debug(f"  {f}")
        return traj_files

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

        # Collect all traj files
        trajectories = self._collect_trajectories(trajectory)

        # Try to auto-detect wrapped trajectory
        traj_dir = os.path.dirname(trajectory) or "."
        traj_base = os.path.splitext(os.path.basename(trajectory))[0]

        wrapped_pattern = os.path.join(traj_dir, f"{traj_base}*wrapped*.xtc")

        wrapped_files = [
            os.path.join(traj_dir, f)
            for f in os.listdir(traj_dir)
            if "wrapped" in f and f.endswith(".xtc")
        ]

        wrapped_trajectories = natsorted(wrapped_files) if wrapped_files else None

        if wrapped_trajectories:
            logger.info(f"Detected wrapped trajectories: {len(wrapped_trajectories)} file(s)")
        else:
            logger.warning("No wrapped trajectories found — falling back to raw DCD (may break contact analyses)")

        logger.info(f"Topology: {topology}")
        logger.info(f"Trajectories: {len(trajectories)} file(s)")
        if len(trajectories) == 1:
            logger.info(f"Using trajectory: {trajectories[0]}")
        else:
            logger.info("Using multiple trajectories in sorted order.")

        results = {}

        # ==============================================================
        # Time-series analyses (RMSD, RMSF, interactions, etc.)
        # ==============================================================
        traj_analysis_cfg = post_cfg.get("trajectory_analysis", {})
        time_series = traj_analysis_cfg.get("time_series", [])

        if "rmsd" in time_series:
            logger.info("Running RMSD analysis...")
            rmsd_outdir = os.path.join(output_dir, "rmsd")
            os.makedirs(rmsd_outdir, exist_ok=True)
            rmsd_task = RMSDAnalysisTask(
                topology=topology,
                trajectory=trajectories,
                start=start,
                stop=stop,
                step=step,
                output_dir=rmsd_outdir,
                ligand_resname=ligand_resname,
            )
            results["rmsd"] = rmsd_task.run()

        if "rmsf" in time_series:
            logger.info("Running RMSF analysis...")
            rmsf_outdir = os.path.join(output_dir, "rmsf")
            os.makedirs(rmsf_outdir, exist_ok=True)
            rmsf_task = RMSFAnalysisTask(
                topology=topology,
                trajectory=trajectories,
                start=start,
                stop=stop,
                step=step,
                output_dir=rmsf_outdir,
            )
            results["rmsf"] = rmsf_task.run()

        if "interactions" in time_series:
            logger.info("Running interaction analysis (ProLIF)...")
            interactions_outdir = os.path.join(output_dir, "interactions")
            os.makedirs(interactions_outdir, exist_ok=True)
            interactions_task = InteractionFingerprintTask(
                topology=topology,
                trajectory=wrapped_trajectories or trajectories,
                ligand_selection=f"resname {ligand_resname}",
                protein_selection="(protein or resname WAT) and byres around 20.0 group ligand",
                start=start,
                stop=stop,
                step=step,
                output_dir=interactions_outdir,
            )
            results["interactions"] = interactions_task.run()

        # ==============================================================
        # Graph / Network Analyses
        # ==============================================================
        graph_analyses = traj_analysis_cfg.get("graph_analyses", [])

        if "protein_ligand_communities" in graph_analyses:
            logger.info("Running protein–ligand community analysis...")
            task = ProteinLigandCommunityTask(
                topology=topology,
                trajectory=wrapped_trajectories or trajectories,
                ligand_resname=ligand_resname,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "protein_ligand_communities"),
            )
            results["protein_ligand_communities"] = task.run()

        if "hydration_site_energy" in graph_analyses:
            logger.info("Running hydration site energy analysis...")
            task = HydrationSiteEnergyTask(
                topology=topology,
                trajectory=wrapped_trajectories or trajectories,
                ligand_resname=ligand_resname,
                water_resname=water_resname,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "hydration_site_energy"),
            )
            results["hydration_site_energy"] = task.run()

        if "temporal_motif_persistence" in graph_analyses:
            logger.info("Running temporal motif persistence analysis...")
            task = TemporalMotifPersistenceTask(
                topology=topology,
                trajectory=trajectories,
                ligand_resname=ligand_resname,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "temporal_motif_persistence"),
            )
            results["temporal_motif_persistence"] = task.run()

        if "network_embedding_analysis" in graph_analyses:
            logger.info("Running network embedding analysis...")
            task = NetworkEmbeddingAnalysisTask(
                topology=topology,
                trajectory=trajectories,
                ligand_resname=ligand_resname,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "network_embedding_analysis"),
            )
            results["network_embedding_analysis"] = task.run()

        if "protein_protein_network_embedding" in graph_analyses:
            logger.info("Running protein-protein network embedding analysis...")
            task = ProteinProteinNetworkEmbeddingTask(
                topology=topology,
                trajectory=trajectories,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "protein_protein_network_embedding_analysis"),
            )
            results["protein_protein_network_embedding_analysis"] = task.run()

        # ==============================================================
        # Solvent-mediated hydrogen bond analysis
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
                trajectory=wrapped_trajectories or trajectories,
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
            logger.info("Solvent analysis disabled - skipping.")

        return results
