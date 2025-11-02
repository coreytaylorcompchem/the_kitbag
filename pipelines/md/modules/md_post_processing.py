import os
from modules.automated_analyses import (
    HydrogenBondAnalysisTask,
    RMSDAnalysisTask,
    RMSFAnalysisTask,
    InteractionFingerprintTask,
    ProteinLigandCommunityTask,        # new
    HydrationSiteEnergyTask,           # new
    TemporalMotifPersistenceTask,      # new
    NetworkEmbeddingAnalysisTask,      # new
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
        
        # Prolif
        if "interactions" in time_series:
            logger.info("Running interactions analysis...")
            interactions_outdir = os.path.join(output_dir, "interactions")
            os.makedirs(interactions_outdir, exist_ok=True)
            interactions_task = InteractionFingerprintTask(
                topology=topology,
                trajectory=trajectory,
                ligand_selection=f"resname {ligand_resname}",  # respect YAML
                protein_selection="(protein or resname WAT) and byres around 20.0 group ligand",
                start=start,
                stop=stop,
                step=step,
                output_dir=interactions_outdir,
            )

            # Select ligand from YAML
            ligand = interactions_task.u.select_atoms(interactions_task.ligand_selection)
            if len(ligand) == 0:
                raise ValueError(f"No atoms found for ligand selection: {interactions_task.ligand_selection}")

            # Now select protein around the ligand AtomGroup
            protein = interactions_task.u.select_atoms(
                f"(protein or resname WAT) and byres around 20.0 group atomgroup",
                atomgroup=ligand  # pass the AtomGroup explicitly
            )
            if len(protein) == 0:
                raise ValueError("No protein atoms found around the ligand")

            results["interactions"] = interactions_task.run()

        # ==============================================================
        # 2️⃣  Graph / Network Analyses (advanced)
        # ==============================================================
        graph_analyses = traj_analysis_cfg.get("graph_analyses", [])

        # --- Protein–Ligand Community Detection ---
        if "protein_ligand_communities" in graph_analyses:
            logger.info("Running protein–ligand community analysis...")
            plin_task = ProteinLigandCommunityTask(
                topology=topology,
                trajectory=trajectory,
                ligand_resname=ligand_resname,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "protein_ligand_communities"),
            )
            results["protein_ligand_communities"] = plin_task.run()

        # --- Hydration Site Energetics ---
        if "hydration_site_energy" in graph_analyses:
            logger.info("Running hydration site energy analysis...")
            hse_task = HydrationSiteEnergyTask(
                topology=topology,
                trajectory=trajectory,
                ligand_resname=ligand_resname,
                water_resname=water_resname,
                start=start,
                stop=stop,
                step=step,
                # binding_site_cutoff=post_cfg.get("solvent_analyses", {}).get("binding_site_cutoff", 5.0),
                output_dir=os.path.join(output_dir, "hydration_site_energy"),
            )
            results["hydration_site_energy"] = hse_task.run()

        # --- Temporal Motif Persistence ---
        if "temporal_motif_persistence" in graph_analyses:
            logger.info("Running temporal motif persistence analysis...")
            tmp_task = TemporalMotifPersistenceTask(
                topology=topology,
                trajectory=trajectory,
                ligand_resname=ligand_resname,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "temporal_motif_persistence"),
            )
            results["temporal_motif_persistence"] = tmp_task.run()

        # --- Network Embedding Analysis ---
        if "network_embedding_analysis" in graph_analyses:
            logger.info("Running network embedding analysis...")
            ne_task = NetworkEmbeddingAnalysisTask(
                topology=topology,
                trajectory=trajectory,
                ligand_resname=ligand_resname,
                start=start,
                stop=stop,
                step=step,
                output_dir=os.path.join(output_dir, "network_embedding_analysis"),
            )
            results["network_embedding_analysis"] = ne_task.run()

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
