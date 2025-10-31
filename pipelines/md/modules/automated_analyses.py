import os
import warnings
import json

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import networkx as nx
import tqdm as tqdm
import numpy as np


# suppress MDA warnings.

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis")

from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

from typing import List, Optional, Dict, Any
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


@register_task(
    "solvent_hbonds",
    category="Analyses",
    description="Compute direct and water-mediated hydrogen bonds between binding site and ligand."
)
class HydrogenBondAnalysisTask:
    """
    Task to calculate both direct and water-mediated hydrogen bonds
    between binding site and ligand, then plot summary results.
    """

    def __init__(self,
                 topology: str,
                 trajectory: str,
                 binding_site_resids: Optional[List[int]] = None,
                 ligand_resname: str = "LIG",
                 water_resname: str = "HOH",
                 top_n: int = 20,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 1,
                 binding_site_cutoff: float = 5.0,
                 output_dir: str = "output",
                 context: Optional[Dict[str, Any]] = None):
        self.topology = topology
        self.trajectory = trajectory
        self.binding_site_resids = binding_site_resids
        self.ligand_resname = ligand_resname
        self.water_resname = water_resname
        self.start = start
        self.stop = stop
        self.step = step
        self.binding_site_cutoff = binding_site_cutoff
        self.output_dir = output_dir
        self.context = context or {}
        self.top_n = top_n

        # Load MDAnalysis Universe
        self.u = mda.Universe(self.topology, self.trajectory)

        # Automatically detect binding site if not provided
        if not self.binding_site_resids:
            self.u.trajectory[0]
            sel = self.u.select_atoms(f"protein and around {self.binding_site_cutoff} resname {self.ligand_resname}")
            self.binding_site_resids = list(sel.resids)
            logger.info(f"Detected {len(self.binding_site_resids)} binding site residues.")

        self.direct_paths = []
        self.water_paths = []

    def _compute_hbonds(self, water_mediated: bool = False):
        paths = []
        binding_site_atoms = self.u.select_atoms(
            "protein and resid " + " ".join(str(r) for r in self.binding_site_resids)
        )

        for ts in tqdm.tqdm(self.u.trajectory[self.start:self.stop:self.step],
                            desc="Direct HBonds" if not water_mediated else "Water-mediated HBonds"):
            if water_mediated:
                donors_sel = f"prop mass > 13 and resname {self.ligand_resname} " \
                             f"or (prop mass > 13 and resname {self.water_resname} and around 10 resname {self.ligand_resname}) " \
                             f"or (prop mass > 13 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                hydrogens_sel = f"prop mass < 2 and resname {self.ligand_resname} " \
                                f"or (prop mass < 2 and resname {self.water_resname} and around 10 resname {self.ligand_resname}) " \
                                f"or (prop mass < 2 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                acceptors_sel = donors_sel
                between = [f"protein or resname {self.ligand_resname}", f"resname {self.water_resname}"]
                path_length = 3
            else:
                donors_sel = f"prop mass > 13 and resname {self.ligand_resname} " \
                             f"or (prop mass > 13 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                hydrogens_sel = f"prop mass < 2 and resname {self.ligand_resname} " \
                                f"or (prop mass < 2 and protein and index {' '.join(map(str, binding_site_atoms.indices))})"
                acceptors_sel = donors_sel
                between = ["protein", f"resname {self.ligand_resname}"]
                path_length = 2

            hba = HBA(
                universe=self.u,
                donors_sel=donors_sel,
                hydrogens_sel=hydrogens_sel,
                acceptors_sel=acceptors_sel,
                between=between,
                d_h_cutoff=1.2,
                d_h_a_angle_cutoff=120,
                d_a_cutoff=3.0
            )

            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r"No hydrogen bonds were found given angle of",
                    category=UserWarning,
                    module="MDAnalysis.analysis.hydrogenbonds.hbond_analysis"
                )
            hba.run(start=ts.frame, stop=ts.frame + 1)
            unique_hbs = [list(i) for i in hba.hbonds[:, 1:-2] if i[0] != i[1]]

            G = nx.Graph()
            for edge in unique_hbs:
                donor = f"{self.u.atoms[int(edge[0])].resname}{self.u.atoms[int(edge[0])].resid}"
                acceptor = f"{self.u.atoms[int(edge[2])].resname}{self.u.atoms[int(edge[2])].resid}"
                if self.u.atoms[int(edge[0])].resname != self.water_resname:
                    donor += f": {self.u.atoms[int(edge[0])].name}"
                if self.u.atoms[int(edge[2])].resname != self.water_resname:
                    acceptor += f": {self.u.atoms[int(edge[2])].name}"
                G.add_edge(donor, acceptor)

            sources = [n for n in G.nodes if self.ligand_resname not in n and self.water_resname not in n]
            targets = [n for n in G.nodes if self.ligand_resname in n]

            for target in targets:
                if G.has_node(target):
                    for s in sources:
                        if nx.has_path(G, source=s, target=target):
                            paths_iter = nx.all_shortest_paths(G, source=s, target=target)
                            for p in paths_iter:
                                if len(p) == path_length:
                                    paths.append([list(p), ts.frame])
        return paths

    def _plot_summary(self, results):
        """Generate and save summary plot of water-mediated HBond probabilities."""
        if not results["water_mediated"]:
            logger.warning("No water-mediated HBonds found — skipping plot.")
            return None

        logger.info(f"Number of water-mediated paths: {len(results['water_mediated'])}")

        df = pd.DataFrame({
            "path": [" -- ".join(p[0]) for p in results["water_mediated"]],
            "frame": [p[1] for p in results["water_mediated"]],
        })

        frame_count = len(self.u.trajectory[self.start:self.stop:self.step])
        counts = df["path"].value_counts(normalize=True) * 100
        df_plot = counts.reset_index()
        df_plot.columns = ["Water bridge", "Probability (%)"]

        df_plot = df_plot.head(self.top_n)

        # Plot water-mediated interactions

        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(8, max(4, len(df_plot) * 0.3)))
        sns.barplot(data=df_plot, y="Water bridge", x="Probability (%)", color="#5DADE2")
        plt.title("Water-mediated Hydrogen Bonds")
        plt.xlabel("Occurrence probability (%)")
        plt.ylabel("Residue ↔ Ligand Water Bridge")
        plt.tight_layout()

        out_plot = os.path.join(self.output_dir, "solvent_hbond_summary.png")
        plt.savefig(out_plot, bbox_inches="tight", dpi=300)
        plt.close()
        logger.info(f"Saved solvent hydrogen bond summary plot to {out_plot}")

        
        # Plot direct h-bonds
        
        if results["direct"]:
            df_direct = pd.DataFrame({
                "path": [" -- ".join(p[0]) for p in results["direct"]],
                "frame": [p[1] for p in results["direct"]],
            })
            counts = df_direct["path"].value_counts(normalize=True) * 100
            df_direct_plot = counts.reset_index()
            df_direct_plot.columns = ["Direct H-bond", "Probability (%)"]

            df_direct_plot = df_direct_plot[df_direct_plot["Probability (%)"] > 1.0]
            df_direct_plot = df_direct_plot.head(20)

            plt.figure(figsize=(8, max(4, len(df_direct_plot) * 0.3)))
            sns.barplot(data=df_direct_plot, y="Direct H-bond", x="Probability (%)", color="#58D68D")
            plt.title("Direct Hydrogen Bonds (Protein ↔ Ligand)")
            plt.xlabel("Occurrence probability (%)")
            plt.ylabel("Residue ↔ Ligand H-Bond")
            plt.tight_layout()

            out_plot_direct = os.path.join(self.output_dir, "direct_hbond_summary.png")
            plt.savefig(out_plot_direct, bbox_inches="tight", dpi=300)
            plt.close()
            logger.info(f"Saved direct hydrogen bond summary plot to {out_plot_direct}")
        else:
            logger.warning("No direct hydrogen bonds found — skipping direct plot.")

        return out_plot


    def run(self):
        logger.info("Calculating direct hydrogen bonds...")
        self.direct_paths = self._compute_hbonds(water_mediated=False)

        logger.info("Calculating water-mediated hydrogen bonds...")
        self.water_paths = self._compute_hbonds(water_mediated=True)

        results = {
            "direct": self.direct_paths,
            "water_mediated": self.water_paths,
            "binding_site_resids": self.binding_site_resids
        }

        plot_path = self._plot_summary(results)
        results["plot_file"] = plot_path

        # Dump JSON of plotted data for later use (all probs will be dumped here).

        out_json = os.path.join(self.output_dir, "solvent_hbonds.json")
        def _json_safe(obj):
            if isinstance(obj, (np.integer, np.int64, np.int32)):
                return int(obj)
            if isinstance(obj, (np.floating, np.float64, np.float32)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

        with open(out_json, "w") as f:
            json.dump(results, f, indent=2, default=_json_safe)
        logger.info(f"Saved solvent hydrogen bond data to {out_json}")

        return results

@register_task(
    "rmsd_analysis",
    category="Analyses",
    description="Compute RMSD of protein backbone and ligand separately over time."
)
class RMSDAnalysisTask:
    """
    Compute RMSD for protein backbone and ligand (if present), relative to the first frame.
    """

    def __init__(self,
                 topology: str,
                 trajectory: str,
                 ligand_resname: str = "LIG",
                 start: int = 0,
                 stop: int = -1,
                 step: int = 1,
                 output_dir: str = "output_rmsd",
                 context: Optional[Dict[str, Any]] = None):
        self.topology = topology
        self.trajectory = trajectory
        self.ligand_resname = ligand_resname
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.context = context or {}

        os.makedirs(self.output_dir, exist_ok=True)
        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.info("Starting RMSD analysis...")

        protein = self.u.select_atoms("protein and backbone")
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")

        if len(protein) == 0:
            raise ValueError("No protein atoms found in topology for RMSD calculation.")

        if len(ligand) == 0:
            logger.warning(f"No ligand atoms found for resname {self.ligand_resname}. Skipping ligand RMSD.")

        # Reference (first frame)
        ref = mda.Universe(self.topology)
        ref_protein = ref.select_atoms("protein and backbone")
        ref_ligand = ref.select_atoms(f"resname {self.ligand_resname}") if len(ligand) > 0 else None

        # Compute RMSD for protein backbone
        rmsd_protein = rms.RMSD(self.u, ref, select="protein and backbone",
                                start=self.start, stop=self.stop, step=self.step).run()

        df_rmsd = pd.DataFrame({
            "Time (ps)": rmsd_protein.rmsd[:, 1],
            "RMSD (Å)": rmsd_protein.rmsd[:, 2],
            "Component": "Protein (backbone)"
        })

        # Compute RMSD for ligand (optional)
        if len(ligand) > 0 and ref_ligand is not None:
            rmsd_ligand = rms.RMSD(self.u, ref, select=f"resname {self.ligand_resname}",
                                   start=self.start, stop=self.stop, step=self.step).run()
            df_ligand = pd.DataFrame({
                "Time (ps)": rmsd_ligand.rmsd[:, 1],
                "RMSD (Å)": rmsd_ligand.rmsd[:, 2],
                "Component": f"Ligand ({self.ligand_resname})"
            })
            df_rmsd = pd.concat([df_rmsd, df_ligand], ignore_index=True)

        # Save JSON data
        out_json = os.path.join(self.output_dir, "rmsd_data.json")
        df_rmsd.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved RMSD data to {out_json}")

        # Plot
        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df_rmsd, x="Time (ps)", y="RMSD (Å)", hue="Component", lw=2)
        plt.title("RMSD Over Time")
        plt.tight_layout()

        out_plot = os.path.join(self.output_dir, "rmsd_plot.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved RMSD plot to {out_plot}")

        return {"rmsd_data": out_json, "rmsd_plot": out_plot}

@register_task(
    "rmsf_analysis",
    category="Analyses",
    description="Compute per-residue RMSF for protein (all atoms and Cα only)."
)
class RMSFAnalysisTask:
    """
    Compute the Root Mean Square Fluctuation (RMSF) per residue for:
      - all protein atoms (averaged by residue)
      - alpha carbons (Cα only)
    
    Trajectory is aligned to the protein backbone reference before RMSF
    to remove overall rotational and translational motion.
    """

    def __init__(self,
                 topology: str,
                 trajectory: str,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 1,
                 output_dir: str = "output_rmsf",
                 context: Optional[Dict[str, Any]] = None):
        self.topology = topology
        self.trajectory = trajectory
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.context = context or {}

        os.makedirs(self.output_dir, exist_ok=True)
        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.info("Starting RMSF analysis...")

        protein = self.u.select_atoms("protein")
        calphas = self.u.select_atoms("protein and name CA")

        if len(protein) == 0:
            raise ValueError("No protein atoms found for RMSF calculation.")

        # Align trajectory to backbone reference
        logger.info("Aligning trajectory to protein backbone reference...")
        align.AlignTraj(self.u, self.u, select="protein and backbone",
                        in_memory=True).run()

        # ------------------------------------------------------------------
        # RMSF for all protein atoms (average per residue)
        # ------------------------------------------------------------------
        rmsf_all = rms.RMSF(protein).run()
        df_all_atoms = pd.DataFrame({
            "Residue": [atom.resid for atom in protein.atoms],
            "RMSF (Å)": rmsf_all.rmsf,
        })

        df_all = df_all_atoms.groupby("Residue", as_index=False).mean()
        df_all["Component"] = "Protein (all atoms)"

        # ------------------------------------------------------------------
        # RMSF for alpha carbons (Cα only)
        # ------------------------------------------------------------------
        rmsf_ca = rms.RMSF(calphas).run()
        df_ca = pd.DataFrame({
            "Residue": [atom.resid for atom in calphas.atoms],
            "RMSF (Å)": rmsf_ca.rmsf,
            "Component": "Cα only"
        })

        # Combine
        df_rmsf = pd.concat([df_all, df_ca], ignore_index=True)
        logger.info(f"Computed RMSF for {len(df_all)} residues (all atoms) "
                    f"and {len(df_ca)} alpha carbons.")

        # ------------------------------------------------------------------
        # Save results
        # ------------------------------------------------------------------
        out_json = os.path.join(self.output_dir, "rmsf_data.json")
        df_rmsf.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved RMSF data to {out_json}")

        # ------------------------------------------------------------------
        # Plot per-residue RMSF
        # ------------------------------------------------------------------
        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(10, 5))
        sns.lineplot(data=df_rmsf, x="Residue", y="RMSF (Å)", hue="Component", lw=2)
        plt.title("Per-Residue RMSF")
        plt.tight_layout()

        out_plot = os.path.join(self.output_dir, "rmsf_plot.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved RMSF plot to {out_plot}")

        return {"rmsf_data": out_json, "rmsf_plot": out_plot}

