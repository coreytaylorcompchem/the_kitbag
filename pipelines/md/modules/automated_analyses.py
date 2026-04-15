import os
import warnings
import json
import itertools
import h5py
import re

import networkx as nx

from tqdm import tqdm
import numpy as np
import pandas as pd

from collections import Counter

from typing import List, Optional, Dict, Any

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from modules.utils.mda_utils import _bit_to_color_value, _get_inv_color_mapper, _get_color_mapper
from modules.utils.component_detection import ComponentDetector

import prolif as plf
from prolif.plotting.utils import separated_interaction_colors

# suppress MDA warnings.

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis")

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.transformations import unwrap, center_in_box, wrap

from sklearn.decomposition import IncrementalPCA, PCA
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN

from node2vec import Node2Vec
import community as community_louvain

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


@register_task(
    "solvent_hbonds",
    category="Post-proc; graph analyses",
    description="Compute direct and water-mediated hbonds with ligand."
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
                 context: Optional[Dict[str, Any]] = None
                 ):
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

        self.u = mda.Universe(self.topology, self.trajectory)

        # Automatically detect binding site if not provided
        if not self.binding_site_resids:
            self.u.trajectory[0]
            sel = self.u.select_atoms(f"protein and around {self.binding_site_cutoff} resname {self.ligand_resname}")
            self.binding_site_resids = list(sel.resids)
            logger.info(f"Detected {len(self.binding_site_resids)} binding site residues.")

        self.direct_paths = []
        self.water_paths = []
        self.detector = ComponentDetector(
            self.u,
            ligand_resname=self.ligand_resname if hasattr(self, "ligand_resname") else "UNK"
        )

        self.components = self.context.get("components")

        if self.components is None:
            self.components = self.detector.detect()
            self.context["components"] = self.components

    def _compute_hbonds(self, water_mediated: bool = False):
        paths = []
        binding_site_atoms = self.u.select_atoms(
            "protein and resid " + " ".join(str(r) for r in self.binding_site_resids)
        )

        for ts in tqdm(self.u.trajectory[self.start:self.stop:self.step],
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
            logger.warning("No water-mediated hbonds found — skipping plot.")
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
    category="Post-proc; traj analyses",
    description="Compute RMSD of protein bb and ligand."
)
class RMSDAnalysisTask:
    """
    Compute RMSD for protein backbone and ligand (if present) relative to the first frame.
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
        # Show trajectory loading progress (non-intrusive)
        if isinstance(self.trajectory, list):
            logger.info(f"Preparing to load {len(self.trajectory)} trajectory files...")

            for traj_file in tqdm(self.trajectory, desc="Loading trajectories", unit="file"):
                logger.debug(f"Queued: {traj_file}")

        # Now load as usual (no behavior change)
        self.u = mda.Universe(self.topology, self.trajectory)

        protein = self.u.select_atoms("protein")

        self.u.trajectory.add_transformations(
            unwrap(self.u.atoms),              # fix broken systems
            center_in_box(protein, wrap=True)  # keep protein centered
        )
        
        # Detect system type automatically
        self.detector = ComponentDetector(
            self.u,
            ligand_resname=self.ligand_resname if hasattr(self, "ligand_resname") else "UNK"
        )

        self.components = self.context.get("components")

        if self.components is None:
            self.components = self.detector.detect()
            self.context["components"] = self.components
    
    def run(self):

        # More robusterer selection strings from segments

        def segid_backbone_selection(atomgroup):
            if atomgroup is None or len(atomgroup) == 0:
                return None
            segids = sorted(set(atomgroup.segids))
            seg_str = " or ".join([f"segid {seg}" for seg in segids])
            return f"({seg_str}) and backbone"

        def segid_all_selection(atomgroup):
            if atomgroup is None or len(atomgroup) == 0:
                return None
            segids = sorted(set(atomgroup.segids))
            seg_str = " or ".join([f"segid {seg}" for seg in segids])
            return f"({seg_str})"

        logger.debug("Starting RMSD analysis...")

        # Use detected components
        comp = self.components
        receptor = comp["receptor"]
        partner = comp["partner"]
        ligand = comp["ligand"]

        protein = self.u.select_atoms("protein and backbone")

        if len(protein) == 0:
            raise ValueError("No protein atoms found in topology for RMSD calculation.")

        # Final selection strings (MDAnalysis weirdness)
        rec_sel = segid_backbone_selection(receptor)
        partner_sel = segid_backbone_selection(partner) if comp["has_partner"] else None
        lig_sel = f"resname {self.ligand_resname}" if comp["has_ligand"] else None

        # Safer alignment selection string from receptor segments

        align_selection = rec_sel if rec_sel is not None else "protein and backbone"

        logger.info("Starting trajectory alignment (this may take a few minutes)...")

        # --- Reference frame (frame 0) ---
        self.u.trajectory[0]
        ref_atoms = self.u.select_atoms(align_selection)
        ref_pos = ref_atoms.positions.copy()

        # Preselect mobile atoms ONCE (important for speed)
        mobile_atoms = self.u.select_atoms(align_selection)

        n_frames = len(self.u.trajectory)

        for ts in tqdm(self.u.trajectory, total=n_frames, desc="Aligning trajectory", unit="frame"):

            # Compute optimal rotation + translation
            R, rmsd_val = align.rotation_matrix(
                mobile_atoms.positions - mobile_atoms.center_of_mass(),
                ref_pos - ref_atoms.center_of_mass()
            )

            # Apply transformation to ALL atoms (not just selection!)
            self.u.atoms.positions = np.dot(
                self.u.atoms.positions - mobile_atoms.center_of_mass(), R.T
            ) + ref_atoms.center_of_mass()

        logger.info("Alignment complete.")

        # Ligand detection
        if ligand is None:
            logger.warning("No ligand detected. Skipping ligand RMSD.")

        # Reference (frame0)
        ref = mda.Universe(self.topology, self.trajectory)
        ref.trajectory[0]

        # Align reference once
        logger.info("Aligning reference structure...")

        align.AlignTraj(ref, ref, select=align_selection, in_memory=True).run()

        logger.info("Reference alignment complete.")

        # Extract reference positions
        ref_atoms = {
            "receptor": ref.select_atoms(rec_sel) if rec_sel else None,
            "partner": ref.select_atoms(partner_sel) if partner_sel else None,
            "ligand": ref.select_atoms(lig_sel) if lig_sel else None,
        }

        dfs = []

        # Receptor RMSD by AtomGroups
        times = []
        rmsd_data = {
            "Receptor": [],
            "Partner (VHH)": [],
            f"Ligand ({self.ligand_resname})": []
        }

        # Frame indices
        traj = self.u.trajectory[self.start:self.stop:self.step]
        n_frames = len(traj)

        logger.info(f"Processing {n_frames} frames for RMSD...")

        for ts in tqdm(
                traj,
                desc="RMSD (all trajectories)",
                unit="frame",
                mininterval=1.0  # avoids too frequent updates
            ):
            time_ns = ts.time / 1000.0
            times.append(time_ns)

            # Receptor RMSD
            if rec_sel:
                mob = self.u.select_atoms(rec_sel)
                ref_sel = ref_atoms["receptor"]
                val = rms.rmsd(mob.positions, ref_sel.positions, center=True, superposition=True)
                rmsd_data["Receptor"].append(val)

            # Partner RMSD
            if partner_sel:
                mob = self.u.select_atoms(partner_sel)
                ref_sel = ref_atoms["partner"]
                val = rms.rmsd(mob.positions, ref_sel.positions, center=True, superposition=True)
                rmsd_data["Partner (VHH)"].append(val)

            # Ligand RMSD
            if lig_sel:
                mob = self.u.select_atoms(lig_sel)
                ref_sel = ref_atoms["ligand"]
                val = rms.rmsd(mob.positions, ref_sel.positions, center=True, superposition=True)
                rmsd_data[f"Ligand ({self.ligand_resname})"].append(val)

            for comp, values in rmsd_data.items():
                if len(values) == 0:
                    continue
                dfs.append(pd.DataFrame({
                    "Time (ns)": times,
                    "RMSD (Å)": values,
                    "Component": comp
                }))

            df_rmsd = pd.concat(dfs, ignore_index=True)

        # Dump plot data to JSON
        out_json = os.path.join(self.output_dir, "rmsd_data.json")
        df_rmsd.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved RMSD data to {out_json}")

        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df_rmsd, x="Time (ns)", y="RMSD (Å)", hue="Component", lw=2)
        plt.title("RMSD over Time")
        plt.xlabel("Simulation Time (ns)")

        out_plot = os.path.join(self.output_dir, "rmsd_plot.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved RMSD plot to {out_plot}")

        return {"rmsd_data": out_json, "rmsd_plot": out_plot}

@register_task(
    "rmsf_analysis",
    category="Post-proc; traj analyses",
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
        
        # Detect system type automatically
        self.detector = ComponentDetector(
            self.u,
            ligand_resname=self.ligand_resname if hasattr(self, "ligand_resname") else "UNK"
        )

        self.components = self.context.get("components")

        if self.components is None:
            self.components = self.detector.detect()
            self.context["components"] = self.components

    def run(self):
        logger.debug("Starting RMSF analysis...")

        comp = self.components
        receptor = comp["receptor"]
        partner = comp["partner"]

        protein = self.u.select_atoms("protein")
        calphas = self.u.select_atoms("protein and name CA")

        if len(protein) == 0:
            raise ValueError("No protein atoms found for RMSF calculation.")

        # Align trajectory to bb reference
        logger.debug("Aligning trajectory to protein backbone reference...")
        align.AlignTraj(self.u, self.u, select="protein and backbone",
                        in_memory=True).run()

        # Per-residue RMSF
        dfs = []

        # helper; compute RMSF safely
        def compute_rmsf(atomgroup, label):
            atoms = atomgroup.atoms
            rmsf_vals = rms.RMSF(atoms).run()

            if len(atoms) != len(rmsf_vals.rmsf):
                raise ValueError(
                    f"Mismatch: {len(atoms)} atoms vs {len(rmsf_vals.rmsf)} RMSF values"
                )

            df_atoms = pd.DataFrame({
                "Residue": [int(atom.resid) for atom in atoms],
                "Chain": [str(atom.segid) for atom in atoms],
                "RMSF (Å)": rmsf_vals.rmsf,
            })

            df_res = df_atoms.groupby(["Chain", "Residue"], as_index=False).mean()
            df_res["Component"] = label

            return df_res
    
        # Cα RMSF per component
        def compute_rmsf_ca(atomgroup, label):
            if atomgroup is None or len(atomgroup) == 0:
                return None
            atoms = atomgroup.select_atoms("name CA")
            if len(atoms) == 0:
                return None

            rmsf_vals = rms.RMSF(atoms).run()

            df_atoms = pd.DataFrame({
                "Residue": [int(atom.resid) for atom in atoms],
                "Chain": [str(atom.segid) for atom in atoms],
                "RMSF (Å)": rmsf_vals.rmsf,
                "Component": label
            })

            # Group only the RMSF values by Chain+Residue
            df_res = df_atoms.groupby(["Chain", "Residue"], as_index=False)[["RMSF (Å)"]].mean()
            df_res["Component"] = label  # re-add the label after aggregation

            return df_res

        # Receptor RMSF
        if receptor is not None:
            dfs.append(compute_rmsf(receptor, "Receptor (all atoms)"))
            dfs.append(compute_rmsf_ca(receptor, "Receptor (Cα only)"))

        # Partner RMSF
        if comp["has_partner"]:
            dfs.append(compute_rmsf(partner, "Partner (VHH, all atoms)"))
            dfs.append(compute_rmsf_ca(partner, "Partner (VHH, Cα only)"))

        # # Cα RMSF
        # atoms_ca = calphas.atoms
        # rmsf_ca = rms.RMSF(atoms_ca).run()

        # df_ca = pd.DataFrame({
        #     "Residue": [int(atom.resid) for atom in atoms_ca],
        #     "Chain": [str(atom.segid) for atom in atoms_ca],
        #     "RMSF (Å)": rmsf_ca.rmsf,
        #     "Component": "Cα only (global)"
        # })

        # dfs.append(df_ca)

        df_rmsf = pd.concat(dfs, ignore_index=True)

        # ============================================================
        # Combine residue indexing across chains
        # ============================================================

        df_rmsf = df_rmsf.sort_values(["Chain", "Residue"]).reset_index(drop=True)

        chain_offsets = {}
        offset = 0

        for chain in df_rmsf["Chain"].unique():
            chain_df = df_rmsf[df_rmsf["Chain"] == chain]

            chain_min = chain_df["Residue"].min()
            chain_max = chain_df["Residue"].max()

            chain_offsets[chain] = offset - chain_min + 1
            offset += (chain_max - chain_min + 1)

        df_rmsf["Residue_cont"] = df_rmsf.apply(
            lambda row: row["Residue"] + chain_offsets[row["Chain"]],
            axis=1
        )

        logger.debug(f"Computed RMSF for {len(df_rmsf)} total entries.")

        # Save data
        out_json = os.path.join(self.output_dir, "rmsf_data.json")
        df_rmsf.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved RMSF data to {out_json}")

        # ============================================================
        # Plotting
        # ============================================================

        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(12, 5))

        # Plot RMSF lines
        sns.lineplot(
            data=df_rmsf,
            x="Residue_cont",
            y="RMSF (Å)",
            hue="Component",
            lw=2
        )

        # Build chain:resid labels for x-axis
        df_rmsf["Label"] = df_rmsf["Chain"] + ":" + df_rmsf["Residue"].astype(str)

        # Use all protein components for ticks
        protein_components = ["Receptor (all atoms)", "Partner (VHH, all atoms)"]
        df_ticks = df_rmsf[df_rmsf["Component"].isin(protein_components)]
        df_ticks = df_ticks.drop_duplicates(subset=["Residue_cont"]).sort_values("Residue_cont")

        # Subsample ticks to avoid overcrowding
        n_ticks = 20
        step = max(1, len(df_ticks) // n_ticks)
        tick_positions = df_ticks["Residue_cont"].iloc[::step]
        tick_labels = df_ticks["Label"].iloc[::step]

        plt.xticks(tick_positions, tick_labels, rotation=45, ha="right", fontsize=8)

        # Show chain boundaries
        for chain, off in chain_offsets.items():
            plt.axvline(off, linestyle="--", alpha=0.3)

        plt.grid(True, axis="y")
        plt.grid(False, axis="x")
        plt.title("Per-Residue RMSF")
        plt.xlabel("Residue (chain:resid)")
        plt.tight_layout()

        out_plot = os.path.join(self.output_dir, "rmsf_plot.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved RMSF plot to {out_plot}")

        return {"rmsf_data": out_json, "rmsf_plot": out_plot}

def _bit_to_color_value(series):
    return series.astype(int)

@register_task(
    "interaction_fingerprint",
    category="Post-proc; traj analyses",
    description="Compute IFP using ProLIF, generate barcode plot."
)
class InteractionFingerprintTask:
    """
    Generate a time-series interaction 'barcode' plot using ProLIF.
    """

    def __init__(
        self,
        topology: str,
        trajectory: str,
        ligand_selection: str = "resname UNK",
        protein_selection: str = "(protein or resname WAT) and byres around 20.0 group ligand",
        start: int = 0,
        stop: int = -1,
        step: int = 1,
        output_dir: str = "output_fp",
        n_frame_ticks: int = 10,
        residues_tick_location: str = "top",
        figsize: tuple = (8, 10),
        dpi: int = 100,
        context: Optional[Dict[str, Any]] = None,
    ):
        self.topology = topology
        self.trajectory = trajectory
        self.ligand_selection = ligand_selection
        self.protein_selection = protein_selection
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.n_frame_ticks = n_frame_ticks
        self.residues_tick_location = residues_tick_location
        self.figsize = figsize
        self.dpi = dpi
        self.context = context or {}

        os.makedirs(self.output_dir, exist_ok=True)
        self.u = mda.Universe(self.topology, self.trajectory)

        self.detector = ComponentDetector(self.u)

        self.components = self.context.get("components")

        if self.components is None:
            self.components = self.detector.detect()
            self.context["components"] = self.components

    def run(self):
        logger.debug("Starting ProLIF interaction fingerprint analysis...")

        comp = self.components

        if comp["has_ligand"]:
            ligand = comp["ligand"]
            protein = self.u.select_atoms(
                "(protein or resname WAT) and byres around 20.0 group ligand",
                ligand=ligand
            )
            mode = "ligand"
        else:
            if not comp["has_partner"]:
                raise ValueError("No ligand or partner detected for interaction analysis.")

            # Treat VHH as ligand
            ligand = comp["partner"]
            protein = comp["receptor"]
            mode = "protein-protein"

        if len(protein) == 0:
            raise ValueError("No atoms found for interaction partner selection.")

        # protein = self.u.select_atoms(self.protein_selection, ligand=ligand)
        # if len(protein) == 0:
        #     raise ValueError(f"No atoms found for protein selection: {self.protein_selection}")

        # Prolif FP calculation
        logger.debug("Running ProLIF fingerprint calculation...")
        
        fp = plf.Fingerprint()

        ligand.guess_bonds()
        protein.guess_bonds()

        fp.run(self.u.trajectory[self.start:self.stop:self.step], ligand, protein)
        fp_df = fp.to_dataframe()
        
        logger.debug(fp_df.columns.get_level_values(2).unique()) #debug to check all interactions are being found

        out_data_path = os.path.join(self.output_dir, "fp_data.pkl")
        fp_df.to_pickle(out_data_path)
        logger.info(f"Saved fingerprint DataFrame to {out_data_path}")

        # Generate barcode plot
        logger.info("Generating interaction barcode plot...")

        # As Prolif output doesn't return time in ps, we need to reconstruct it from the trajectories
        # And then convert to ns for plotting.
        logger.debug("Reconstructing continuous simulation time across trajectory segments...")

        times_ns = []
        offset_ns = 0.0

        times_ns = np.array([
            ts.time / 1000.0  # ps → ns
            for ts in self.u.trajectory[self.start:self.stop:self.step]
        ])

        fp_transposed = fp_df.astype(np.uint8).T.apply(_bit_to_color_value, axis=1)

        color_mapper = _get_color_mapper()  
        inv_color_mapper = _get_inv_color_mapper()

        # Convert MultiIndex df to integer array for plotting
        plot_data = fp_transposed.copy()

        for idx in plot_data.index:
            interaction_type = idx[2]  # i.e. third level of MultiIndex
            color_mapper.get(interaction_type, 0)
            plot_data.loc[idx] = plot_data.loc[idx] * color_mapper.get(interaction_type, 0)

        plot_data_array = plot_data.values.astype(int)

        interaction_order = [None] + list(separated_interaction_colors.keys())

        cmap = ListedColormap(
            ["white"] + [separated_interaction_colors[name] for name in separated_interaction_colors]
        )

        sns.set(style="whitegrid", context="talk")
        fig, ax = plt.subplots(1, 1, figsize=self.figsize, dpi=self.dpi)

        im = ax.imshow(
            plot_data_array,
            aspect="auto",
            interpolation="none",
            cmap=cmap,
            vmin=0,
            vmax=max(color_mapper.values()),
        )

        # Format X-axis: nicer frame ticks
        frames = fp_transposed.columns.astype(int)
        num_ticks = min(self.n_frame_ticks, len(frames)) 
        tick_indices = np.round(np.linspace(0, len(frames) - 1, num_ticks)).astype(int)

        times_ns_int = np.round(times_ns).astype(int)

        ax.set_xticks(tick_indices)
        ax.set_xticklabels(times_ns_int[tick_indices])
        ax.set_xlabel("Simulation Time (ns)")

        # Format Y-axis
        residues = [f"{prot}" for lig, prot, inter in fp_transposed.index]
        ax.set_yticks(np.arange(len(residues)))
        ax.set_yticklabels(residues)

        # Legend
        unique_values = np.unique(plot_data_array)
        unique_values = [v for v in unique_values if v != 0]
        legend_colors = {inv_color_mapper[v]: im.cmap(v) for v in unique_values}
        patches = [Patch(color=color, label=interaction) for interaction, color in legend_colors.items()]
        ax.legend(handles=patches, loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=False, ncol=3)
        
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave 5% at bottom for legend

        out_plot_path = os.path.join(self.output_dir, "interaction_barcode.png")

        ligands = sorted({lig[:3] for lig, prot, inter in fp_transposed.index})
        ligand_str = ", ".join(ligands)

        if mode == "ligand":
            ax.set_title(f"Ligand interaction fingerprint for ligand {ligand_str}.")
        else:
            ax.set_title("Protein–protein interaction fingerprint (receptor vs partner).")
        plt.tight_layout()
        plt.savefig(out_plot_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        logger.info(f"Saved interaction barcode plot to {out_plot_path}")

        return {"fp_data": out_data_path, "fp_plot": out_plot_path}

@register_task(
    "cluster_analysis",
    category="Post-proc; traj analyses",
    description="Cluster trajectory based on RMSD and output centroid structures."
)
class ClusterAnalysisTask:
    """
    Perform RMSD-based clustering on a trajectory and output:
      - cluster populations
      - centroid structures (PDB)
    """

    def __init__(self,
                 topology: str,
                 trajectory: str,
                 selection: str = "protein and backbone",
                 n_clusters: int = 5,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 1,
                 output_dir: str = "output_clusters",
                 context: Optional[Dict[str, Any]] = None):

        self.topology = topology
        self.trajectory = trajectory
        self.selection = selection
        self.n_clusters = n_clusters
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.context = context or {}

        os.makedirs(self.output_dir, exist_ok=True)

        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.info("Starting cluster analysis...")

        # --------------------------------------------------
        # Align trajectory
        # --------------------------------------------------
        logger.debug("Aligning trajectory...")
        align.AlignTraj(
            self.u, self.u,
            select=self.selection,
            in_memory=True
        ).run()

        atoms = self.u.select_atoms(self.selection)

        # --------------------------------------------------
        # Collect coordinates
        # --------------------------------------------------
        coords = []
        frame_indices = []

        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            coords.append(atoms.positions.copy())
            frame_indices.append(ts.frame)

        coords = np.array(coords)  # shape: (n_frames, n_atoms, 3)

        n_frames = coords.shape[0]
        logger.info(f"Collected {n_frames} frames for clustering.")

        # --------------------------------------------------
        # Compute RMSD distance matrix
        # --------------------------------------------------
        logger.debug("Computing RMSD distance matrix...")

        dist_matrix = np.zeros((n_frames, n_frames))

        for i in range(n_frames):
            for j in range(i + 1, n_frames):
                d = rms.rmsd(coords[i], coords[j], superposition=False)
                dist_matrix[i, j] = d
                dist_matrix[j, i] = d

        # --------------------------------------------------
        # Clustering (Agglomerative)
        # --------------------------------------------------
        logger.debug("Running agglomerative clustering...")

        from sklearn.cluster import AgglomerativeClustering

        clustering = AgglomerativeClustering(
            n_clusters=self.n_clusters,
            metric="precomputed",
            linkage="average"
        )

        labels = clustering.fit_predict(dist_matrix)

        # --------------------------------------------------
        # Cluster statistics
        # --------------------------------------------------
        df = pd.DataFrame({
            "Frame": frame_indices,
            "Cluster": labels
        })

        cluster_counts = df["Cluster"].value_counts().sort_index()
        cluster_props = cluster_counts / cluster_counts.sum()

        df_summary = pd.DataFrame({
            "Cluster": cluster_counts.index,
            "Count": cluster_counts.values,
            "Proportion": cluster_props.values
        })

        out_json = os.path.join(self.output_dir, "cluster_summary.json")
        df_summary.to_json(out_json, orient="records", indent=2)

        logger.info(f"Saved cluster summary to {out_json}")

        # --------------------------------------------------
        # Find centroids
        # --------------------------------------------------
        logger.debug("Extracting cluster centroids...")

        centroids = {}

        for clust_id in sorted(df["Cluster"].unique()):
            cluster_frames = np.where(labels == clust_id)[0]

            submatrix = dist_matrix[np.ix_(cluster_frames, cluster_frames)]

            # centroid = frame with minimal average distance
            avg_dist = submatrix.mean(axis=1)
            centroid_idx_local = np.argmin(avg_dist)
            centroid_idx_global = cluster_frames[centroid_idx_local]

            centroids[clust_id] = centroid_idx_global

        # --------------------------------------------------
        # Write centroid PDBs
        # --------------------------------------------------
        logger.debug("Writing centroid structures...")

        for clust_id, frame_idx in centroids.items():
            self.u.trajectory[frame_idx]

            out_pdb = os.path.join(
                self.output_dir,
                f"cluster_{clust_id}_centroid.pdb"
            )

            with mda.Writer(out_pdb, self.u.atoms.n_atoms) as W:
                W.write(self.u.atoms)

            logger.info(f"Saved centroid for cluster {clust_id} → {out_pdb}")

        return {
            "cluster_summary": out_json,
            "labels": labels.tolist()
        }

@register_task(
    "free_energy_landscape",
    category="Post-proc; traj analyses",
    description="Compute 2D free energy landscape using PCA projection."
)
class FreeEnergyLandscapeTask:
    """
    Compute a 2D Free Energy Landscape (FEL) using PCA.

    Outputs:
      - FEL grid (numpy)
      - raw projection data (JSON)
      - contour plot (PNG)
    """

    def __init__(self,
                 topology: str,
                 trajectory: str,
                 selection: str = "protein and backbone",
                 n_bins: int = 50,
                 temperature: float = 298.0,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 1,
                 output_dir: str = "output_fel",
                 context: Optional[Dict[str, Any]] = None):

        self.topology = topology
        self.trajectory = trajectory
        self.selection = selection
        self.n_bins = n_bins
        self.temperature = temperature
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.context = context or {}

        os.makedirs(self.output_dir, exist_ok=True)

        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.info("Starting free energy landscape analysis...")
        cluster_labels = self.context.get("cluster_labels", None)

        # --------------------------------------------------
        # Align trajectory
        # --------------------------------------------------
        logger.debug("Aligning trajectory...")
        align.AlignTraj(
            self.u, self.u,
            select=self.selection,
            in_memory=True
        ).run()

        atoms = self.u.select_atoms(self.selection)

        # --------------------------------------------------
        # Collect coordinates
        # --------------------------------------------------
        coords = []

        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            coords.append(atoms.positions.flatten())

        coords = np.array(coords)
        logger.info(f"Collected {coords.shape[0]} frames.")

        # --------------------------------------------------
        # PCA
        # --------------------------------------------------
        logger.debug("Running PCA...")
        from sklearn.decomposition import PCA

        pca = PCA(n_components=2)
        proj = pca.fit_transform(coords)

        x = proj[:, 0]
        y = proj[:, 1]

        # --------------------------------------------------
        # Histogram → probability
        # --------------------------------------------------
        logger.debug("Computing probability density...")
        H, xedges, yedges = np.histogram2d(
            x, y,
            bins=self.n_bins,
            density=True
        )

        # --------------------------------------------------
        # Free energy
        # --------------------------------------------------
        kB = 0.0019872041  # kcal/mol/K
        kT = kB * self.temperature

        F = -kT * np.log(H + 1e-12)
        F -= np.nanmin(F)

        # --------------------------------------------------
        # Save data
        # --------------------------------------------------
        np.save(os.path.join(self.output_dir, "fel_grid.npy"), F)

        df = pd.DataFrame({
            "PC1": x,
            "PC2": y
        })

        out_json = os.path.join(self.output_dir, "fel_projection.json")
        df.to_json(out_json, orient="records", indent=2)

        logger.info(f"Saved FEL projection to {out_json}")

        # ==================================================
        # Cluster thermodynamics
        # ==================================================
        if cluster_labels is not None:
            logger.debug("Computing cluster free energies...")

            cluster_labels = np.array(cluster_labels)

            counts = pd.Series(cluster_labels).value_counts().sort_index()
            probs = counts / counts.sum()

            free_energy = -kT * np.log(probs)

            df_clusters = pd.DataFrame({
                "Cluster": counts.index,
                "Population": probs.values,
                "FreeEnergy (kcal/mol)": free_energy.values
            })

            out_clusters = os.path.join(self.output_dir, "cluster_fel_summary.json")
            df_clusters.to_json(out_clusters, orient="records", indent=2)

            logger.info(f"Saved cluster FEL summary to {out_clusters}")

        # --------------------------------------------------
        # Plot
        # --------------------------------------------------
        logger.debug("Plotting FEL...")

        X, Y = np.meshgrid(xedges[:-1], yedges[:-1])

        plt.figure(figsize=(7, 6))

        contour = plt.contourf(
            X, Y, F.T,
            levels=50
        )

        plt.colorbar(contour, label="Free Energy (kcal/mol)")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("Free Energy Landscape")

        # ==================================================
        # Overlay clusters
        # ==================================================
        if cluster_labels is not None:
            logger.debug("Overlaying cluster labels on FEL...")

            plt.scatter(
                x, y,
                c=cluster_labels,
                s=5,
                alpha=0.5
            )

            # ----------------------------------------------
            # Cluster centroids
            # ----------------------------------------------

            centroids = {}

            for clust_id in np.unique(cluster_labels):
                mask = cluster_labels == clust_id
                cx = x[mask].mean()
                cy = y[mask].mean()

                centroids[clust_id] = (cx, cy)

                plt.text(
                    cx, cy,
                    str(clust_id),
                    fontsize=12,
                    weight="bold",
                    ha="center"
                )

            # ==================================================
            # Transition arrows between clusters
            # ==================================================
            logger.debug("Computing cluster transitions...")

            # Count transitions
            transition_counts = {}

            for i in range(len(cluster_labels) - 1):
                a = cluster_labels[i]
                b = cluster_labels[i + 1]

                if a == b:
                    continue  # skip self-transitions

                key = (a, b)
                transition_counts[key] = transition_counts.get(key, 0) + 1

            if transition_counts:
                max_count = max(transition_counts.values())

                for (a, b), count in transition_counts.items():
                    x1, y1 = centroids[a]
                    x2, y2 = centroids[b]

                    # Normalise arrow thickness
                    width = 1.0 + 4.0 * (count / max_count)

                    plt.arrow(
                        x1, y1,
                        x2 - x1,
                        y2 - y1,
                        alpha=0.4,
                        width=0.0,
                        linewidth=width,
                        length_includes_head=True,
                        head_width=0.1 * np.std(x),
                        head_length=0.1 * np.std(y),
                        color="black"
                    )

        out_plot = os.path.join(self.output_dir, "fel_plot.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved FEL plot to {out_plot}")

        return {
            "fel_projection": out_json,
            "fel_plot": out_plot,
            "fel_grid": os.path.join(self.output_dir, "fel_grid.npy")
        }

from deeptime.decomposition import TICA
from deeptime.clustering import KMeans
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.markov import TransitionCountEstimator

@register_task(
    "msm_analysis",
    category="Post-proc; kinetics",
    description="Rigorous MSM using TICA + clustering + MLE."
)
class MSMAnalysisTask:

    def __init__(self,
                 topology,
                 trajectory,
                 selection="protein and backbone",
                 lagtime=10,
                 n_clusters=100,
                 output_dir="output_msm",
                 context=None):

        self.topology = topology
        self.trajectory = trajectory
        self.selection = selection
        self.lagtime = lagtime
        self.n_clusters = n_clusters
        self.output_dir = output_dir
        self.context = context or {}

        os.makedirs(self.output_dir, exist_ok=True)
        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.info("Starting MSM analysis (deeptime)...")

        atoms = self.u.select_atoms(self.selection)

        # --------------------------------------------------
        # Alignment before featurisation
        # --------------------------------------------------

        align.AlignTraj(
            self.u, self.u,
            select=self.selection,
            in_memory=True
        ).run()

        # --------------------------------------------------
        # Featurisation (flattened coords for now)
        # --------------------------------------------------
        X = []
        for ts in self.u.trajectory:
            X.append(atoms.positions.flatten())
        X = np.array(X)

        # --------------------------------------------------
        # TICA
        # --------------------------------------------------
        tica = TICA(lagtime=self.lagtime, dim=5)
        X_tica = tica.fit_transform(X)

        # --------------------------------------------------
        # Clustering
        # --------------------------------------------------
        kmeans = KMeans(n_clusters=self.n_clusters)
        dtrajs = kmeans.fit_fetch(X_tica)

        # --------------------------------------------------
        # Count matrix
        # --------------------------------------------------
        counts_estimator = TransitionCountEstimator(lagtime=self.lagtime)
        counts = counts_estimator.fit_fetch(dtrajs)

        # --------------------------------------------------
        # MSM
        # --------------------------------------------------
        msm = MaximumLikelihoodMSM(reversible=True)
        model = msm.fit_fetch(counts)

        # --------------------------------------------------
        # Outputs
        # --------------------------------------------------
        np.save(os.path.join(self.output_dir, "transition_matrix.npy"),
                model.transition_matrix)

        np.save(os.path.join(self.output_dir, "stationary_distribution.npy"),
                model.stationary_distribution)

        np.save(os.path.join(self.output_dir, "timescales.npy"),
                model.timescales())

        logger.info("MSM analysis complete.")

        return {
            "transition_matrix": "transition_matrix.npy",
            "stationary_distribution": "stationary_distribution.npy",
            "timescales": "timescales.npy"
        }

@register_task(
    "protein_ligand_communities",
    category="Post-proc; graph analyses",
    description="Detect cooperative residue clusters (communities) in the protein–ligand interaction network."
)
class ProteinLigandCommunityTask:
    def __init__(self,
                 topology: str,
                 trajectory: str,
                 ligand_resname: str = "UNK",
                 start: int = 0,
                 stop: int = -1,
                 step: int = 10,
                 binding_site_cutoff: float = 5.0,
                 output_dir: str = "output_plin_communities",
                 context: Optional[Dict[str, Any]] = None):
        self.topology = topology
        self.trajectory = trajectory
        self.ligand_resname = ligand_resname
        self.start = start
        self.stop = stop
        self.step = step
        self.binding_site_cutoff = binding_site_cutoff
        self.output_dir = output_dir
        self.context = context or {}
        os.makedirs(self.output_dir, exist_ok=True)
        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.info("Building time-averaged residue–ligand interaction network...")
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        if len(ligand) == 0:
            raise ValueError(f"No ligand found for resname {self.ligand_resname}")
        protein = self.u.select_atoms("protein")

        contact_matrix = {}  # resid -> list of 0/1 contacts over time
        frames = list(self.u.trajectory[self.start:self.stop:self.step])

        for res in protein.residues:
            contact_matrix[res.resid] = []

        for ts in tqdm(frames, desc="Contacts"):
            for res in protein.residues:
                d = mda.lib.distances.distance_array(
                    res.atoms.positions, ligand.positions
                ).min()

                # Smooth contact score instead of binary
                contact_matrix[res.resid].append(d)

        contact_freq = {
            resid: np.mean(vals)
            for resid, vals in contact_matrix.items()
        }
        
        # ============================================================
        # Contact persistence
        # ============================================================

        logger.info("Computing residue–ligand contact persistence...")

        contact_counts = {}
        n_frames = 0

        for ts in tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Contacts"):
            n_frames += 1

            for res in protein.residues:
                d = mda.lib.distances.distance_array(
                    res.atoms.positions, ligand.positions
                ).min()

                if d < self.binding_site_cutoff:
                    key = (res.segid, res.resid, res.resname)
                    contact_counts[key] = contact_counts.get(key, 0) + 1

        data = []

        for (chain, resid, resname), count in contact_counts.items():
            persistence = count / n_frames

            data.append({
                "Chain": chain,
                "Residue": resid,
                "Resname": resname,
                "Label": f"{resname}{resid}:{chain}",
                "Persistence": persistence
            })

        df = pd.DataFrame(data)

        if df.empty:
            raise ValueError("No contacts detected. Check ligand selection or cutoff.")

        # Sort by importance
        df = df.sort_values("Persistence", ascending=False)

        # Save JSON

        out_json = os.path.join(self.output_dir, "contact_persistence.json")
        df.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved contact persistence data to {out_json}")


        # Plot

        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(10, 6))

        # Only show meaningful residues (adjust threshold if needed)
        df_plot = df[df["Persistence"] > 0.1]

        # Fallback: if nothing passes threshold, show top 20
        if df_plot.empty:
            df_plot = df.head(20)

        sns.barplot(
            data=df_plot,
            x="Label",
            y="Persistence",
            hue="Chain"
        )

        plt.xticks(rotation=90)
        plt.ylabel("Contact Persistence")
        plt.xlabel("Residue")
        plt.title("Ligand Binding Hotspots")

        plt.tight_layout()

        out_plot = os.path.join(self.output_dir, "contact_persistence_plot.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved contact persistence plot to {out_plot}")


        return {
            "contacts_data": out_json,
            "contacts_plot": out_plot
        }

@register_task(
    "hydration_site_energy",
    category="Post-proc; graph analyses",
    description="Identify hydration sites and rank them by approximate free energy."
)
class HydrationSiteEnergyTask:
    def __init__(self,
                 topology: str,
                 trajectory: str,
                 water_resname: str = "HOH",
                 ligand_resname: str = "UNK",
                 cutoff: float = 5.0,
                 grid_spacing: float = 0.5,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 1,
                 output_dir: str = "output_hydration_energy",
                 context: Optional[Dict[str, Any]] = None):
        self.topology = topology
        self.trajectory = trajectory
        self.water_resname = water_resname
        self.ligand_resname = ligand_resname
        self.cutoff = cutoff
        self.grid_spacing = grid_spacing
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.context = context or {}
        os.makedirs(self.output_dir, exist_ok=True)
        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        def wrap_like_vmd(ts):
            box = ts.dimensions[:3]

            # center on protein
            com = protein.center_of_mass()
            shift = box / 2.0 - com
            self.u.atoms.positions += shift

            # wrap residues consistently
            for res in self.u.residues:
                res_com = res.atoms.center_of_mass()
                delta = (res_com // box) * box * -1
                res.atoms.positions += delta

            return ts

        logger.info("Detecting hydration sites via water oxygen clustering...")

        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        protein = self.u.select_atoms("protein")
        solvent = self.u.select_atoms(f"resname {self.water_resname}")

        self.u.trajectory.add_transformations(
            unwrap(self.u.atoms),
            wrap_like_vmd
        )
        
        # align once
        logger.info("Aligning trajectory...")
        align.AlignTraj(
            self.u,
            self.u,
            select="protein and backbone",
            in_memory=True
        ).run()

        # Extract reference from aligned universe
        ref_pdb_path = os.path.join(self.output_dir, "reference_structure.pdb")

        self.u.trajectory[0]

        ref = self.u.select_atoms(f"protein or resname {self.ligand_resname}").copy()

        with mda.Writer(ref_pdb_path, ref.n_atoms) as W:
            W.write(ref)

        coords = []

        for ts in self.u.trajectory[self.start:self.stop:self.step]:

            water_oxygens = self.u.select_atoms(
                f"resname {self.water_resname} and name O and around {self.cutoff} group ligand",
                ligand=ligand
            )

            coords.append(water_oxygens.positions.copy())
        
        logger.info(f"Water selection atoms: {len(water_oxygens)}")
        logger.info(f"Ligand atoms: {len(ligand)}")

        coords = np.concatenate(coords, axis=0)

        db = DBSCAN(eps=1.5, min_samples=5).fit(coords)
        labels = db.labels_
        clusters = {}

        for label in set(labels):
            if label == -1:
                continue  # noise

            cluster_points = coords[labels == label]
            centroid = cluster_points.mean(axis=0)

            clusters[label] = {
                "centroid": centroid,
                "count": len(cluster_points)
            }

        df = pd.DataFrame(coords, columns=["x", "y", "z"])
        df["Cluster"] = labels
        site_counts = df["Cluster"].value_counts()
        site_probs = site_counts / site_counts.sum()
        R, T = 8.314, 300.0
        total_points = sum(c["count"] for c in clusters.values())

        data = []

        for label, info in clusters.items():
            occupancy = info["count"] / total_points
            dG = -R * T * np.log(occupancy)

            data.append({
                "Site": label,
                "x": info["centroid"][0],
                "y": info["centroid"][1],
                "z": info["centroid"][2],
                "Occupancy": occupancy,
                "DeltaG_kJmol": dG
            })

        df_energy = pd.DataFrame(data).sort_values("Occupancy", ascending=False)

        pdb_lines = []
        for i, row in df_energy.iterrows():
            pdb_lines.append(
                "HETATM{:5d}  O   HOH A{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}           O\n".format(
                    int(i),
                    int(row["Site"]),
                    row["x"], row["y"], row["z"],
                    row["Occupancy"] * 100
                )
            )

        pdb_path = os.path.join(self.output_dir, "hydration_sites.pdb")
        with open(pdb_path, "w") as f:
            f.writelines(pdb_lines)

        logger.info(f"Saved hydration sites PDB to {pdb_path}")

        out_json = os.path.join(self.output_dir, "hydration_sites.json")
        df_energy.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved hydration site free energies to {out_json}")

        plt.figure(figsize=(8, 5))

        sns.scatterplot(
            data=df_energy,
            x="Occupancy",
            y="DeltaG_kJmol",
            size="Occupancy",
            sizes=(50, 300)
        )

        for _, row in df_energy.head(10).iterrows():
            plt.text(row["Occupancy"], row["DeltaG_kJmol"], f"{int(row['Site'])}")

        plt.xlabel("Occupancy")
        plt.ylabel("ΔG (kJ/mol)")
        plt.title("Hydration Site Stability")

        plt.tight_layout()
        out_plot = os.path.join(self.output_dir, "hydration_sites_energy.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved hydration site energy plot to {out_plot}")

        return {"hydration_data": out_json, "hydration_plot": out_plot}

@register_task(
    "temporal_motif_persistence",
    category="Post-proc; graph analyses",
    description="Quantify persistence of small recurring motifs (e.g., water-mediated interactions) in the solvent network."
)
class TemporalMotifPersistenceTask:
    def __init__(self,
                 topology: str,
                 trajectory: str,
                 ligand_resname: str = "UNK",
                 water_resname: str = "HOH",
                 cutoff: float = 3.5,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 10,
                 output_dir: str = "output_motif_persistence",
                 top_n_motifs: int = 100,
                 smoothing_window: int = 20,
                 context: Optional[Dict[str, Any]] = None):
        self.topology = topology
        self.trajectory = trajectory
        self.ligand_resname = ligand_resname
        self.water_resname = water_resname
        self.cutoff = cutoff
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.top_n_motifs = top_n_motifs
        self.smoothing_window = smoothing_window
        self.context = context or {}
        os.makedirs(self.output_dir, exist_ok=True)
        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.info("Analysing protein-ligand temporal motif persistence...")
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        waters = self.u.select_atoms(f"resname {self.water_resname} and name O")
        protein = self.u.select_atoms("protein")

        motif_history = {}

        for ts in tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Frames"):
            ligand_atoms = ligand.atoms.positions
            water_atoms = waters.positions
            prot_atoms = protein.positions
            
            lw_dist = distance_array(ligand_atoms, water_atoms)
            pw_dist = distance_array(prot_atoms, water_atoms)

            ligand_near_wat = np.where(lw_dist < self.cutoff)
            prot_near_wat = np.where(pw_dist < self.cutoff)

            for l_idx, w_idx in zip(*ligand_near_wat):
                for p_idx, w2_idx in zip(*prot_near_wat):
                    if w_idx == w2_idx:
                        motif = (int(l_idx), int(w_idx), int(p_idx))
                        motif_history.setdefault(motif, []).append(ts.frame)

        lifetimes = []
        filtered_motifs = {}
        for motif, frames in motif_history.items():
            if len(frames) > 1:
                lifetime = frames[-1] - frames[0]
                lifetimes.append(lifetime)
                filtered_motifs[motif] = frames

        if not lifetimes:
            logger.warning("No motifs with persistence > 1 frame found.")
            return {}

        # Group motifs by core identity: (LigandRes, WaterIndex, ProteinResNum, ProteinResName)
        grouped_motifs = {}

        for motif in filtered_motifs:
            l_idx, w_idx, p_idx = motif
            ligand_res = ligand.atoms[l_idx].resname
            water_index = w_idx
            prot_atom = protein.atoms[p_idx]
            prot_resnum = prot_atom.resid
            prot_resname = prot_atom.resname
            
            key = (ligand_res, water_index, prot_resnum, prot_resname)
            frames = filtered_motifs[motif]

            if key not in grouped_motifs:
                grouped_motifs[key] = set()
            grouped_motifs[key].update(frames)

        for key in grouped_motifs:
            grouped_motifs[key] = sorted(grouped_motifs[key])

        # Limit to top_n_motifs by lifetime (max span of frames)
        top_grouped_keys = sorted(
            grouped_motifs.keys(),
            key=lambda k: grouped_motifs[k][-1] - grouped_motifs[k][0],
            reverse=True
        )[:self.top_n_motifs]

        top_motifs_data = []
        for key in top_grouped_keys:
            ligand_res, water_index, prot_resnum, prot_resname = key
            frames = grouped_motifs[key]
            start_frame = frames[0]
            end_frame = frames[-1]
            lifetime = end_frame - start_frame + 1

            top_motifs_data.append({
                "LigandRes": ligand_res,
                "WaterIndex": water_index,
                "ProteinResNum": prot_resnum,
                "ProteinResName": prot_resname,
                "StartFrame": start_frame,
                "EndFrame": end_frame,
                "LifetimeFrames": lifetime,
            })

        top_motifs_df = pd.DataFrame(top_motifs_data)

        # Save summary CSV and JSON
        out_csv = os.path.join(self.output_dir, "motif_persistence.csv")
        out_json = os.path.join(self.output_dir, "motif_persistence.json")
        top_motifs_df.to_csv(out_csv, index=False)
        top_motifs_df.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved grouped motif persistence data to {out_csv} and {out_json}")

        # Prepare heatmap presence matrix
        # Collect all unique frames across grouped motifs
        all_frames_set = set()
        for key in top_grouped_keys:
            all_frames_set.update(grouped_motifs[key])
        all_frames_sorted = sorted(all_frames_set)
        frame_to_col = {frame: idx for idx, frame in enumerate(all_frames_sorted)}

        presence_matrix = np.zeros((len(top_grouped_keys), len(all_frames_sorted)), dtype=np.float32)

        for i, key in enumerate(top_grouped_keys):
            frames_present = grouped_motifs[key]
            cols = [frame_to_col[f] for f in frames_present if f in frame_to_col]
            presence_matrix[i, cols] = 1.0

        # Smooth presence matrix
        window = self.smoothing_window
        smoothed_matrix = np.zeros_like(presence_matrix)
        for i in range(presence_matrix.shape[0]):
            row = presence_matrix[i]
            window = min(self.smoothing_window, len(row))
            
            if window < 2:
                smoothed_matrix[i] = row
                continue

            kernel = np.ones(window) / window
            smoothed = np.convolve(row, kernel, mode='same')

            smoothed_matrix[i] = smoothed[:len(row)]

        motif_labels = [f"L:{k[0]} W:{k[1]} P:{k[3]}{k[2]}" for k in top_grouped_keys]

        plt.figure(figsize=(14, 12))
        sns.heatmap(
            smoothed_matrix,
            cmap="Blues",
            cbar_kws={"label": "Fraction of active motif"},
            yticklabels=motif_labels,
            vmin=0, vmax=1
        )

        frame_time_ns = 2.0 * 100 * self.step / 1e6
        num_frames = smoothed_matrix.shape[1]
        tick_positions = list(range(0, num_frames, max(1, num_frames // 10)))
        tick_labels = [f"{pos * frame_time_ns:.2f}" for pos in tick_positions]
        plt.xticks(ticks=tick_positions, labels=tick_labels, rotation=45)

        plt.ylabel(f"Top {self.top_n_motifs} Unique Motifs by Lifetime")
        plt.xlabel("Simulation Time (ns)")
        plt.title("Top Motif Persistence Heatmap (Smoothed & Grouped)")
        plt.tight_layout()
        out_plot = os.path.join(self.output_dir, "motif_persistence_heatmap.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved grouped motif persistence heatmap to {out_plot}")

        plt.figure(figsize=(6, 4))
        sns.histplot(pd.Series(lifetimes), bins=20, color="#58D68D", kde=True)
        plt.xlabel("Lifetime (frames)")
        plt.ylabel("Count")
        plt.title("Temporal Motif Persistence Distribution")
        plt.tight_layout()
        out_hist = os.path.join(self.output_dir, "motif_persistence_hist.png")
        plt.savefig(out_hist, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved motif persistence histogram to {out_hist}")

        return {
            "motif_csv": out_csv,
            "motif_json": out_json,
            "motif_heatmap": out_plot,
            "motif_histogram": out_hist,
        }

@register_task(
    "network_embedding_analysis",
    category="Post-proc; graph analyses",
    description="Convert trajectory frames into residue–ligand contact graphs, then perform Node2Vec + t-SNE embedding to visualise network evolution."
)
class NetworkEmbeddingAnalysisTask:
    def __init__(self,
                 topology: str,
                 trajectory: str,
                 ligand_resname: str = "UNK",
                 distance_cutoff: float = 4.0,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 10,
                 output_dir: str = "output_network_embedding",
                 node2vec_dim: int = 64, # reasonable default - add to YAML?
                 node2vec_walk_length: int = 10,
                 node2vec_num_walks: int = 50,
                 node2vec_p: float = 1.0,
                 node2vec_q: float = 1.0,
                 perplexity: float = 30.0,
                 random_state: int = 42,
                 context: Optional[Dict[str, Any]] = None):

        os.makedirs(output_dir, exist_ok=True)
        self.topology = topology
        self.trajectory = trajectory
        self.ligand_resname = ligand_resname
        self.distance_cutoff = distance_cutoff
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.node2vec_dim = node2vec_dim
        self.node2vec_walk_length = node2vec_walk_length
        self.node2vec_num_walks = node2vec_num_walks
        self.node2vec_p = node2vec_p
        self.node2vec_q = node2vec_q
        self.perplexity = perplexity
        self.random_state = random_state
        self.context = context or {}
        self.u = mda.Universe(topology, trajectory)

    def _frame_to_graph(self, ts):
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        protein = self.u.select_atoms("protein")  # only protein

        # Debug: selections
        if len(ligand) == 0:
            logger.debug(f"Warning: No ligand atoms found for resname {self.ligand_resname} at frame {ts.frame}")
        else:
            logger.debug(f"Found {len(ligand)} ligand atoms at frame {ts.frame}")

        if len(protein) == 0:
            logger.debug(f"Warning: No protein atoms found at frame {ts.frame}")
        
        logger.debug(f"Ligand selection (frame {ts.frame}): {len(ligand)} atoms")
        logger.debug(f"Protein selection (frame {ts.frame}): {len(protein)} atoms")

        lig_center = ligand.center_of_geometry()
        logger.debug(f"Ligand center of geometry (frame {ts.frame}): {lig_center}")
        logger.debug(f"Frame {ts.frame}: Ligand center of geometry: {lig_center}, Protein selection size: {len(protein)}")

        if len(ligand) == 0 or len(protein) == 0:
            return nx.Graph()  # Return empty graph

        G = nx.Graph(frame=ts.frame)

        min_dist = np.min(np.linalg.norm(ligand.positions[:, None, :] - protein.positions[None, :, :], axis=-1))
        logger.debug(f"Frame {ts.frame}: minimum ligand-protein atom distance = {min_dist:.2f} Å")

        G.add_node(self.ligand_resname, type="ligand")

        # Only protein residues within distance_cutoff
        nearby = self.u.select_atoms(f"protein and around {self.distance_cutoff} group ligand", ligand=ligand)
        logger.debug(f"Found {len(nearby.residues)} protein residues within {self.distance_cutoff} Å of ligand at frame {ts.frame}")

        if len(nearby.residues) == 0:
            logger.info(f"Warning: No residues found within {self.distance_cutoff} Å of the ligand at frame {ts.frame}")

        for res in nearby.residues:
            resid_label = f"{res.resname}{res.resid}"
            G.add_node(resid_label, type="residue")
            G.add_edge(resid_label, self.ligand_resname, weight=1.0)

        # Residue–residue edges
        res_list = [f"{res.resname}{res.resid}" for res in nearby.residues]
        for i, r1 in enumerate(res_list):
            for r2 in res_list[i + 1:]:
                G.add_edge(r1, r2, weight=0.5)

        logger.debug(f"Graph at frame {ts.frame} has {len(G.nodes)} nodes and {len(G.edges)} edges")
        return G

    def _generate_graphs(self):
        graphs = []
        for ts in tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Building contact graphs"):
            G = self._frame_to_graph(ts)
            graphs.append(G)
        return graphs

    def _compute_node2vec_embeddings(self, graphs):

        embeddings_file = os.path.join(self.output_dir, "node2vec_embeddings.h5")

        with h5py.File(embeddings_file, 'w') as f:
            emb_dataset = f.create_dataset("embeddings", 
                                        shape=(len(graphs), self.node2vec_dim), 
                                        dtype=np.float32)

            for i, G in tqdm(enumerate(graphs), total=len(graphs), desc="Node2Vec embeddings"):
                if G.number_of_nodes() < 2:
                    emb_dataset[i] = np.zeros(self.node2vec_dim, dtype=np.float32)
                    logger.info(f"Skipping frame {i} because the graph has fewer than 2 nodes")
                    continue

                logger.debug(f"Processing frame {i} with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

                n2v = Node2Vec(G, dimensions=self.node2vec_dim, walk_length=self.node2vec_walk_length, 
                            num_walks=self.node2vec_num_walks, p=self.node2vec_p, q=self.node2vec_q, 
                            workers=1, seed=self.random_state, quiet=True)

                # Debugging: Check graph info before fitting
                logger.debug(f"Graph {i}: number of nodes = {len(G.nodes)}, number of edges = {len(G.edges)}")

                model = n2v.fit(window=5, min_count=1, batch_words=4)
                
                # Compute the mean of node embeddings
                emb = np.mean([model.wv[node] for node in G.nodes if node in model.wv], axis=0)

                if np.isnan(emb).any():
                    logger.debug(f"Warning: NaN values encountered in embedding for frame {i}")

                emb_dataset[i] = emb

                # Debug: Check embedding
                logger.debug(f"Embedding saved for frame {i}: {emb[:5]}...")

            logger.info(f"Embeddings saved to {embeddings_file}")
        return embeddings_file

    def _load_node2vec_embeddings_in_chunks(self, embeddings_file, chunk_size=100):
        with h5py.File(embeddings_file, 'r') as f:
            total_frames = len(f['embeddings'])
            embeddings = []

            for start_idx in range(0, total_frames, chunk_size):
                end_idx = min(start_idx + chunk_size, total_frames)
                chunk = f['embeddings'][start_idx:end_idx]
                chunk = np.nan_to_num(chunk, nan=0.0, posinf=0.0, neginf=0.0)  # Clean NaN/Inf
                embeddings.append(chunk)
        
        return np.concatenate(embeddings, axis=0)

    def _incremental_pca(self, embeddings, chunk_size=100):
        ipca = IncrementalPCA(n_components=50)
        for start_idx in range(0, len(embeddings), chunk_size):
            end_idx = min(start_idx + chunk_size, len(embeddings))
            ipca.partial_fit(embeddings[start_idx:end_idx])
        
        reduced_embeddings = ipca.transform(embeddings)
        return reduced_embeddings

    def _incremental_tsne(self, embeddings, chunk_size=100):
        tsne = TSNE(n_components=2, perplexity=self.perplexity, random_state=self.random_state, init="pca", learning_rate="auto")
        n_frames = len(embeddings)
        tsne_embeddings = []

        for start_idx in range(0, n_frames, chunk_size):
            end_idx = min(start_idx + chunk_size, n_frames)
            chunk = embeddings[start_idx:end_idx]
            tsne_embeddings.append(tsne.fit_transform(chunk))
        
        return np.concatenate(tsne_embeddings, axis=0)

    def _cluster_embeddings(self, emb_2d):
        # TODO: DBSCAN can be tuned; eps may need adjustment depending on t-SNE scale
        db = DBSCAN(eps=2.0, min_samples=5)
        cluster_labels = db.fit_predict(emb_2d)
        return cluster_labels

    def _compute_cluster_contact_frequencies(self, graphs, cluster_labels):
        """
        Count how often each protein residue appears in each cluster.
        Only residues (type="residue") are counted.
        """

        cluster_res_counts = {}
        for i, G in enumerate(graphs):
            cluster = cluster_labels[i]
            if cluster not in cluster_res_counts:
                cluster_res_counts[cluster] = Counter()

            if self.ligand_resname not in G:
                logger.debug(f"Frame {i}: ligand not in graph, skipping")
                continue

            for neighbor in G.neighbors(self.ligand_resname):
                if G.nodes[neighbor].get("type") == "residue":
                    cluster_res_counts[cluster][neighbor] += 1

        all_residues = sorted({res for c in cluster_res_counts for res in cluster_res_counts[c]})
        df = pd.DataFrame(0, index=sorted(cluster_res_counts.keys()), columns=all_residues)

        for cluster, counter in cluster_res_counts.items():
            for res, count in counter.items():
                df.loc[cluster, res] = count
                logger.debug(f"Cluster {cluster}: Residue {res} count = {count}")

        return df

    def run(self):
        logger.info("Starting network embedding analysis...")
        graphs = self._generate_graphs()

        embeddings_file = self._compute_node2vec_embeddings(graphs)
        frame_embeddings = self._load_node2vec_embeddings_in_chunks(embeddings_file, chunk_size=100)

        # Apply Incremental PCA
        reduced_embeddings = self._incremental_pca(frame_embeddings, chunk_size=100)

        # Apply Incremental t-SNE
        emb_2d = self._incremental_tsne(reduced_embeddings, chunk_size=100)

        df_emb = pd.DataFrame({
            "Frame": np.arange(len(emb_2d)),
            "Dim1": emb_2d[:, 0],
            "Dim2": emb_2d[:, 1],
        })

        # Plot
        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(7, 6))
        sc = plt.scatter(df_emb["Dim1"], df_emb["Dim2"],
                         c=df_emb["Frame"], cmap="viridis",
                         s=60, alpha=0.9, edgecolors="k")
        plt.colorbar(sc, label="Frame index / time progression")
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.title("Network Embedding Evolution (contact graph)")
        plt.tight_layout()

        out_plot = os.path.join(self.output_dir, "network_embedding_tsne.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved t-SNE embedding plot to {out_plot}")

        #  Clustering 
        cluster_labels = self._cluster_embeddings(emb_2d)
        df_emb["Cluster"] = cluster_labels

        out_json = os.path.join(self.output_dir, "network_embedding_tsne.json")
        df_emb.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved t-SNE embedding + cluster data to {out_json}")

        #  Plot t-SNE colored by cluster 
        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(7, 6))
        sc = plt.scatter(df_emb["Dim1"], df_emb["Dim2"],
                         c=df_emb["Cluster"], cmap="tab10",
                         s=60, alpha=0.9, edgecolors="k")
        plt.colorbar(sc, label="Cluster ID")
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.title("Network Embedding Clusters (contact graph)")
        plt.tight_layout()
        out_plot = os.path.join(self.output_dir, "network_embedding_tsne_clusters.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved cluster-colored t-SNE plot to {out_plot}")

        # Compute cluster-specific residue contact frequencies
        cluster_contact_freq = self._compute_cluster_contact_frequencies(graphs, cluster_labels)
        df_heatmap = pd.DataFrame(cluster_contact_freq).fillna(0).T

        # Helper - sort by residue number for plotting.

        def resid_sort_key(res_label):
            """
            Sort residue labels by numeric residue number.
            Example: ARG15 -> 15, GLY2 -> 2
            """
            match = re.search(r'(\d+)$', res_label)
            if match:
                return int(match.group(1))
            else:
                return float('inf')  # in case label has no number, push to end

        df_heatmap = df_heatmap.loc[sorted(df_heatmap.index, key=resid_sort_key)]

        plt.figure(figsize=(10, 8))
        sns.heatmap(df_heatmap,
                    annot=False,        # show counts
                    fmt="d",           # format as integers
                    cmap="YlGnBu",
                    cbar_kws={"label": "Contact frequency"})  # colorbar label

        plt.ylabel("Cluster")
        plt.xlabel("Residue Pair")
        plt.title("Residue-Ligand Contact Frequency per Cluster")
        plt.tight_layout()

        heatmap_plot = os.path.join(self.output_dir, "cluster_contact_frequencies.png")
        plt.savefig(heatmap_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved cluster contact frequency heatmap to {heatmap_plot}")
        
        # Save numeric data to disk
        out_csv = os.path.join(self.output_dir, "cluster_contact_frequencies.csv")
        df_heatmap.to_csv(out_csv)
        logger.info(f"Saved cluster contact frequencies to {out_csv}")

        return {
            "embedding_data": out_json,
            "embedding_plot": out_plot,
            "cluster_contact_heatmap": heatmap_plot,
            "cluster_contact_csv": out_csv
        }

@register_task(
    "protein_protein_network_embedding",
    category="Post-proc; graph analyses",
    description="Convert trajectory frames into residue–residue contact graphs, then perform Node2Vec + t-SNE embedding to visualise protein-protein network evolution."
)
class ProteinProteinNetworkEmbeddingTask:
    def __init__(self,
                 topology: str,
                 trajectory: str,
                 distance_cutoff: float = 4.0,
                 start: int = 0,
                 stop: int = -1,
                 step: int = 10,
                 output_dir: str = "output_protein_protein_network",
                 node2vec_dim: int = 64,
                 node2vec_walk_length: int = 10,
                 node2vec_num_walks: int = 50,
                 node2vec_p: float = 1.0,
                 node2vec_q: float = 1.0,
                 perplexity: float = 30.0,
                 random_state: int = 42,
                 context: Optional[Dict[str, Any]] = None):

        os.makedirs(output_dir, exist_ok=True)
        self.topology = topology
        self.trajectory = trajectory
        self.distance_cutoff = distance_cutoff
        self.start = start
        self.stop = stop
        self.step = step
        self.output_dir = output_dir
        self.node2vec_dim = node2vec_dim
        self.node2vec_walk_length = node2vec_walk_length
        self.node2vec_num_walks = node2vec_num_walks
        self.node2vec_p = node2vec_p
        self.node2vec_q = node2vec_q
        self.perplexity = perplexity
        self.random_state = random_state
        self.context = context or {}
        self.u = mda.Universe(topology, trajectory)

    def _compute_node2vec_embeddings(self, graphs):

        embeddings_file = os.path.join(self.output_dir, "node2vec_embeddings.h5")

        with h5py.File(embeddings_file, 'w') as f:
            emb_dataset = f.create_dataset("embeddings", 
                                        shape=(len(graphs), self.node2vec_dim), 
                                        dtype=np.float32)

            for i, G in tqdm(enumerate(graphs), total=len(graphs), desc="Node2Vec embeddings"):
                if G.number_of_nodes() < 2:
                    emb_dataset[i] = np.zeros(self.node2vec_dim, dtype=np.float32)
                    logger.info(f"Skipping frame {i} because the graph has fewer than 2 nodes")
                    continue

                logger.debug(f"Processing frame {i} with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

                n2v = Node2Vec(G, dimensions=self.node2vec_dim, walk_length=self.node2vec_walk_length, 
                            num_walks=self.node2vec_num_walks, p=self.node2vec_p, q=self.node2vec_q, 
                            workers=1, seed=self.random_state, quiet=True)

                # Debug: Check graph info before fitting
                logger.debug(f"Graph {i}: number of nodes = {len(G.nodes)}, number of edges = {len(G.edges)}")

                model = n2v.fit(window=5, min_count=1, batch_words=4)
                
                # Compute the mean of node embeddings
                emb = np.mean([model.wv[node] for node in G.nodes if node in model.wv], axis=0)

                if np.isnan(emb).any():
                    logger.debug(f"Warning: NaN values encountered in embedding for frame {i}")

                emb_dataset[i] = emb

                # Debug: Check embedding
                logger.debug(f"Embedding saved for frame {i}: {emb[:5]}...")

            logger.info(f"Embeddings saved to {embeddings_file}")
        return embeddings_file
    
    def _load_node2vec_embeddings_in_chunks(self, embeddings_file, chunk_size=100):
        with h5py.File(embeddings_file, 'r') as f:
            total_frames = len(f['embeddings'])
            embeddings = []

            for start_idx in range(0, total_frames, chunk_size):
                end_idx = min(start_idx + chunk_size, total_frames)
                chunk = f['embeddings'][start_idx:end_idx]
                chunk = np.nan_to_num(chunk, nan=0.0, posinf=0.0, neginf=0.0) 
                embeddings.append(chunk)
        
        return np.concatenate(embeddings, axis=0)

    def _incremental_pca(self, embeddings, chunk_size=100):
        ipca = IncrementalPCA(n_components=50)
        for start_idx in range(0, len(embeddings), chunk_size):
            end_idx = min(start_idx + chunk_size, len(embeddings))
            ipca.partial_fit(embeddings[start_idx:end_idx])
        
        reduced_embeddings = ipca.transform(embeddings)
        return reduced_embeddings

    def _tsne(self, embeddings):
        tsne = TSNE(
            n_components=2,
            perplexity=min(self.perplexity, len(embeddings)//3),
            random_state=self.random_state,
            init="pca",
            learning_rate="auto"
        )
        return tsne.fit_transform(embeddings)

    def _cluster_embeddings(self, emb_2d):
        from sklearn.preprocessing import StandardScaler

        emb_scaled = StandardScaler().fit_transform(emb_2d)

        db = DBSCAN(eps=0.5, min_samples=5)
        return db.fit_predict(emb_scaled)


    def _frame_to_graph_pp(self, ts):
        protein = self.u.select_atoms("protein")

        if len(protein) == 0:
            logger.debug(f"Warning: No protein atoms found at frame {ts.frame}")
            return nx.Graph()

        G = nx.Graph(frame=ts.frame)

        residues = protein.residues

        # Add nodes
        for res in residues:
            resid_label = f"{res.resname}{res.resid}"
            G.add_node(resid_label, type="residue")

        # Add edges
        for i, res_i in enumerate(residues):
            label_i = f"{res_i.resname}{res_i.resid}"

            for j, res_j in enumerate(residues[i + 1:], start=i + 1):
                label_j = f"{res_j.resname}{res_j.resid}"

                # Remove trivial local contacts
                seq_sep = abs(res_i.resid - res_j.resid)
                if seq_sep < 5:
                    continue

                # CA–CA distance
                ca_i = res_i.atoms.select_atoms("name CA")
                ca_j = res_j.atoms.select_atoms("name CA")

                if len(ca_i) == 0 or len(ca_j) == 0:
                    continue

                dist = np.linalg.norm(ca_i.positions[0] - ca_j.positions[0])

                if dist <= 8.0:
                    G.add_edge(label_i, label_j, weight=1.0)

        return G

    def _generate_graphs_pp(self):
        graphs = []
        for ts in tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Building protein-protein graphs"):     
            G = self._frame_to_graph_pp(ts)
            graphs.append(G)
        return graphs
    
    def _graphs_to_adjacency_embeddings(self, graphs):
        """
        Convert each graph into a flattened adjacency matrix vector.
        Ensures consistent node ordering across frames.
        """
        # Use node list from first graph (assumes consistent residue set)
        node_list = sorted(graphs[0].nodes())

        embeddings = []
        for G in graphs:
            # Ensure all nodes exist
            for node in node_list:
                if node not in G:
                    G.add_node(node)

            A = nx.to_numpy_array(G, nodelist=node_list)
            vec = A[np.triu_indices(len(node_list), k=1)]  # upper triangle
            embeddings.append(vec)

        return np.array(embeddings, dtype=np.float32)

    def run(self):
        logger.info("Starting protein-protein network embedding analysis...")
        graphs = self._generate_graphs_pp()

        # Compute top contacts across trajectory
        from collections import Counter

        contact_counts = Counter()
        for G in graphs:
            for u, v in G.edges:
                pair = tuple(sorted((u, v)))
                contact_counts[pair] += 1

        df_contacts = pd.DataFrame.from_dict(contact_counts, orient="index", columns=["Count"])
        df_contacts["Frequency"] = df_contacts["Count"] / len(graphs)
        df_contacts = df_contacts.sort_values("Frequency", ascending=False).head(30)

        out_csv_contacts = os.path.join(self.output_dir, "top_contacts.csv")
        df_contacts.to_csv(out_csv_contacts)
        logger.info(f"Saved top contacts to {out_csv_contacts}")

        embeddings = self._graphs_to_adjacency_embeddings(graphs)

        # Reduce dimensionality before t-SNE
        pca = PCA(n_components=min(50, embeddings.shape[1]))
        reduced_embeddings = pca.fit_transform(embeddings)

        # t-SNE (global)
        emb_2d = self._tsne(reduced_embeddings)

        df_emb = pd.DataFrame({
            "Frame": np.arange(len(emb_2d)),
            "Dim1": emb_2d[:, 0],
            "Dim2": emb_2d[:, 1],
        })

        # Plot t-SNE
        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(7, 6))
        sc = plt.scatter(df_emb["Dim1"], df_emb["Dim2"],
                         c=df_emb["Frame"], cmap="viridis",
                         s=60, alpha=0.9, edgecolors="k")
        plt.colorbar(sc, label="Frame index / time progression")
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.title("Protein-Protein Network Embedding Evolution")
        plt.tight_layout()
        out_plot = os.path.join(self.output_dir, "protein_protein_network_tsne.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved t-SNE embedding plot to {out_plot}")

        # Clustering
        cluster_labels = self._cluster_embeddings(emb_2d)
        df_emb["Cluster"] = cluster_labels
        out_json = os.path.join(self.output_dir, "protein_protein_network_tsne.json")
        df_emb.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved t-SNE embedding + cluster data to {out_json}")

        # Cluster-coloured t-SNE
        sns.set(style="whitegrid", context="talk")
        plt.figure(figsize=(7, 6))
        sc = plt.scatter(df_emb["Dim1"], df_emb["Dim2"],
                         c=df_emb["Cluster"], cmap="tab10",
                         s=60, alpha=0.9, edgecolors="k")
        plt.colorbar(sc, label="Cluster ID")
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.title("Protein-Protein Network Embedding Clusters")
        plt.tight_layout()
        out_plot = os.path.join(self.output_dir, "protein_protein_network_tsne_clusters.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved cluster-colored t-SNE plot to {out_plot}")

        # Compute cluster-specific residue contact frequencies
        cluster_res_counts = self._compute_cluster_contact_frequencies_pp(graphs, cluster_labels)
        df_heatmap = pd.DataFrame(cluster_res_counts).fillna(0)

        # Keep only top N most frequent contacts
        top_n = 50
        total_counts = df_heatmap.sum(axis=0)
        top_pairs = total_counts.sort_values(ascending=False).head(top_n).index
        df_heatmap = df_heatmap[top_pairs]

        plt.figure(figsize=(10, 6))
        sns.heatmap(df_heatmap, annot=False, cmap="YlGnBu")
        plt.xlabel("Cluster")
        plt.ylabel("Residue")
        plt.title("Residue-Residue Contact Frequency per Cluster")
        plt.tight_layout()
        heatmap_plot = os.path.join(self.output_dir, "protein_protein_cluster_contact_frequencies.png")
        plt.savefig(heatmap_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved protein-protein cluster contact frequency heatmap to {heatmap_plot}")

        out_csv = os.path.join(self.output_dir, "protein_protein_cluster_contact_frequencies.csv")
        df_heatmap.to_csv(out_csv)
        logger.info(f"Saved protein-protein cluster contact frequencies to {out_csv}")

        return {
            "embedding_data": out_json,
            "embedding_plot": out_plot,
            "cluster_contact_heatmap": heatmap_plot,
            "cluster_contact_csv": out_csv
        }

    def _compute_cluster_contact_frequencies_pp(self, graphs, cluster_labels):
        """
        Count how often each residue appears in contacts within its cluster.
        """
        from collections import Counter
        cluster_res_counts = {}
        for i, G in enumerate(graphs):
            cluster = cluster_labels[i]
            if cluster not in cluster_res_counts:
                cluster_res_counts[cluster] = Counter()
            for u, v in G.edges:
                pair = tuple(sorted((u, v)))
                cluster_res_counts[cluster][pair] += 1

        all_pairs = sorted({
            pair for c in cluster_res_counts for pair in cluster_res_counts[c]
        })
        df = pd.DataFrame(0, index=cluster_res_counts.keys(), columns=all_pairs)

        for cluster, counter in cluster_res_counts.items():
            for res, count in counter.items():
                df.loc[cluster, res] = count
        
        # Normalise by number of frames
        n_frames = len(graphs)
        df = df / n_frames

        return df
