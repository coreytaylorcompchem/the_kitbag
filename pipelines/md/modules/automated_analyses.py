import os
import warnings
import json
import itertools
import h5py

import networkx as nx
import tqdm as tqdm
import numpy as np

from typing import List, Optional, Dict, Any
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import pandas as pd

from modules.utils.mda_utils import _bit_to_color_value, _get_inv_color_mapper, _get_color_mapper

import prolif as plf
from prolif.plotting.utils import separated_interaction_colors

# suppress MDA warnings.

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis")

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

from sklearn.decomposition import IncrementalPCA
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
        self.u = mda.Universe(self.topology, self.trajectory)

    def run(self):
        logger.debug("Starting RMSD analysis...")

        protein = self.u.select_atoms("protein and backbone")
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        align.AlignTraj(self.u, self.u, select="protein and backbone",
                        in_memory=True).run()

        if len(protein) == 0:
            raise ValueError("No protein atoms found in topology for RMSD calculation.")

        if len(ligand) == 0:
            logger.warning(f"No ligand atoms found for resname {self.ligand_resname}. Skipping ligand RMSD.")

        # Reference (frame0)
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

    def run(self):
        logger.debug("Starting RMSF analysis...")

        protein = self.u.select_atoms("protein")
        calphas = self.u.select_atoms("protein and name CA")

        if len(protein) == 0:
            raise ValueError("No protein atoms found for RMSF calculation.")

        # Align trajectory to backbone reference
        logger.debug("Aligning trajectory to protein backbone reference...")
        align.AlignTraj(self.u, self.u, select="protein and backbone",
                        in_memory=True).run()

        # ------------------------------------------------------------------
        # Per-residue RMSF
        # ------------------------------------------------------------------
        rmsf_all = rms.RMSF(protein).run()
        df_all_atoms = pd.DataFrame({
            "Residue": [atom.resid for atom in protein.atoms],
            "RMSF (Å)": rmsf_all.rmsf,
        })

        df_all = df_all_atoms.groupby("Residue", as_index=False).mean()
        df_all["Component"] = "Protein (all atoms)"

        # ------------------------------------------------------------------
        # Cα RMSF
        # ------------------------------------------------------------------
        rmsf_ca = rms.RMSF(calphas).run()
        df_ca = pd.DataFrame({
            "Residue": [atom.resid for atom in calphas.atoms],
            "RMSF (Å)": rmsf_ca.rmsf,
            "Component": "Cα only"
        })

        df_rmsf = pd.concat([df_all, df_ca], ignore_index=True)
        logger.debug(f"Computed RMSF for {len(df_all)} residues (all atoms) "
                    f"and {len(df_ca)} alpha carbons.")

        out_json = os.path.join(self.output_dir, "rmsf_data.json")
        df_rmsf.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved RMSF data to {out_json}")

        # Plot
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

    def run(self):
        logger.debug("Starting ProLIF interaction fingerprint analysis...")

        # Select lig
        ligand = self.u.select_atoms(self.ligand_selection)
        if len(ligand) == 0:
            raise ValueError(f"No atoms found for ligand selection: {self.ligand_selection}")

        # Select prot
        protein = self.u.select_atoms(self.protein_selection, ligand=ligand)
        if len(protein) == 0:
            raise ValueError(f"No atoms found for protein selection: {self.protein_selection}")

        # Prolif FP calculation
        logger.debug("Running ProLIF fingerprint calculation...")
        
        fp = plf.Fingerprint()
        fp.run(self.u.trajectory[self.start:self.stop:self.step], ligand, protein)
        fp_df = fp.to_dataframe()

        out_data_path = os.path.join(self.output_dir, "fp_data.pkl")
        fp_df.to_pickle(out_data_path)
        logger.info(f"Saved fingerprint DataFrame to {out_data_path}")

        # Generate barcode plot
        logger.info("Generating interaction barcode plot...")

        # Calculate time per frame in ns
        times_ps = [ts.time for ts in self.u.trajectory[self.start:self.stop:self.step]]
        times_ns = np.array(times_ps) / 1000  # convert ps to ns

        fp_transposed = fp_df.astype(np.uint8).T.apply(_bit_to_color_value, axis=1)

        color_mapper = _get_color_mapper()  
        inv_color_mapper = _get_inv_color_mapper()

        # Convert MultiIndex df to integer array for plotting
        plot_data = fp_transposed.copy()
        for idx in plot_data.index:
            interaction_type = idx[2]  # i.e. third level of MultiIndex
            plot_data.loc[idx] = plot_data.loc[idx] * color_mapper.get(interaction_type, 0)

        plot_data_array = plot_data.values.astype(int)

        cmap = ListedColormap([separated_interaction_colors.get(name, "white") for name in [None]+list(color_mapper.keys())])

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

        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"{times_ns[i]:.1f}" for i in tick_indices])
        ax.set_xlabel("Simulation Time (ns)")

        # Format Y-axis
        residues = [f"{prot}" for lig, prot, inter in fp_transposed.index]
        ax.set_yticks(np.arange(len(residues)))
        ax.set_yticklabels(residues)

        # Legend
        unique_values = np.unique(plot_data_array)
        unique_values = [v for v in unique_values if v != 0]  # remove 0 (white)
        legend_colors = {inv_color_mapper[v]: im.cmap(v) for v in unique_values}
        patches = [Patch(color=color, label=interaction) for interaction, color in legend_colors.items()]
        ax.legend(handles=patches, loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=False, ncol=3)
        
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave 5% at bottom for legend

        out_plot_path = os.path.join(self.output_dir, "interaction_barcode.png")

        ligands = sorted({lig[:3] for lig, prot, inter in fp_transposed.index})
        ligand_str = ", ".join(ligands)

        ax.set_title(f"Ligand interaction fingerprint for ligand {ligand_str}.")
        plt.tight_layout()
        plt.savefig(out_plot_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        logger.info(f"Saved interaction barcode plot to {out_plot_path}")

        return {"fp_data": out_data_path, "fp_plot": out_plot_path}

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
        protein = self.u.select_atoms(f"protein and around {self.binding_site_cutoff} resname {self.ligand_resname}")

        contact_counts = {}
        for ts in tqdm.tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Contacts"):
            for res in protein.residues:
                res_atoms = res.atoms
                d = mda.lib.distances.distance_array(res_atoms.positions, ligand.positions).min()
                if d < self.binding_site_cutoff:
                    key = res.resid
                    contact_counts[key] = contact_counts.get(key, 0) + 1

        # Build weighted graph
        G = nx.Graph()
        for res in protein.residues:
            G.add_node(res.resid, name=res.resname)
        for i, j in itertools.combinations(contact_counts.keys(), 2):
            G.add_edge(i, j, weight=(contact_counts[i] + contact_counts[j]) / 2)

        # Community detection (Louvain)
        logger.info("Detecting residue communities...")
        
        partition = community_louvain.best_partition(G, weight='weight')
        nx.set_node_attributes(G, partition, "community")

        df_nodes = pd.DataFrame([
            {"Resid": n, "Resname": G.nodes[n]["name"], "Community": G.nodes[n]["community"]}
            for n in G.nodes
        ])
        out_json = os.path.join(self.output_dir, "plin_communities.json")
        df_nodes.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved community assignments to {out_json}")

        # Plot network
        sns.set(style="white", context="talk")
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(G, seed=42, weight='weight')
        communities = [G.nodes[n]["community"] for n in G.nodes]
        cmap = sns.color_palette("husl", len(set(communities)))
        nx.draw(G, pos, node_color=[cmap[c] for c in communities], with_labels=True, font_size=8)
        plt.title("Protein–Ligand Interaction Communities")
        out_plot = os.path.join(self.output_dir, "plin_community_network.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved PLIN community network to {out_plot}")

        return {"communities_data": out_json, "communities_plot": out_plot}

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
        logger.info("Detecting hydration sites via water oxygen clustering...")
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        water_oxygens = self.u.select_atoms(f"resname {self.water_resname} and name O and around {self.cutoff} resname {self.ligand_resname}")
        coords = []

        for ts in tqdm.tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Waters"):
            coords.append(water_oxygens.positions.copy())
        coords = np.concatenate(coords, axis=0)

        db = DBSCAN(eps=1.5, min_samples=5).fit(coords)
        labels = db.labels_

        df = pd.DataFrame(coords, columns=["x", "y", "z"])
        df["Cluster"] = labels
        site_counts = df["Cluster"].value_counts()
        site_probs = site_counts / site_counts.sum()
        R, T = 8.314, 300.0
        df_energy = pd.DataFrame({
            "Site": site_probs.index,
            "Occupancy": site_probs.values,
            "DeltaG_kJmol": -R * T * np.log(site_probs.values)
        })

        out_json = os.path.join(self.output_dir, "hydration_sites.json")
        df_energy.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved hydration site free energies to {out_json}")

        plt.figure(figsize=(6, 4))
        sns.barplot(data=df_energy, x="Site", y="DeltaG_kJmol", color="#5DADE2")
        plt.ylabel("ΔG (kJ/mol)")
        plt.title("Relative Hydration Site Free Energies")
        plt.tight_layout()
        out_plot = os.path.join(self.output_dir, "hydration_sites_energy.png")
        plt.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved hydration site energy plot to {out_plot}")

        return {"hydration_data": out_json, "hydration_plot": out_plot}

@register_task(
    "temporal_motif_persistence",
    category="Post-proc; graph analyses",
    description="Quantify persistence of small recurring motifs (e.g., ligand–water–residue triangles) in the solvent network."
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
        logger.info("Analyzing temporal motif persistence...")
        ligand = self.u.select_atoms(f"resname {self.ligand_resname}")
        waters = self.u.select_atoms(f"resname {self.water_resname} and name O")
        protein = self.u.select_atoms("protein")

        motif_history = {}

        for ts in tqdm.tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Frames"):
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

            # Append all frames for this core motif key
            if key not in grouped_motifs:
                grouped_motifs[key] = set()
            grouped_motifs[key].update(frames)

        # Convert frames sets to sorted lists
        for key in grouped_motifs:
            grouped_motifs[key] = sorted(grouped_motifs[key])

        # Limit to top_n_motifs by lifetime (max span of frames)
        top_grouped_keys = sorted(
            grouped_motifs.keys(),
            key=lambda k: grouped_motifs[k][-1] - grouped_motifs[k][0],
            reverse=True
        )[:self.top_n_motifs]

        # Prepare DataFrame rows with combined start, end, lifetime
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
            smoothed_matrix[i] = np.convolve(presence_matrix[i], np.ones(window) / window, mode='same')

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

        # Plot histogram of lifetimes (from original motifs, optional)
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
                 node2vec_dim: int = 64,
                 node2vec_walk_length: int = 10,
                 node2vec_num_walks: int = 50,
                 node2vec_p: float = 1.0,
                 node2vec_q: float = 1.0,
                 perplexity: float = 30.0,
                 random_state: int = 42,
                 context: Optional[Dict[str, Any]] = None):
        import os
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
        for ts in tqdm.tqdm(self.u.trajectory[self.start:self.stop:self.step], desc="Building contact graphs"):
            G = self._frame_to_graph(ts)
            graphs.append(G)
        return graphs

    def _compute_node2vec_embeddings(self, graphs):

        embeddings_file = os.path.join(self.output_dir, "node2vec_embeddings.h5")

        with h5py.File(embeddings_file, 'w') as f:
            emb_dataset = f.create_dataset("embeddings", 
                                        shape=(len(graphs), self.node2vec_dim), 
                                        dtype=np.float32)

            for i, G in tqdm.tqdm(enumerate(graphs), total=len(graphs), desc="Node2Vec embeddings"):
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
        # DBSCAN can be tuned; eps may need adjustment depending on t-SNE scale
        db = DBSCAN(eps=2.0, min_samples=5)
        cluster_labels = db.fit_predict(emb_2d)
        return cluster_labels

    def _compute_cluster_contact_frequencies(self, graphs, cluster_labels):
        """
        Count how often each protein residue appears in each cluster.
        Only residues (type="residue") are counted.
        """
        from collections import Counter
        import pandas as pd

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

        # Convert to DataFrame
        all_residues = sorted({res for c in cluster_res_counts for res in cluster_res_counts[c]})
        df = pd.DataFrame(0, index=sorted(cluster_res_counts.keys()), columns=all_residues)

        for cluster, counter in cluster_res_counts.items():
            for res, count in counter.items():
                df.loc[cluster, res] = count
                logger.debug(f"Cluster {cluster}: Residue {res} count = {count}")

        return df

    def run(self):
        logger.info("Starting network embedding analysis (from trajectory)...")
        graphs = self._generate_graphs()

        embeddings_file = self._compute_node2vec_embeddings(graphs)
        frame_embeddings = self._load_node2vec_embeddings_in_chunks(embeddings_file, chunk_size=100)

        # Step 1: Apply Incremental PCA
        reduced_embeddings = self._incremental_pca(frame_embeddings, chunk_size=100)

        # Step 2: Apply Incremental t-SNE
        emb_2d = self._incremental_tsne(reduced_embeddings, chunk_size=100)

        df_emb = pd.DataFrame({
            "Frame": np.arange(len(emb_2d)),
            "Dim1": emb_2d[:, 0],
            "Dim2": emb_2d[:, 1],
        })

        # ---- Plot ----
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

        # ----- Step 3: Clustering -----
        cluster_labels = self._cluster_embeddings(emb_2d)
        df_emb["Cluster"] = cluster_labels

        out_json = os.path.join(self.output_dir, "network_embedding_tsne.json")
        df_emb.to_json(out_json, orient="records", indent=2)
        logger.info(f"Saved t-SNE embedding + cluster data to {out_json}")

        # ----- Plot: t-SNE colored by cluster -----
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

        # ----- Step 4–5: Compute cluster-specific residue contact frequencies -----
        cluster_contact_freq = self._compute_cluster_contact_frequencies(graphs, cluster_labels)
        df_heatmap = pd.DataFrame(cluster_contact_freq).fillna(0).T

        import re
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

        # Sort residues by numeric order
        df_heatmap = df_heatmap.loc[sorted(df_heatmap.index, key=resid_sort_key)]

        plt.figure(figsize=(10, 8))
        sns.heatmap(df_heatmap,
                    annot=False,        # show counts
                    fmt="d",           # format as integers
                    cmap="YlGnBu",
                    cbar_kws={"label": "Contact frequency"})  # colorbar label

        plt.xlabel("Cluster")
        plt.ylabel("Residue")
        plt.title("Residue-Ligand Contact Frequency per Cluster")
        plt.tight_layout()

        heatmap_plot = os.path.join(self.output_dir, "cluster_contact_frequencies.png")
        plt.savefig(heatmap_plot, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved cluster contact frequency heatmap to {heatmap_plot}")
        # Also save numeric data
        out_csv = os.path.join(self.output_dir, "cluster_contact_frequencies.csv")
        df_heatmap.to_csv(out_csv)
        logger.info(f"Saved cluster contact frequencies to {out_csv}")

        return {
            "embedding_data": out_json,
            "embedding_plot": out_plot,
            "cluster_contact_heatmap": heatmap_plot,
            "cluster_contact_csv": out_csv
        }

