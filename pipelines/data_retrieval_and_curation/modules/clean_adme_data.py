import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from functools import reduce
import seaborn as sns
from collections import Counter

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from sklearn.manifold import TSNE
from PIL import Image, ImageDraw
import networkx as nx

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from modules.utils.convert_adme import (
    convert_permeability,
    convert_met_stab,
    convert_ppb,
    convert_solubility,
    convert_logp_logd,
    convert_cyp_activity,
    convert_herg,
)
from modules.utils.detect_adme import detect_papp_direction, detect_metstab_system

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

# Assay categories
permeability_assays = ["caco", "mdck", "pampa", "p-gp", "bcrp", "mrp"]
solubility_assays = ["solubility", "logs", "logS"]
metstab_assays = ["microsomal", "microsomes", "hepato", "hepatocyte", "hepatic"]
ppb_assays = ["plasma protein binding", "ppb", "protein binding"]
logp_assays = ["logp", "log p", "partition coefficient"]
logd_assays = ["logd", "distribution coefficient"]
# herg_assays = ["herg", "kcnq1", "ikrv", "ether-a-go-go"]
cyp_assays = ["cyp1a2","cyp2c9","cyp2c19","cyp2d6","cyp3a4","cyp3a5","cyp2b6","cyp2e1", "cyp2c8"]

@register_task(
    "harmonise_units",
    category="ADME",
    description="Harmonise units for ADME and type-B tox endpoints with before/after plots"
)
def harmonise_units(config, enriched=None):
    output_dir = config.get("output", {}).get("directory", "outputs/adme")
    plots_dir = os.path.join(output_dir, "_plots")
    os.makedirs(plots_dir, exist_ok=True)

    cleaned = {}

    # Tox assay names from config
    config_targets = config.get("targets", {})
    tox_assays = {k.strip().lower() for k in config_targets.keys()}

    # Simple species extractor
    def extract_species(text):
        if not text:
            return "unknown"
        t = str(text).lower()
        species_map = {
            "human": "Human", "homo sapiens": "Human",
            "mouse": "Mouse", "mus musculus": "Mouse",
            "rat": "Rat", "rattus": "Rat",
            "monkey": "Monkey", "macaca": "Monkey", "cyno": "Monkey"
        }
        for k, v in species_map.items():
            if k in t:
                return v
        return "unknown"

    for assay_name, records in enriched.items():
        logger.debug(f"Processing assay: {assay_name}")
        records = [r for r in records if isinstance(r, dict)]
        if not records:
            continue

        lname = assay_name.lower()

        units_before = [r.get("standard_units") for r in records if r.get("standard_units")]
        new_records = []

        # -----------------------------
        # ADME assays
        # -----------------------------
        if any(x in lname for x in permeability_assays):
            for r in records:
                r["papp_direction"] = detect_papp_direction(r.get("assay_description"))
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, new_unit = convert_permeability(val, unit)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, new_unit
                    new_records.append(r)

        elif any(x in lname for x in solubility_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                smiles = r.get("smiles")
                new_val = convert_solubility(val, unit, smiles)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, "log10(mol/L)"
                    new_records.append(r)

        elif any(x in lname for x in logp_assays + logd_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val = convert_logp_logd(val, unit)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, "log_unitless"
                    new_records.append(r)

        elif any(x in lname for x in metstab_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, new_unit = convert_met_stab(val, unit)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, new_unit
                    r["species"] = extract_species(r.get("target_organism") or r.get("assay_description"))
                    system = detect_metstab_system(r.get("assay_description"))
                    if system == "in_vitro":
                        r["metstab_system"] = system
                        new_records.append(r)

        elif any(x in lname for x in ppb_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val = convert_ppb(val, unit)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, "fraction_unbound"
                    r["species"] = extract_species(r.get("target_organism") or r.get("assay_description"))
                    new_records.append(r)

        elif any(x in lname for x in cyp_assays):  # CYPs are ADME (type A)
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, endpoint, new_unit = convert_cyp_activity(val, unit)
                if new_val is None:
                    continue
                r["standard_value"], r["standard_units"] = new_val, new_unit
                r["cyp_endpoint"] = endpoint  # e.g., IC50 or inhibition
                if endpoint == "IC50":
                    r["IC50"] = new_val
                    r["IC50_unit"] = new_unit
                elif endpoint == "inhibition":
                    r["inhibition"] = new_val
                    r["inhibition_unit"] = new_unit
                # Also keep for generic pivoting
                r["tox_endpoint"] = endpoint
                new_records.append(r)
        # -----------------------------
        # Type-B tox endpoints (hERG, Nav1.5, etc.) from config
        # -----------------------------
        elif any(lname == t for t in tox_assays):

            for r in records:
                val = r.get("standard_value")
                unit = r.get("standard_units")
                stype = str(r.get("standard_type", "")).strip().upper()

                if val is None:
                    continue

                try:
                    val = float(val)
                except:
                    continue

                endpoint_type = None

                # concentration endpoints
                if stype in {"IC50", "EC50", "KI", "KD"}:
                    endpoint_type = stype

                # percent inhibition
                elif "%" in str(unit).lower() or "INHIBITION" in stype:
                    endpoint_type = "inhibition"

                if endpoint_type is None:
                    continue

                # unit normalisation
                if endpoint_type in {"IC50", "EC50", "KI", "KD"}:
                    u = str(unit).lower()

                    factor = 1
                    if "um" in u or "μm" in u:
                        factor = 1e3
                    elif "mm" in u:
                        factor = 1e6

                    val = val * factor
                    unit = "nM"

                r["standard_value"] = val
                r["standard_units"] = unit
                r["tox_type"] = endpoint_type
                r["endpoint"] = assay_name

                # debug log for verification
                logger.debug({
                    "smiles": r.get("smiles"),
                    "assay_name": assay_name,
                    "tox_type": r["tox_type"],
                    "standard_value": r["standard_value"],
                    "standard_units": r["standard_units"]
                })
                
                new_records.append(r)

        # -----------------------------
        # Fallback
        # -----------------------------
        else:
            for r in records:
                val = r.get("standard_value")
                if val is None:
                    continue
                r["endpoint"] = assay_name
                new_records.append(r)

        # -----------------------------
        # Plot before/after proportions
        # -----------------------------

        # logger.debug(f"{assay_name} endpoint types: {Counter(r['tox_type'] for r in new_records)}")

        units_after = [r.get("standard_units") for r in new_records if r.get("standard_units")]
        if units_before or units_after:
            all_units = sorted(set(units_before) | set(units_after))
            before_counts = Counter(units_before)
            after_counts = Counter(units_after)
            n_before = sum(before_counts.values()) or 1
            n_after = sum(after_counts.values()) or 1
            proportions_before = [before_counts.get(u,0)/n_before for u in all_units]
            proportions_after = [after_counts.get(u,0)/n_after for u in all_units]

            x = range(len(all_units))
            width = 0.35

            plt.figure(figsize=(8,4))
            plt.bar([i - width/2 for i in x], proportions_before, width=width, color='skyblue', label='Before')
            plt.bar([i + width/2 for i in x], proportions_after, width=width, color='lightgreen', label='After')
            plt.xticks(x, all_units, rotation=45, fontsize=6)
            plt.ylabel("Proportion of measurements")
            plt.title(f"{assay_name} unit harmonisation (before/after)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{assay_name}_units_before_after.png"), dpi=300)
            plt.close()
        
        # debug checks of assay data

        if lname in tox_assays:
            logger.debug(f"{assay_name} example records after harmonisation:")
            for r in new_records[:5]:
                logger.debug({
                    "standard_type": r.get("standard_type"),
                    "endpoint": r.get("endpoint"),
                    "standard_units": r.get("standard_units"),
                    "standard_value": r.get("standard_value")
                })

        cleaned[assay_name] = new_records

    return cleaned

@register_task(
    "build_multitask_dataset",
    category="ADME",
    description="Construct multi-task dataset including ADME (type-A) and type-B tox endpoints"
)
def build_multitask_dataset(config, cleaned=None):
    output_cfg = config.get("output", {})
    out_dir = output_cfg.get("directory", "outputs/adme")
    filename = output_cfg.get("filename", "chembl_mtl_dataset.csv")
    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, filename)

    dfs = []
    allowed_species = {"Human", "Mouse", "Rat", "Monkey", "unknown"}

    # Helper to recover tox values (IC50, inhibition, etc.)
    def recover_tox_value(row, endpoint_col):
        val = row["standard_value"]
        unit = str(row.get("standard_units", "")).lower() if row.get("standard_units") else ""
        ep = row[endpoint_col]
        if unit in ["", "none", "no unit"] and ep in {"IC50", "EC50", "Ki", "Kd"}:
            val = val * 1000
            unit = "nM"
        elif "um" in unit:
            val = val * 1000
            unit = "nM"
        elif unit == "m":
            val = val * 1e9
            unit = "nM"
        return val, ep, unit

    for assay_name, records in cleaned.items():
        records = [r for r in records if isinstance(r, dict)]
        if not records:
            continue

        df = pd.DataFrame(records)

        # debug checks of assay data

        logger.debug(f"{assay_name} DataFrame columns after harmonisation: {df.columns.tolist()}")

        logger.debug(f"Processing assay: {assay_name}")
        logger.debug(f"df columns: {df.columns.tolist()}")
        if len(records) > 0:
            logger.debug(f"Sample record: {records[0]}")

        if assay_name.lower() in {"herg", "nav1.5"}:
            logger.debug(f"{assay_name} columns in DataFrame: {df.columns.tolist()}")
            logger.debug(df[["standard_type","endpoint","standard_units"]].head())
    
        df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
        df = df.dropna(subset=["standard_value"])
        if df.empty:
            continue

        # -----------------------------
        # Type-A assays (ADME, CYPs)
        # -----------------------------
        group_cols = []
        if "species" in df.columns:
            group_cols.append("species")
        if "papp_direction" in df.columns:
            group_cols.append("papp_direction")
        if "cyp_endpoint" in df.columns:
            group_cols.append("cyp_endpoint")

        if group_cols:
            df[group_cols] = df[group_cols].fillna("unknown")
            df_grouped = df.groupby(["smiles"] + group_cols)["standard_value"].agg(["mean", "std", "count"]).reset_index()
            df_pivot = df_grouped.pivot(index="smiles", columns=group_cols)
            # Properly join MultiIndex columns without splitting letters
            df_pivot.columns = [
                f"{assay_name}_{'_'.join(map(str, col))}_{stat}" if isinstance(col, tuple) else f"{assay_name}_{col}_{stat}"
                for stat, col in df_pivot.columns
            ]
            dfs.append(df_pivot.reset_index())
            continue

        # -----------------------------
        # Type-B assays (tox endpoints: hERG, Nav, etc.)
        # -----------------------------
        if "tox_type" in df.columns:

            logger.debug(f"Building pivot for tox assay: {assay_name}")
            logger.debug(
                f"DataFrame head before grouping:\n"
                f"{df[['smiles','tox_type','standard_value','endpoint']].head(10)}"
            )

            df_grouped = (
                df.groupby(["smiles", "tox_type"])["standard_value"]
                .agg(["mean", "std", "count"])
                .reset_index()
            )

            logger.debug(f"Grouped DataFrame:\n{df_grouped.head(10)}")

            df_pivot = df_grouped.pivot(index="smiles", columns="tox_type")

            # flatten MultiIndex columns
            df_pivot.columns = [
                f"{assay_name}_{tox}_{stat}"
                for stat, tox in df_pivot.columns
            ]

            dfs.append(df_pivot.reset_index())
            continue

        # -----------------------------
        # Fallback for remaining assays
        # -----------------------------
        df_subset = df.groupby("smiles")["standard_value"].mean().reset_index()
        df_subset.rename(columns={"standard_value": assay_name}, inplace=True)
        dfs.append(df_subset)

    if not dfs:
        return pd.DataFrame()

    # Merge all assay tables
    
    mtl_df = reduce(lambda l, r: pd.merge(l, r, on="smiles", how="outer"), dfs)
    mtl_df.to_csv(output_path, index=False)
    logger.debug(f"Saved multi-task dataset CSV to: {output_path}")

    return mtl_df

@register_task(
    "generate_diagnostics",
    category="ADME",
    description="Generate diagnostic plots and tables for ADME dataset"
)
def generate_diagnostics(config, mtl_df=None):
    """
    Generates:
    1. Unit distributions per assay
    2. Pairwise assay overlap table and heatmap
    3. Network plot to show task connectivity
    4. t-SNE of Morgan fingerprints
    5. Scaffold frequency histogram and top scaffold grid
    
    Saves outputs in <output_dir>_plots
    """

    if mtl_df is None or mtl_df.empty:
        logger.warning("mtl_df is None or empty - skipping diagnostics")
        return {"plots_dir": None}
    
    output_dir = config.get("output", {}).get("directory", "outputs/adme")
    plots_dir = os.path.join(output_dir, "_plots")
    os.makedirs(plots_dir, exist_ok=True)

    # # Pairwise overlap
    # # Only include original assay columns, exclude computed statistics
    # stat_suffixes = ["_std", "_min", "_max", "_count"]
    # task_cols = [
    #     c for c in mtl_df.columns
    #     if c != 'smiles' and not any(c.endswith(suffix) for suffix in stat_suffixes)
    # ]
    # mask = mtl_df[task_cols].notna()

    # overlap_counts = pd.DataFrame(index=task_cols, columns=task_cols, dtype=int)
    # overlap_frac = pd.DataFrame(index=task_cols, columns=task_cols, dtype=float)

    # for t1 in task_cols:
    #     for t2 in task_cols:
    #         both = (mask[t1] & mask[t2]).sum()
    #         overlap_counts.loc[t1, t2] = both
    #         min_count = min(mask[t1].sum(), mask[t2].sum())
    #         overlap_frac.loc[t1, t2] = both / min_count if min_count > 0 else np.nan

    # overlap_counts.to_csv(os.path.join(plots_dir, "pairwise_overlap_counts.csv"))
    # overlap_frac.to_csv(os.path.join(plots_dir, "pairwise_overlap_fraction.csv"))

    # # Clean labels for plotting (remove "_mean")
    # plot_labels = [c.replace("_mean", "") for c in task_cols]

    # overlap_plot = overlap_frac.copy()
    # overlap_plot.index = plot_labels
    # overlap_plot.columns = plot_labels

    # plt.figure(figsize=(10,8))
    # sns.heatmap(overlap_plot.astype(float), annot=True, fmt=".2f", cmap="viridis")
    # plt.title("Pairwise Assay Overlap Fraction")
    # plt.tight_layout()
    # plt.savefig(os.path.join(plots_dir, "pairwise_overlap_heatmap.png"), dpi=300)
    # plt.close()

    # Pairwise overlap
    # Only include mean columns (ignore other stats)
    task_cols = [c for c in mtl_df.columns if c.endswith("_mean")]
    mask = mtl_df[task_cols].notna()

    overlap_counts = pd.DataFrame(index=task_cols, columns=task_cols, dtype=int)
    overlap_frac = pd.DataFrame(index=task_cols, columns=task_cols, dtype=float)

    for t1 in task_cols:
        for t2 in task_cols:
            both = (mask[t1] & mask[t2]).sum()
            overlap_counts.loc[t1, t2] = both
            min_count = min(mask[t1].sum(), mask[t2].sum())
            overlap_frac.loc[t1, t2] = both / min_count if min_count > 0 else np.nan

    # Save raw counts and fraction tables
    overlap_counts.to_csv(os.path.join(plots_dir, "pairwise_overlap_counts.csv"))
    overlap_frac.to_csv(os.path.join(plots_dir, "pairwise_overlap_fraction.csv"))

    # Clean labels for plotting (remove "_mean")
    plot_labels = [c.replace("_mean", "") for c in task_cols]
    overlap_plot = overlap_frac.copy()
    overlap_plot.index = plot_labels
    overlap_plot.columns = plot_labels

    g = sns.clustermap(
        overlap_plot.astype(float),
        cmap="viridis",
        annot=True,
        fmt=".2f",
        figsize=(15,15)
    )

    g.fig.suptitle("Pairwise Assay Overlap Fraction (clustered)", y=1.02)
    plt.savefig(os.path.join(plots_dir, "pairwise_overlap_heatmap.png"), dpi=300)
    plt.close()

    # Task–Task Correlation Heatmap
    # Only use numeric mean columns
    numeric_cols = [c for c in mtl_df.columns if c.endswith("_mean")]
    corr_df = mtl_df[numeric_cols].corr()

    plt.figure(figsize=(12,10))
    sns.heatmap(corr_df, cmap="coolwarm", center=0, annot=False)
    plt.title("Task–Task Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "task_correlation_heatmap.png"), dpi=300)
    plt.close()

    # Correlation vs Overlap Scatter
    # Compute max overlap for each task pair
    max_overlap = overlap_frac.where(~np.eye(len(overlap_frac),dtype=bool)).max(axis=0)

    cor_vs_overlap = pd.DataFrame({
        "task": numeric_cols,
        "max_overlap": [max_overlap.get(c, np.nan) for c in numeric_cols],
        "mean_correlation": [corr_df[c].drop(c).mean() for c in numeric_cols]
    })

    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=cor_vs_overlap,
        x="max_overlap",
        y="mean_correlation"
    )
    plt.xlabel("Maximum Fraction Overlap with Any Task")
    plt.ylabel("Mean Correlation with Other Tasks")
    plt.title("Correlation vs Overlap per Task")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "correlation_vs_overlap.png"), dpi=300)
    plt.close()

    # Task connectivity graph
    G = nx.Graph()

    # Add nodes
    for task in task_cols:
        G.add_node(task)

    # Add edges based on overlap
    for t1 in task_cols:
        for t2 in task_cols:
            if t1 == t2:
                continue
            weight = overlap_frac.loc[t1, t2]
            if pd.notna(weight) and weight > 0.05:  # threshold to reduce clutter
                G.add_edge(t1, t2, weight=weight)

    # Draw graph
    plt.figure(figsize=(10,8))

    pos = nx.spring_layout(G, seed=42)

    edges = G.edges(data=True)
    weights = [d["weight"]*5 for (_,_,d) in edges]  # scale for visibility

    nx.draw_networkx_nodes(G, pos, node_size=1500, node_color="lightblue")
    nx.draw_networkx_labels(G, pos, font_size=9)
    nx.draw_networkx_edges(G, pos, width=weights, alpha=0.6)

    plt.title("Task Connectivity Graph (shared compound overlap)")
    plt.axis("off")
    plt.tight_layout()

    plt.savefig(os.path.join(plots_dir, "task_connectivity_graph.png"), dpi=300)
    plt.close()

    # TSNE of all fingerprints
    mols = [Chem.MolFromSmiles(sm) for sm in mtl_df.smiles]
    fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius=2, nBits=1024) for m in mols]

    def fp_to_array(fp):
        arr = np.zeros((fp.GetNumBits(),), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    fp_arrays = np.array([fp_to_array(fp) for fp in fps if fp is not None])
    tsne = TSNE(n_components=2, random_state=42, perplexity=5)
    fp_tsne = tsne.fit_transform(fp_arrays)

    plt.figure(figsize=(8,6))
    sns.scatterplot(x=fp_tsne[:,0], y=fp_tsne[:,1])
    plt.title("t-SNE of Molecular Fingerprints")
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "tsne_fingerprints.png"), dpi=300)
    plt.close()

    # t-SNE grid with highlighted endpoints
    # Group columns by endpoint type, include raw and _mean as needed
    def get_group_cols(mtl_df, keywords):
        return [c for c in mtl_df.columns if any(k.lower() in c.lower() for k in keywords)]

    endpoint_groups = {
        "PhysChem": get_group_cols(mtl_df, ["LogP", "LogD", "logS", "ppb"]),
        "Metabolism": get_group_cols(mtl_df, ["microsomal", "hepatocyte", "hepato"]),
        "Permeability": get_group_cols(mtl_df, ["caco", "mdck", "pampa", "p-gp", "bcrp", "mrp"]),
        "Tox": get_group_cols(mtl_df, ["herg", "cyp"])
    }

    fig, axes = plt.subplots(2, 2, figsize=(16,14))
    axes = axes.flatten()

    for idx, (group_name, cols) in enumerate(endpoint_groups.items()):
        ax = axes[idx]
        
        if not cols:
            ax.scatter(fp_tsne[:,0], fp_tsne[:,1], color="lightgrey", s=30, alpha=0.5)
            ax.set_title(f"t-SNE colored by {group_name} (no endpoints found)")
            ax.set_xlabel("t-SNE 1")
            ax.set_ylabel("t-SNE 2")
            continue

        # Determine which compounds have at least one value present in this group
        present_mask = mtl_df[cols].notna().any(axis=1)
        
        # Plot: overlay data points with all data to highlight where the data is on the tsne
        ax.scatter(
            fp_tsne[present_mask,0], 
            fp_tsne[present_mask,1], 
            color="dodgerblue", 
            s=20, 
            alpha=0.35,
            label="Present"
        )

        ax.scatter(
            fp_tsne[~present_mask,0], 
            fp_tsne[~present_mask,1], 
            color="lightgrey", 
            s=15, 
            alpha=0.1, 
            label="Absent"
        )

        ax.set_title(f"t-SNE highlighting {group_name}")
        ax.set_xlabel("t-SNE 1")
        ax.set_ylabel("t-SNE 2")
        ax.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "tsne_grid_highlighted.png"), dpi=300)
    plt.close()

    # Histogram of scaffold frequency

    scaffolds = [MurckoScaffold.GetScaffoldForMol(m) for m in mols if m is not None]
    scaffold_smiles = [Chem.MolToSmiles(s) for s in scaffolds if s is not None]
    scaffold_counts = Counter(scaffold_smiles)

    plt.figure(figsize=(8,5))
    plt.hist(list(scaffold_counts.values()), bins=30, log=True)
    plt.xlabel("Scaffold Frequency")
    plt.ylabel("Number of Scaffolds (log scale)")
    plt.title("Distribution of Scaffold Frequencies")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "scaffold_frequency.png"), dpi=300)
    plt.close()

    # Top 10 scaffolds grid
    
    top_n = 10
    top_scaffolds = scaffold_counts.most_common(top_n)
    scaffold_mols = [Chem.MolFromSmiles(smi) for smi, _ in top_scaffolds]
    scaffold_labels = [f"{count} molecules" for _, count in top_scaffolds]

    mol_size = (200, 200)
    images = [Draw.MolToImage(mol, size=mol_size) for mol in scaffold_mols]

    def add_label(img, label, size=(200, 220)):
        new_img = Image.new("RGB", size, "white")
        new_img.paste(img, (0, 0))
        draw = ImageDraw.Draw(new_img)
        draw.text((10, 200), label, fill="black")
        return new_img

    labeled_images = [add_label(img, label) for img, label in zip(images, scaffold_labels)]
    cols = 5
    rows = (len(labeled_images) + cols - 1) // cols
    grid_width = mol_size[0] * cols
    grid_height = 220 * rows
    grid_img = Image.new("RGB", (grid_width, grid_height), "white")

    for idx, img in enumerate(labeled_images):
        x = (idx % cols) * mol_size[0]
        y = (idx // cols) * 220
        grid_img.paste(img, (x, y))

    grid_img.save(os.path.join(plots_dir, "top_10_scaffolds.png"))

    logger.info(f"Plots saved to: {plots_dir}")

    return {"df": mtl_df, "plots_dir": plots_dir}