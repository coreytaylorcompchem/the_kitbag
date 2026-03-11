import re
import os

from functools import reduce
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from sklearn.manifold import TSNE
from PIL import Image, ImageDraw
import networkx as nx

from modules.utils.convert_adme import (
    convert_permeability, 
    convert_met_stab,
    convert_herg, 
    convert_cyp_activity, 
    convert_logp_logd, 
    convert_ppb,
    convert_solubility

)

from modules.utils.detect_adme import (
    detect_metstab_system,
    detect_papp_direction,
    build_dual_endpoint_dataset
)

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

# TODO: re-factor this to store in a dict or something.

permeability_assays = ["caco", "mdck", "pampa", "p-gp", "bcrp", "mrp"]
solubility_assays = ["solubility", "logS", "logs"]
metstab_assays = ["microsomal", "microsomes", "hepato", "hepatocyte", "hepatic"]
ppb_assays = ["plasma protein binding", "ppb", "protein binding"]
logp_assays = ["logp", "log p", "partition coefficient"]
logd_assays = ["logd", "distribution coefficient"]
herg_assays = ["herg", "kcnq1", "ikrv", "ether-a-go-go"]
cyp_assays = ["cyp1a2", "cyp2c9", "cyp2c19", "cyp2d6", "cyp3a4", "cyp3a5", "cyp2b6", "cyp2e1"]

@register_task(
    "harmonise_units",
    category="ADME",
    description="Harmonise units for permeability, solubility, microsomal/hepatocyte stability, and PPB with before/after proportion plots"
)
def harmonise_units(config, enriched=None):
    """
    Converts ADME units to canonical form and plots before/after unit proportions:
      - Permeability → 10^-6 cm/s, adds AB/BA direction
      - Solubility → log10(mol/L)
      - Microsomal/Hepatocyte stability → canonical clearance units
      - PPB → fraction unbound (0-1)
    Extracts species from assay_description if target_organism is missing.
    """
    output_dir = config.get("output", {}).get("directory", "outputs/adme")
    plots_dir = os.path.join(output_dir, "_plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    cleaned = {}

    # Helper; extract species from description
    def extract_species(text):
        if text is None:
            return "unknown"

        t = str(text).lower()

        species_map = {
            "human": "Human",
            "homo sapiens": "Human",

            "mouse": "Mouse",
            "mus musculus": "Mouse",

            "rat": "Rat",
            "rattus": "Rat",

            "monkey": "Monkey",
            "macaca": "Monkey",
            "cyno": "Monkey",
        }

        for k, v in species_map.items():
            if k in t:
                return v

        return "unknown"

    for assay_name, records in enriched.items():
        records = [r for r in records if isinstance(r, dict)]
        if not records:
            continue

        units_before = [r.get("standard_units") for r in records if r.get("standard_units")]
        new_records = []
        lname = assay_name.lower()

        if any(x in lname for x in permeability_assays):
            for r in records:
                r["papp_direction"] = detect_papp_direction(r.get("assay_description"))
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, new_unit = convert_permeability(val, unit)
                if new_val is None:
                    continue
                r["standard_value"] = new_val
                r["standard_units"] = new_unit
                new_records.append(r)

        elif any(x in lname for x in solubility_assays):
            for r in records:
                val = r.get("standard_value")
                unit = r.get("standard_units")
                smiles = r.get("smiles")
                new_val = convert_solubility(val, unit, smiles)

                if new_val is None:
                    continue
                r["standard_value"] = new_val
                r["standard_units"] = "log10(mol/L)"

                new_records.append(r)
        elif any(x in lname for x in logp_assays + logd_assays):
            for r in records:
                val = r.get("standard_value")
                unit = r.get("standard_units")
                new_val = convert_logp_logd(val, unit)

                if new_val is None:
                    continue
                r["standard_value"] = new_val
                r["standard_units"] = "log_unitless"

                new_records.append(r)

        elif any(x in lname for x in metstab_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, new_unit = convert_met_stab(val, unit)

                if new_val is None:
                    continue
                r["standard_value"] = new_val
                r["standard_units"] = new_unit
                # Use target_organism if present, else extract from description
                species_text = r.get("target_organism") or r.get("assay_description")
                r["species"] = extract_species(species_text)

                system = detect_metstab_system(r.get("assay_description"))

                # keep only in vitro
                if system != "in_vitro":
                    continue

                r["metstab_system"] = system

                new_records.append(r)
        # PPB
        elif any(x in lname for x in ppb_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val = convert_ppb(val, unit)
                if new_val is None:
                    continue
                r["standard_value"] = new_val
                r["standard_units"] = "fraction_unbound"
                # Use target_organism if present, else extract from description
                species_text = r.get("target_organism") or r.get("assay_description")
                r["species"] = extract_species(species_text)
                new_records.append(r)
        # hERG
        elif any(x in lname for x in herg_assays):

            for r in records:

                val, unit = r.get("standard_value"), r.get("standard_units")

                new_val, endpoint = convert_herg(val, unit)

                if new_val is None:
                    continue

                r["standard_value"] = new_val

                if endpoint == "IC50_nM":
                    r["standard_units"] = "nM"
                    r["herg_endpoint"] = "IC50"

                elif endpoint == "inhibition_pct":
                    r["standard_units"] = "%"
                    r["herg_endpoint"] = "inhibition"

                new_records.append(r)
        #CYP
        elif any(x in lname for x in cyp_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, endpoint, new_unit = convert_cyp_activity(val, unit)
                if new_val is None:
                    continue

                if endpoint == "IC50":
                    r["cyp_endpoint"] = "IC50"
                    r["IC50"] = new_val
                    r["IC50_unit"] = new_unit
                elif endpoint == "inhibition":
                    r["cyp_endpoint"] = "inhibition"
                    r["inhibition"] = new_val
                    r["inhibition_unit"] = new_unit
                
                r["standard_value"] = new_val
                r["standard_units"] = new_unit

                new_records.append(r)
        else:
            new_records = records

        units_after = [r.get("standard_units") for r in new_records if r.get("standard_units")]

        # Plot before vs after proportions
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
            plt.xticks(x, all_units, rotation=45, fontsize=8)
            plt.ylabel("Proportion of measurements")
            plt.title(f"{assay_name} unit harmonisation (before/after)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{assay_name}_units_before_after.png"), dpi=300)
            plt.close()

        cleaned[assay_name] = new_records

    return cleaned

@register_task(
    "build_multitask_dataset",
    category="ADME",
    description="Construct multi-task learning dataset with AB/BA/unknown splitting for permeability assays, including replicate statistics"
)
def build_multitask_dataset(config, cleaned=None):

    output_cfg = config.get("output", {})
    out_dir = output_cfg.get("directory", "outputs/adme")
    filename = output_cfg.get("filename", "chembl_mtl_dataset.csv")
    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, filename)

    dfs = []

    allowed_species = {"Human", "Mouse", "Rat", "Monkey", "unknown"}

    for assay_name, records in cleaned.items():

        records = [r for r in records if isinstance(r, dict)]
        if not records:
            continue

        df = pd.DataFrame(records)

        df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
        df = df.dropna(subset=["standard_value"])
        if df.empty:
            continue

        assay_lower = assay_name.lower()

        # Normalise solubility endpoint
        if "solubility" in assay_lower or "logs" in assay_lower:
            assay_name = "logS"

        # -----------------------------
        # Permeability assays
        # -----------------------------
        if any(a in assay_lower for a in permeability_assays):

            if "papp_direction" not in df.columns:
                df["papp_direction"] = df["assay_description"].apply(
                    lambda x: detect_papp_direction(x)
                )

            df["papp_direction"] = df["papp_direction"].fillna("unknown")

            df_grouped = (
                df.groupby(["smiles", "papp_direction"])["standard_value"]
                .agg(["mean", "std", "count"])
                .reset_index()
            )

            df_pivot = df_grouped.pivot(index="smiles", columns="papp_direction")

            df_pivot.columns = [
                f"{assay_name}_{direction}_{stat}" for stat, direction in df_pivot.columns
            ]

            ab_col = f"{assay_name}_AB_mean"
            ba_col = f"{assay_name}_BA_mean"

            if ab_col in df_pivot.columns and ba_col in df_pivot.columns:

                df_pivot[f"{assay_name}_efflux_ratio"] = (
                    df_pivot[ba_col] / df_pivot[ab_col].replace(0, np.nan)
                )

                df_pivot[f"{assay_name}_log_efflux_ratio"] = np.log10(
                    df_pivot[f"{assay_name}_efflux_ratio"]
                )

            df_pivot = df_pivot.reset_index()

            dfs.append(df_pivot)

        # -----------------------------
        # Metabolic stability
        # -----------------------------
        elif any(a in assay_lower for a in metstab_assays):

            df["species"] = df.get("species", "unknown")
            df["species"] = df["species"].fillna("unknown")

            df = df[df["species"].isin(allowed_species)]
            if df.empty:
                continue

            df_grouped = (
                df.groupby(["smiles", "species"])["standard_value"]
                .agg(["mean", "std", "count"])
                .reset_index()
            )

            df_pivot = df_grouped.pivot(index="smiles", columns="species")

            df_pivot.columns = [
                f"{assay_name}_{species}_{stat}" for stat, species in df_pivot.columns
            ]

            df_pivot = df_pivot.reset_index()

            dfs.append(df_pivot)

        # -----------------------------
        # hERG special handling TODO: Refactor this similar to CYP code
        # -----------------------------
        elif any(a in assay_lower for a in herg_assays):

            # Recover incorrectly labelled units
            def recover_herg_units(row):

                val = row["standard_value"]
                unit = str(row.get("standard_units", "")).lower()

                # percent inhibition
                if "%" in unit:
                    return val, "inhibition"

                # nM
                if "nm" in unit:
                    return val, "IC50"

                # µM or mislabeled values
                if "um" in unit or unit in ["", "none", "no unit"]:
                    return val * 1000, "IC50"

                # M
                if "m" == unit:
                    return val * 1e9, "IC50"

                return None, None

            df[["herg_value", "herg_type"]] = df.apply(
                lambda r: pd.Series(recover_herg_units(r)), axis=1
            )

            df = df.dropna(subset=["herg_value"])

            df_grouped = (
                df.groupby(["smiles", "herg_type"])["herg_value"]
                .agg(["mean", "std", "count"])
                .reset_index()
            )

            df_pivot = df_grouped.pivot(index="smiles", columns="herg_type")

            df_pivot.columns = [
                f"hERG_{etype}_{stat}" for stat, etype in df_pivot.columns
            ]

            df_pivot = df_pivot.reset_index()

            dfs.append(df_pivot)
        
        # hERG assays
        elif any(a in assay_lower for a in herg_assays):

            df_pivot = build_dual_endpoint_dataset(
                df,
                assay_name,
                endpoint_col="herg_endpoint",
                value_cols={
                    "IC50": "IC50",
                    "inhibition": "inhibition"
                }
            )

            dfs.append(df_pivot)

        # -----------------------------
        # Other assays
        # -----------------------------
        elif any(a in assay_lower for a in cyp_assays):

            df_pivot = build_dual_endpoint_dataset(
                df,
                assay_name,
                endpoint_col="cyp_endpoint",
                value_cols={
                    "IC50": "IC50",
                    "inhibition": "inhibition"
                }
            )

            dfs.append(df_pivot)
        else:

            if "species" in df.columns:

                df["species"] = df["species"].fillna("unknown")

                df = df[df["species"].isin(allowed_species)]
                if df.empty:
                    continue

                df_grouped = (
                    df.groupby(["smiles", "species"])["standard_value"]
                    .agg(["mean", "std", "count"])
                    .reset_index()
                )

                df_pivot = df_grouped.pivot(index="smiles", columns="species")

                df_pivot.columns = [
                    f"{assay_name}_{species}_{stat}" for stat, species in df_pivot.columns
                ]

                df_pivot = df_pivot.reset_index()

                dfs.append(df_pivot)

            else:

                df_subset = df.groupby("smiles")["standard_value"].mean().reset_index()

                df_subset = df_subset.rename(columns={"standard_value": assay_name})

                dfs.append(df_subset)

    if not dfs:
        return pd.DataFrame()

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