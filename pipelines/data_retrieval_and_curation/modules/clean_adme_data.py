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
from rdkit.Chem import Descriptors

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from sklearn.manifold import TSNE
from PIL import Image, ImageDraw
import networkx as nx

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

from pipeline.task_registry import register_task

def detect_papp_direction(description):
        if not description:
            return None
        desc = description.lower().replace("→", "->").replace("–", "-").replace("—", "-")
        ab_patterns = [r"a\s*->\s*b", r"apical\s*(?:to|->|-)\s*basolateral", r"\bab\b"]
        ba_patterns = [r"b\s*->\s*a", r"basolateral\s*(?:to|->|-)\s*apical", r"\bba\b"]
        for p in ab_patterns:
            if re.search(p, desc):
                return "AB"
        for p in ba_patterns:
            if re.search(p, desc):
                return "BA"
        return None

def convert_permeability(value, unit):
    """
    Convert various permeability units to canonical 10^-6 cm/s.
    Returns (float_value, canonical_unit) or (None, None) if unconvertible.
    """
    if value is None or unit is None:
        return None, None

    u = unit.lower().replace(" ", "").replace("'", "").replace("*", "").replace("^", "-").strip()

    # Canonical unit
    canonical = "10^-6 cm/s"

    # Map known variants directly
    direct_map = {
        "10-6cm/s": canonical,
        "10^-6cm/s": canonical,
        "10e-6cm/s": canonical,
        "10-6cms": canonical,
        "10^-6cms": canonical,
        "10e-6cms": canonical,
        "ucm/s": canonical,
        "papp10e6cm/s": canonical,
        "10^6cm/s": canonical,
        "10-7cm/s": canonical,   # convert to 0.1 * 10^-6 cm/s
        "10^-7cm/s": canonical,  # same here
        "10-8cm/s": canonical,   # 0.01 * 10^-6 cm/s
        "10^-8cm/s": canonical,
    }

    if u in direct_map:
        factor = 1.0
        if u in ["10-7cm/s", "10^-7cm/s"]:
            factor = 0.1
        elif u in ["10-8cm/s", "10^-8cm/s"]:
            factor = 0.01
        return float(value) * factor, canonical

    # Convert cm/s, nm/s, um/s
    if "cms" in u:
        return float(value) * 1e6, canonical
    if "nm/s" in u:
        return float(value) * 1e-1, canonical  # 1 nm/s = 0.1 * 10^-6 cm/s
    if "um/s" in u or "μm/s" in u:
        return float(value) * 1e-2, canonical
    if "cmhr-1" in u:
        return float(value) * 2.778e-6, canonical  # cm/hr → cm/s × 10^-6

    # Reject obviously incompatible units
    reject_patterns = ["%", "pmol", "nmol", "min", "sec", "umol", "ug", "a.u.", "uL/hr/cm2"]
    if any(r in u for r in reject_patterns):
        return None, None

    # Fallback: unhandled
    return None, None

def convert_solubility(value, unit, smiles):
    if value is None or unit is None:
        return None

    try:
        value = float(value)
    except:
        return None

    unit = unit.lower().strip()

    # nM -> mol/L
    if unit == "nm":
        mol_l = value * 1e-9

    # ug/mL -> mol/L (requires MW)
    elif unit in ["ug/ml", "ug.mL-1", "µg/ml"]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mw = Descriptors.MolWt(mol)
        mol_l = (value / mw) * 1e-3

    else:
        return None

    # Convert to logS
    if mol_l <= 0:
        return None

    return np.log10(mol_l)

def convert_met_stab(value, unit):
    """
    Convert microsomal/hepatocyte clearance units to canonical mL/min/g liver or uL/min/10^6 cells.
    Returns (float_value, canonical_unit) or (None, None) if unconvertible.
    """
    if value is None or unit is None:
        return None, None

    u = str(unit).lower().replace(" ", "").replace("*","").replace("^","-").replace("(","").replace(")","").strip()

    # Canonical units
    if "cell" in u:
        canonical = "uL/min/10^6 cells"
    else:
        canonical = "mL/min/g liver"

    try:
        val = float(value)
    except:
        return None, None

    # Mapping known variants
    direct_map = {
        "ul/min": 1.0 if canonical=="uL/min/10^6 cells" else None,
        "ul.min-1.(10^6cells)-1":1.0,
        "uL/min/1e6cells":1.0,
        "ml.min-1.kg-1": 1.0,   # assume mL/min/g conversion will be handled later
        "ml.min-1.g-1": 1.0,
        "hr": None,  # needs half-life conversion, skip for now
        "min": 1.0,
    }

    for k,v in direct_map.items():
        if k in u and v is not None:
            return val * v, canonical

    # Reject obviously incompatible units
    reject_patterns = ["%", "nmol", "nM", "ug", "mM", "pmol"]
    if any(r in u for r in reject_patterns):
        return None, None

    return val, canonical

def convert_ppb(value, unit):
    """
    Convert PPB to fraction unbound (0-1)
    """
    if value is None or unit is None:
        return None

    try:
        val = float(value)
    except:
        return None

    u = str(unit).lower().replace(" ", "").strip()

    if u in ["%", "percent", "%bound"]:
        return val / 100.0
    elif u in ["fu","fractionunbound","fraction_unbound"]:
        return val
    else:
        # units like ug/mL, pmol, etc are ignored
        return None

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
    
    permeability_assays = ["caco", "mdck", "pampa", "p-gp", "bcrp", "mrp"]
    solubility_assays = ["solubility", "logS", "logs"]
    metstab_assays = ["microsomal", "microsomes", "hepato", "hepatocyte", "hepatic"]
    ppb_assays = ["plasma protein binding", "ppb", "protein binding"]

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
                new_records.append(r)

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
    """
    Build a multi-task dataset from cleaned/enriched activity records.
    Captures mean, std, and replicate count for permeability and metabolic stability assays.
    Species extracted from assay_description are used for pivoting.
    """
    output_cfg = config.get("output", {})
    out_dir = output_cfg.get("directory", "outputs/adme")
    filename = output_cfg.get("filename", "chembl_mtl_dataset.csv")
    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, filename)

    dfs = []

    permeability_assays = ["caco", "mdck", "pampa", "p-gp", "bcrp", "mrp"]
    metstab_assays = ["microsomal", "microsomes", "hepato", "hepatocyte", "hepatic"]

    # Species we care about
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

        # Permeability assays 
        if any(a in assay_lower for a in permeability_assays):
            if "papp_direction" not in df.columns:
                df["papp_direction"] = df["assay_description"].apply(lambda x: detect_papp_direction(x))
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
            df_pivot = df_pivot.reset_index()
            dfs.append(df_pivot)

        # Metabolic stability assays 
        elif any(a in assay_lower for a in metstab_assays):
            df["species"] = df.get("species", "unknown")
            df["species"] = df["species"].fillna("unknown")

            # Filter unwanted species
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

        # Other assays (PPB, logP/D, etc) 
        else:
            if "species" in df.columns:
                df["species"] = df["species"].fillna("unknown")

                # Filter unwanted species
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
    logger.info(f"Saved multi-task dataset CSV to: {output_path}")

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

    # Clustered heatmap for better visualization
    g = sns.clustermap(
        overlap_plot.astype(float),
        cmap="viridis",
        annot=True,
        fmt=".2f",
        figsize=(12,12)
    )

    g.fig.suptitle("Pairwise Assay Overlap Fraction (clustered)", y=1.02)
    plt.savefig(os.path.join(plots_dir, "pairwise_overlap_heatmap.png"), dpi=300)
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

    # TSNE of fingerprints
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

    # Scaffold analysis
    scaffolds = [MurckoScaffold.GetScaffoldForMol(m) for m in mols if m is not None]
    scaffold_smiles = [Chem.MolToSmiles(s) for s in scaffolds if s is not None]
    scaffold_counts = Counter(scaffold_smiles)

    # Histogram of scaffold frequency
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