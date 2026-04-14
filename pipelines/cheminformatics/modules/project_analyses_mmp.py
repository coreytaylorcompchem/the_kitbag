import csv
import logging
import subprocess
import json
from pathlib import Path
from tqdm import tqdm
import hashlib
import io

import numpy as np
from PIL import Image
import seaborn as sns
import networkx as nx
import pandas as pd
from collections import Counter

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, Descriptors, Crippen, rdMolDescriptors, QED, rdRGroupDecomposition, Scaffolds
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from mpl_toolkits.axes_grid1 import ImageGrid

from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap

from scipy.stats import linregress, ttest_ind, wasserstein_distance
from scipy.spatial.distance import mahalanobis

import matplotlib
matplotlib.use('Agg')
import matplotlib.offsetbox as offsetbox
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

from concurrent.futures import ProcessPoolExecutor, as_completed
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def write_mmpdb_inputs(df, output_prefix, smiles_col="smiles", id_col="id", props_cols=None):
    output_prefix = Path(output_prefix)

    if id_col not in df.columns:
        df[id_col] = df.index.astype(str)

    smi_path = output_prefix.with_suffix(".smi")
    with open(smi_path, "w") as f:
        for _, row in df.iterrows():
            smi = row.get(smiles_col)
            mol_id = row.get(id_col)
            if pd.notna(smi) and pd.notna(mol_id):
                f.write(f"{smi}\t{mol_id}\n")
    logger.info(f"SMILES file written: {smi_path}")

    prop_path = None
    if props_cols:
        missing_cols = [col for col in props_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing property columns in DataFrame: {missing_cols}")

        props_df = df[[id_col] + props_cols].copy()
        props_df.rename(columns={id_col: "id"}, inplace=True)
        prop_path = output_prefix.with_suffix(".props")
        props_df.to_csv(prop_path, index=False, sep="\t", encoding="utf-8", lineterminator="\n")
        logger.info(f"Properties file written: {prop_path}")
    else:
        logger.info("[ℹ] No props_cols provided — skipping props file.")

    return smi_path, prop_path

def run_transform(mmpdb_file, smiles, property_name):
    cmd = [
        "mmpdb", "transform",
        str(mmpdb_file),
        "--smiles", smiles,
        "--property", property_name
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_lines = result.stdout.strip().splitlines()

        # Add property name to each line to identify property in merged output
        # mmpdb is painful to use here so we need to call the program on each SMILES and aggregate.

        return property_name, smiles, output_lines
    except subprocess.CalledProcessError as e:
        logger.warning(f"Transform failed for SMILES {smiles} with property {property_name}")
        logger.debug(e.stderr)
        return property_name, smiles, None

@register_task("mmp_analysis", category="Project-based analyses", description="Matched Molecular Pairs (mmpdb).")
def mmp_analysis(config, data=None):
    input_file = config.get("input_file")
    cfg = config.get("mmp_analysis", {})
    activity_col = cfg.get("activity_col", "pActivity")
    print(activity_col)
    output_dir = Path(config.get("output", {}).get("directory", "outputs/mmp"))
    output_dir.mkdir(parents=True, exist_ok=True)
    out_filename = config.get("output", {}).get("filename", Path(input_file).stem + "_mmp.tsv")

    df = pd.read_csv(input_file)
    df = df.dropna(subset=["smiles", activity_col])

    # Recode toxic and reactive flags to numeric 0/1 if present
    if 'toxic_flag' in df.columns:
        df['toxic_flag'] = df['toxic_flag'].map({'N': 0, 'Y': 1})
    if 'reactive_flag' in df.columns:
        df['reactive_flag'] = df['reactive_flag'].map({'N': 0, 'Y': 1})

    # Rename activity_col to "property" as required by mmpdb
    df = df.rename(columns={activity_col: "property"})

    # Define property columns to process
    props_cols = ["property", "mw", "logp", "hbd", "hba", "rotatable_bonds", "tpsa", "qed", "stereocenters"] # TODO: add to yaml
    # props_cols = ["property", "mw"] # Keep for testing
    if 'toxic_flag' in df.columns:
        props_cols.append("toxic_flag")
    if 'reactive_flag' in df.columns:
        props_cols.append("reactive_flag")

    smi_path, prop_csv_path = write_mmpdb_inputs(
        df,
        output_prefix=output_dir / "temp_input",
        smiles_col="smiles",
        id_col="id",
        props_cols=props_cols
    )

    fragdb = output_dir / "temp.fragdb"
    mmpdb_file = output_dir / "temp.mmpdb"
    transform_out = output_dir / out_filename

    # Run mmpdb fragment/index in separate commands (!)
    cmd1 = ["mmpdb", "fragment", str(smi_path), "-o", str(fragdb)]
    cmd2 = ["mmpdb", "index", str(fragdb), "--properties", str(prop_csv_path), "-o", str(mmpdb_file)]

    logger.info(f"[mmp_analysis] Running: {' '.join(cmd1)}")
    subprocess.run(cmd1, check=True)

    logger.info(f"[mmp_analysis] Running: {' '.join(cmd2)}")
    subprocess.run(cmd2, check=True)

    smiles_list = df["smiles"].tolist()
    all_results = []
    header_line = None

    # Mapping of renamed property back to original column name
    prop_name_map = {"property": activity_col, "mw": "mw", "toxic_flag": "toxic_flag", "reactive_flag": "reactive_flag"}

    with ProcessPoolExecutor() as executor:
        futures = []
        for prop in props_cols:
            for smi in smiles_list:
                futures.append(executor.submit(run_transform, mmpdb_file, smi, prop))

        for future in as_completed(futures):
            prop, smi, output_lines = future.result()
            if output_lines:
                if header_line is None:
                    header_line = output_lines[0]

                # Get the real/original property name (hack)
                real_prop_name = prop_name_map.get(prop, prop)

                for line in output_lines[1:]:
                    if line.strip():
                        all_results.append({
                            "property": real_prop_name,
                            "smiles": smi,
                            "data": line
                        })
            else:
                logger.warning(f"[!] No output for SMILES {smi} with property {prop}")

    if header_line is None:
        raise RuntimeError("No MMP transform results generated.")

    columns = ["property", "smiles"] + header_line.strip().split("\t")

    rows = []
    for res in all_results:
        row_fields = res["data"].strip().split("\t")
        row = [res["property"], res["smiles"]] + row_fields
        rows.append(row)

    with open(transform_out, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)

    result_df = pd.read_csv(transform_out, sep="\t")
    logger.info(f"MMP: got {len(result_df)} transform rows")

    return (Path(input_file).stem, result_df)

############### Plotting helper for mmp_report

def mol_to_img(smiles, img_size=(100, 100)):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=img_size)
    return img

@register_task("mmp_report", category="Project-based analyses", description="Generate MMP transform summary plots and tables.")
def mmp_report(config, data=None):

    def mol_to_img(smiles, size=(200, 100)):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Draw.MolToImage(mol, size=size)
        return None

    # 1. Setup 
    input_file = config.get("input_file")
    output_dir = Path(config.get("output", {}).get("directory", "outputs/mmp_report"))
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_file, sep='\t')

    min_count = config.get("min_count", 10)
    significance_level = config.get("significance_level", 0.05)

    df_filtered = df[df['property_count'] >= min_count].copy()
    df_filtered['transform_label'] = df_filtered['property_from_smiles'] + " → " + df_filtered['property_to_smiles']

    #  2. Structure image caching 
    mol_images_from = {smi: mol_to_img(smi) for smi in pd.unique(df_filtered['property_from_smiles'])}
    mol_images_to = {smi: mol_to_img(smi) for smi in pd.unique(df_filtered['property_to_smiles'])}

    #  3. Multi-Property Heatmap from Long Format 
    # 1. Create transform plot labels
    df_filtered['transform_label'] = df_filtered['property_from_smiles'] + " → " + df_filtered['property_to_smiles']

    # 2. Aggregate per transform + all properties as set in the MMP analysis
    agg_df = df_filtered.groupby(['transform_label', 'property']).agg({
        'property_avg': 'mean',
        'property_count': 'sum'
    }).reset_index()

    # 3. Filter by min_count as set in the yaml
    agg_df = agg_df[agg_df['property_count'] >= min_count]

    # 4. Pivot to wide
    pivot_df = agg_df.pivot(index='property', columns='transform_label', values='property_avg')

    # 5. Normalise each row (property) to z-score for comparability
    pivot_df = pivot_df.apply(lambda row: (row - row.mean()) / row.std(ddof=0) if row.std(ddof=0) > 0 else row - row.mean(), axis=1)

    # 6. Plotting
    n_properties = len(pivot_df)
    n_transforms = pivot_df.shape[1]

    fig_width = max(10, n_transforms * 0.4)
    fig_height = max(4, n_properties * 0.8)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 4], wspace=0.05)

    ax_imgs = fig.add_subplot(gs[0])
    ax_heatmap = fig.add_subplot(gs[1])

    sns.heatmap(
        pivot_df,
        cmap='RdYlGn',
        center=0,
        cbar_kws={"label": "Avg Property Change (z-score)"},
        yticklabels=True,
        xticklabels=False,  # We'll set custom ticks below
        ax=ax_heatmap,
        vmin=-2, vmax=2
    )

    # Add integer index labels to x-axis (another hack)
    ax_heatmap.set_xticks(np.arange(n_transforms) + 0.5)
    ax_heatmap.set_xticklabels([str(i) for i in range(n_transforms)], rotation=90, fontsize=8)

    ax_heatmap.set_xlabel("Transform Index")
    ax_heatmap.set_ylabel("")

    # 8. Add molecule images with arrows
    transform_labels = pivot_df.columns.tolist()
    transform_smiles = [label.split(" → ") for label in transform_labels]
    for i, (from_smiles, to_smiles) in enumerate(transform_smiles):
        from_img = mol_images_from.get(from_smiles)
        to_img = mol_images_to.get(to_smiles)
        y_pos = n_transforms - i - 1

        if from_img:
            ax_imgs.imshow(from_img, extent=(0, 0.4, y_pos + 0.1, y_pos + 0.9), aspect='auto')
        if to_img:
            ax_imgs.imshow(to_img, extent=(0.6, 1, y_pos + 0.1, y_pos + 0.9), aspect='auto')
        ax_imgs.annotate('→', xy=(0.5, y_pos + 0.5), ha='center', va='center', fontsize=14)

    ax_imgs.set_xlim(0, 1)
    ax_imgs.set_ylim(0, n_transforms)
    ax_imgs.axis('off')

    plt.suptitle(f"Multi-Property Change per Transform (Count ≥ {min_count})", fontsize=14)
    heatmap_file = output_dir / "mmp_transform_multi_property_heatmap.png"
    plt.savefig(heatmap_file, bbox_inches='tight', dpi=300)
    plt.close()

    #  4. Network plot with images (needs work)
    sig_df = df_filtered[df_filtered['property_p_value'] < significance_level]

    G = nx.DiGraph()
    for _, row in sig_df.iterrows():
        G.add_edge(row['property_from_smiles'], row['property_to_smiles'], weight=row['property_avg'])

    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, k=0.6, iterations=100)
    nx.draw_networkx_nodes(G, pos, node_size=300, node_color='lightblue')

    edges = G.edges(data=True)
    weights = [e[2]['weight'] for e in edges]
    edge_colors = ['green' if w > 0 else 'red' for w in weights]
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, arrows=True, arrowsize=15)

    ax = plt.gca()
    for node, (x, y) in pos.items():
        img = mol_images_from.get(node) or mol_images_to.get(node)
        if img is not None:
            ab = AnnotationBbox(OffsetImage(img, zoom=0.4), (x, y), frameon=False)
            ax.add_artist(ab)

    plt.title(f"Significant Transform Network (p < {significance_level}, Count ≥ {min_count})")
    plt.axis('off')
    network_file = output_dir / "mmp_transform_network_structures.png"
    plt.savefig(network_file, bbox_inches='tight')
    plt.close()

    #  5. Normalised histograms per property (needs work)
    prop_cols = [
        col for col in df_filtered.columns
        if (col.startswith('property_delta') or col == 'property_avg')
        and col not in ['property_count', 'property_p_value']
    ]

    hist_df = df_filtered.copy()
    for col in prop_cols:
        mean = hist_df[col].mean()
        std = hist_df[col].std()
        hist_df[col + '_norm'] = (hist_df[col] - mean) / std if std > 0 else hist_df[col]

    fig, axes = plt.subplots(len(prop_cols), 1, figsize=(8, 4 * len(prop_cols)))
    if len(prop_cols) == 1:
        axes = [axes]

    for i, col in enumerate(prop_cols):
        ax = axes[i]
        sns.histplot(hist_df[col + '_norm'], bins=50, kde=True, ax=ax)
        ax.axvline(0, color='black', linestyle='--')
        ax.set_title(f"Distribution of Normalized {col} (z-score)")
        ax.set_xlabel("Normalized Value")
        ax.set_ylabel("Number of Transforms")

    plt.tight_layout()
    hist_file = output_dir / "mmp_property_histograms_normalized.png"
    plt.savefig(hist_file)
    plt.close()

    #  6. Top transform tables 
    top_improving = df_filtered.sort_values('property_avg', ascending=False).head(10)
    top_detrimental = df_filtered.sort_values('property_avg', ascending=True).head(10)

    top_improving_file = output_dir / "top_10_improving_transforms_summary.csv"
    top_detrimental_file = output_dir / "top_10_detrimental_transforms_summary.csv"

    top_improving.to_csv(top_improving_file, index=False)
    top_detrimental.to_csv(top_detrimental_file, index=False)

    logger.debug("Top 10 Improving Transforms:")
    logger.debug(top_improving[['property_from_smiles', 'property_to_smiles', 'property_avg']])
    logger.debug("\nTop 10 Detrimental Transforms:")
    logger.debug(top_detrimental[['property_from_smiles', 'property_to_smiles', 'property_avg']])

    return {
        "heatmap_png": str(heatmap_file),
        "network_png": str(network_file),
        "histogram_png": str(hist_file),
        "top_improving_csv": str(top_improving_file),
        "top_detrimental_csv": str(top_detrimental_file),
    }

def compute_fingerprint(smiles, fp_type="morgan", radius=2, nbits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if fp_type.lower() == "morgan":
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    elif fp_type.lower() == "rdkit":
        return Chem.RDKFingerprint(mol)
    else:
        raise ValueError(f"Unsupported fingerprint type: {fp_type}")

def pairwise_similarity(args):
    i, smi_i, fp_i, id_i, act_i, series_i, data, cutoff, delta_threshold, activity_col = args

    cliffs = []

    for j, (smi_j, fp_j, act_j, id_j, series_j) in enumerate(
        zip(
            data["smiles"],
            data["fps"],
            data[activity_col],
            data["FORMATTED_ID"],
            data["Chemical series"]
        )
    ):
        if i >= j:
            continue

        if fp_i is None or fp_j is None:
            continue

        sim = DataStructs.TanimotoSimilarity(fp_i, fp_j)

        if sim >= cutoff:
            act_i = data.loc[i, activity_col]
            delta = abs(act_i - act_j)

            if delta >= delta_threshold:
                cliffs.append({
                    "smiles_1": smi_i,
                    "smiles_2": smi_j,
                    "formatted_id_1": id_i,
                    "formatted_id_2": id_j,
                    "series_1": series_i,
                    "series_2": series_j,
                    "similarity": sim,
                    "delta_activity": delta,
                    "activity_1": act_i,
                    "activity_2": act_j
                })

    return cliffs

@register_task("sar_cliff_analysis", category="Project-based analyses",
               description="Identify and visualise activity cliffs.")
def sar_cliff_analysis(config, data=None):

    #  1. Load configs
    input_file = config.get("input_file")

    sar_cliff_cfg = config.get("sar_cliff_analysis", {})    
    activity_col = sar_cliff_cfg.get("activity_col", "pActivity")
    similarity_cutoff = sar_cliff_cfg.get("similarity_cutoff", 0.85)
    delta_threshold = sar_cliff_cfg.get("delta_threshold", 1.0)
    fp_type = sar_cliff_cfg.get("fp_type", "morgan")
    fp_radius = sar_cliff_cfg.get("fp_radius", 2)
    fp_nbits = sar_cliff_cfg.get("fp_nbits", 2048)
    max_pairs = sar_cliff_cfg.get("max_pairs", 500000)
    output_dir = Path(sar_cliff_cfg.get("output_dir", "outputs/sar_cliffs"))
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.debug(f"[SAR Cliff] Loading input file: {input_file}")
    df = pd.read_csv(input_file)
    # ensure numeric
    df[activity_col] = pd.to_numeric(df[activity_col], errors="coerce")

    # check if series in in the df

    if "Chemical series" not in df.columns:
        df["Chemical series"] = "Unknown"

    # remove NaN and inf
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["smiles", activity_col])

    df = df.reset_index(drop=True)

    #  2. Fingerprint computation 
    logger.debug("[SAR Cliff] Computing fingerprints...")
    df["fps"] = [compute_fingerprint(smi, fp_type, fp_radius, fp_nbits) for smi in df["smiles"]]
    df = df[df["fps"].notna()]

    #  3. Pairwise similarity in parallel 
    logger.debug(f"[SAR Cliff] Computing pairwise similarities (cutoff={similarity_cutoff})...")

    tasks = [
        (
            i,
            row["smiles"],
            row["fps"],
            row["FORMATTED_ID"],
            row[activity_col],
            row["Chemical series"],
            df,
            similarity_cutoff,
            delta_threshold,
            activity_col
        )
        for i, row in df.iterrows()
    ]

    all_cliffs = []
    with ProcessPoolExecutor() as executor:
        for future in as_completed(executor.submit(pairwise_similarity, t) for t in tasks):
            all_cliffs.extend(future.result())
            if len(all_cliffs) > max_pairs:
                break

    cliff_df = pd.DataFrame(all_cliffs)
    cliff_df.sort_values("delta_activity", ascending=False, inplace=True)
    logger.info(f"[SAR Cliff] Found {len(cliff_df)} cliff pairs (Δ≥{delta_threshold}, sim≥{similarity_cutoff}).")

    #  4. Save
    output_csv = output_dir / "sar_cliff_pairs.csv"
    cliff_df.to_csv(output_csv, index=False)
    logger.info(f"[SAR Cliff] Results saved to: {output_csv}")

    #  5. Visualise
    if not cliff_df.empty:
        # Top 10 cliffs - with molecule images

        def delta_bin(delta, step=0.5):
            lower = np.floor(delta / step) * step
            upper = lower + step
            return f"delta_act_{lower:.1f}-{upper:.1f}"

        # where to save images

        img_root = output_dir / "cliff_images"
        img_root.mkdir(exist_ok=True)

        # add colums to df
        
        cliff_df["delta_bin"] = cliff_df["delta_activity"].apply(delta_bin)
        cliff_df = cliff_df[cliff_df["series_1"] == cliff_df["series_2"]]
        cliff_df["series"] = cliff_df["series_1"]

        plt.figure(figsize=(9, 7))

        sns.scatterplot(
            data=cliff_df,
            x="similarity",
            y="delta_activity",
            hue="series",
            palette="tab10",
            s=40,
            alpha=0.7
        )

        plt.axvline(similarity_cutoff, color='red', linestyle='--')
        plt.axhline(delta_threshold, color='orange', linestyle='--')

        plt.xlabel("Tanimoto Similarity")
        plt.ylabel(f"|Δ {activity_col}|")
        plt.title("Activity Cliff Landscape (by Chemical Series)")

        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  # avoids clutter

        scatter_file = output_dir / "sar_cliff_scatter.png"
        plt.savefig(scatter_file, bbox_inches='tight', dpi=300)
        plt.close()

        max_imgs = sar_cliff_cfg.get("max_images_per_bin", None)

        for series_name, series_group in cliff_df.groupby("series"):

            # sanitise folder name
            safe_series = str(series_name).replace(" ", "_").replace("/", "_")
            series_dir = img_root / safe_series
            series_dir.mkdir(exist_ok=True)

            for bin_name, bin_group in series_group.groupby("delta_bin"):

                if max_imgs:
                    bin_group = bin_group.head(max_imgs)

                bin_dir = series_dir / bin_name
                bin_dir.mkdir(parents=True, exist_ok=True)

                for idx, row in bin_group.iterrows():
                    mol1 = Chem.MolFromSmiles(row["smiles_1"])
                    mol2 = Chem.MolFromSmiles(row["smiles_2"])

                    img = Draw.MolsToGridImage(
                        [mol1, mol2],
                        legends=[
                            f"{row['formatted_id_1']}\nAct={row['activity_1']:.2f}\nSim={row['similarity']:.2f}",
                            f"{row['formatted_id_2']}\nAct={row['activity_2']:.2f}\nSim={row['similarity']:.2f}",
                        ],
                        molsPerRow=2,
                        subImgSize=(250, 250)
                    )

                    img_path = bin_dir / f"cliff_{idx}.png"
                    img.save(img_path)

        return {
            "cliff_csv": str(output_csv),
            "scatter_png": str(scatter_file),
        }

    else:
        logger.warning("[SAR Cliff] No activity cliffs found.")
        return None

try:
    import umap
except ImportError:
    umap = None
    logger.warning("UMAP not installed. Falling back to t-SNE for dimensionality reduction.")


@register_task(
    "chemical_space_drift",
    category="Project-based analyses",
    description="Track week-to-week evolution of chemical space via (UMAP / t-SNE)."
)
def chemical_space_drift(config, data=None):
    """Visualise chemical space evolution colored by time or activity."""

    # 1. Load configs
    task_conf = config.get("chemical_space_drift", {})
    input_file = config.get("input_file")
    output_dir = Path(task_conf.get("output_dir", "outputs/chem_space"))
    output_dir.mkdir(parents=True, exist_ok=True)

    fp_radius = task_conf.get("fp_radius", 2)
    fp_nbits = task_conf.get("fp_nbits", 1024)
    reducer_type = task_conf.get("reducer", "umap").lower()
    activity_col = task_conf.get("activity_col", "pActivity")
    date_col = task_conf.get("date_col", "date")
    color_by = task_conf.get("color_by", "week")

    logger.debug("[ChemicalSpaceDrift] Starting task...")

    # 2. Load inputs
    df = None
    if data is not None:
        # Case 1: Tuple
        if isinstance(data, tuple) and len(data) == 2 and isinstance(data[1], pd.DataFrame):
            df = data[1]
        # Case 2: Dict of task outputs
        elif isinstance(data, dict):
            for v in data.values():
                if isinstance(v, pd.DataFrame):
                    df = v
                    break
        # Case 3: Already a DataFrame
        elif isinstance(data, pd.DataFrame):
            df = data

    # Fallback to CSV
    if df is None and input_file and Path(input_file).exists():
        df = pd.read_csv(input_file)

    if df is None or not isinstance(df, pd.DataFrame):
        raise TypeError(f"[ChemicalSpaceDrift] Expected DataFrame input, got {type(df)}")

    logger.info(f"[ChemicalSpaceDrift] Data loaded with {len(df)} rows and {len(df.columns)} columns.")

    #  3. Which columns is SMILES? Let's give options and hope I don't accidentally mess this up later... 
    smiles_col = next((c for c in ["smiles", "SMILES", "canonical_smiles"] if c in df.columns), None)
    if smiles_col is None:
        raise KeyError(f"[ChemicalSpaceDrift] No 'smiles' column found. Available: {list(df.columns)}")

    #  4. Generate fingerprints 
    logger.debug(f"[ChemicalSpaceDrift] Generating Morgan fingerprints (radius={fp_radius}, nBits={fp_nbits})...")
    fps = []
    for smi in tqdm(df[smiles_col], desc="Computing fingerprints"):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            fps.append(np.zeros(fp_nbits, dtype=int))
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=fp_radius, nBits=fp_nbits)
        arr = np.zeros((fp_nbits,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fps.append(arr)
    fps = np.array(fps)

    #  5. Dimensionality reduction 
    X = StandardScaler().fit_transform(fps)
    if reducer_type == "tsne" or umap is None:
        reducer = TSNE(n_components=2, random_state=42, perplexity=30)
        method = "t-SNE"
    else:
        reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
        method = "UMAP"

    logger.debug(f"[ChemicalSpaceDrift] Performing {method} reduction...")
    embedding = reducer.fit_transform(X)
    df["x"], df["y"] = embedding[:, 0], embedding[:, 1]

    #  6. Handle time grouping 
    if date_col in df.columns:
        df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
        df["week"] = df[date_col].dt.strftime("%Y-%U")

        # Specify to most recent N weeks from yaml
        recent_weeks = task_conf.get("recent_weeks", 12)
        unique_weeks = sorted(df["week"].dropna().unique())
        if len(unique_weeks) > recent_weeks:
            recent_weeks_set = set(unique_weeks[-recent_weeks:])
            df = df[df["week"].isin(recent_weeks_set)]
            logger.info(f"[ChemicalSpaceDrift] Filtering to last {recent_weeks} weeks "
                        f"({unique_weeks[-recent_weeks]} → {unique_weeks[-1]}).")
    else:
        df["week"] = "unknown"

    #  7. Visualisation 
    color_label = color_by if color_by in df.columns else (
        activity_col if activity_col in df.columns else "week"
    )

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="x", y="y",
        hue=df[color_label],
        palette="viridis", s=40, alpha=0.8
    )
    plt.title(f"{method} projection of chemical space – colored by {color_label}")
    plt.tight_layout()

    plot_path = output_dir / f"chemical_space_{method.lower()}_weekly.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()

    emb_csv = output_dir / "chemical_space_embedding.csv"
    df.to_csv(emb_csv, index=False)
    logger.info(f"[ChemicalSpaceDrift] Results saved: {emb_csv}")

    logger.info(f"[ChemicalSpaceDrift] Plot saved: {plot_path}")

    return {
        "embedding_csv": str(emb_csv),
        "projection_plot": str(plot_path),
        "n_molecules": len(df),
        "method": method
    }

@register_task(
    "scaffold_enrichment_trends",
    category="Project-based analyses",
    description="Track scaffold enrichment trends over time."
)
def scaffold_enrichment_trends(config, data=None):
    """Compute weekly scaffold enrichment trends."""

    # 1. Load config
    task_conf = config.get("scaffold_enrichment_trends", {})
    input_file = config.get("input_file")
    output_dir = Path(task_conf.get("output_dir", "outputs/scaffold_enrichment"))
    output_dir.mkdir(parents=True, exist_ok=True)

    activity_col = task_conf.get("activity_col", "pActivity")
    date_col = task_conf.get("date_col", "date")

    logger.info("[ScaffoldEnrichment] Starting task...")

    df = None
    if data is not None:
        if isinstance(data, tuple) and len(data) == 2 and isinstance(data[1], pd.DataFrame):
            df = data[1]
        elif isinstance(data, dict):
            for v in data.values():
                if isinstance(v, pd.DataFrame):
                    df = v
                    break
        elif isinstance(data, pd.DataFrame):
            df = data

    if df is None and input_file and Path(input_file).exists():
        df = pd.read_csv(input_file)

    if df is None or not isinstance(df, pd.DataFrame):
        raise TypeError(f"[ScaffoldEnrichment] Expected DataFrame input, got {type(df)}")

    logger.info(f"[ScaffoldEnrichment] Data loaded with {len(df)} rows and {len(df.columns)} columns.")

    smiles_col = next((c for c in ["smiles", "SMILES", "canonical_smiles"] if c in df.columns), None)
    if smiles_col is None:
        raise KeyError(f"[ScaffoldEnrichment] No 'smiles' column found. Available: {list(df.columns)}")

    if date_col not in df.columns:
        logger.warning(f"[ScaffoldEnrichment] No '{date_col}' column found — using 'week'='unknown'.")
        df["week"] = "unknown"
    else:
        df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
        df["week_period"] = df[date_col].dt.to_period("W")
        df["week"] = df["week_period"].astype(str)  # for plotting

    # 2. Bemis–Murcko scaffolds
    logger.debug("[ScaffoldEnrichment] Computing Murcko scaffolds...")
    method = task_conf.get("scaffold_method", "murcko")
    min_size = task_conf.get("min_scaffold_size", 0)

    scaffolds = []
    for smi in tqdm(df[smiles_col], desc="Extracting scaffolds"):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            scaffolds.append(None)
            continue

        # if min_size and mol.GetNumHeavyAtoms() < min_size:
        #     scaffolds.append(None)
        #     continue

        try:
            scaf = None

            if method == "murcko":
                scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)

            elif method == "generic":
                scaf_mol = MurckoScaffold.GetScaffoldForMol(mol)
                scaf_mol = MurckoScaffold.MakeScaffoldGeneric(scaf_mol)
                scaf = Chem.MolToSmiles(scaf_mol)

            elif method == "brics":
                frags = BRICS.BRICSDecompose(mol)
                scaf = max(frags, key=len) if frags else None

            # canonicalize if valid
            if scaf:
                scaf_mol = Chem.MolFromSmiles(scaf)
                if min_size and scaf_mol and scaf_mol.GetNumHeavyAtoms() < min_size:
                    scaf = None

            scaffolds.append(scaf)

        except Exception:
            scaffolds.append(None)

    df["scaffold"] = scaffolds

    df_valid = df.dropna(subset=["scaffold"]).copy()
    
    if df_valid.empty:
        logger.warning("[ScaffoldEnrichment] No valid scaffolds found.")
        return None

    #  3. Weekly scaffold stats 
    logger.debug("[ScaffoldEnrichment] Aggregating weekly scaffold statistics...")
    agg_df = (
        df_valid.groupby(["week_period", "scaffold"])
        .agg(
            count=("scaffold", "size"),
            mean_activity=(activity_col, "mean")
        )
        .reset_index()
    )
    
    agg_df = agg_df.sort_values("week_period")
    agg_df["week_num"] = agg_df["week_period"].rank(method="dense").astype(int)

    agg_df = agg_df.sort_values(["scaffold", "week_period"])

    agg_df["prev_mean"] = agg_df.groupby("scaffold")["mean_activity"].shift(1)
    agg_df["delta_mean"] = agg_df["mean_activity"] - agg_df["prev_mean"]

    #  4. Identify improving / declining scaffolds

    # Use last N weeks instead of single latest week
    window_weeks = task_conf.get("trend_window_weeks", 4)
    print(window_weeks)

    recent_weeks = sorted(agg_df["week_period"].dropna().unique())[-window_weeks:]

    recent_data = agg_df[agg_df["week_period"].isin(recent_weeks)].copy()

    # Keep only valid deltas
    recent_data = recent_data.dropna(subset=["delta_mean"])

    # Aggregate trend per scaffold (mean delta over window)
    trend_df = (
        recent_data.groupby("scaffold")
        .agg(
            mean_delta=("delta_mean", "mean"),
            n_points=("delta_mean", "count"),
            latest_activity=("mean_activity", "last")
        )
        .reset_index()
    )

    # Optional robustness filter
    trend_df = trend_df[trend_df["n_points"] >= 2]

    # Fallback if too strict
    if trend_df.empty:
        logger.warning("[ScaffoldEnrichment] No valid multi-week trends - relaxing constraints.")
        trend_df = (
            recent_data.groupby("scaffold")
            .agg(mean_delta=("delta_mean", "mean"))
            .reset_index()
        )

    top_improving = trend_df.nlargest(5, "mean_delta")
    top_declining = trend_df.nsmallest(5, "mean_delta")

    # 5. Visualisation
    image_zoom = 0.5  # hardcoding size of 2D structures for the plots. Might need to adjust.

    top_n = task_conf.get("top_n_scaffolds", 8)
    top_scaffolds = (
        agg_df.groupby("scaffold")["count"].sum()
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )
    agg_df["week"] = agg_df["week_period"].astype(str)
    plot_df = agg_df[agg_df["scaffold"].isin(top_scaffolds)]

    week_labels = (
        agg_df[["week_num", "week"]]
        .drop_duplicates()
        .sort_values("week_num")
    )

    tick_step = task_conf.get("xtick_step", 4)

    week_labels_sub = week_labels.iloc[::tick_step]

    n_ticks = len(week_labels_sub)
    width = max(10, min(20, n_ticks * 1.2))
    plt.figure(figsize=(width, 6))
    
    ax = plt.gca()

    from PIL import Image, ImageDraw

    def make_placeholder_image(text="NA", size=(100, 100)):
        img = Image.new("RGB", size, color="white")
        draw = ImageDraw.Draw(img)
        draw.text((10, 40), text, fill="black")
        return img

    scaffold_imgs = {}

    for smi in top_scaffolds:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                img = Draw.MolToImage(mol, size=(100, 100))
            else:
                raise ValueError("Invalid SMILES")

        except Exception:
            logger.warning(f"[ScaffoldEnrichment] Failed to render scaffold: {smi}")
            img = make_placeholder_image("Invalid")

        scaffold_imgs[smi] = img

    # Helper for mapping week labels to x-axis positions
    xtick_labels = [label.get_text() for label in ax.get_xticklabels()]
    xtick_positions = ax.get_xticks()
    def get_x_position(x_val):
        if isinstance(x_val, (int, float)):
            return x_val
        try:
            idx = xtick_labels.index(str(x_val))
            return xtick_positions[idx]
        except ValueError:
            return xtick_positions[-1]

    # Split scaffolds into multi-point and single-point sets, depending on how many data points there are.
    multi_point_scaffolds = []
    single_point_scaffolds = []
    for scaffold in top_scaffolds:
        count_points = len(plot_df[plot_df["scaffold"] == scaffold])
        if count_points > 1:
            multi_point_scaffolds.append(scaffold)
        else:
            single_point_scaffolds.append(scaffold)

    # Plot lines for scaffolds with multiple points
    palette = sns.color_palette("tab10", n_colors=len(multi_point_scaffolds))
    color_map = dict(zip(multi_point_scaffolds, palette))

    for scaffold in multi_point_scaffolds:
        scaffold_data = plot_df[plot_df["scaffold"] == scaffold]

        sns.lineplot(
            data=scaffold_data,
            x="week_num",
            y="mean_activity",
            color=color_map[scaffold],
            lw=2,
            ax=ax,
        )
    
    ax.set_xticks(week_labels_sub["week_num"])
    ax.set_xticklabels(week_labels_sub["week"], rotation=45)

    # Remove default legend if present (hack)
    if ax.legend_:
        ax.legend_.remove()
    
    logger.debug(f"Multi_point_scaffolds: {len(multi_point_scaffolds)}")

    for scaffold in multi_point_scaffolds:
        scaffold_data = plot_df[plot_df["scaffold"] == scaffold]
        logger.debug(
            f"Scaffold={scaffold[:30]}... "
            f"points={len(scaffold_data)} "
            f"has_img={scaffold in scaffold_imgs}"
        )

    # Add images for line plots at the far right of the plot 
    y_min, y_max = ax.get_ylim()
    y_spacing = (y_max - y_min) / (len(multi_point_scaffolds) + 1)

    # --- prevent overlapping endpoint images ---

    # 1. Collect positions ONCE
    scaffold_positions = []

    for scaffold in multi_point_scaffolds:
        scaffold_data = plot_df[plot_df["scaffold"] == scaffold]
        if scaffold_data.empty:
            continue

        last_row = (
            scaffold_data
            .sort_values(["week_period", "mean_activity"])
            .drop_duplicates(subset=["week_period"], keep="last")
            .iloc[-1]
        )

        scaffold_positions.append({
            "scaffold": scaffold,
            "x": float(last_row["week_num"]),
            "y": float(last_row["mean_activity"])
        })

    # 2. Sort by y (stable ordering)
    scaffold_positions = sorted(scaffold_positions, key=lambda d: d["y"])

    # 3. Place with alternating offsets (prevents stacking)
    placed_y = []
    min_sep = (y_max - y_min) * 0.05

    for i, item in enumerate(scaffold_positions):
        y = item["y"]

        # alternate direction: up, down, up, down...
        direction = 1 if i % 2 == 0 else -1

        # resolve collisions iteratively
        for _ in range(10):  # small bounded loop = stable
            collision = False
            for py in placed_y:
                if abs(y - py) < min_sep:
                    y += direction * min_sep
                    direction *= -1  # flip direction
                    collision = True
            if not collision:
                break

        # keep INSIDE bounds with margin (prevents edge stacking)
        margin = min_sep
        y = max(y_min + margin, min(y, y_max - margin))

        placed_y.append(y)

        scaffold = item["scaffold"]
        x = min(item["x"], ax.get_xlim()[1] - 0.5)

        img = scaffold_imgs.get(scaffold)
        if img:
            imagebox = offsetbox.OffsetImage(img, zoom=image_zoom)

            ab = offsetbox.AnnotationBbox(
                imagebox,
                (x, y),
                xybox=(25, 0),
                xycoords='data',
                boxcoords="offset points",
                frameon=True,
                bboxprops=dict(
                    edgecolor=color_map[scaffold],
                    linewidth=2
                ),
                annotation_clip=False
            )
            ax.add_artist(ab)

    # ax.set_xlim(x_left, x_right + x_offset * 2)

    plt.title(f"Top {top_n} Scaffold Enrichment Trends (mean {activity_col})")
    plt.xticks(rotation=45)
    plt.xlabel("Week")
    plt.ylabel(f"Mean {activity_col}")
    plt.tight_layout()

    trend_plot = output_dir / "scaffold_trends_plot.png"
    plt.savefig(trend_plot, dpi=300)
    plt.close()

    trends_csv = output_dir / "scaffold_trends.csv"
    agg_df.to_csv(trends_csv, index=False)

    logger.info(f"[ScaffoldEnrichment] Trends saved to: {trends_csv}")
    logger.info(f"[ScaffoldEnrichment] Plot saved to: {trend_plot}")

    improving_csv = output_dir / "top5_improving_scaffolds.csv" 
    declining_csv = output_dir / "top5_declining_scaffolds.csv" 
    top_improving.to_csv(improving_csv, index=False) 
    top_declining.to_csv(declining_csv, index=False)

    return {
        "trends_csv": str(trends_csv),
        "top_improving_csv": str(improving_csv),
        "top_declining_csv": str(declining_csv),
        "trend_plot": str(trend_plot),
        "n_valid_scaffolds": len(df_valid)
    }

@register_task("physchem_property_drift", category="Project-based analyses",
               description="Monitor property balance and drift over time (MW, clogP, TPSA, etc.).")
def physchem_property_drift(config, data=None):
    """
    Track trends in key physicochemical properties across weekly design cycles.
    """

    def penalty(val, low, high):
        if pd.isna(val):
            return np.nan
        if low <= val <= high:
            return 0
        return min(abs(val - low), abs(val - high))

    def compute_slope(series):
        y = series.values
        x = np.arange(len(y))
        if len(y) < 3 or np.all(np.isnan(y)):
            return np.nan
        return linregress(x, y).slope

    def direction(val, threshold=0.2):
        if pd.isna(val):
            return "NA"
        if val > threshold:
            return "↑"
        elif val < -threshold:
            return "↓"
        return "→"

    def window_shift_test(prev, curr):
        if len(prev) < 2 or len(curr) < 2:
            return np.nan
        return ttest_ind(prev, curr, nan_policy='omit').pvalue
        
    PROPERTY_TARGETS = {
        "MW": (250, 500),
        "clogP": (1.0, 3.5),
        "TPSA": (40, 100),
        "HBD": (0, 3),
        "HBA": (2, 8),
        "RotBonds": (0, 5),
    }

    PROPERTY_WEIGHTS = {k: 1.0 for k in PROPERTY_TARGETS}

    DESCRIPTOR_FUNCS = {
        "MW": Descriptors.MolWt,
        "clogP": Crippen.MolLogP,
        "TPSA": rdMolDescriptors.CalcTPSA,
        "HBD": rdMolDescriptors.CalcNumHBD,
        "HBA": rdMolDescriptors.CalcNumHBA,
        "Fsp3": rdMolDescriptors.CalcFractionCSP3,
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds,
        "RingCount": rdMolDescriptors.CalcNumRings,
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings,
        "HeavyAtoms": Descriptors.HeavyAtomCount,
    }

    # 1. Load configs
    input_file = config.get("input_file")
    drift_cfg = config.get("physchem_property_drift", {})

    date_col = drift_cfg.get("date_col", "publication_date")
    output_dir = Path(drift_cfg.get("output_dir", "outputs/physchem_drift"))
    rolling_window = drift_cfg.get("rolling_window", 4)  # e.g. 4 weeks
    output_dir.mkdir(parents=True, exist_ok=True)
    
    trend_dir = output_dir / "trends"
    penalty_dir = output_dir / "penalties"
    multi_dir = output_dir / "multivariate"
    slope_dir = output_dir / "slopes"
    direction_dir = output_dir / "directional"
    dist_dir = output_dir / "distributions"
    summary_dir = output_dir / "summaries"

    for d in [trend_dir, penalty_dir, multi_dir, slope_dir, direction_dir, dist_dir, summary_dir]:
        d.mkdir(parents=True, exist_ok=True)

    logger.info(f"[PhysChem Drift] Loading input file: {input_file}")
    df = pd.read_csv(input_file)
    df = df.dropna(subset=["smiles", date_col])
    df[date_col] = pd.to_datetime(df[date_col])

    #  2. Calculate properties 
    logger.debug("[PhysChem Drift] Calculating basic properties...")
    mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]

    # Manual list (DESCRIPTOR_FUNCS) - perhaps we should automate or select in the yaml

    for name, func in DESCRIPTOR_FUNCS.items():
        df[name] = [func(m) if m else np.nan for m in mols]

    # 3. Weekly aggregation
    df["week"] = df[date_col].dt.to_period("W").apply(lambda r: r.start_time)
    props = list(DESCRIPTOR_FUNCS.keys())
    series_col = drift_cfg.get("series_col", None)

    if series_col and series_col in df.columns:
        grouped = df.groupby([series_col, "week"])[props].mean()
        agg = (
            grouped
            .groupby(level=0)
            .rolling(window=rolling_window, min_periods=1)
            .mean()
            .droplevel(0) 
        )
    else:
        grouped = df.groupby("week")[props].mean()
        agg = grouped.rolling(window=rolling_window, min_periods=1).mean()
    
    numeric_cols = agg.select_dtypes(include=[np.number]).columns

    # 3a property penalty

    penalty_df = agg.copy()

    for p in PROPERTY_TARGETS:
        if p in penalty_df.columns:
            low, high = PROPERTY_TARGETS[p]
            penalty_df[p] = penalty_df[p].apply(lambda x: penalty(x, low, high))

    penalty_df["total_penalty"] = sum(
        PROPERTY_WEIGHTS[p] * penalty_df[p]
        for p in PROPERTY_TARGETS if p in penalty_df.columns
    )

    # 3b multivariate shift

    valid = agg[numeric_cols].dropna()

    if len(valid) > 5:
        cov = np.cov(valid.T)
        inv_cov = np.linalg.pinv(cov)
        mean_vec = valid.mean()

        def mahal(row):
            if row.isna().any():
                return np.nan
            return mahalanobis(row, mean_vec, inv_cov)

        agg["multivariate_drift"] = agg[numeric_cols].apply(mahal, axis=1)
    else:
        agg["multivariate_drift"] = np.nan

    # 3c Slopes

    slopes = {}

    if series_col and series_col in df.columns:
        for series in df[series_col].dropna().unique():
            if series not in agg.index.get_level_values(0):
                continue
            sub = agg.xs(series, level=0).sort_index()
            slopes[series] = {
                col: compute_slope(sub[col]) for col in numeric_cols
            }
    else:
        slopes = {col: compute_slope(agg[col]) for col in numeric_cols}

    # 3d Statistical shift

    shift_pvals = {}

    def split_windows(series, window):
        if len(series) < window * 2:
            return None, None
        return series.iloc[-2*window:-window], series.iloc[-window:]

    if series_col and series_col in df.columns:
        for series in df[series_col].dropna().unique():
            if series not in agg.index.get_level_values(0):
                continue
            sub = agg.xs(series, level=0).sort_index()

            shift_pvals[series] = {}
            for col in numeric_cols:
                prev, curr = split_windows(sub[col].dropna(), rolling_window)
                if prev is not None:
                    shift_pvals[series][col] = window_shift_test(prev, curr)
    else:
        shift_pvals = {}
        for col in numeric_cols:
            prev, curr = split_windows(agg[col].dropna(), rolling_window)
            if prev is not None:
                shift_pvals[col] = window_shift_test(prev, curr)
    
    # 3e Distribution shift

    dist_shift = {}

    if series_col and series_col in df.columns:
        for series in df[series_col].dropna().unique():
            sub_df = df[df[series_col] == series]
            sub_df = sub_df.sort_values("week")

            dist_shift[series] = {}

            for col in numeric_cols:
                prev = sub_df.iloc[:-rolling_window][col].dropna()
                curr = sub_df.iloc[-rolling_window:][col].dropna()

                if len(prev) > 1 and len(curr) > 1:
                    dist_shift[series][col] = wasserstein_distance(prev, curr)
    else:
        for col in numeric_cols:
            prev = df.iloc[:-rolling_window][col].dropna()
            curr = df.iloc[-rolling_window:][col].dropna()

            if len(prev) > 1 and len(curr) > 1:
                dist_shift[col] = wasserstein_distance(prev, curr)
    
    # 3f Directional fingerprint

    directional = {}

    if series_col and series_col in df.columns:
        for series in df[series_col].dropna().unique():
            if series not in agg.index.get_level_values(0):
                continue

            sub = agg.xs(series, level=0).sort_index()

            if len(sub) > 1:
                delta = sub.iloc[-1] - sub.iloc[-2]
                directional[series] = {
                    col: direction(delta[col]) for col in numeric_cols
                }
    else:
        if len(agg) > 1:
            delta = agg.iloc[-1] - agg.iloc[-2]
            directional = {
                col: direction(delta[col]) for col in numeric_cols
            }
    
    # 3g alerts

    alerts = {}

    def generate_alerts(latest, slope_dict, pvals, multi_drift):
        a = []

        if "clogP" in slope_dict and slope_dict["clogP"] > 0.2:
            a.append("Rapid lipophilicity increase")

        if "TPSA" in latest and latest["TPSA"] < 50:
            a.append("Low polarity risk")

        if "MW" in latest and latest["MW"] > 500:
            a.append("High molecular weight risk")

        if multi_drift is not None and multi_drift > 3:
            a.append("Significant multivariate drift")

        for k, v in (pvals or {}).items():
            if v is not None and v < 0.05:
                a.append(f"Significant shift in {k}")

        return list(set(a))

    if series_col and series_col in df.columns:
        for series in df[series_col].dropna().unique():
            if series not in agg.index.get_level_values(0):
                continue

            sub = agg.xs(series, level=0).sort_index()
            latest = sub.iloc[-1]

            alerts[series] = generate_alerts(
                latest,
                slopes.get(series, {}),
                shift_pvals.get(series, {}),
                latest.get("multivariate_drift", None)
            )
    else:
        latest = agg.iloc[-1]
        alerts = generate_alerts(
            latest,
            slopes,
            shift_pvals,
            latest.get("multivariate_drift", None)
        )
    
    penalty_df.to_csv(summary_dir / "property_penalty.csv")

    pd.DataFrame.from_dict(slopes, orient="index").to_csv(trend_dir / "trend_slopes.csv")
    pd.DataFrame.from_dict(shift_pvals, orient="index").to_csv(summary_dir / "shift_pvalues.csv")
    pd.DataFrame.from_dict(dist_shift, orient="index").to_csv(dist_dir / "distribution_shift.csv")
    pd.DataFrame.from_dict(directional, orient="index").to_csv(direction_dir / "directional_fingerprint.csv", encoding="utf-8-sig")

    with open(output_dir / "alerts.json", "w") as f:
        json.dump(alerts, f, indent=2)
    
    #### plot the above quantities

    logger.info("[PhysChem Drift] Plotting advanced diagnostics...")

    # A. Penalty trend (stacked by property)
    plt.figure(figsize=(10,6))

    if series_col and series_col in df.columns:
        for series in penalty_df.index.get_level_values(0).unique():
            sub = penalty_df.xs(series, level=0).sort_index()

            props = [p for p in PROPERTY_TARGETS if p in sub.columns]
            stacked = sub[props]

            stacked.plot.area(alpha=0.6)

            plt.title(f"Property Penalty Decomposition — {series}")
            plt.ylabel("Penalty Contribution")
            plt.xlabel("Week")

            plt.savefig(penalty_dir / f"penalty_breakdown_{series}.png", dpi=300)
            plt.close()

    else:
        props = [p for p in PROPERTY_TARGETS if p in penalty_df.columns]
        penalty_df[props].plot.area(alpha=0.6)

        plt.title("Property Penalty Decomposition")
        plt.ylabel("Penalty Contribution")
        plt.xlabel("Week")

        plt.savefig(penalty_dir / "penalty_breakdown.png", dpi=300)
        plt.close()


    # B1. Multivariate drift heatmap over time

    try:
        if series_col and series_col in df.columns:

            # reshape to Series x Week
            heat_df = (
                agg["multivariate_drift"]
                .unstack(level=0)   # columns = series
                .T                 # rows = series
            )

            # deal with dates
            
            heat_df.columns = pd.to_datetime(heat_df.columns).strftime("%Y-%m-%d")

            def series_key(x):
                # split numeric prefix + suffix (handles 4b, 4c)
                import re
                m = re.match(r"(\d+)", str(x))
                return int(m.group(1)) if m else float("inf")

            heat_df = heat_df.loc[sorted(heat_df.index, key=series_key)]

            fig, ax = plt.subplots(figsize=(12, max(4, len(heat_df) * 0.4)))

            sns.heatmap(
                heat_df,
                cmap="plasma_r",
                cbar_kws={"label": "Mahalanobis Drift"},
                linewidths=0.5,
                ax=ax
            )

            ax.set_title("Multivariate Drift Over Time (All Series)")
            ax.set_xlabel("Week")
            ax.set_ylabel("Series")

            # --- FIX: ensure full x-axis labels are visible ---
            ax.set_xticklabels(
                ax.get_xticklabels(),
                rotation=45,
                ha="right"
            )

            plt.tight_layout()

            # extra safety margin (prevents clipping on savefig)
            plt.subplots_adjust(bottom=0.25)

            plt.savefig(
                multi_dir / "multivariate_timeseries_heatmap.png",
                dpi=300,
                bbox_inches="tight"
            )
            plt.close()

    except Exception as e:
        logger.warning(f"Multivariate time heatmap failed: {e}")

    # --- B2. Multivariate drift heatmap (latest snapshot) ---
    try:
        if series_col and series_col in df.columns:
            latest_vals = {}

            for series in agg.index.get_level_values(0).unique():
                sub = agg.xs(series, level=0).sort_index()
                if len(sub) > 0:
                    latest_vals[series] = sub["multivariate_drift"].iloc[-1]

            if latest_vals:
                heat_df = pd.DataFrame.from_dict(latest_vals, orient="index", columns=["Drift"])

                plt.figure(figsize=(4, max(4, len(heat_df)*0.4)))
                sns.heatmap(
                    heat_df,
                    cmap="Reds",
                    annot=True,
                    fmt=".2f"
                )

                plt.title("Latest Multivariate Drift (by Series)")

                plt.savefig(multi_dir / "multivariate_heatmap.png", dpi=300)
                plt.close()

    except Exception as e:
        logger.warning(f"Multivariate heatmap failed: {e}")

    # --- C. Directional heatmap (robust) ---
    try:
        mapping = {"↑": 1, "↓": -1, "→": 0, "NA": 0}

        dir_df = pd.DataFrame.from_dict(directional, orient="index")

        if dir_df.empty:
            logger.warning("Directional heatmap skipped: no data")
        else:
            dir_numeric = dir_df.replace(mapping)
            dir_numeric = dir_numeric.apply(pd.to_numeric, errors="coerce").fillna(0).astype(float)

            plt.figure(figsize=(10,4))
            sns.heatmap(
                dir_numeric,
                cmap="coolwarm",
                center=0,
                annot=dir_df,
                fmt=""
            )

            plt.title("Directional Drift Fingerprint")
            plt.xlabel("Property")
            plt.ylabel("Series")

            plt.savefig(direction_dir / "directional_heatmap.png", dpi=300)
            plt.close()

    except Exception as e:
        logger.warning(f"Directional heatmap failed: {e}")

    # --- D. Slope bar charts (per series) ---
    if series_col and isinstance(slopes, dict):
        for series, vals in slopes.items():
            if not isinstance(vals, dict):
                continue

            plt.figure(figsize=(8,4))
            pd.Series(vals).sort_values().plot(kind="barh")

            plt.title(f"Drift Velocity — {series}")
            plt.xlabel("Slope")

            plt.savefig(slope_dir / f"slope_{series}.png", dpi=300)
            plt.close()


    # --- E. Alerts text summary (human readable) ---
    with open(summary_dir / "alerts.txt", "w") as f:
        for series, msgs in alerts.items():
            f.write(f"{series}:\n")
            for m in msgs:
                f.write(f"  - {m}\n")
            f.write("\n")

    # --- F. Distribution shift (per series, FIXED) ---
    key_props = ["MW", "clogP", "TPSA"]

    if series_col and series_col in df.columns:
        for series in df[series_col].dropna().unique():
            sub_df = df[df[series_col] == series].sort_values("week")

            for col in key_props:
                if col not in sub_df.columns:
                    continue

                prev = sub_df.iloc[:-rolling_window][col].dropna()
                curr = sub_df.iloc[-rolling_window:][col].dropna()

                if len(prev) >= 3 and len(curr) >= 3:
                    try:
                        plt.figure(figsize=(6,4))
                        sns.kdeplot(prev, label="Previous")
                        sns.kdeplot(curr, label="Current")

                        plt.title(f"{series} — Distribution Shift: {col}")
                        plt.legend()

                        plt.savefig(dist_dir / f"{series}_dist_shift_{col}.png", dpi=300)
                        plt.close()

                    except Exception as e:
                        logger.warning(f"Dist plot failed for {series}-{col}: {e}")
    else:
        for col in key_props:
            if col not in df.columns:
                continue

            prev = df.iloc[:-rolling_window][col].dropna()
            curr = df.iloc[-rolling_window:][col].dropna()

            if len(prev) >= 3 and len(curr) >= 3:
                try:
                    plt.figure(figsize=(6,4))
                    sns.kdeplot(prev, label="Previous")
                    sns.kdeplot(curr, label="Current")

                    plt.title(f"Distribution Shift: {col}")
                    plt.legend()

                    plt.savefig(dist_dir / f"dist_shift_{col}.png", dpi=300)
                    plt.close()

                except Exception as e:
                    logger.warning(f"Dist plot failed for {col}: {e}")

    #  4. Trend visualisation 
    logger.info("[PhysChem Drift] Plotting property trends...")

    trend_file = None
    numeric_cols = agg.select_dtypes(include=[np.number]).columns
    deltas = {}

    # --- Z-score normalization ---
    agg_z = agg.copy()

    for col in numeric_cols:
        if agg[col].std() == 0 or agg[col].isna().all():
            agg_z[col] = 0
        else:
            agg_z[col] = (agg[col] - agg[col].mean()) / agg[col].std()

    # --- Define logical property groups + palettes ---
    property_groups = {
        "Size & Lipophilicity": ["MW", "clogP", "TPSA"],
        "HBD/HBA": ["HBD", "HBA"],
        "3D Character": ["Fsp3"],
        "Topology": ["RotBonds", "RingCount", "AromaticRings", "HeavyAtoms"],
    }

    group_palettes = {
        "Size & Lipophilicity": plt.cm.Blues,
        "HBD/HBA": plt.cm.Greens,
        "3D Character": plt.cm.Purples,
        "Topology": plt.cm.Oranges,
    }

    if series_col and series_col in df.columns:
        series_values = df[series_col].dropna().unique()

        for series in series_values:
            if series not in agg.index.get_level_values(0):
                continue

            sub = agg_z.xs(series, level=0).sort_index()
            sub_raw = agg.xs(series, level=0).sort_index()

            n_groups = len(property_groups)
            fig, axes = plt.subplots(n_groups, 1, figsize=(10, 3 * n_groups), sharex=True)

            if n_groups == 1:
                axes = [axes]

            for ax, (group_name, group_props) in zip(axes, property_groups.items()):
                palette = group_palettes[group_name]
                valid_props = [p for p in group_props if p in sub.columns]

                colors = palette(np.linspace(0.4, 0.9, len(valid_props)))

                for p, c in zip(valid_props, colors):
                    ax.plot(sub.index, sub[p], marker='o', label=p, color=c)

                ax.axhline(0, linestyle="--", linewidth=1)
                ax.set_title(group_name)
                ax.grid(True)
                ax.legend(fontsize=8)

            axes[-1].set_xlabel("Week")
            fig.suptitle(f"Property Trends (Z-score) — {series}", fontsize=14)

            plt.tight_layout(rect=[0, 0, 1, 0.97])
            plt.savefig(trend_dir / f"trend_{series}.png", dpi=300)
            plt.close()

            # --- deltas
            if len(sub_raw) > 1:
                latest = sub_raw.iloc[-1][numeric_cols]
                prev = sub_raw.iloc[-2][numeric_cols]
                deltas[series] = (latest - prev).to_dict()

    else:
        trend_file = output_dir / "property_trends.png"

        n_groups = len(property_groups)
        fig, axes = plt.subplots(n_groups, 1, figsize=(10, 3 * n_groups), sharex=True)

        if n_groups == 1:
            axes = [axes]

        for ax, (group_name, group_props) in zip(axes, property_groups.items()):
            palette = group_palettes[group_name]
            valid_props = [p for p in group_props if p in agg_z.columns]

            colors = palette(np.linspace(0.4, 0.9, len(valid_props)))

            for p, c in zip(valid_props, colors):
                ax.plot(agg_z.index, agg_z[p], marker='o', label=p, color=c)

            ax.axhline(0, linestyle="--", linewidth=1)
            ax.set_title(group_name)
            ax.grid(True)
            ax.legend(fontsize=8)

        axes[-1].set_xlabel("Week")
        fig.suptitle("Rolling Property Trends (Z-score)", fontsize=14)

        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(trend_file, dpi=300)
        plt.close()

        # --- deltas (RAW values) ---
        if len(agg) > 1:
            latest = agg.iloc[-1][numeric_cols]
            prev = agg.iloc[-2][numeric_cols]
            deltas = (latest - prev).to_dict()

    # --- 6. Save summary ---
    summary_file = output_dir / "property_drift_summary.csv"
    agg.to_csv(summary_file)
    logger.info(f"[PhysChem Drift] Summary saved to: {summary_file}")

    return {
        "trend_plot": str(trend_file) if trend_file else None,
        "summary_csv": str(summary_file),
        "property_changes": deltas,
        "penalty_csv": str(output_dir / "property_penalty.csv"),
        "trend_slopes": slopes,
        "shift_pvalues": shift_pvals,
        "distribution_shift": dist_shift,
        "directional": directional,
        "alerts": alerts,
    }

@register_task("druglikeness_indices_trend", category="Project-based analyses",
               description="Calculate QED, MPO, CNS MPO trends and alert on drops in drug-likeness.")
def druglikeness_indices_trend(config, data=None):
    """
    Calculate and monitor drug-likeness indices (QED, MPO, CNS MPO) over time.
    """
    # 1. Load config
    input_file = config.get("input_file")
    dl_cfg = config.get("druglikeness_indices_trend", {})

    date_col = dl_cfg.get("date_col", "publication_date")
    output_dir = Path(dl_cfg.get("output_dir", "outputs/druglikeness_trends"))
    rolling_window = dl_cfg.get("rolling_window", 4)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.debug(f"[Drug-likeness Trend] Loading input file: {input_file}")
    df = pd.read_csv(input_file)
    df = df.dropna(subset=["smiles", date_col])
    df[date_col] = pd.to_datetime(df[date_col])

    #  2. Calculate QED and MPO-like scores 
    logger.debug("[Drug-likeness Trend] Calculating indices (QED, MPO, CNS MPO)...")
    mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
    df["QED"] = [QED.qed(m) if m else np.nan for m in mols]

    # Basic properties (required for MPO-like score) - will need to adjust this when we have real data
    df["MW"] = [Descriptors.MolWt(m) if m else np.nan for m in mols]
    df["clogP"] = [Crippen.MolLogP(m) if m else np.nan for m in mols]
    df["TPSA"] = [rdMolDescriptors.CalcTPSA(m) if m else np.nan for m in mols]

    # Simplified MPO-like metrics
    df["MPO"] = 0.5 * ((6 - np.clip(df["clogP"], 0, 6)) / 6 + (140 - np.clip(df["TPSA"], 0, 140)) / 140)
    df["CNS_MPO"] = df["MPO"] - 0.2 * (df["MW"] > 450)  # toy CNS MPO estimate

    # 3. Weekly aggregation
    df["week"] = df[date_col].dt.to_period("W").apply(lambda r: r.start_time)
    agg = df.groupby("week")[["QED", "MPO", "CNS_MPO"]].median().rolling(rolling_window, min_periods=1).mean()

    # 4. Visualisation
    plt.figure(figsize=(10, 6))
    for metric in ["QED", "MPO", "CNS_MPO"]:
        plt.plot(agg.index, agg[metric], marker='o', label=metric)
    plt.legend()
    plt.title("Drug-likeness Index Trends (Rolling Medians)")
    plt.xlabel("Week")
    plt.ylabel("Score (0–1)")
    plt.grid(True)
    trend_file = output_dir / "druglikeness_trends.png"
    plt.savefig(trend_file, dpi=300, bbox_inches="tight")
    plt.close()

    if len(agg) > 1:
        latest = agg.iloc[-1]
        if latest["CNS_MPO"] < 4.0:
            logger.warning("[Drug-likeness Trend] CNS MPO score dropped below 4.0.")
        if latest["MPO"] > 0.8:
            logger.info("[Drug-likeness Trend] > 80% of designs within oral drug-like region.")

    # 6. Save results 
    output_csv = output_dir / "druglikeness_indices_summary.csv"
    agg.to_csv(output_csv)
    logger.info(f"[Drug-likeness Trend] Summary saved to: {output_csv}")

    return {
        "trend_plot": str(trend_file),
        "summary_csv": str(output_csv)
    }

@register_task("rgroup_frequency_tracking",
               category="Project-based analyses",
               description="Track introduction and activity impact of new R-groups week-to-week.")
def rgroup_frequency_tracking(config, data=None):
    cfg = config.get("rgroup_frequency_tracking", {})
    input_file = config.get("input_file")
    activity_col = cfg.get("activity_col", "pActivity")
    date_col = cfg.get("date_col", "publication_date")
    top_n = cfg.get("top_n", 10)
    window_weeks = cfg.get("window_weeks", 1)
    output_dir = Path(cfg.get("output_dir", "outputs/mmp/rgroup_frequency"))
    output_dir.mkdir(parents=True, exist_ok=True)

    # Better way to handle flexible inputs ---
    if data is None:
        df = pd.read_csv(input_file)
    elif isinstance(data, pd.DataFrame):
        df = data
    elif isinstance(data, dict):
        # If dict contains a path to a CSV output from previous step, then...
        possible_path = data.get("filtered_csv") or data.get("output_csv") or data.get("summary_csv")
        if possible_path and Path(possible_path).exists():
            df = pd.read_csv(possible_path)
        else:
            # fallback: reload original filtered file
            logger.warning("[RGroup Tracking] Got dict input; loading from input_file instead.")
            df = pd.read_csv(input_file)
    else:
        raise TypeError(f"Unexpected data type: {type(data)}")
    
    df = df.dropna(subset=["smiles", activity_col, date_col])
    df[date_col] = pd.to_datetime(df[date_col])
    df["week"] = df[date_col].dt.to_period("W").apply(lambda r: r.start_time)

    logger.debug("[RGroup Tracking] Performing R-group decomposition...")
    mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
    scaffolds = [Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol(m) for m in mols]
    
    core_counts = Counter(Chem.MolToSmiles(s) for s in scaffolds if s)

    unique_cores = [
        core for core, count in core_counts.items()
        if count >= 5   # configurable threshold
    ]

    rgroups_all = []
    for core_smiles in unique_cores:
        core = Chem.MolFromSmiles(core_smiles)
        matches = [(i, m) for i, m in enumerate(mols) if m and core and m.HasSubstructMatch(core)]
        
        if len(matches) < 3:
            continue
        
        # Split before calling RDKit
        match_indices = [i for i, _ in matches]
        match_mols = [m for _, m in matches]

        groups, _ = rdRGroupDecomposition.RGroupDecompose([core], match_mols, asRows=True)
        
        if not groups:
            continue

        for mol_idx, row in enumerate(groups):
            original_idx = match_indices[mol_idx]  # Map back to df

            for label, frag in row.items():
                if label == "Core" or frag is None:
                    continue
                if frag:
                    rgroups_all.append({
                        "core": core_smiles,
                        "series": df.iloc[original_idx]["Chemical series"],
                        "rgroup_label": label,
                        "rgroup_smiles": Chem.MolToSmiles(frag),
                        activity_col: df.iloc[original_idx][activity_col],
                        "week": df.iloc[original_idx]["week"]             
                    })
    rgroups_df = pd.DataFrame(rgroups_all)
    if rgroups_df.empty:
        logger.warning("[RGroup Tracking] No R-groups extracted.")
        return None

    grouped = (
        rgroups_df.groupby(["series", "core", "week", "rgroup_smiles"])
        .agg(
            count=(activity_col, "size"),
            mean_activity=(activity_col, "mean")
        )
        .reset_index()
    )

    # Find top N new R-groups introduced by week
    grouped = grouped.sort_values(["week", "count"], ascending=[True, False])
    all_weeks = grouped["week"].sort_values().unique()
    if len(all_weeks) < 2:
        logger.warning("[RGroup Tracking] Not enough weeks for comparison.")
        return None

    latest, prev = all_weeks[-1], all_weeks[-2]
    current_week = grouped[grouped["week"] == latest]
    previous_week = grouped[grouped["week"] == prev]

    new_rgroups = set(current_week["rgroup_smiles"]) - set(previous_week["rgroup_smiles"])
    new_df = current_week[current_week["rgroup_smiles"].isin(new_rgroups)]

    if new_df.empty:
        logger.warning("[RGroup Tracking] No new R-groups found; falling back to most frequent.")
        top_new = current_week.nlargest(top_n, "count")
    else:
        top_new = new_df.nlargest(top_n, "count")

    summary_csv = output_dir / f"top_rgroups_{latest.date()}.csv"
    top_new.to_csv(summary_csv, index=False)
    
    plt.figure(figsize=(8,5))
    plt.barh(top_new["rgroup_smiles"], top_new["mean_activity"], color="teal")
    plt.xlabel("Average pActivity")
    plt.ylabel("R-group (SMILES)")
    plt.title(f"Top {top_n} New R-groups Introduced - {latest.date()}")
    plt.tight_layout()
    plt.savefig(output_dir / f"rgroup_top_{latest.date()}.png", dpi=300)
    plt.close()

    logger.info(f"[RGroup Tracking] Saved summary: {summary_csv}")
    return {"summary_csv": str(summary_csv)}

def draw_core_with_large_atom_numbers(mol, size=(300, 300), font_size=40):
    """
    Draw the core molecule with large atom numbers for attachment points.
    Returns a PIL Image.
    """
    if mol is None:
        return None

    # Create 2D coords if not present
    if not mol.GetNumConformers():
        rdMolDraw2D.PrepareMolForDrawing(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.baseFontSize = font_size / 100

    # Label dummy atoms (atomic number == 0) or atoms with R-group-like mapping
    atom_labels = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            # try atom map numbers, else index+1
            amap = atom.GetAtomMapNum()
            label = f"R{amap}" if amap else str(atom.GetIdx() + 1)
            atom_labels[atom.GetIdx()] = label

    for idx, label in atom_labels.items():
        opts.atomLabels[idx] = label

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()

    png_data = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_data))

@register_task("rgroup_sar_tree",
               category="Project-based analyses",
               description="Build hierarchical SAR trees (core → R-group → mean potency).")
def rgroup_sar_tree(config, data=None):
    cfg = config.get("rgroup_sar_tree", {})
    input_file = config.get("input_file")
    output_dir = Path(cfg.get("output_dir", "outputs/mmp/rgroup_sar_tree"))
    output_dir.mkdir(parents=True, exist_ok=True)

    activity_col = cfg.get("activity_col", "pActivity")
    date_col = cfg.get("date_col", "publication_date")
    core_method = cfg.get("core_method", "murcko")
    min_variants = cfg.get("min_variants", 3)

    # Better handling of multiple input types
    if data is None:
        df = pd.read_csv(input_file)
    elif isinstance(data, pd.DataFrame):
        df = data
    elif isinstance(data, dict):
        possible_paths = [
            data.get("filtered_csv"),
            data.get("output_csv"),
            data.get("summary_csv"),
        ]
        loaded = False
        for path in possible_paths:
            if path and Path(path).exists():
                df_try = pd.read_csv(path)
                if {"smiles", activity_col}.issubset(df_try.columns):
                    df = df_try
                    loaded = True
                    break
        if not loaded:
            logger.warning("[RGroup SAR Tree] Input dict lacks molecule data; reloading from input_file.")
            df = pd.read_csv(input_file)
    else:
        raise TypeError(f"Unexpected data type: {type(data)}")

    # Quick sanity check
    if not {"smiles", activity_col}.issubset(df.columns):
        raise KeyError(f"Input dataset for SAR tree must contain ['smiles', '{activity_col}'], found: {df.columns.tolist()}")

    df = df.dropna(subset=["smiles", activity_col])
    mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]

    # Get core scaffolds
    logger.info(f"[RGroup SAR Tree] Using core method: {core_method}")
    if core_method == "murcko":
        cores = [Scaffolds.MurckoScaffold.GetScaffoldForMol(m) for m in mols]
    else:
        raise NotImplementedError(f"Core method {core_method} not implemented yet.")
    df["core"] = [Chem.MolToSmiles(c) if c else None for c in cores]

    # Group by core and build decomp
    all_records = []
    for core_smiles, subset in df.groupby("core"):
        if not core_smiles or len(subset) < min_variants:
            continue
        core = Chem.MolFromSmiles(core_smiles)
        match_indices = subset.index.tolist()
        match_mols = [Chem.MolFromSmiles(s) for s in subset["smiles"]]

        groups, _ = rdRGroupDecomposition.RGroupDecompose([core], match_mols)
        for row, (_, subrow) in zip(groups, subset.iterrows()):
            for label, frag in row.items():
                if frag:
                    all_records.append({
                        "core": core_smiles,
                        "series": subrow["Chemical series"],
                        "rgroup_label": label,
                        "rgroup_smiles": Chem.MolToSmiles(frag),
                        activity_col: subrow[activity_col],
                    })

    rg_df = pd.DataFrame(all_records)
    if rg_df.empty:
        logger.warning("[RGroup SAR Tree] No decompositions generated.")
        return None

    sar_summary = (
        rg_df.groupby(["series", "core", "rgroup_label", "rgroup_smiles"])
        .agg(mean_pActivity=(activity_col, "mean"),
             std_pActivity=(activity_col, "std"),
             count=(activity_col, "size"))
        .reset_index()
    )

    output_csv = output_dir / "rgroup_sar_tree.csv"
    sar_summary.to_csv(output_csv, index=False)
    logger.info(f"[RGroup SAR Tree] SAR table saved: {output_csv}")

    # Outputting a hierarchical json for later use (interactive html or something else silly)
    tree_json = {}
    for (series, core), subdf in sar_summary.groupby(["series", "core"]):
        tree_json.setdefault(series, {})[core] = (
            subdf[["rgroup_label", "rgroup_smiles", "mean_pActivity", "std_pActivity"]]
            .sort_values("mean_pActivity", ascending=False)
            .to_dict(orient="records")
        )

    json_path = output_dir / "rgroup_sar_tree.json"
    with open(json_path, "w") as fout:
        json.dump(tree_json, fout, indent=2)
    
    # Begin per-core entered SAR viz
    logger.debug("[RGroup SAR Tree] Generating centered per-core SAR visualizations (per-core, with true core fragments)...")
    try:

        cmap = cm.get_cmap("viridis")

        for series, cores in tree_json.items():

            # create per-series directory
            safe_series = str(series).replace(" ", "_").replace("/", "_")
            series_dir = output_dir / safe_series
            series_dir.mkdir(exist_ok=True)

            for core, rgroups in cores.items():
                if not rgroups:
                    continue

                # Try to find the true decomposition core
                core_frag = None
                for rg in rgroups:
                    if str(rg["rgroup_label"]).lower() == "core":
                        core_frag = rg["rgroup_smiles"]
                        break
                if not core_frag:
                    core_frag = core  # fallback

                # Build graph
                G = nx.DiGraph()
                core_label = "CORE"
                G.add_node(core_label, kind="core", smiles=core_frag, pAct=None)

                for rg in rgroups:
                    if str(rg["rgroup_label"]).lower() == "core":
                        continue

                    r_smiles = rg["rgroup_smiles"]
                    pAct = rg["mean_pActivity"]
                    label = rg["rgroup_label"]

                    # ✅ FIX: unique node IDs (prevents overwriting)
                    node_id = f"{label}_{hash(r_smiles)}"

                    G.add_node(node_id, kind="rgroup", smiles=r_smiles, pAct=pAct, label=label)
                    G.add_edge(core_label, node_id)

                # Layout (adaptive)
                n = len(G.nodes)

                # dynamic radius (key fix)
                radius = max(1.5, 0.4 * n)

                # dynamic image scaling
                img_size = max(80, int(140 - n * 2))      # shrink with more nodes
                zoom = max(0.25, 0.5 - n * 0.01)

                # reduce label clutter if many nodes
                show_labels = n <= 12

                pos = {core_label: np.array([0.0, 0.0])}
                angle_step = 2 * np.pi / max(1, n - 1)

                i = 0
                for node in G.nodes:
                    if node == core_label:
                        continue
                    angle = i * angle_step
                    pos[node] = np.array([
                        radius * np.cos(angle),
                        radius * np.sin(angle)
                    ])
                    i += 1

                # Normalize activity
                pActs = [d.get("pAct") for _, d in G.nodes(data=True) if d["pAct"] is not None]
                if not pActs:
                    logger.warning(f"[RGroup SAR Tree] No activity values for core: {core}")
                    continue

                norm = plt.Normalize(min(pActs), max(pActs))

                fig, ax = plt.subplots(figsize=(7, 7))
                ax.axis("off")

                # Edges
                nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#666666", width=1.2, arrows=False)

                # Core
                core_mol = Chem.MolFromSmiles(core_frag)
                if core_mol is None:
                    core_mol = Chem.MolFromSmiles(core)

                if core_mol:
                    core_img = draw_core_with_large_atom_numbers(core_mol, size=(300, 300), font_size=50)
                    if core_img:
                        im = OffsetImage(core_img, zoom=0.5)
                        ab = AnnotationBbox(im, (0, 0), frameon=True,
                                            bboxprops=dict(facecolor="#ffcc00", edgecolor="#333333", lw=1.5))
                        ax.add_artist(ab)

                # R-groups
                for node, data in G.nodes(data=True):
                    if data["kind"] != "rgroup":
                        continue

                    x, y = pos[node]
                    mol = Chem.MolFromSmiles(data["smiles"])

                    if mol is None:
                        continue

                    img = Draw.MolToImage(mol, size=(120, 120))
                    im = OffsetImage(img, zoom=0.45)
                    ab = AnnotationBbox(im, (x, y), frameon=False)
                    ax.add_artist(ab)
                    ax.set_title(f"{series} | Core SAR", fontsize=10)

                    # label (use stored label, not node_id)
                    ax.text(x, y + 0.25, data["label"],
                            ha="center", va="bottom",
                            fontsize=8, fontweight="bold", color="#333333")

                    # activity
                    if data["pAct"] is not None:
                        ax.text(x, y - 0.25, f"{data['pAct']:.2f}",
                                ha="center", va="top", fontsize=8,
                                color=cm.viridis(norm(data["pAct"])))

                # Colorbar
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = fig.colorbar(sm, ax=ax, shrink=0.75, pad=0.03)
                cbar.set_label("Mean pActivity")

                # Save
                core_hash = hashlib.sha1(core_frag.encode()).hexdigest()[:8]
                plot_path = series_dir / f"rgroup_tree_{core_hash}.png"
                fig.savefig(plot_path, dpi=300, bbox_inches="tight")
                plt.close(fig)

                logger.info(f"[RGroup SAR Tree] Saved: {plot_path}")

    except Exception as e:
        import traceback
        traceback.print_exc()
        logger.warning(f"[RGroup SAR Tree] Plot generation failed: {e}")

    return {"summary_csv": str(output_csv), "hierarchy_json": str(json_path)}