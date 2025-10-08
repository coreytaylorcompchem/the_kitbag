import csv
import logging
import subprocess
from pathlib import Path
from tqdm import tqdm

import matplotlib
matplotlib.use('Agg')
import matplotlib.offsetbox as offsetbox
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np
from PIL import Image
import seaborn as sns
import networkx as nx
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from mpl_toolkits.axes_grid1 import ImageGrid

from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import umap

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
    logger.info(f"[✔] SMILES file written: {smi_path}")

    prop_path = None
    if props_cols:
        missing_cols = [col for col in props_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing property columns in DataFrame: {missing_cols}")

        props_df = df[[id_col] + props_cols].copy()
        props_df.rename(columns={id_col: "id"}, inplace=True)
        prop_path = output_prefix.with_suffix(".props")
        props_df.to_csv(prop_path, index=False, sep="\t", encoding="utf-8", lineterminator="\n")
        logger.info(f"[✔] Properties file written: {prop_path}")
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
        logger.warning(f"[!] Transform failed for SMILES {smiles} with property {property_name}")
        logger.debug(e.stderr)
        return property_name, smiles, None

@register_task("mmp_analysis", category="Project-based analyses", description="Matched Molecular Pairs (mmpdb).")
def mmp_analysis(config, data=None):
    input_file = config.get("input_file")
    activity_col = config.get("activity_col", "pActivity")
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
    # props_cols = ["property", "mw"]
    props_cols = ["property", "mw", "logp", "hbd", "hba", "rotatable_bonds", "tpsa", "qed", "stereocenters"]
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

    # Run mmpdb fragment/index
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

############### Plotting helper

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

    # --- 1. Setup ---
    input_file = config.get("input_file")
    output_dir = Path(config.get("output", {}).get("directory", "outputs/mmp_report"))
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_file, sep='\t')

    min_count = config.get("min_count", 10)
    significance_level = config.get("significance_level", 0.05)

    df_filtered = df[df['property_count'] >= min_count].copy()
    df_filtered['transform_label'] = df_filtered['property_from_smiles'] + " → " + df_filtered['property_to_smiles']

    # --- 2. Structure image caching ---
    mol_images_from = {smi: mol_to_img(smi) for smi in pd.unique(df_filtered['property_from_smiles'])}
    mol_images_to = {smi: mol_to_img(smi) for smi in pd.unique(df_filtered['property_to_smiles'])}

    # --- 3. Multi-Property Heatmap from Long Format ---

    # 1. Create transform label
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

    # 6. Plotting - dynamically resize
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

    # Add integer index labels to x-axis
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

    # 9. Save
    plt.suptitle(f"Multi-Property Change per Transform (Count ≥ {min_count})", fontsize=14)
    heatmap_file = output_dir / "mmp_transform_multi_property_heatmap.png"
    plt.savefig(heatmap_file, bbox_inches='tight', dpi=300)
    plt.close()

    # --- 4. Network plot with images --- (needs work)
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

    # --- 5. Normalised histograms per property --- (needs work)
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

    # --- 6. Top transform tables ---
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
    i, smi_i, fp_i, data, cutoff, delta_threshold, activity_col = args
    cliffs = []
    for j, (smi_j, fp_j, act_j) in enumerate(zip(data["smiles"], data["fps"], data[activity_col])):
        if i >= j:
            continue
        if fp_i is None or fp_j is None:
            continue
        sim = DataStructs.TanimotoSimilarity(fp_i, fp_j)
        if sim >= cutoff:
            delta = abs(data.loc[i, activity_col] - act_j)
            if delta >= delta_threshold:
                cliffs.append({
                    "smiles_1": smi_i,
                    "smiles_2": smi_j,
                    "similarity": sim,
                    "delta_activity": delta,
                    "activity_1": data.loc[i, activity_col],
                    "activity_2": act_j
                })
    return cliffs

@register_task("sar_cliff_analysis", category="Project-based analyses",
               description="Identify and visualise activity cliffs.")
def sar_cliff_analysis(config, data=None):

    # --- 1. Configuration ---
    input_file = config.get("input_file")

    sar_cliff_cfg = config.get("sar_cliff", {})    
    activity_col = sar_cliff_cfg.get("activity_col", "pActivity")
    similarity_cutoff = sar_cliff_cfg.get("similarity_cutoff", 0.85)
    delta_threshold = sar_cliff_cfg.get("delta_threshold", 1.0)
    fp_type = sar_cliff_cfg.get("fp_type", "morgan")
    fp_radius = sar_cliff_cfg.get("fp_radius", 2)
    fp_nbits = sar_cliff_cfg.get("fp_nbits", 2048)
    max_pairs = sar_cliff_cfg.get("max_pairs", 500000)
    output_dir = Path(sar_cliff_cfg.get("output_dir", "outputs/sar_cliffs"))
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"[SAR Cliff] Loading input file: {input_file}")
    df = pd.read_csv(input_file)
    df = df.dropna(subset=["smiles", activity_col])
    df = df.reset_index(drop=True)

    # --- 2. Fingerprint computation ---
    logger.info("[SAR Cliff] Computing fingerprints...")
    df["fps"] = [compute_fingerprint(smi, fp_type, fp_radius, fp_nbits) for smi in df["smiles"]]
    df = df[df["fps"].notna()]

    # --- 3. Pairwise similarity in parallel ---
    logger.info(f"[SAR Cliff] Computing pairwise similarities (cutoff={similarity_cutoff})...")

    tasks = [
        (i, row["smiles"], row["fps"], df, similarity_cutoff, delta_threshold, activity_col)
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

    # --- 4. Save raw results ---
    output_csv = output_dir / "sar_cliff_pairs.csv"
    cliff_df.to_csv(output_csv, index=False)
    logger.info(f"[SAR Cliff] Results saved to: {output_csv}")

    # --- 5. Visualization ---
    if not cliff_df.empty:
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            data=cliff_df,
            x="similarity",
            y="delta_activity",
            s=40,
            alpha=0.7
        )
        plt.axvline(similarity_cutoff, color='red', linestyle='--', label=f"sim ≥ {similarity_cutoff}")
        plt.axhline(delta_threshold, color='orange', linestyle='--', label=f"Δact ≥ {delta_threshold}")
        plt.xlabel("Tanimoto Similarity")
        plt.ylabel(f"|Δ {activity_col}|")
        plt.title("Activity Cliff Landscape")
        plt.legend()
        scatter_file = output_dir / "sar_cliff_scatter.png"
        plt.savefig(scatter_file, bbox_inches='tight', dpi=300)
        plt.close()

        # Top 10 cliffs — with molecule images
        top_cliffs = cliff_df.head(10)
        imgs = []
        for _, row in top_cliffs.iterrows():
            mol1 = Chem.MolFromSmiles(row["smiles_1"])
            mol2 = Chem.MolFromSmiles(row["smiles_2"])
            img = Draw.MolsToGridImage(
                [mol1, mol2],
                legends=[
                    f"Act={row['activity_1']:.2f}",
                    f"Act={row['activity_2']:.2f}"
                ],
                molsPerRow=2,
                subImgSize=(200, 200)
            )
            imgs.append(img)

        for i, img in enumerate(imgs):
            img_path = output_dir / f"top_cliff_{i+1}.png"
            img.save(img_path)

        # Summary CSV
        top_csv = output_dir / "top_10_sar_cliffs.csv"
        top_cliffs.to_csv(top_csv, index=False)
        logger.info(f"[SAR Cliff] Top 10 cliffs written to: {top_csv}")

        return {
            "cliff_csv": str(output_csv),
            "scatter_png": str(scatter_file),
            "top_cliffs_csv": str(top_csv),
        }

    else:
        logger.warning("[SAR Cliff] No activity cliffs found.")
        return None

try:
    import umap
except ImportError:
    umap = None
    warnings.warn("UMAP not installed. Falling back to t-SNE for dimensionality reduction.")


@register_task(
    "chemical_space_drift",
    category="Project-based analyses",
    description="Track week-to-week evolution of chemical space via (UMAP / t-SNE)."
)
def chemical_space_drift(config, data=None):
    """Visualize chemical space evolution colored by time or activity."""

    # --- 1. Configuration ---
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

    logger.info("[ChemicalSpaceDrift] Starting task...")

    # --- 2. Load input ---
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

    # Fallback to CSV input
    if df is None and input_file and Path(input_file).exists():
        df = pd.read_csv(input_file)

    if df is None or not isinstance(df, pd.DataFrame):
        raise TypeError(f"[ChemicalSpaceDrift] Expected DataFrame input, got {type(df)}")

    logger.info(f"[ChemicalSpaceDrift] Data loaded with {len(df)} rows and {len(df.columns)} columns.")

    # --- 3. Locate SMILES column ---
    smiles_col = next((c for c in ["smiles", "SMILES", "canonical_smiles"] if c in df.columns), None)
    if smiles_col is None:
        raise KeyError(f"[ChemicalSpaceDrift] No 'smiles' column found. Available: {list(df.columns)}")

    # --- 4. Generate fingerprints ---
    logger.info(f"[ChemicalSpaceDrift] Generating Morgan fingerprints (radius={fp_radius}, nBits={fp_nbits})...")
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

    # --- 5. Dimensionality reduction ---
    X = StandardScaler().fit_transform(fps)
    if reducer_type == "tsne" or umap is None:
        reducer = TSNE(n_components=2, random_state=42, perplexity=30)
        method = "t-SNE"
    else:
        reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
        method = "UMAP"

    logger.info(f"[ChemicalSpaceDrift] Performing {method} reduction...")
    embedding = reducer.fit_transform(X)
    df["x"], df["y"] = embedding[:, 0], embedding[:, 1]

    # --- 6. Handle time grouping ---
    if date_col in df.columns:
        df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
        df["week"] = df[date_col].dt.strftime("%Y-%U")

        # Optional trimming to most recent N weeks
        recent_weeks = task_conf.get("recent_weeks", 12)
        unique_weeks = sorted(df["week"].dropna().unique())
        if len(unique_weeks) > recent_weeks:
            recent_weeks_set = set(unique_weeks[-recent_weeks:])
            df = df[df["week"].isin(recent_weeks_set)]
            logger.info(f"[ChemicalSpaceDrift] Filtering to last {recent_weeks} weeks "
                        f"({unique_weeks[-recent_weeks]} → {unique_weeks[-1]}).")
    else:
        df["week"] = "unknown"

    # --- 7. Visualization ---
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

    # --- 8. Export data ---
    emb_csv = output_dir / "chemical_space_embedding.csv"
    df.to_csv(emb_csv, index=False)
    logger.info(f"[ChemicalSpaceDrift] Results saved: {emb_csv}")

    logger.info(f"[ChemicalSpaceDrift] Plot saved: {plot_path}")

    # --- 9. Return dictionary of artifacts (like sar_cliff_analysis) ---
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

    # --- 1. Config ---
    task_conf = config.get("scaffold_enrichment_trends", {})
    input_file = config.get("input_file")
    output_dir = Path(task_conf.get("output_dir", "outputs/scaffold_enrichment"))
    output_dir.mkdir(parents=True, exist_ok=True)

    activity_col = task_conf.get("activity_col", "pActivity")
    date_col = task_conf.get("date_col", "date")

    logger.info("[ScaffoldEnrichment] Starting task...")

    # --- 2. Load input ---
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

    # --- 3. Locate SMILES column ---
    smiles_col = next((c for c in ["smiles", "SMILES", "canonical_smiles"] if c in df.columns), None)
    if smiles_col is None:
        raise KeyError(f"[ScaffoldEnrichment] No 'smiles' column found. Available: {list(df.columns)}")

    if date_col not in df.columns:
        logger.warning(f"[ScaffoldEnrichment] No '{date_col}' column found — using 'week'='unknown'.")
        df["week"] = "unknown"
    else:
        df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
        df["week"] = df[date_col].dt.strftime("%Y-%U")

    # --- 4. Compute Bemis–Murcko scaffolds ---
    logger.info("[ScaffoldEnrichment] Computing Murcko scaffolds...")
    scaffolds = []
    for smi in tqdm(df[smiles_col], desc="Extracting scaffolds"):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            scaffolds.append(None)
            continue
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            scaffolds.append(scaf)
        except Exception:
            scaffolds.append(None)
    df["scaffold"] = scaffolds

    df_valid = df.dropna(subset=["scaffold"]).copy()
    if df_valid.empty:
        logger.warning("[ScaffoldEnrichment] No valid scaffolds found.")
        return None

    # --- 5. Weekly scaffold stats ---
    logger.info("[ScaffoldEnrichment] Aggregating weekly scaffold statistics...")
    agg_df = (
        df_valid.groupby(["week", "scaffold"])
        .agg(
            count=("scaffold", "size"),
            mean_activity=(activity_col, "mean")
        )
        .reset_index()
    )

    # Identify change vs previous week
    agg_df["prev_mean"] = agg_df.groupby("scaffold")["mean_activity"].shift(1)
    agg_df["delta_mean"] = agg_df["mean_activity"] - agg_df["prev_mean"]

    # --- 6. Identify improving / declining scaffolds ---
    latest_week = agg_df["week"].dropna().unique()
    latest_week = sorted(latest_week)[-1] if len(latest_week) else None

    if latest_week:
        latest_data = agg_df[agg_df["week"] == latest_week]
        top_improving = latest_data.nlargest(5, "delta_mean", keep="all")
        top_declining = latest_data.nsmallest(5, "delta_mean", keep="all")
    else:
        top_improving = top_declining = pd.DataFrame()

    # --- 7. Save outputs ---
    trends_csv = output_dir / "scaffold_trends.csv"
    agg_df.to_csv(trends_csv, index=False)

    improving_csv = output_dir / "top5_improving_scaffolds.csv"
    declining_csv = output_dir / "top5_declining_scaffolds.csv"
    top_improving.to_csv(improving_csv, index=False)
    top_declining.to_csv(declining_csv, index=False)

    # --- 8. Visualization ---
    image_zoom = 0.4  # hardcoding size of 2D structures for the plots 

    top_n = task_conf.get("top_n_scaffolds", 8)
    top_scaffolds = (
        agg_df.groupby("scaffold")["count"].sum()
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )
    plot_df = agg_df[agg_df["scaffold"].isin(top_scaffolds)]

    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    # Prepare scaffold images
    scaffold_imgs = {}
    for smi in top_scaffolds:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            img = Draw.MolToImage(mol, size=(100, 100))
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

    # Split scaffolds into multi-point and single-point sets
    multi_point_scaffolds = []
    single_point_scaffolds = []
    for scaffold in top_scaffolds:
        count_points = len(plot_df[plot_df["scaffold"] == scaffold])
        if count_points > 1:
            multi_point_scaffolds.append(scaffold)
        else:
            single_point_scaffolds.append(scaffold)

    # Plot lines for scaffolds with multiple points
    if multi_point_scaffolds:
        sns.lineplot(
            data=plot_df[plot_df["scaffold"].isin(multi_point_scaffolds)],
            x="week",
            y="mean_activity",
            hue="scaffold",
            palette="tab10",
            lw=2,
            ax=ax,
        )

    # Remove default legend if present
    if ax.legend_:
        ax.legend_.remove()

    # --- Add images for line plots at the far right of the plot ---
    x_left, x_right = ax.get_xlim()
    x_offset = (x_right - x_left) * 0.1  # move 10% beyond plot limit
    y_min, y_max = ax.get_ylim()
    y_spacing = (y_max - y_min) / (len(multi_point_scaffolds) + 1)

    for i, scaffold in enumerate(multi_point_scaffolds):
        scaffold_data = plot_df[plot_df["scaffold"] == scaffold]
        if scaffold_data.empty:
            continue
        # new position for the scaffold image (further right)
        y_pos = y_max - (i + 1) * y_spacing
        x_img = x_right + x_offset

        # place scaffold image
        img = scaffold_imgs.get(scaffold)
        if img:
            imagebox = offsetbox.OffsetImage(img, zoom=image_zoom)
            ab = offsetbox.AnnotationBbox(
                imagebox, (x_img, y_pos),
                frameon=False, box_alignment=(0, 0.5)
            )
            ax.add_artist(ab)

    # --- Plot single data points with pastel scatter colors ---
    if single_point_scaffolds:
        pastel_palette = sns.color_palette("Pastel1", n_colors=len(single_point_scaffolds))
        for color, scaffold in zip(pastel_palette, single_point_scaffolds):
            scaffold_data = plot_df[plot_df["scaffold"] == scaffold]
            if scaffold_data.empty:
                continue

            x_raw = scaffold_data["week"].iloc[0]
            y = scaffold_data["mean_activity"].iloc[0]
            x = get_x_position(x_raw)

            # Plot single point
            ax.scatter(x, y, s=60, color=color, edgecolor='black', zorder=5)

            # Add scaffold image slightly above the point
            img = scaffold_imgs.get(scaffold)
            if img:
                imagebox = offsetbox.OffsetImage(img, zoom=image_zoom)
                ab = offsetbox.AnnotationBbox(
                    imagebox, (x, y + 0.02),
                    frameon=False, box_alignment=(0.5, 0)
                )
                ax.add_artist(ab)

    # Expand limits to fit right-hand images fully
    ax.set_xlim(x_left, x_right + x_offset * 2)

    plt.title(f"Top {top_n} Scaffold Enrichment Trends (mean {activity_col})")
    plt.xticks(rotation=45)
    plt.xlabel("Week")
    plt.ylabel(f"Mean {activity_col}")
    plt.tight_layout()

    trend_plot = output_dir / "scaffold_trends_plot.png"
    plt.savefig(trend_plot, dpi=300)
    plt.close()


    logger.info(f"[ScaffoldEnrichment] Trends saved to: {trends_csv}")
    logger.info(f"[ScaffoldEnrichment] Plot saved to: {trend_plot}")

    # --- 9. Return dictionary of outputs ---
    return {
        "trends_csv": str(trends_csv),
        "top_improving_csv": str(improving_csv),
        "top_declining_csv": str(declining_csv),
        "trend_plot": str(trend_plot),
        "n_valid_scaffolds": len(df_valid)
    }