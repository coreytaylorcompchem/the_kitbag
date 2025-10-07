import csv
import logging
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np
from PIL import Image
import seaborn as sns
import networkx as nx
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from mpl_toolkits.axes_grid1 import ImageGrid

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
        # mmpdb is painful to use hereso we need to call the program on each SMILES and aggregate.
        return property_name, smiles, output_lines
    except subprocess.CalledProcessError as e:
        logger.warning(f"[!] Transform failed for SMILES {smiles} with property {property_name}")
        logger.debug(e.stderr)
        return property_name, smiles, None

@register_task("mmp_analysis", category="Project-based analyses", description="Matched Molecular Pair (MMP) analysis.")
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
    props_cols = ["property", "mw"]
    # props_cols = ["property", "mw", "logp", "hbd", "hba", "rotatable_bonds", "tpsa", "qed", "stereocenters"]
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

    # --- 4. Network plot with images ---
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



