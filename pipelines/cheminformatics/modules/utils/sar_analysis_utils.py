import re
from pathlib import Path
import hashlib
import io

import numpy as np
from PIL import Image
import seaborn as sns

import pandas as pd


from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

import matplotlib
matplotlib.use('Agg')
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import combinations

from PIL import Image

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def safe_name(x, max_len=80):
    x = str(x)
    x = re.sub(r"[^\w\-.]+", "_", x)
    return x[:max_len].strip("_") or "unknown"

def smiles_to_pil(smiles, size=(160, 160)):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    try:
        return Draw.MolToImage(mol, size=size)
    except Exception:
        return None

def plot_ranked_rgroup_grid(
    df,
    output_path,
    activity_col="mean_pActivity",
    std_col="std_pActivity",
    count_col="count",
    rgroup_label="R",
    title=None,
    img_size=160,
    n_cols=4,
    dpi=300,
    descending=True,
    show_std=True,
    show_count=True,
    color_by_activity=True,
    show_colorbar=True,
):
    """
    Plot R-groups as a ranked structure grid.

    Required columns:
      - rgroup_smiles
      - activity_col

    Optional columns:
      - std_col
      - count_col
    """

    if df is None or df.empty:
        return None

    required_cols = {"rgroup_smiles", activity_col}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(
            f"plot_ranked_rgroup_grid missing required columns: {missing}. "
            f"Available columns: {df.columns.tolist()}"
        )

    df = df.copy()
    df[activity_col] = pd.to_numeric(df[activity_col], errors="coerce")
    df = df.sort_values(activity_col, ascending=not descending).reset_index(drop=True)

    n = len(df)
    n_cols = min(n_cols, max(1, n))
    n_rows = int(np.ceil(n / n_cols))

    fig_w = max(8, n_cols * 2.6)
    fig_h = max(3.5, n_rows * 2.8 + 0.7)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h))
    axes = np.array(axes).reshape(-1)

    vals = df[activity_col].dropna()

    has_activity_range = (
        color_by_activity
        and len(vals) > 1
        and vals.max() > vals.min()
    )

    if has_activity_range:
        norm = plt.Normalize(vals.min(), vals.max())
    else:
        norm = plt.Normalize(0, 1)

    cmap = plt.cm.viridis

    # Prepare all axes, including unused ones
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_facecolor("white")

    for i, row in df.iterrows():
        ax = axes[i]

        smi = row["rgroup_smiles"]
        img = smiles_to_pil(smi, size=(img_size, img_size))

        p = row.get(activity_col, np.nan)

        if has_activity_range and pd.notna(p):
            tile_color = cmap(norm(p))
        else:
            tile_color = "#444444"

        if img is not None:
            ax.imshow(img)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(True)

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(2.5)
            spine.set_edgecolor(tile_color)

        label_lines = [f"{rgroup_label} rank {i + 1}"]

        if pd.notna(p):
            label_lines.append(f"mean={p:.2f}")

        if show_std and std_col in df.columns:
            std = pd.to_numeric(row.get(std_col, np.nan), errors="coerce")
            if pd.notna(std):
                label_lines.append(f"sd={std:.2f}")

        if show_count and count_col in df.columns:
            count = pd.to_numeric(row.get(count_col, np.nan), errors="coerce")
            if pd.notna(count):
                label_lines.append(f"n={int(count)}")

        ax.set_title("\n".join(label_lines), fontsize=8)

    # Leave room at bottom for horizontal colourbar
    if title:
        fig.suptitle(title, fontsize=12)

    fig.tight_layout(rect=[0, 0.10, 1, 0.95])

    if has_activity_range and show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        cbar = fig.colorbar(
            sm,
            ax=axes.tolist(),
            orientation="horizontal",
            fraction=0.035,
            pad=0.04,
            aspect=40,
        )
        cbar.set_label(activity_col)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path

def make_top_bottom_panel(df, top_n=10, bottom_n=10, activity_col="mean_pActivity", descending=True):
    df = df.copy().dropna(subset=[activity_col])

    if df.empty:
        return df

    df = df.sort_values(activity_col, ascending=not descending)

    top = df.head(top_n).copy()
    top["sar_bucket"] = "top"

    bottom = df.tail(bottom_n).copy()
    bottom["sar_bucket"] = "bottom"

    out = pd.concat([top, bottom], ignore_index=True)
    out = out.drop_duplicates(subset=["rgroup_smiles", "rgroup_label"])

    return out

def iter_pages(df, page_size=24):
    for start in range(0, len(df), page_size):
        yield start // page_size + 1, df.iloc[start:start + page_size].copy()

def compute_pairwise_rgroup_combinations(
    wide_df,
    activity_col,
    min_count=3,
):
    r_cols = get_rgroup_columns(wide_df)

    records = []

    for r1, r2 in combinations(r_cols, 2):
        sub = wide_df.dropna(subset=[r1, r2, activity_col]).copy()

        if sub.empty:
            continue

        grouped = (
            sub.groupby([r1, r2])
            .agg(
                mean_activity=(activity_col, "mean"),
                std_activity=(activity_col, "std"),
                count=(activity_col, "size"),
            )
            .reset_index()
        )

        grouped = grouped[grouped["count"] >= min_count]

        for _, row in grouped.iterrows():
            records.append({
                "rgroup_label_1": r1,
                "rgroup_value_1": row[r1],
                "rgroup_label_2": r2,
                "rgroup_value_2": row[r2],
                "mean_activity": row["mean_activity"],
                "std_activity": row["std_activity"],
                "count": row["count"],
            })

    return pd.DataFrame(records)

# def compute_pairwise_rgroup_residuals(
#     wide_df,
#     activity_col,
#     min_count=3,
# ):
#     r_cols = [c for c in wide_df.columns if str(c).startswith("R")]
#     overall_mean = wide_df[activity_col].mean()

#     single_effects = {}

#     for r in r_cols:
#         tmp = (
#             wide_df.dropna(subset=[r, activity_col])
#             .groupby(r)[activity_col]
#             .agg(["mean", "count"])
#             .reset_index()
#         )
#         tmp = tmp[tmp["count"] >= min_count]

#         for _, row in tmp.iterrows():
#             single_effects[(r, row[r])] = row["mean"] - overall_mean

#     records = []

#     for r1, r2 in combinations(r_cols, 2):
#         sub = wide_df.dropna(subset=[r1, r2, activity_col]).copy()

#         grouped = (
#             sub.groupby([r1, r2])[activity_col]
#             .agg(["mean", "std", "count"])
#             .reset_index()
#         )

#         grouped = grouped[grouped["count"] >= min_count]

#         for _, row in grouped.iterrows():
#             v1 = row[r1]
#             v2 = row[r2]

#             e1 = single_effects.get((r1, v1), np.nan)
#             e2 = single_effects.get((r2, v2), np.nan)

#             if pd.isna(e1) or pd.isna(e2):
#                 continue

#             expected = overall_mean + e1 + e2
#             observed = row["mean"]
#             residual = observed - expected

#             records.append({
#                 "rgroup_label_1": r1,
#                 "rgroup_value_1": v1,
#                 "rgroup_label_2": r2,
#                 "rgroup_value_2": v2,
#                 "observed_mean_activity": observed,
#                 "expected_additive_activity": expected,
#                 "residual_activity": residual,
#                 "std_activity": row["std"],
#                 "count": row["count"],
#             })

#     out = pd.DataFrame(records)

#     if not out.empty:
#         out = out.sort_values("residual_activity", ascending=False)

#     return out

def plot_pairwise_combo_heatmap(
    combo_df,
    r1,
    r2,
    value_col,
    output_path,
    title=None,
    max_levels=25,
    dpi=300,
):
    sub = combo_df[
        (combo_df["rgroup_label_1"] == r1) &
        (combo_df["rgroup_label_2"] == r2)
    ].copy()

    if sub.empty:
        return None

    # Keep top populated/value entries to avoid enormous plots
    sub = sub.sort_values("count", ascending=False).head(max_levels * max_levels)

    pivot = sub.pivot_table(
        index="rgroup_value_1",
        columns="rgroup_value_2",
        values=value_col,
        aggfunc="mean"
    )

    if pivot.empty:
        return None

    fig_w = max(8, min(18, pivot.shape[1] * 0.45))
    fig_h = max(6, min(18, pivot.shape[0] * 0.45))

    plt.figure(figsize=(fig_w, fig_h))

    center = 0 if "residual" in value_col else None

    sns.heatmap(
        pivot,
        cmap="coolwarm" if center == 0 else "viridis",
        center=center,
        linewidths=0.3,
        cbar_kws={"label": value_col}
    )

    plt.title(title or f"{r1} × {r2}: {value_col}")
    plt.xlabel(r2)
    plt.ylabel(r1)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

    return output_path

def get_rgroup_columns(wide_df):
    return [
        c for c in wide_df.columns
        if isinstance(c, str)
        and c.startswith("R")
        and c[1:].isdigit()
    ]

def plot_pairwise_combo_structure_heatmap(
    combo_df,
    r1,
    r2,
    value_col,
    output_path,
    count_col="count",
    title=None,
    max_levels=20,
    img_size=120,
    dpi=300,
):
    """
    Plot pairwise R-group combination heatmap with molecule images
    instead of SMILES tick labels.

    Expects combo_df columns:
      rgroup_label_1
      rgroup_value_1
      rgroup_label_2
      rgroup_value_2
      value_col
      count_col
    """

    if combo_df is None or combo_df.empty:
        return None

    sub = combo_df[
        (combo_df["rgroup_label_1"] == r1) &
        (combo_df["rgroup_label_2"] == r2)
    ].copy()

    if sub.empty:
        return None

    sub[value_col] = pd.to_numeric(sub[value_col], errors="coerce")
    sub[count_col] = pd.to_numeric(sub[count_col], errors="coerce")

    sub = sub.dropna(subset=[value_col])
    if sub.empty:
        return None

    # Keep most populated combinations first to avoid huge unreadable heatmaps
    top_rows = (
        sub.groupby("rgroup_value_1")[count_col]
        .sum()
        .sort_values(ascending=False)
        .head(max_levels)
        .index
    )

    top_cols = (
        sub.groupby("rgroup_value_2")[count_col]
        .sum()
        .sort_values(ascending=False)
        .head(max_levels)
        .index
    )

    sub = sub[
        sub["rgroup_value_1"].isin(top_rows) &
        sub["rgroup_value_2"].isin(top_cols)
    ].copy()

    if sub.empty:
        return None

    value_pivot = sub.pivot_table(
        index="rgroup_value_1",
        columns="rgroup_value_2",
        values=value_col,
        aggfunc="mean"
    )

    count_pivot = sub.pivot_table(
        index="rgroup_value_1",
        columns="rgroup_value_2",
        values=count_col,
        aggfunc="sum"
    )

    if value_pivot.empty:
        return None

    row_smiles = list(value_pivot.index)
    col_smiles = list(value_pivot.columns)

    n_rows = len(row_smiles)
    n_cols = len(col_smiles)

    fig_w = max(8, min(22, 2.0 + n_cols * 0.8))
    fig_h = max(7, min(24, 2.0 + n_rows * 0.8))

    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = gridspec.GridSpec(
        nrows=2,
        ncols=2,
        width_ratios=[1.4, max(4, n_cols)],
        height_ratios=[1.4, max(4, n_rows)],
        wspace=0.02,
        hspace=0.02,
    )

    ax_blank = fig.add_subplot(gs[0, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_heat = fig.add_subplot(gs[1, 1])

    ax_blank.axis("off")
    ax_top.axis("off")
    ax_left.axis("off")

    # Main heatmap
    center = 0 if "residual" in value_col or "nonadd" in value_col.lower() else None

    sns.heatmap(
        value_pivot,
        cmap="coolwarm" if center == 0 else "viridis",
        center=center,
        ax=ax_heat,
        cbar_kws={
            "label": value_col,
            "orientation": "horizontal",
            "pad": 0.12,
            "fraction": 0.05,
        },
        linewidths=0.4,
        linecolor="white",
        xticklabels=False,
        yticklabels=False,
    )

    # Annotate cells with value and count
    for i, row_smi in enumerate(row_smiles):
        for j, col_smi in enumerate(col_smiles):
            val = value_pivot.loc[row_smi, col_smi]
            if pd.isna(val):
                continue

            cnt = count_pivot.loc[row_smi, col_smi] if (
                row_smi in count_pivot.index and col_smi in count_pivot.columns
            ) else np.nan

            if pd.notna(cnt):
                txt = f"{val:.2f}\nn={int(cnt)}"
            else:
                txt = f"{val:.2f}"

            ax_heat.text(
                j + 0.5,
                i + 0.5,
                txt,
                ha="center",
                va="center",
                fontsize=7,
                color="black",
            )

    ax_heat.set_xlabel(r2)
    ax_heat.set_ylabel(r1)

    # Top molecule images for column R-groups
    ax_top.set_xlim(0, n_cols)
    ax_top.set_ylim(0, 1)

    for j, smi in enumerate(col_smiles):
        img = smiles_to_pil(smi, size=(img_size, img_size))
        if img is None:
            continue

        imagebox = OffsetImage(img, zoom=0.45)
        ab = AnnotationBbox(
            imagebox,
            (j + 0.5, 0.5),
            frameon=False,
            xycoords="data",
        )
        ax_top.add_artist(ab)

    # Left molecule images for row R-groups
    ax_left.set_xlim(0, 1)
    ax_left.set_ylim(n_rows, 0)

    for i, smi in enumerate(row_smiles):
        img = smiles_to_pil(smi, size=(img_size, img_size))
        if img is None:
            continue

        imagebox = OffsetImage(img, zoom=0.45)
        ab = AnnotationBbox(
            imagebox,
            (0.5, i + 0.5),
            frameon=False,
            xycoords="data",
        )
        ax_left.add_artist(ab)

    if title:
        fig.suptitle(title, fontsize=12)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path

def compute_pairwise_rgroup_residuals(
    wide_df,
    activity_col,
    min_count=2,
):
    """
    Estimate pairwise non-additivity.

    residual = observed mean activity for R_i/R_j combination
               - expected additive activity from individual R effects

    Positive residual:
      combination is better than expected.

    Negative residual:
      combination is worse than expected.
    """

    if wide_df is None or wide_df.empty:
        return pd.DataFrame()

    r_cols = get_rgroup_columns(wide_df)

    if len(r_cols) < 2:
        return pd.DataFrame()

    work = wide_df.copy()
    work[activity_col] = pd.to_numeric(work[activity_col], errors="coerce")
    work = work.dropna(subset=[activity_col])

    if work.empty:
        return pd.DataFrame()

    overall_mean = work[activity_col].mean()

    single_effects = {}

    for r in r_cols:
        tmp = (
            work.dropna(subset=[r, activity_col])
            .groupby(r)[activity_col]
            .agg(["mean", "count"])
            .reset_index()
        )

        tmp = tmp[tmp["count"] >= min_count]

        for _, row in tmp.iterrows():
            single_effects[(r, row[r])] = row["mean"] - overall_mean

    records = []

    for r1, r2 in combinations(r_cols, 2):
        sub = work.dropna(subset=[r1, r2, activity_col]).copy()

        if sub.empty:
            continue

        grouped = (
            sub.groupby([r1, r2])[activity_col]
            .agg(["mean", "std", "count"])
            .reset_index()
        )

        grouped = grouped[grouped["count"] >= min_count]

        for _, row in grouped.iterrows():
            v1 = row[r1]
            v2 = row[r2]

            e1 = single_effects.get((r1, v1), np.nan)
            e2 = single_effects.get((r2, v2), np.nan)

            if pd.isna(e1) or pd.isna(e2):
                continue

            expected = overall_mean + e1 + e2
            observed = row["mean"]
            residual = observed - expected

            records.append({
                "rgroup_label_1": r1,
                "rgroup_value_1": v1,
                "rgroup_label_2": r2,
                "rgroup_value_2": v2,
                "observed_mean_activity": observed,
                "expected_additive_activity": expected,
                "residual_activity": residual,
                "std_activity": row["std"],
                "count": row["count"],
            })

    out = pd.DataFrame(records)

    if not out.empty:
        out = out.sort_values("residual_activity", ascending=False)

    return out

def summarise_top_pairwise_combinations(
    residual_combo_df,
    top_n=20,
    min_count=2,
    score_with_count=True,
):
    """
    Extract best and worst pairwise R-group combinations.

    Primary ranking:
      residual_activity

    Tie-break / support:
      combo_score = residual_activity * log1p(count)

    For pIC50/pKi/pEC50-like activities:
      positive residual = better than additive expectation
      negative residual = worse than additive expectation
    """

    if residual_combo_df is None or residual_combo_df.empty:
        return pd.DataFrame()

    df = residual_combo_df.copy()

    required = {
        "rgroup_label_1",
        "rgroup_value_1",
        "rgroup_label_2",
        "rgroup_value_2",
        "observed_mean_activity",
        "expected_additive_activity",
        "residual_activity",
        "count",
    }

    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"summarise_top_pairwise_combinations missing columns: {missing}. "
            f"Available columns: {df.columns.tolist()}"
        )

    df["count"] = pd.to_numeric(df["count"], errors="coerce")
    df["residual_activity"] = pd.to_numeric(df["residual_activity"], errors="coerce")
    df["observed_mean_activity"] = pd.to_numeric(df["observed_mean_activity"], errors="coerce")
    df["expected_additive_activity"] = pd.to_numeric(df["expected_additive_activity"], errors="coerce")

    df = df.dropna(subset=["residual_activity", "count"])
    df = df[df["count"] >= min_count].copy()

    if df.empty:
        return pd.DataFrame()

    if score_with_count:
        df["combo_score"] = df["residual_activity"] * np.log1p(df["count"])
    else:
        df["combo_score"] = df["residual_activity"]

    best = (
        df.sort_values(
            ["residual_activity", "combo_score"],
            ascending=[False, False]
        )
        .head(top_n)
        .copy()
    )
    best["combo_class"] = "best_positive_residual"

    worst = (
        df.sort_values(
            ["residual_activity", "combo_score"],
            ascending=[True, True]
        )
        .head(top_n)
        .copy()
    )
    worst["combo_class"] = "worst_negative_residual"

    out = pd.concat([best, worst], ignore_index=True)

    out["combo_label"] = (
        out["rgroup_label_1"].astype(str)
        + " + "
        + out["rgroup_label_2"].astype(str)
    )

    # Final plotting order: positive residuals first, then negative residuals.
    out["combo_class_order"] = out["combo_class"].map({
        "best_positive_residual": 0,
        "worst_negative_residual": 1,
    })

    out = out.sort_values(
        ["combo_class_order", "residual_activity", "combo_score"],
        ascending=[True, False, False]
    ).drop(columns=["combo_class_order"])

    return out.reset_index(drop=True)

def plot_top_pairwise_combinations(
    combo_summary_df,
    output_path,
    title=None,
    top_n_each=20,
    img_size=140,
    dpi=300,
):
    """
    Plot best/worst pairwise R-group combinations as structure rows.

    Visual encoding:
      - coloured horizontal bar = residual_activity
      - dashed vertical marker = combo_score

    Ranking:
      - primary: residual_activity
      - secondary: combo_score

    Intended to be called per series/core.
    """

    if combo_summary_df is None or combo_summary_df.empty:
        return None

    df = combo_summary_df.copy()

    required = {
        "combo_class",
        "rgroup_label_1",
        "rgroup_value_1",
        "rgroup_label_2",
        "rgroup_value_2",
        "observed_mean_activity",
        "expected_additive_activity",
        "residual_activity",
        "count",
        "combo_score",
    }

    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"plot_top_pairwise_combinations missing columns: {missing}. "
            f"Available columns: {df.columns.tolist()}"
        )

    numeric_cols = [
        "observed_mean_activity",
        "expected_additive_activity",
        "residual_activity",
        "count",
        "combo_score",
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["residual_activity", "combo_score", "count"]).copy()

    if df.empty:
        return None

    best = (
        df[df["combo_class"].isin(["best_positive_residual", "best_positive_nonadditivity"])]
        .sort_values(["residual_activity", "combo_score"], ascending=[False, False])
        .head(top_n_each)
    )

    worst = (
        df[df["combo_class"].isin(["worst_negative_residual", "worst_negative_nonadditivity"])]
        .sort_values(["residual_activity", "combo_score"], ascending=[True, True])
        .head(top_n_each)
    )

    plot_df = pd.concat([best, worst], ignore_index=True)

    if plot_df.empty:
        return None

    n = len(plot_df)

    residual_vals = plot_df["residual_activity"].dropna()
    score_vals = plot_df["combo_score"].dropna()

    max_abs = max(
        abs(residual_vals.min()) if not residual_vals.empty else 0,
        abs(residual_vals.max()) if not residual_vals.empty else 0,
        abs(score_vals.min()) if not score_vals.empty else 0,
        abs(score_vals.max()) if not score_vals.empty else 0,
    )

    if max_abs == 0 or pd.isna(max_abs):
        max_abs = 1.0

    max_abs = max_abs * 1.10

    # Colour bars by residual value.
    residual_norm = plt.Normalize(-max_abs, max_abs)
    residual_cmap = plt.cm.coolwarm_r

    fig_w = 14.0
    row_h = 1.20
    legend_h = 0.45
    title_h = 0.55 if title else 0.10
    fig_h = max(4, n * row_h + legend_h + title_h)

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Add one header row above the data rows.
    gs = gridspec.GridSpec(
        nrows=n + 1,
        ncols=4,
        height_ratios=[0.35] + [1.0] * n,
        width_ratios=[1.15, 1.15, 3.0, 4.6],
        hspace=0.38,
        wspace=0.35,
        figure=fig,
    )

    # Global legend/header row
    ax_header = fig.add_subplot(gs[0, :])
    ax_header.axis("off")
    ax_header.text(
        0.0,
        0.5,
        "Coloured bar = residual activity   |   dashed vertical marker = support-weighted score",
        ha="left",
        va="center",
        fontsize=9,
        color="#333333",
        transform=ax_header.transAxes,
    )

    for row_idx, (_, row) in enumerate(plot_df.iterrows(), start=1):
        ax1 = fig.add_subplot(gs[row_idx, 0])
        ax2 = fig.add_subplot(gs[row_idx, 1])
        ax_metric = fig.add_subplot(gs[row_idx, 2])
        ax_text = fig.add_subplot(gs[row_idx, 3])

        for ax in [ax1, ax2, ax_text]:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_frame_on(False)

        smi1 = row["rgroup_value_1"]
        smi2 = row["rgroup_value_2"]

        img1 = smiles_to_pil(smi1, size=(img_size, img_size))
        img2 = smiles_to_pil(smi2, size=(img_size, img_size))

        if img1 is not None:
            ax1.imshow(img1)

        if img2 is not None:
            ax2.imshow(img2)

        ax1.set_title(str(row["rgroup_label_1"]), fontsize=9)
        ax2.set_title(str(row["rgroup_label_2"]), fontsize=9)

        residual = row["residual_activity"]
        score = row["combo_score"]

        # Metric axis: residual bar + score marker.
        ax_metric.set_xlim(-max_abs, max_abs)
        ax_metric.set_ylim(0, 1)
        ax_metric.set_yticks([])
        ax_metric.set_frame_on(False)

        # Zero line.
        ax_metric.axvline(0, color="black", linewidth=0.8)

        # Residual as coloured main bar.
        residual_color = residual_cmap(residual_norm(residual))

        ax_metric.barh(
            y=0.5,
            width=residual,
            height=0.38,
            color=residual_color,
            edgecolor="#333333",
            linewidth=0.6,
        )

        # Score as dashed vertical marker.
        ax_metric.axvline(
            score,
            ymin=0.18,
            ymax=0.82,
            color="#111111",
            linestyle="--",
            linewidth=1.5,
        )

        # Show x-axis only on the last data row.
        if row_idx == n:
            ax_metric.set_xticks([-max_abs, 0, max_abs])
            ax_metric.set_xticklabels(
                [f"{-max_abs:.2f}", "0", f"{max_abs:.2f}"],
                fontsize=8,
            )
            ax_metric.set_xlabel("Residual activity and support-weighted score", fontsize=8)
        else:
            ax_metric.set_xticks([])

        combo_class = str(row["combo_class"]).replace("_", " ")

        txt = (
            f"{combo_class}\n"
            f"observed={row['observed_mean_activity']:.2f} | "
            f"expected={row['expected_additive_activity']:.2f}\n"
            f"residual={row['residual_activity']:+.2f} | "
            f"n={int(row['count'])} | "
            f"score={row['combo_score']:+.2f}"
        )

        ax_text.text(
            0.0,
            0.5,
            txt,
            va="center",
            ha="left",
            fontsize=9,
        )

    if title:
        fig.suptitle(title, fontsize=13, y=0.995)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path

def smiles_hash(smiles, n=10):
    return hashlib.sha1(str(smiles).encode()).hexdigest()[:n]


def save_smiles_image(smiles, output_path, size=(180, 180)):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None

    try:
        img = Draw.MolToImage(mol, size=size)
        img.save(output_path)
        return output_path
    except Exception:
        return None


def dump_combo_structure_images(
    combo_summary_df,
    output_dir,
    labelled_core_smiles=None,
    core_img_size=260,
    rgroup_img_size=180,
):
    """
    Save core and R-group structure images used in top pairwise combo summaries.

    Creates:
      structures/core_labelled.png
      structures/rgroups/<Rlabel>_<hash>.png
      structures/structure_manifest.csv
    """

    if combo_summary_df is None or combo_summary_df.empty:
        return None

    structure_dir = Path(output_dir) / "structures"
    rgroup_dir = structure_dir / "rgroups"

    structure_dir.mkdir(parents=True, exist_ok=True)
    rgroup_dir.mkdir(parents=True, exist_ok=True)

    manifest = []

    # Labelled core
    if labelled_core_smiles:
        core_mol = Chem.MolFromSmiles(labelled_core_smiles)
        if core_mol is not None:
            core_img = draw_core_with_rgroup_labels(
                core_mol,
                size=(core_img_size, core_img_size),
                font_size=50,
            )
            if core_img is not None:
                core_path = structure_dir / "core_labelled.png"
                core_img.save(core_path)

                manifest.append({
                    "kind": "core",
                    "rgroup_label": None,
                    "smiles": labelled_core_smiles,
                    "image_path": str(core_path),
                })

    # Unique R-group structures from both sides of the pairwise combos
    records = []

    for _, row in combo_summary_df.iterrows():
        records.append({
            "rgroup_label": row["rgroup_label_1"],
            "smiles": row["rgroup_value_1"],
        })
        records.append({
            "rgroup_label": row["rgroup_label_2"],
            "smiles": row["rgroup_value_2"],
        })

    rg_img_df = pd.DataFrame(records).drop_duplicates()

    for _, row in rg_img_df.iterrows():
        rlabel = str(row["rgroup_label"])
        smi = str(row["smiles"])
        h = smiles_hash(smi)

        img_path = rgroup_dir / f"{safe_name(rlabel)}_{h}.png"

        saved = save_smiles_image(
            smi,
            img_path,
            size=(rgroup_img_size, rgroup_img_size),
        )

        if saved is not None:
            manifest.append({
                "kind": "rgroup",
                "rgroup_label": rlabel,
                "smiles": smi,
                "image_path": str(img_path),
            })

    manifest_df = pd.DataFrame(manifest)
    manifest_csv = structure_dir / "structure_manifest.csv"
    manifest_df.to_csv(manifest_csv, index=False)

    return {
        "structure_dir": str(structure_dir),
        "manifest_csv": str(manifest_csv),
    }

def is_hydrogen_rgroup_smiles(smiles, treat_dummy_only_as_h=True):
    """
    Return True if an R-group fragment represents Hydrogen / no real substituent.

    RDKit R-group decomposition may encode H-like substituents as e.g.
      [H][*:1]
      [*:1][H]
      [1*][H]
      [*:1]

    We treat fragments with only H atoms and dummy atoms as Hydrogen.
    Optionally, dummy-only fragments are also treated as Hydrogen.
    """

    if smiles is None or pd.isna(smiles):
        return False

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return False

    real_atoms = [
        atom for atom in mol.GetAtoms()
        if atom.GetAtomicNum() != 0
    ]

    # Some RDKit R-group outputs may effectively be only the attachment dummy.
    if len(real_atoms) == 0:
        return bool(treat_dummy_only_as_h)

    return all(atom.GetAtomicNum() == 1 for atom in real_atoms)


def filter_non_hydrogen_pairwise_combos(
    combo_df,
    value_col_1="rgroup_value_1",
    value_col_2="rgroup_value_2",
    treat_dummy_only_as_h=True,
):
    """
    Keep only pairwise combinations where both R-group values are non-Hydrogen.
    """

    if combo_df is None or combo_df.empty:
        return pd.DataFrame()

    df = combo_df.copy()

    required = {value_col_1, value_col_2}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"filter_non_hydrogen_pairwise_combos missing columns: {missing}. "
            f"Available columns: {df.columns.tolist()}"
        )

    is_h_1 = df[value_col_1].apply(
        lambda x: is_hydrogen_rgroup_smiles(
            x,
            treat_dummy_only_as_h=treat_dummy_only_as_h,
        )
    )

    is_h_2 = df[value_col_2].apply(
        lambda x: is_hydrogen_rgroup_smiles(
            x,
            treat_dummy_only_as_h=treat_dummy_only_as_h,
        )
    )

    return df[~is_h_1 & ~is_h_2].copy()

def draw_core_with_rgroup_labels(mol, size=(300, 300), font_size=40):
    """
    Draw a decomposition core with visible R-group labels.

    Handles dummy atoms labelled by:
      - atom map number, e.g. [*:1]
      - isotope, e.g. [1*]
      - fallback atom index
    """

    if mol is None:
        return None

    try:
        mol = Chem.Mol(mol)
        rdMolDraw2D.PrepareMolForDrawing(mol)
    except Exception:
        pass

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.baseFontSize = font_size / 100

    atom_labels = {}

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            amap = atom.GetAtomMapNum()
            isotope = atom.GetIsotope()

            if amap:
                label = f"R{amap}"
            elif isotope:
                label = f"R{isotope}"
            else:
                label = f"R?{atom.GetIdx()}"

            atom_labels[atom.GetIdx()] = label

    for idx, label in atom_labels.items():
        opts.atomLabels[idx] = label

    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()

    png_data = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_data))