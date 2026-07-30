import hashlib
import re

from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import json

from rdkit import Chem
from itertools import product
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def get_core_attachment_labels(labelled_core_smiles):
    """
    Return attachment labels actually present on a labelled core.

    Example:
      c1cc([*:1])ccc([*:2])c1 -> ["1", "2"]

    Supports atom-map dummy notation:
      [*:1], [*:2], ...
    """

    if labelled_core_smiles is None or pd.isna(labelled_core_smiles):
        return []

    return sorted(
        set(re.findall(r"\[\*:(\d+)\]", str(labelled_core_smiles))),
        key=lambda x: int(x)
    )


def get_core_attachment_rcols(labelled_core_smiles, available_r_cols):
    """
    Convert labels present on the labelled core into R-column names.

    Example:
      labels ["1", "2", "8"] -> ["R1", "R2", "R8"]

    Only returns columns that are actually available in the matrix.
    """

    labels = get_core_attachment_labels(labelled_core_smiles)

    active_r_cols = [
        f"R{label}"
        for label in labels
        if f"R{label}" in available_r_cols
    ]

    return active_r_cols


def get_representative_labelled_core(work_df):
    """
    Pick the most common labelled_core for a series/core group.
    """

    if "labelled_core" not in work_df.columns:
        return None

    labelled_core_values = (
        work_df["labelled_core"]
        .dropna()
        .astype(str)
        .mode()
    )

    if len(labelled_core_values) == 0:
        return None

    return labelled_core_values.iloc[0]

def summarise_fw_prediction_distributions(
    observed_df,
    virtual_df,
    output_path,
    pred_col="fw_predicted_activity",
    actual_col=None,
):
    records = []

    def add_summary(label, values):
        values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()

        if values.empty:
            return

        records.append({
            "set": label,
            "n": len(values),
            "mean": values.mean(),
            "median": values.median(),
            "p10": values.quantile(0.10),
            "p25": values.quantile(0.25),
            "p75": values.quantile(0.75),
            "p90": values.quantile(0.90),
            "p95": values.quantile(0.95),
            "max": values.max(),
        })

    if observed_df is not None and pred_col in observed_df.columns:
        add_summary("observed_fw_predictions", observed_df[pred_col])

    if observed_df is not None and actual_col and actual_col in observed_df.columns:
        add_summary("observed_actual", observed_df[actual_col])

    if virtual_df is not None and pred_col in virtual_df.columns:
        virt = virtual_df.copy()

        if "reconstruction_success" in virt.columns:
            virt = virt[virt["reconstruction_success"] == True]

        if "has_unresolved_dummy_atoms" in virt.columns:
            virt = virt[virt["has_unresolved_dummy_atoms"] == False]

        add_summary("virtual_fw_predictions_valid_reconstructed", virt[pred_col])

    summary_df = pd.DataFrame(records)
    summary_df.to_csv(output_path, index=False)

    return output_path

def plot_virtual_reconstruction_failure_diagnostics(
    candidate_df,
    reconstruction_r_cols,
    output_path,
    series=None,
    core_hash=None,
    dpi=300,
):
    """
    Per-core diagnostic plot for virtual reconstruction.

    Panels:
      1. Reconstruction status counts.
      2. Unresolved attachment label counts.
      3. Missing token counts by reconstruction R-column.
    """

    if candidate_df is None or candidate_df.empty:
        return None

    df = candidate_df.copy()

    required_cols = {
        "reconstruction_success",
        "has_unresolved_dummy_atoms",
        "reconstructed_smiles",
    }

    missing = required_cols - set(df.columns)
    if missing:
        logger.warning(
            f"[Free-Wilson] Cannot plot reconstruction diagnostics, missing: {missing}"
        )
        return None

    if "remaining_attachment_labels" not in df.columns:
        df["remaining_attachment_labels"] = df["reconstructed_smiles"].apply(
            extract_attachment_labels_from_smiles
        )

    if "n_remaining_attachment_labels" not in df.columns:
        df["n_remaining_attachment_labels"] = df["remaining_attachment_labels"].apply(len)

    if "reconstruction_status" not in df.columns:
        df["reconstruction_status"] = np.select(
            [
                df["reconstruction_success"] & ~df["has_unresolved_dummy_atoms"],
                ~df["reconstruction_success"],
                df["has_unresolved_dummy_atoms"],
            ],
            [
                "final_smiles_ok",
                "reconstruction_null",
                "unresolved_dummy_atoms",
            ],
            default="unknown",
        )

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    # Panel 1: status counts
    status_order = [
        "final_smiles_ok",
        "unresolved_dummy_atoms",
        "reconstruction_null",
        "unknown",
    ]

    status_counts = (
        df["reconstruction_status"]
        .value_counts()
        .reindex(status_order)
        .dropna()
        .reset_index()
    )
    status_counts.columns = ["status", "count"]

    sns.barplot(
        data=status_counts,
        x="count",
        y="status",
        ax=axes[0],
        color="#4C78A8",
    )
    axes[0].set_title("Reconstruction status")
    axes[0].set_xlabel("Count")
    axes[0].set_ylabel("")

    total = len(df)
    for i, row in status_counts.iterrows():
        frac = row["count"] / total if total else 0
        axes[0].text(
            row["count"],
            i,
            f" {row['count']} ({frac:.1%})",
            va="center",
            fontsize=8,
        )

    # Panel 2: unresolved attachment labels
    label_records = []

    failed_df = df[df["reconstruction_status"] != "final_smiles_ok"].copy()

    for labels in failed_df["remaining_attachment_labels"]:
        if labels:
            for lab in labels:
                label_records.append({"label": f"R{lab}"})
        else:
            label_records.append({"label": "none/parse_fail"})

    if label_records:
        label_df = (
            pd.DataFrame(label_records)
            .value_counts("label")
            .reset_index(name="count")
            .sort_values("count", ascending=False)
        )
    else:
        label_df = pd.DataFrame(columns=["label", "count"])

    if not label_df.empty:
        sns.barplot(
            data=label_df,
            x="count",
            y="label",
            ax=axes[1],
            color="#F58518",
        )
    axes[1].set_title("Unresolved labels in failures")
    axes[1].set_xlabel("Failure count")
    axes[1].set_ylabel("")

    # Panel 3: missing token / NA by reconstruction column
    missing_records = []

    for r in reconstruction_r_cols:
        recon_col = f"recon_{r}"
        source_col = recon_col if recon_col in df.columns else r

        if source_col not in df.columns:
            missing_records.append({
                "rgroup_col": r,
                "metric": "column_missing",
                "count": len(df),
            })
            continue

        vals = df[source_col].astype("object")
        vals_str = vals.astype(str)

        missing_records.append({
            "rgroup_col": r,
            "metric": "__MISSING__",
            "count": int((vals_str == "__MISSING__").sum()),
        })
        missing_records.append({
            "rgroup_col": r,
            "metric": "NaN",
            "count": int(vals.isna().sum()),
        })

    missing_df = pd.DataFrame(missing_records)

    if not missing_df.empty:
        pivot = missing_df.pivot_table(
            index="rgroup_col",
            columns="metric",
            values="count",
            aggfunc="sum",
            fill_value=0,
        )

        sns.heatmap(
            pivot,
            annot=True,
            fmt=".0f",
            cmap="Reds",
            linewidths=0.5,
            ax=axes[2],
            cbar_kws={"label": "Count"},
        )

    axes[2].set_title("Missing values by R-position")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("")

    title_bits = ["Virtual reconstruction diagnostics"]
    if series is not None:
        title_bits.append(str(series))
    if core_hash is not None:
        title_bits.append(f"core {core_hash}")

    fig.suptitle(" | ".join(title_bits), fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path

def plot_overall_virtual_reconstruction_diagnostics(
    virtual_df,
    output_path,
    dpi=300,
    top_n_cores=30,
):
    """
    Overall diagnostic plot across all FW virtual candidates.

    Panels:
      1. Overall reconstruction status.
      2. Failure rate by series.
      3. Top failing cores.
      4. Unresolved attachment labels overall.
    """

    if virtual_df is None or virtual_df.empty:
        logger.warning("[Free-Wilson] Overall reconstruction plot skipped: no virtual data.")
        return None

    df = virtual_df.copy()

    if "reconstruction_status" not in df.columns:
        df["reconstruction_status"] = np.select(
            [
                df["reconstruction_success"] & ~df["has_unresolved_dummy_atoms"],
                ~df["reconstruction_success"],
                df["has_unresolved_dummy_atoms"],
            ],
            [
                "final_smiles_ok",
                "reconstruction_null",
                "unresolved_dummy_atoms",
            ],
            default="unknown",
        )

    if "remaining_attachment_labels" not in df.columns:
        df["remaining_attachment_labels"] = df["reconstructed_smiles"].apply(
            extract_attachment_labels_from_smiles
        )

    df["is_failure"] = df["reconstruction_status"] != "final_smiles_ok"

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # Panel 1: overall status
    status_counts = (
        df["reconstruction_status"]
        .value_counts()
        .reset_index()
    )
    status_counts.columns = ["status", "count"]

    sns.barplot(
        data=status_counts,
        x="count",
        y="status",
        ax=axes[0, 0],
        color="#4C78A8",
    )
    axes[0, 0].set_title("Overall reconstruction status")
    axes[0, 0].set_xlabel("Count")
    axes[0, 0].set_ylabel("")

    total = len(df)
    for i, row in status_counts.iterrows():
        frac = row["count"] / total if total else 0
        axes[0, 0].text(
            row["count"],
            i,
            f" {row['count']} ({frac:.1%})",
            va="center",
            fontsize=8,
        )

    # Panel 2: failure rate by series
    if "series" in df.columns:
        series_summary = (
            df.groupby("series")
            .agg(
                n_candidates=("is_failure", "size"),
                n_failures=("is_failure", "sum"),
            )
            .reset_index()
        )
        series_summary["failure_rate"] = (
            series_summary["n_failures"] / series_summary["n_candidates"]
        )
        series_summary = series_summary.sort_values("failure_rate", ascending=False)

        sns.barplot(
            data=series_summary,
            x="failure_rate",
            y="series",
            ax=axes[0, 1],
            color="#E45756",
        )
        axes[0, 1].set_title("Failure rate by series")
        axes[0, 1].set_xlabel("Failure rate")
        axes[0, 1].set_ylabel("")
        axes[0, 1].set_xlim(0, 1)
    else:
        axes[0, 1].axis("off")

    # Panel 3: top failing cores
    if {"series", "core_hash"}.issubset(df.columns):
        core_summary = (
            df.groupby(["series", "core_hash"])
            .agg(
                n_candidates=("is_failure", "size"),
                n_failures=("is_failure", "sum"),
            )
            .reset_index()
        )
        core_summary["failure_rate"] = (
            core_summary["n_failures"] / core_summary["n_candidates"]
        )
        core_summary["series_core"] = (
            core_summary["series"].astype(str)
            + " | "
            + core_summary["core_hash"].astype(str)
        )

        plot_core_summary = (
            core_summary
            .sort_values(["failure_rate", "n_failures"], ascending=[False, False])
            .head(top_n_cores)
        )

        sns.barplot(
            data=plot_core_summary,
            x="failure_rate",
            y="series_core",
            ax=axes[1, 0],
            color="#F58518",
        )
        axes[1, 0].set_title(f"Top {top_n_cores} failing cores")
        axes[1, 0].set_xlabel("Failure rate")
        axes[1, 0].set_ylabel("")
        axes[1, 0].set_xlim(0, 1)
    else:
        axes[1, 0].axis("off")

    # Panel 4: unresolved labels overall
    failed_df = df[df["is_failure"]].copy()

    label_records = []
    for labels in failed_df["remaining_attachment_labels"]:
        if labels:
            for lab in labels:
                label_records.append({"label": f"R{lab}"})
        else:
            label_records.append({"label": "none/parse_fail"})

    if label_records:
        label_df = (
            pd.DataFrame(label_records)
            .value_counts("label")
            .reset_index(name="count")
            .sort_values("count", ascending=False)
        )

        sns.barplot(
            data=label_df,
            x="count",
            y="label",
            ax=axes[1, 1],
            color="#72B7B2",
        )
        axes[1, 1].set_title("Unresolved attachment labels overall")
        axes[1, 1].set_xlabel("Failure count")
        axes[1, 1].set_ylabel("")
    else:
        axes[1, 1].axis("off")

    fig.suptitle("Free-Wilson virtual reconstruction diagnostics", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    logger.info(f"[Free-Wilson] Overall reconstruction diagnostics saved: {output_path}")

    return output_path

def extract_attachment_labels_from_smiles(smiles):
    """
    Extract mapped dummy attachment labels from a SMILES string.

    Examples:
      N([*:1])C -> ["1"]
      c1cc([*:2])ccc1 -> ["2"]
    """
    if smiles is None or pd.isna(smiles):
        return []

    return re.findall(r"\[\*:(\d+)\]", str(smiles))


def count_attachment_labels(smiles):
    """
    Count unresolved mapped dummy attachment labels in a SMILES string.
    """
    return len(extract_attachment_labels_from_smiles(smiles))


def summarise_virtual_reconstruction_failures(
    candidate_df,
    reconstruction_r_cols,
    output_dir,
    series,
    core,
):
    """
    Write detailed diagnostics for virtual reconstruction failures.

    Outputs:
      fw_virtual_reconstruction_failures.csv
      fw_virtual_reconstruction_failure_summary.csv
      fw_virtual_missing_rgroup_diagnostics.csv
      fw_virtual_failure_rgroup_value_counts.csv
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    core_hash = hashlib.sha1(str(core).encode()).hexdigest()[:8]

    if candidate_df is None or candidate_df.empty:
        return {
            "failure_csv": None,
            "failure_summary_csv": None,
            "missing_diag_csv": None,
            "value_counts_csv": None,
            "n_fail": 0,
            "n_total": 0,
        }

    df = candidate_df.copy()

    # Make sure expected diagnostic columns exist.
    if "reconstruction_success" not in df.columns:
        df["reconstruction_success"] = df["reconstructed_smiles"].notna()

    if "has_unresolved_dummy_atoms" not in df.columns:
        df["has_unresolved_dummy_atoms"] = df["reconstructed_smiles"].apply(
            smiles_has_dummy_atoms
        )

    df["remaining_attachment_labels"] = df["reconstructed_smiles"].apply(
        extract_attachment_labels_from_smiles
    )

    df["n_remaining_attachment_labels"] = df["remaining_attachment_labels"].apply(len)

    df["remaining_attachment_labels_str"] = df["remaining_attachment_labels"].apply(
        lambda xs: ",".join(xs) if xs else ""
    )

    failure_mask = (
        ~df["reconstruction_success"]
        | df["has_unresolved_dummy_atoms"]
        | (df["n_remaining_attachment_labels"] > 0)
    )

    failed_df = df[failure_mask].copy()

    failure_csv = None
    failure_summary_csv = None
    missing_diag_csv = None
    value_counts_csv = None

    if not failed_df.empty:
        failure_csv = output_dir / "fw_virtual_reconstruction_failures.csv"
        failed_df.to_csv(failure_csv, index=False)

    # Summary by unresolved attachment labels
    if not failed_df.empty:
        label_records = []

        for labels in failed_df["remaining_attachment_labels"]:
            if not labels:
                label_records.append({"remaining_attachment_label": "NONE_OR_PARSE_FAILURE"})
            else:
                for lab in labels:
                    label_records.append({"remaining_attachment_label": f"R{lab}"})

        label_summary_df = (
            pd.DataFrame(label_records)
            .value_counts("remaining_attachment_label")
            .reset_index(name="count")
        )

        label_summary_df["series"] = series
        label_summary_df["core"] = core
        label_summary_df["core_hash"] = core_hash
        label_summary_df["fraction_of_failures"] = (
            label_summary_df["count"] / len(failed_df)
        )

        failure_summary_csv = output_dir / "fw_virtual_reconstruction_failure_summary.csv"
        label_summary_df.to_csv(failure_summary_csv, index=False)

    # Missing-value diagnostics for all candidate rows
    missing_diag_records = []

    # for r in reconstruction_r_cols:
    #     if r not in df.columns:
    #         missing_diag_records.append({
    #             "series": series,
    #             "core": core,
    #             "core_hash": core_hash,
    #             "rgroup_col": r,
    #             "n_candidates": len(df),
    #             "column_present": False,
    #             "n_na": np.nan,
    #             "n_missing_token": np.nan,
    #             "n_unique_values": np.nan,
    #         })
    #         continue

    #     vals = df[r].astype("object")
    #     vals_str = vals.astype(str)

    #     missing_diag_records.append({
    #         "series": series,
    #         "core": core,
    #         "core_hash": core_hash,
    #         "rgroup_col": r,
    #         "n_candidates": len(df),
    #         "column_present": True,
    #         "n_na": int(vals.isna().sum()),
    #         "n_missing_token": int((vals_str == "__MISSING__").sum()),
    #         "n_unique_values": int(vals_str.nunique()),
    #     })

    # missing_diag_df = pd.DataFrame(missing_diag_records)
    # missing_diag_csv = output_dir / "fw_virtual_missing_rgroup_diagnostics.csv"
    # missing_diag_df.to_csv(missing_diag_csv, index=False)

    for r in reconstruction_r_cols:
        recon_col = f"recon_{r}"
        source_col = recon_col if recon_col in df.columns else r

        if source_col not in df.columns:
            missing_diag_records.append({
                "series": series,
                "core": core,
                "core_hash": core_hash,
                "rgroup_col": r,
                "diagnostic_source_col": source_col,
                "n_candidates": len(df),
                "column_present": False,
                "n_na": np.nan,
                "n_missing_token": np.nan,
                "n_unique_values": np.nan,
            })
            continue

        vals = df[source_col].astype("object")
        vals_str = vals.astype(str)

        missing_diag_records.append({
            "series": series,
            "core": core,
            "core_hash": core_hash,
            "rgroup_col": r,
            "diagnostic_source_col": source_col,
            "n_candidates": len(df),
            "column_present": True,
            "n_na": int(vals.isna().sum()),
            "n_missing_token": int((vals_str == "__MISSING__").sum()),
            "n_unique_values": int(vals_str.nunique()),
        })

    # Value counts in failures by R-column, useful for spotting a bad default or missing position.
    value_count_records = []

    if not failed_df.empty:
        for r in reconstruction_r_cols:
            recon_col = f"recon_{r}"
            source_col = recon_col if recon_col in failed_df.columns else r

            if source_col not in failed_df.columns:
                continue

            vc = (
                failed_df[source_col]
                .astype(str)
                .value_counts(dropna=False)
                .head(30)
            )

            for value, count in vc.items():
                value_count_records.append({
                    "series": series,
                    "core": core,
                    "core_hash": core_hash,
                    "rgroup_col": r,
                    "diagnostic_source_col": source_col,
                    "rgroup_value": value,
                    "failure_count": int(count),
                    "fraction_of_failures": float(count / len(failed_df)),
                })

    value_counts_df = pd.DataFrame(value_count_records)
    value_counts_csv = output_dir / "fw_virtual_failure_rgroup_value_counts.csv"
    value_counts_df.to_csv(value_counts_csv, index=False)

    logger.debug(
        f"[Free-Wilson] Virtual reconstruction failure diagnostics for "
        f"series={series}, core={core_hash}: "
        f"failures={len(failed_df)}/{len(df)}"
    )

    if failure_summary_csv:
        logger.debug(
            f"[Free-Wilson] Failure summary saved: {failure_summary_csv}"
        )

    return {
        "failure_csv": str(failure_csv) if failure_csv else None,
        "failure_summary_csv": str(failure_summary_csv) if failure_summary_csv else None,
        "missing_diag_csv": str(missing_diag_csv) if missing_diag_csv else None,
        "value_counts_csv": str(value_counts_csv) if value_counts_csv else None,
        "n_fail": len(failed_df),
        "n_total": len(df),
    }

def mol_has_dummy_atoms(mol):
    if mol is None:
        return True

    return any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())

def smiles_has_dummy_atoms(smiles):
    """
    Return True if a SMILES contains unresolved dummy atoms.
    """
    if smiles is None or pd.isna(smiles):
        return True

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return True

    return mol_has_dummy_atoms(mol)

def get_group_default_rgroups(work_df, reconstruction_cols):
    """
    Get one representative value for every active R column in the group.

    Important:
      - Do not use "__MISSING__" as a default.
      - Only reconstruction columns that correspond to actual attachment points
        should be passed into this function.
    """

    defaults = {}

    for r in reconstruction_cols:
        if r not in work_df.columns:
            continue

        vals = work_df[r].dropna().astype(str)
        vals = vals[vals != "__MISSING__"]

        if vals.empty:
            continue

        mode_vals = vals.mode()
        defaults[r] = mode_vals.iloc[0] if len(mode_vals) > 0 else vals.iloc[0]

    return defaults

def get_rgroup_columns(wide_df):
    return [
        c for c in wide_df.columns
        if isinstance(c, str)
        and c.startswith("R")
        and c[1:].isdigit()
    ]

def make_fw_design_matrix(df, r_cols):
    """
    One-hot encode R-group columns for Free-Wilson modelling.

    Missing R-groups are encoded explicitly as '__MISSING__'
    so that matrix shape is stable.
    """
    x = df[r_cols].copy()

    for c in r_cols:
        x[c] = x[c].astype("object").where(x[c].notna(), "__MISSING__")
        x[c] = x[c].astype(str)

    return pd.get_dummies(x, columns=r_cols, prefix=r_cols)


def fit_free_wilson_model(x, y, model_type="ridge", ridge_alpha=1.0):
    """
    Fit vanilla additive Free-Wilson model.

    model_type:
      - ridge: regularized FW, recommended for sparse SAR matrices
      - linear: ordinary least squares
    """
    if model_type == "linear":
        model = LinearRegression()
    elif model_type == "ridge":
        model = Ridge(alpha=ridge_alpha)
    else:
        raise ValueError(f"Unsupported Free-Wilson model_type: {model_type}")

    model.fit(x, y)
    return model


def calculate_regression_metrics(y_true, y_pred):
    y_true = pd.to_numeric(pd.Series(y_true), errors="coerce")
    y_pred = pd.to_numeric(pd.Series(y_pred), errors="coerce")

    valid = y_true.notna() & y_pred.notna()
    y_true = y_true[valid]
    y_pred = y_pred[valid]

    if len(y_true) < 2:
        return {
            "n": len(y_true),
            "r2": np.nan,
            "rmse": np.nan,
            "mae": np.nan,
        }

    mse = mean_squared_error(y_true, y_pred)

    return {
        "n": len(y_true),
        "r2": r2_score(y_true, y_pred),
        "rmse": float(np.sqrt(mse)),
        "mae": mean_absolute_error(y_true, y_pred),
    }


def generate_free_wilson_candidates(
    core_df,
    r_cols,
    id_col="compound_id",
    max_candidates=50000,
    max_changed_positions=None,
):
    """
    Enumerate unsynthesised R-group recombinations using observed R-group values.

    If max_changed_positions is None:
      generate the full Cartesian product of observed R-group levels.

    If max_changed_positions is an integer:
      keep candidates that differ from at least one observed molecule
      by <= max_changed_positions.
      This is useful for conservative enumeration.
    """

    if core_df.empty or not r_cols:
        return pd.DataFrame()

    work = core_df.copy()

    for r in r_cols:
        work[r] = work[r].astype("object").where(work[r].notna(), "__MISSING__")
        work[r] = work[r].astype(str)

    rg_values = {
        r: sorted(work[r].dropna().unique().tolist())
        for r in r_cols
    }

    # If any R column has no usable values, stop.
    if any(len(v) == 0 for v in rg_values.values()):
        return pd.DataFrame()

    total_possible = 1
    for r in r_cols:
        total_possible *= len(rg_values[r])

    if total_possible > max_candidates:
        logger.warning(
            f"[Free-Wilson] Candidate space too large: {total_possible} combinations. "
            f"Generating only first {max_candidates} combinations."
        )

    records = []
    for i, combo in enumerate(product(*[rg_values[r] for r in r_cols])):
        if i >= max_candidates:
            break
        records.append(dict(zip(r_cols, combo)))

    cand_df = pd.DataFrame(records)

    if cand_df.empty:
        return cand_df

    observed_keys = work[r_cols].drop_duplicates().copy()
    observed_keys["_observed"] = True

    cand_df = cand_df.merge(
        observed_keys,
        on=r_cols,
        how="left"
    )

    cand_df = cand_df[cand_df["_observed"].isna()].drop(columns=["_observed"])

    if cand_df.empty:
        return cand_df

    # Optional conservative filter: keep only candidates close to existing rows
    if max_changed_positions is not None:
        observed_tuples = [
            tuple(row)
            for row in work[r_cols].drop_duplicates().itertuples(index=False, name=None)
        ]

        def min_changes_to_observed(row):
            candidate_tuple = tuple(row[r] for r in r_cols)
            return min(
                sum(a != b for a, b in zip(candidate_tuple, obs))
                for obs in observed_tuples
            )

        cand_df["min_changed_positions_vs_observed"] = cand_df.apply(
            min_changes_to_observed,
            axis=1
        )

        cand_df = cand_df[
            cand_df["min_changed_positions_vs_observed"] <= max_changed_positions
        ].copy()
    else:
        cand_df["min_changed_positions_vs_observed"] = np.nan

    return cand_df.reset_index(drop=True)


def plot_fw_predicted_vs_actual(df, activity_col, pred_col, output_path, title=None, dpi=300):
    if df is None or df.empty:
        return None

    plot_df = df.copy()
    plot_df[activity_col] = pd.to_numeric(plot_df[activity_col], errors="coerce")
    plot_df[pred_col] = pd.to_numeric(plot_df[pred_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[activity_col, pred_col])

    if plot_df.empty:
        return None

    metrics = calculate_regression_metrics(plot_df[activity_col], plot_df[pred_col])

    plt.figure(figsize=(6, 6))

    sns.scatterplot(
        data=plot_df,
        x=activity_col,
        y=pred_col,
        s=45,
        alpha=0.75
    )

    min_val = min(plot_df[activity_col].min(), plot_df[pred_col].min())
    max_val = max(plot_df[activity_col].max(), plot_df[pred_col].max())

    pad = (max_val - min_val) * 0.05 if max_val > min_val else 0.5
    min_val -= pad
    max_val += pad

    plt.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="black", linewidth=1)

    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)

    plt.xlabel(f"Actual {activity_col}")
    plt.ylabel(f"Predicted {activity_col}")
    plt.title(title or "Free-Wilson predicted vs actual")

    txt = (
        f"n={metrics['n']}\n"
        f"R²={metrics['r2']:.2f}\n"
        f"RMSE={metrics['rmse']:.2f}\n"
        f"MAE={metrics['mae']:.2f}"
    )

    plt.text(
        0.05,
        0.95,
        txt,
        transform=plt.gca().transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

    return output_path


def plot_fw_residuals(df, pred_col, residual_col, output_path, title=None, dpi=300):
    if df is None or df.empty:
        return None

    plot_df = df.copy()
    plot_df[pred_col] = pd.to_numeric(plot_df[pred_col], errors="coerce")
    plot_df[residual_col] = pd.to_numeric(plot_df[residual_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[pred_col, residual_col])

    if plot_df.empty:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    sns.scatterplot(
        data=plot_df,
        x=pred_col,
        y=residual_col,
        s=40,
        alpha=0.75,
        ax=axes[0]
    )
    axes[0].axhline(0, linestyle="--", color="black", linewidth=1)
    axes[0].set_xlabel(f"Predicted {pred_col}")
    axes[0].set_ylabel("Residual")
    axes[0].set_title("Residuals vs predicted")

    sns.histplot(
        plot_df[residual_col],
        bins=25,
        kde=True,
        ax=axes[1]
    )
    axes[1].axvline(0, linestyle="--", color="black", linewidth=1)
    axes[1].set_xlabel("Residual")
    axes[1].set_title("Residual distribution")

    fig.suptitle(title or "Free-Wilson residual diagnostics")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

    return output_path


def plot_fw_coefficients(coeff_df, output_path, top_n=25, dpi=300):
    if coeff_df is None or coeff_df.empty:
        return None

    df = coeff_df.copy()
    df["coefficient"] = pd.to_numeric(df["coefficient"], errors="coerce")
    df = df.dropna(subset=["coefficient"])

    if df.empty:
        return None

    top_pos = df.nlargest(top_n, "coefficient")
    top_neg = df.nsmallest(top_n, "coefficient")

    plot_df = pd.concat([top_pos, top_neg], ignore_index=True)
    plot_df = plot_df.drop_duplicates(subset=["feature"])
    plot_df = plot_df.sort_values("coefficient")

    plt.figure(figsize=(9, max(5, len(plot_df) * 0.28)))

    colors = ["#2ca25f" if v > 0 else "#de2d26" for v in plot_df["coefficient"]]

    plt.barh(plot_df["feature"], plot_df["coefficient"], color=colors)
    plt.axvline(0, color="black", linewidth=1)
    plt.xlabel("Free-Wilson coefficient")
    plt.ylabel("R-group level")
    plt.title("Largest positive and negative Free-Wilson contributions")
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()

    return output_path

def plot_fw_prediction_distributions(
    observed_df,
    virtual_df,
    output_path,
    pred_col="fw_predicted_activity",
    actual_col=None,
    dpi=300,
):
    """
    Compare distribution of FW predictions for observed vs virtual molecules.

    Produces:
      left: layered density histogram
      right: ECDF

    If actual_col is provided and present, observed actual activity is optionally
    shown as a third distribution.
    """

    if observed_df is None:
        observed_df = pd.DataFrame()

    if virtual_df is None:
        virtual_df = pd.DataFrame()

    obs = observed_df.copy()
    virt = virtual_df.copy()

    if pred_col not in obs.columns and pred_col not in virt.columns:
        logger.warning(
            f"[Free-Wilson] Cannot plot prediction distributions: missing '{pred_col}'."
        )
        return None

    plot_records = []

    if pred_col in obs.columns:
        vals = pd.to_numeric(obs[pred_col], errors="coerce").dropna()
        plot_records.extend([
            {"set": "Observed FW predictions", "value": v}
            for v in vals
        ])

    if pred_col in virt.columns:
        # For virtuals, use only successfully reconstructed final molecules if column exists.
        virt_plot = virt.copy()

        if "reconstruction_success" in virt_plot.columns:
            virt_plot = virt_plot[virt_plot["reconstruction_success"] == True]

        if "has_unresolved_dummy_atoms" in virt_plot.columns:
            virt_plot = virt_plot[virt_plot["has_unresolved_dummy_atoms"] == False]

        vals = pd.to_numeric(virt_plot[pred_col], errors="coerce").dropna()
        plot_records.extend([
            {"set": "Virtual FW predictions", "value": v}
            for v in vals
        ])

    if actual_col and actual_col in obs.columns:
        vals = pd.to_numeric(obs[actual_col], errors="coerce").dropna()
        plot_records.extend([
            {"set": "Observed actual", "value": v}
            for v in vals
        ])

    plot_df = pd.DataFrame(plot_records)

    if plot_df.empty:
        logger.warning("[Free-Wilson] Prediction distribution plot skipped: no data.")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Layered histogram / density
    sns.histplot(
        data=plot_df,
        x="value",
        hue="set",
        stat="density",
        common_norm=False,
        element="step",
        fill=True,
        alpha=0.25,
        bins=35,
        kde=True,
        ax=axes[0],
    )

    axes[0].set_xlabel("Activity")
    axes[0].set_ylabel("Density")
    axes[0].set_title("Free-Wilson prediction distributions")

    # ECDF
    sns.ecdfplot(
        data=plot_df,
        x="value",
        hue="set",
        ax=axes[1],
    )

    axes[1].set_xlabel("Activity")
    axes[1].set_ylabel("Cumulative fraction")
    axes[1].set_title("Cumulative distribution")

    # Add useful quantile summary as text
    summary_lines = []

    for label, sub in plot_df.groupby("set"):
        vals = sub["value"].dropna()

        if vals.empty:
            continue

        summary_lines.append(
            f"{label}: n={len(vals)}, "
            f"median={vals.median():.2f}, "
            f"P90={vals.quantile(0.90):.2f}, "
            f"P95={vals.quantile(0.95):.2f}, "
            f"max={vals.max():.2f}"
        )

    axes[1].text(
        0.02,
        0.02,
        "\n".join(summary_lines),
        transform=axes[1].transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    logger.info(f"[Free-Wilson] Prediction distribution plot saved: {output_path}")

    return output_path

def run_free_wilson_for_group(
    group_df,
    r_cols,
    activity_col,
    cfg,
    output_dir,
    series,
    core,
    reconstruction_r_cols=None,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    id_col = cfg.get("id_col", "compound_id")
    model_type = cfg.get("model_type", "ridge")
    ridge_alpha = cfg.get("ridge_alpha", 1.0)

    diagnostics_cfg = cfg.get("diagnostics", {})
    make_plots = diagnostics_cfg.get("make_plots", True)
    make_residual_plots = diagnostics_cfg.get("make_residual_plots", True)
    make_coefficient_plots = diagnostics_cfg.get("make_coefficient_plots", True)
    top_n_coefficients = diagnostics_cfg.get("top_n_coefficients", 25)

    candidate_cfg = cfg.get("candidate_generation", {})
    generate_candidates = candidate_cfg.get("enabled", True)
    max_candidates = candidate_cfg.get("max_candidates_per_group", 50000)
    max_changed_positions = candidate_cfg.get("max_changed_positions", None)

    work = group_df.copy()

    all_matrix_r_cols = get_rgroup_columns(work)

    # Pick the actual labelled core for this series/core group.
    labelled_core_smiles = get_representative_labelled_core(work)

    # Infer which R columns actually exist as attachment points on this core.
    core_attachment_r_cols = get_core_attachment_rcols(
        labelled_core_smiles,
        available_r_cols=all_matrix_r_cols,
    )

    if reconstruction_r_cols is None:
        reconstruction_r_cols = core_attachment_r_cols
    else:
        # Intersect caller-provided reconstruction columns with actual core attachment columns.
        reconstruction_r_cols = [
            r for r in reconstruction_r_cols
            if r in core_attachment_r_cols
        ]

    logger.debug(
        f"[Free-Wilson] Core attachment columns for series={series}, "
        f"core={hashlib.sha1(str(core).encode()).hexdigest()[:8]}: "
        f"{reconstruction_r_cols}"
    )

    work[activity_col] = pd.to_numeric(work[activity_col], errors="coerce")
    work = work.dropna(subset=[activity_col])

    if work.empty:
        return None

    # Keep FW model columns restricted to actual core attachment positions.
    # This prevents global matrix columns from entering groups where
    # those attachment points do not exist.
    r_cols = [
        r for r in r_cols
        if r in reconstruction_r_cols
    ]

    if not r_cols:
        logger.info(
            f"[Free-Wilson] Skipping series={series}, "
            f"core={hashlib.sha1(str(core).encode()).hexdigest()[:8]}: "
            f"no FW model columns overlap actual core attachment points."
        )
        return None

    # Ensure categorical R columns are string-like and non-null.
    for r in r_cols:
        work[r] = work[r].astype("object").where(work[r].notna(), "__MISSING__")
        work[r] = work[r].astype(str)

    # Normalise reconstruction R columns too.
    for r in reconstruction_r_cols:
        if r in work.columns:
            work[r] = work[r].astype("object").where(work[r].notna(), "__MISSING__")
            work[r] = work[r].astype(str)

    x_fw = make_fw_design_matrix(work, r_cols)
    y = work[activity_col].values

    model = fit_free_wilson_model(
        x_fw,
        y,
        model_type=model_type,
        ridge_alpha=ridge_alpha,
    )

    work["fw_predicted_activity"] = model.predict(x_fw)
    work["fw_residual"] = work[activity_col] - work["fw_predicted_activity"]
    work["fw_abs_error"] = work["fw_residual"].abs()

    if "labelled_core" in work.columns:
        work["reconstructed_smiles_from_rgroups"] = work.apply(
            lambda row: reconstruct_molecule_from_labelled_core(
                row.get("labelled_core"),
                {
                    r: row.get(r)
                    for r in reconstruction_r_cols
                },
            ),
            axis=1,
        )

        if "smiles" in work.columns:
            work["canonical_input_smiles"] = work["smiles"].apply(
                lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(str(s)), canonical=True)
                if pd.notna(s) and Chem.MolFromSmiles(str(s)) is not None
                else None
            )

            work["reconstruction_matches_input"] = (
                work["canonical_input_smiles"] == work["reconstructed_smiles_from_rgroups"]
            )

    metrics = calculate_regression_metrics(
        work[activity_col],
        work["fw_predicted_activity"],
    )

    observed_csv = output_dir / "fw_observed_predictions.csv"
    work.to_csv(observed_csv, index=False)

    # Coefficients
    coeff_df = pd.DataFrame({
        "feature": x_fw.columns,
        "coefficient": model.coef_,
    })

    coeff_df["series"] = series
    coeff_df["core"] = core
    coeff_df["model_intercept"] = model.intercept_

    coeff_csv = output_dir / "fw_coefficients.csv"
    coeff_df.to_csv(coeff_csv, index=False)

    # Candidate generation
    candidate_csv = None
    candidate_df = pd.DataFrame()

    n_recon = 0
    n_final_ok = 0
    n_dummy = 0

    failure_diag = {
        "failure_csv": None,
        "failure_summary_csv": None,
        "missing_diag_csv": None,
        "value_counts_csv": None,
        "n_fail": 0,
        "n_total": 0,
    }

    reconstruction_diag_plot = None

    if generate_candidates:
        candidate_df = generate_free_wilson_candidates(
            work,
            r_cols=r_cols,
            id_col=id_col,
            max_candidates=max_candidates,
            max_changed_positions=max_changed_positions,
        )

        if not candidate_df.empty:
            x_cand = make_fw_design_matrix(candidate_df, r_cols)
            x_cand = x_cand.reindex(columns=x_fw.columns, fill_value=0)

            candidate_df["fw_predicted_activity"] = model.predict(x_cand)

            candidate_df["series"] = series
            candidate_df["core"] = core
            candidate_df["core_hash"] = hashlib.sha1(str(core).encode()).hexdigest()[:8]

            # Use the most common labelled core for this series/core group.
            candidate_df["labelled_core"] = labelled_core_smiles
            candidate_df["active_reconstruction_r_cols"] = ",".join(reconstruction_r_cols)

            candidate_df["candidate_id"] = [
                f"FWC_{hashlib.sha1('|'.join(map(str, row)).encode()).hexdigest()[:12]}"
                for row in candidate_df[r_cols].itertuples(index=False, name=None)
            ]

            default_rgroups = get_group_default_rgroups(
                work,
                reconstruction_r_cols,
            )

            def build_virtual_rgroup_assignment(row):
                """
                Combine:
                - default/core-wide R-groups, including constant positions
                - candidate variable R-groups from FW enumeration
                """
                assignment = dict(default_rgroups)

                for r in r_cols:
                    if r in row.index and pd.notna(row.get(r)):
                        assignment[r] = row.get(r)

                return assignment


            # Build complete reconstruction assignments, including constant/non-modelled R positions.
            candidate_df["_rgroup_assignment_dict"] = candidate_df.apply(
                build_virtual_rgroup_assignment,
                axis=1,
            )

            candidate_df["rgroup_assignment_for_reconstruction"] = candidate_df["_rgroup_assignment_dict"].apply(
                lambda d: json.dumps(d, sort_keys=True)
            )

            # Also expose the actual reconstruction assignment as columns for diagnostics.
            # These are the values actually used to rebuild the molecule.
            for r in reconstruction_r_cols:
                candidate_df[f"recon_{r}"] = candidate_df["_rgroup_assignment_dict"].apply(
                    lambda d: d.get(r, np.nan)
                )

            candidate_df["reconstructed_smiles"] = candidate_df["_rgroup_assignment_dict"].apply(
                lambda assignment: reconstruct_molecule_from_labelled_core(
                    labelled_core_smiles,
                    assignment,
                )
            )

            candidate_df = candidate_df.drop(columns=["_rgroup_assignment_dict"])

            candidate_df["reconstruction_success"] = candidate_df["reconstructed_smiles"].notna()

            candidate_df["has_unresolved_dummy_atoms"] = candidate_df["reconstructed_smiles"].apply(
                smiles_has_dummy_atoms
            )

            candidate_df["final_smiles"] = candidate_df["reconstructed_smiles"]

            candidate_df["remaining_attachment_labels"] = candidate_df["reconstructed_smiles"].apply(
                extract_attachment_labels_from_smiles
            )

            candidate_df["n_remaining_attachment_labels"] = candidate_df["remaining_attachment_labels"].apply(len)

            candidate_df["remaining_attachment_labels_str"] = candidate_df["remaining_attachment_labels"].apply(
                lambda xs: ",".join([f"R{x}" for x in xs]) if xs else ""
            )

            candidate_df["reconstruction_status"] = np.select(
                [
                    candidate_df["reconstruction_success"] & ~candidate_df["has_unresolved_dummy_atoms"],
                    ~candidate_df["reconstruction_success"],
                    candidate_df["has_unresolved_dummy_atoms"],
                ],
                [
                    "final_smiles_ok",
                    "reconstruction_null",
                    "unresolved_dummy_atoms",
                ],
                default="unknown",
            )

            n_total = len(candidate_df)
            n_recon = int(candidate_df["reconstruction_success"].sum())
            n_dummy = int(candidate_df["has_unresolved_dummy_atoms"].sum())
            n_final_ok = int(
                (
                    candidate_df["reconstruction_success"]
                    & ~candidate_df["has_unresolved_dummy_atoms"]
                ).sum()
            )

            core_hash = hashlib.sha1(str(core).encode()).hexdigest()[:8]

            logger.debug(
                f"[Free-Wilson] Virtual reconstruction summary for series={series}, core={core_hash}: "
                f"final_ok={n_final_ok}/{n_total}, "
                f"reconstructed_non_null={n_recon}/{n_total}, "
                f"unresolved_dummy_atoms={n_dummy}/{n_total}"
            )

            failure_diag = summarise_virtual_reconstruction_failures(
                candidate_df=candidate_df,
                reconstruction_r_cols=reconstruction_r_cols,
                output_dir=output_dir,
                series=series,
                core=core,
            )

            reconstruction_diag_plot = output_dir / "fw_virtual_reconstruction_diagnostics.png"

            plot_virtual_reconstruction_failure_diagnostics(
                candidate_df=candidate_df,
                reconstruction_r_cols=reconstruction_r_cols,
                output_path=reconstruction_diag_plot,
                series=series,
                core_hash=core_hash,
            )

            candidate_df = candidate_df.sort_values(
                "fw_predicted_activity",
                ascending=False
            ).reset_index(drop=True)

        candidate_csv = output_dir / "fw_virtual_predictions.csv"
        candidate_df.to_csv(candidate_csv, index=False)

    # Diagnostics
    pred_plot = None
    residual_plot = None
    coeff_plot = None

    if make_plots:
        pred_plot = output_dir / "fw_predicted_vs_actual.png"
        plot_fw_predicted_vs_actual(
            work,
            activity_col=activity_col,
            pred_col="fw_predicted_activity",
            output_path=pred_plot,
            title=f"Free-Wilson predicted vs actual | {series}",
        )

    if make_residual_plots:
        residual_plot = output_dir / "fw_residual_diagnostics.png"
        plot_fw_residuals(
            work,
            pred_col="fw_predicted_activity",
            residual_col="fw_residual",
            output_path=residual_plot,
            title=f"Free-Wilson residual diagnostics | {series}",
        )

    if make_coefficient_plots:
        coeff_plot = output_dir / "fw_top_coefficients.png"
        plot_fw_coefficients(
            coeff_df,
            output_path=coeff_plot,
            top_n=top_n_coefficients,
        )

    summary = {
        "series": series,
        "core": core,
        "core_hash": hashlib.sha1(str(core).encode()).hexdigest()[:8],
        "n_observed": len(work),
        "n_rgroup_columns": len(r_cols),
        "n_core_attachment_points": len(reconstruction_r_cols),
        "core_attachment_r_cols": ",".join(reconstruction_r_cols),
        "n_fw_features": x_fw.shape[1],
        "n_virtual_candidates": len(candidate_df) if candidate_df is not None else 0,
        "n_virtual_reconstructed_non_null": n_recon,
        "n_virtual_final_smiles_ok": n_final_ok,
        "n_virtual_unresolved_dummy_atoms": n_dummy,
        "n_virtual_reconstruction_failures": failure_diag.get("n_fail", 0),
        "virtual_reconstruction_failure_csv": failure_diag.get("failure_csv"),
        "virtual_reconstruction_failure_summary_csv": failure_diag.get("failure_summary_csv"),
        "virtual_missing_rgroup_diagnostics_csv": failure_diag.get("missing_diag_csv"),
        "virtual_failure_rgroup_value_counts_csv": failure_diag.get("value_counts_csv"),
        "virtual_reconstruction_diagnostics_png": str(reconstruction_diag_plot) if reconstruction_diag_plot else None,
        "model_type": model_type,
        "ridge_alpha": ridge_alpha if model_type == "ridge" else np.nan,
        "intercept": model.intercept_,
        "r2": metrics["r2"],
        "rmse": metrics["rmse"],
        "mae": metrics["mae"],
        "observed_predictions_csv": str(observed_csv),
        "virtual_predictions_csv": str(candidate_csv) if candidate_csv else None,
        "coefficients_csv": str(coeff_csv),
        "predicted_vs_actual_png": str(pred_plot) if pred_plot else None,
        "residual_diagnostics_png": str(residual_plot) if residual_plot else None,
        "coefficient_plot_png": str(coeff_plot) if coeff_plot else None,
    }

    return {
        "observed_predictions": work,
        "virtual_predictions": candidate_df,
        "coefficients": coeff_df,
        "summary": summary,
    }

def get_dummy_attachment_label(atom):
    """
    Return attachment label for a dummy atom.

    Supports:
      [*:1] via atom map number
      [1*] via isotope
    """
    if atom.GetAtomicNum() != 0:
        return None

    amap = atom.GetAtomMapNum()
    if amap:
        return int(amap)

    isotope = atom.GetIsotope()
    if isotope:
        return int(isotope)

    return None


def mol_to_clean_smiles(mol, require_no_dummy_atoms=True):
    if mol is None:
        return None

    try:
        Chem.SanitizeMol(mol)

        if require_no_dummy_atoms and mol_has_dummy_atoms(mol):
            return None

        return Chem.MolToSmiles(mol, canonical=True)

    except Exception:
        return None

def remove_dummy_atoms_and_connect(mol):
    """
    Remove mapped dummy atoms by connecting their neighbours.

    This assumes each attachment label appears exactly twice:
      once on the labelled core
      once on the R-group fragment.

    Example:
      core-[*:1] + R-[*:1] -> core-R
    """

    if mol is None:
        return None

    rw = Chem.RWMol(mol)

    label_to_dummy_idxs = {}

    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue

        label = get_dummy_attachment_label(atom)
        if label is None:
            continue

        label_to_dummy_idxs.setdefault(label, []).append(atom.GetIdx())

    # Work out which real atoms should be connected for each label
    bonds_to_add = []
    dummy_idxs_to_remove = []

    for label, dummy_idxs in label_to_dummy_idxs.items():
        if len(dummy_idxs) != 2:
            # Not enough information to connect safely
            continue

        neighbour_idxs = []

        for dummy_idx in dummy_idxs:
            dummy_atom = rw.GetAtomWithIdx(dummy_idx)
            neighbours = [n.GetIdx() for n in dummy_atom.GetNeighbors()]

            # We expect one neighbour per dummy atom
            if len(neighbours) != 1:
                neighbour_idxs = []
                break

            neighbour_idxs.append(neighbours[0])

        if len(neighbour_idxs) != 2:
            continue

        a1, a2 = neighbour_idxs

        if a1 != a2 and rw.GetBondBetweenAtoms(a1, a2) is None:
            bonds_to_add.append((a1, a2))

        dummy_idxs_to_remove.extend(dummy_idxs)

    for a1, a2 in bonds_to_add:
        rw.AddBond(a1, a2, Chem.BondType.SINGLE)

    # Remove dummy atoms from highest index to lowest
    for idx in sorted(set(dummy_idxs_to_remove), reverse=True):
        rw.RemoveAtom(idx)

    out = rw.GetMol()

    try:
        Chem.SanitizeMol(out)

        # If any dummy atoms remain, reconstruction is incomplete.
        if mol_has_dummy_atoms(out):
            return None

        return out

    except Exception:
        return None


def reconstruct_molecule_from_labelled_core(
    labelled_core_smiles,
    rgroup_assignment,
):
    """
    Reconstruct full molecule from labelled core and R-group fragments.

    Parameters
    ----------
    labelled_core_smiles : str
        Core with mapped dummy atoms, e.g. c1cc([*:1])ccc([*:2])c1

    rgroup_assignment : dict
        Example:
          {
            "R1": "CO[*:1]",
            "R2": "C[*:2]",
            "R3": "[H][*:3]"
          }

    Returns
    -------
    canonical_smiles : str or None
    """

    if labelled_core_smiles is None or pd.isna(labelled_core_smiles):
        return None

    mols = []

    core_mol = Chem.MolFromSmiles(str(labelled_core_smiles))
    if core_mol is None:
        return None

    mols.append(core_mol)

    for _, frag_smiles in sorted(rgroup_assignment.items()):
        if frag_smiles is None or pd.isna(frag_smiles):
            continue

        frag_smiles = str(frag_smiles)

        if frag_smiles == "__MISSING__":
            continue

        frag_mol = Chem.MolFromSmiles(frag_smiles)
        if frag_mol is None:
            continue

        mols.append(frag_mol)

    if len(mols) == 1:
        return mol_to_clean_smiles(core_mol)

    combo = mols[0]
    for frag in mols[1:]:
        combo = Chem.CombineMols(combo, frag)

    reconstructed = remove_dummy_atoms_and_connect(combo)

    return mol_to_clean_smiles(reconstructed)

def validate_observed_reconstruction(
    matrix_df,
    activity_col,
    output_dir,
    smiles_col="smiles",
):
    """
    Reconstruct observed molecules from labelled_core + R-groups
    and compare against original SMILES.

    Outputs:
      reconstruction_validation.csv
      reconstruction_summary.csv
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if smiles_col not in matrix_df.columns:
        logger.warning(
            f"[Free-Wilson] Cannot validate reconstruction: missing '{smiles_col}' column."
        )
        return None, None

    r_cols = get_rgroup_columns(matrix_df)

    records = []

    for _, row in matrix_df.iterrows():
        labelled_core = row.get("labelled_core")

        rgroup_assignment = {
            r: row.get(r)
            for r in r_cols
            if r in matrix_df.columns
        }

        reconstructed_smiles = reconstruct_molecule_from_labelled_core(
            labelled_core,
            rgroup_assignment,
        )

        original_smiles = row.get(smiles_col)

        original_can = None
        if original_smiles is not None and pd.notna(original_smiles):
            mol = Chem.MolFromSmiles(str(original_smiles))
            if mol is not None:
                original_can = Chem.MolToSmiles(mol, canonical=True)

        match = (
            original_can is not None
            and reconstructed_smiles is not None
            and original_can == reconstructed_smiles
        )

        records.append({
            "compound_id": row.get("compound_id"),
            "series": row.get("series"),
            "core": row.get("core"),
            "original_smiles": original_smiles,
            "original_canonical_smiles": original_can,
            "reconstructed_smiles": reconstructed_smiles,
            "reconstruction_success": reconstructed_smiles is not None,
            "canonical_match": match,
            activity_col: row.get(activity_col),
        })

    validation_df = pd.DataFrame(records)

    validation_csv = output_dir / "reconstruction_validation.csv"
    validation_df.to_csv(validation_csv, index=False)

    n_total = len(validation_df)
    n_success = int(validation_df["reconstruction_success"].sum()) if n_total else 0
    n_match = int(validation_df["canonical_match"].sum()) if n_total else 0

    summary_df = pd.DataFrame([{
        "n_total": n_total,
        "n_reconstruction_success": n_success,
        "n_canonical_match": n_match,
        "reconstruction_success_rate": n_success / n_total if n_total else np.nan,
        "canonical_match_rate": n_match / n_total if n_total else np.nan,
    }])

    summary_csv = output_dir / "reconstruction_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    logger.info(
        f"[Free-Wilson] Reconstruction validation: "
        f"success={n_success}/{n_total}, canonical_match={n_match}/{n_total}"
    )

    return validation_csv, summary_csv
