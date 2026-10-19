import re
import math

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import kendalltau, pearsonr, spearmanr

from schrodinger import structure

from modules.utils.glideem_helpers import (
    glide_smiles_by_lignum,
    glide_score,
    glide_score_components,
    property_value,
    compound_id,
    smiles_value,
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def analyse_correlations(
    complex_jobs,
    poseviewers,
    ligand_file,
    source_ligand_file,
    activity_column,
    output_dir,
    min_pairs=5,
):
    activities = {}
    ligand_metadata = {}

    for st in structure.StructureReader(
        str(ligand_file)
    ):
        cid = compound_id(st)

        if not cid:
            continue

        try:
            value = float(
                property_value(
                    st,
                    activity_column,
                )
            )
        except (TypeError, ValueError):
            continue

        if not math.isfinite(value):
            continue

        activities.setdefault(
            cid,
            value,
        )

        ligand_metadata.setdefault(
            cid,
            {
                "smiles": smiles_value(st),
            },
        )
    if not activities:
        raise RuntimeError(f"No numeric values found for activity property '{activity_column}'")

    docked_ids = set()
    scored_ids = set()
    matched_ids = set()

    by_name = {
        j["name"]: j
        for j in complex_jobs
    }

    details = []

    docked_ids = set()
    scored_ids = set()
    matched_ids = set()

    for name, pv in poseviewers.items():

        best = {}

        smiles_by_lignum = glide_smiles_by_lignum(
            pv
        )

        for st in structure.StructureReader(
            str(pv)
        ):
            cid = compound_id(st)
            score = glide_score(st)

            lignum = st.property.get(
                "i_i_glide_lignum"
            )

            try:
                lignum = int(lignum)
            except (TypeError, ValueError):
                lignum = None

            if cid:
                docked_ids.add(cid)

            if cid and score is not None:
                scored_ids.add(cid)

            if (
                cid in activities
                and score is not None
                and (
                    cid not in best
                    or score < best[cid]["score"]
                )
            ):
                best[cid] = {
                    "score": score,
                    "structure": st,
                    "lignum": lignum,
                }

                matched_ids.add(cid)

        logger.debug(
            "[%s] Docked IDs with scores: %d; "
            "matched to activity: %d",
            name,
            len(best),
            sum(
                cid in activities
                for cid in best
            ),
        )

        for cid, result in best.items():

            score = result["score"]
            docked_st = result["structure"]

            lignum = result["lignum"]

            smiles = (
                smiles_by_lignum.get(lignum)
                if lignum is not None
                else None
            )

            components = glide_score_components(
                docked_st
            )

            metadata = ligand_metadata.get(
                cid,
                {}
            )

            details.append(
                {
                    "complex_pose": name,
                    "source_name": by_name[name]["source_name"],
                    "pose_number": by_name[name]["pose_number"],

                    "compound_id": cid,
                    "smiles": smiles,

                    "experimental_potency": activities[cid],
                    "activity_column": activity_column,

                    "glide_score": score,

                    **components,

                    "poseviewer_path": str(
                        Path(pv).resolve()
                    ),

                    "reference_ligand_path": str(
                        by_name[name][
                            "reference"
                        ].resolve()
                    ),

                    "grid_path": str(
                        by_name[name][
                            "grid"
                        ].resolve()
                    ),
                }
            )


    if not details:
        raise RuntimeError(
            "No matched docking score/activity pairs. "
            f"Activity IDs={len(activities)}, "
            f"docked IDs={len(docked_ids)}, "
            f"scored IDs={len(scored_ids)}, "
            f"matched IDs={len(matched_ids)}."
        )
    
    # if not details:
    #     raise RuntimeError("No matched docking score/activity pairs")

    output_dir.mkdir(parents=True, exist_ok=True)

    detail_df = pd.DataFrame(details)
    detail_df.to_csv(output_dir / "matched_scores_and_activities.csv", index=False)
    summaries = []
    for name, group in detail_df.groupby("complex_pose"):
        if len(group) < min_pairs:
            logger.warning("%s has only %d matched pairs; excluded", name, len(group))
            continue
        x = group["glide_score"].to_numpy(float)
        y = group["experimental_potency"].to_numpy(float)
        p, s, k = pearsonr(x, y), spearmanr(x, y), kendalltau(x, y)
        # pActivity should normally anticorrelate with docking score, hence larger -(rP+rS) is better.
        summaries.append({"complex_pose": name, "n": len(group), "pearson_r": p.statistic, "pearson_p": p.pvalue,
                          "spearman_rho": s.statistic, "spearman_p": s.pvalue, "kendall_tau": k.statistic,
                          "kendall_p": k.pvalue, "rank_score": -(p.statistic + s.statistic)})
    summary = pd.DataFrame(summaries)
    if summary.empty:
        raise RuntimeError("No complex had enough matched pairs for correlation analysis")
    summary = summary.sort_values(["rank_score", "spearman_rho", "pearson_r"], ascending=[False, True, True]).reset_index(drop=True)
    summary.insert(0, "rank", np.arange(1, len(summary) + 1))
    summary.to_csv(output_dir / "complex_pose_correlations.csv", index=False)

    for name, group in detail_df.groupby("complex_pose"):

        if len(group) < min_pairs:
            logger.warning(
                "%s has only %d matched pairs; "
                "correlation plot skipped",
                name,
                len(group),
            )
            continue

        source_name = group[
                "source_name"
            ].iloc[0]

        pose_number = int(
            group[
                "pose_number"
            ].iloc[0]
        )

        plot_dir = (
            output_dir /
            source_name
        )

        plot_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        pose_table = group[
            [
                "compound_id",
                "smiles",
                "experimental_potency",
                "glide_score",
                "docking_score",
                "glide_emodel",
                "glide_energy",
                "glide_evdw",
                "glide_ecoul",
                "glide_lipo",
                "glide_hbond",
                "glide_erotb",
                "glide_esite",
                "glide_einternal",
            ]
        ].copy()

        pose_table = pose_table.sort_values(
            "experimental_potency",
            ascending=False,
        )

        pose_table.to_csv(
            plot_dir /
            (
                f"pose{pose_number:02d}"
                "_docking_vs_affinity.csv"
            ),
            index=False,
        )

        # Conventional experimental vs predicted orientation:
        # experimental values on X, predicted values on Y.
        x = group[
            "experimental_potency"
        ].to_numpy(float)

        y = group[
            "glide_score"
        ].to_numpy(float)

        # ---------------------------------------------
        # Correlation statistics
        # ---------------------------------------------

        pearson = pearsonr(x, y)
        spearman = spearmanr(x, y)
        kendall = kendalltau(x, y)

        # ---------------------------------------------
        # Linear regression / best-fit line
        # ---------------------------------------------

        slope, intercept = np.polyfit(
            x,
            y,
            1,
        )

        x_fit = np.linspace(
            x.min(),
            x.max(),
            200,
        )

        y_fit = (
            slope * x_fit
            + intercept
        )

        # ---------------------------------------------
        # Plot
        # ---------------------------------------------

        fig, ax = plt.subplots(
            figsize=(7, 6)
        )

        ax.scatter(
            x,
            y,
            alpha=0.8,
            label="Compounds",
        )

        ax.plot(
            x_fit,
            y_fit,
            linewidth=2,
            label="Linear fit",
        )

        ax.set_title(name)

        ax.set_xlabel(
            activity_column
        )

        ax.set_ylabel(
            "Glide score"
        )

        # ---------------------------------------------
        # Statistics annotation
        # ---------------------------------------------

        stats_text = (
            f"N = {len(group)}\n"
            f"Pearson r = {pearson.statistic:.3f} "
            f"(p = {pearson.pvalue:.2g})\n"
            f"Spearman ρ = {spearman.statistic:.3f} "
            f"(p = {spearman.pvalue:.2g})\n"
            f"Kendall τ = {kendall.statistic:.3f} "
            f"(p = {kendall.pvalue:.2g})"
        )

        ax.text(
            0.03,
            0.97,
            stats_text,
            transform=ax.transAxes,
            verticalalignment="top",
            horizontalalignment="left",
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.85,
            },
        )

        ax.legend()

        fig.tight_layout()

        fig.savefig(
            plot_dir /
            (
                f"pose{pose_number:02d}"
                "_experimental_vs_predicted.png"
            ),
            dpi=200,
        )

        plt.close(fig)
    best = summary.iloc[0]
    (output_dir / "best_binding_pose.txt").write_text(
        f"Best complex pose: {best.complex_pose}\nN: {best.n}\nPearson r: {best.pearson_r:.6f}\n"
        f"Spearman rho: {best.spearman_rho:.6f}\nKendall tau: {best.kendall_tau:.6f}\n"
    )
    return summary