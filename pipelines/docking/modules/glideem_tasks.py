from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr
import matplotlib.pyplot as plt
from schrodinger import structure

from pipeline.logger import setup_logger
from pipeline.task_registry import register_task

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","HID","HIE","HIP","ASH","GLH","LYN","CYX"}
IGNORE = {"HOH","WAT","SOL","DOD","NA","K","CL","CA","MG","ZN","MN","FE","SO4","PO4","GOL","EDO","PEG"}
ID_FIELDS = ("s_user_compound_id", "Compound ID", "Compound_ID", "ID", "Name", "s_m_title")


def norm(value):
    return re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")

def property_value(st, requested):
    target = norm(requested)
    for key, value in st.property.items():
        nk = norm(key)
        if nk == target or nk.endswith("_" + target) or target in nk:
            return value
    return None


def compound_id(st):
    for field in ID_FIELDS:
        value = property_value(st, field)
        if value is not None and str(value).strip():
            return str(value).strip()
    return st.title.strip()


def glide_score(st):
    props = {norm(k): v for k, v in st.property.items()}
    preferred = ("R_I_GLIDE_GSCORE", "R_I_DOCKING_SCORE", "R_GLIDE_GSCORE", "R_I_GLIDE_DOCKING_SCORE")
    for key in preferred:
        try:
            if key in props:
                return float(props[key])
        except (TypeError, ValueError):
            pass
    for key, value in props.items():
        if "GLIDE" in key and ("GSCORE" in key or "DOCKING_SCORE" in key):
            try:
                return float(value)
            except (TypeError, ValueError):
                pass
    return None


def is_protein_atom(atom):
    try:
        return bool(atom.isProtein())
    except Exception:
        return (atom.pdbres or "").strip().upper() in AA


def heavy_atoms(st, indices):
    return sum(st.atom[i].atomic_number > 1 for i in indices)


@register_task("prepare_glideem_receptor", category="GlideEM", description="Prepare a receptor using Protein Preparation Wizard")
def prepare_glideem_receptor(backend, ligand, config, **kwargs):
    job = kwargs["job"]
    work = Path(kwargs["work_dir"])
    source = Path(job["receptor"]).expanduser().resolve()
    output = work / f"{job['job_name']}_prepared.maegz"
    if output.exists() and not config.get("force", False):
        return output
    prep = config["glideem"].get("prepare_receptor", {})
    cmd = [backend.exe("prepwizard"), source, output, "-WAIT", "-HOST", backend.host,
           "-fillsidechains", "-disulfides", "-rehtreat", "-propka_pH", str(prep.get("pH", 7.4)),
           "-rmsd", str(prep.get("rmsd", 0.3))]
    backend.run(cmd, work, work / f"{job['job_name']}_prepwizard.log")
    if not output.exists():
        raise RuntimeError(f"PrepWizard did not create {output}")
    return output


@register_task("run_glideem", category="GlideEM", description="Run configured GlideEM command")
def run_glideem(backend, ligand, config, **kwargs):
    job = kwargs["job"]
    work = Path(kwargs["work_dir"])
    receptor = Path(kwargs["prepared_receptor"])
    output = work / f"{job['job_name']}_glideem_complexes.maegz"
    values = {
        "schrodinger": backend.schrodinger_root,
        "host": backend.host,
        "job_name": job["job_name"],
        "receptor": receptor,
        "density_map": Path(job["density_map"]).expanduser().resolve(),
        "ligand": Path(job["reference_ligand"]).expanduser().resolve(),
        "output": output,
    }
    if output.exists() and not config.get("force", False):
        return output
    return backend.run_glideem_template(values, config["glideem"], work)


def extract_complexes(archive: Path, out_root: Path, source_name: str, force: bool = False):
    jobs = []
    pose = 0
    for ct_number, st in enumerate(structure.StructureReader(str(archive)), 1):
        protein = [a.index for a in st.atom if is_protein_atom(a)]
        if not protein:
            logger.warning("%s CT %d has no protein; skipped", archive.name, ct_number)
            continue
        residues = defaultdict(list)
        for atom in st.atom:
            if is_protein_atom(atom):
                continue
            rn = (atom.pdbres or "").strip().upper()
            if rn not in IGNORE:
                residues[((atom.chain or "").strip(), int(atom.resnum), rn)].append(atom.index)
        candidates = [(heavy_atoms(st, idx), key, idx) for key, idx in residues.items() if heavy_atoms(st, idx) >= 6]
        if not candidates:
            logger.warning("%s CT %d has no ligand; skipped", archive.name, ct_number)
            continue
        pose += 1
        name = f"{source_name}_pose{pose:02d}"
        directory = out_root / name
        directory.mkdir(parents=True, exist_ok=True)
        receptor = directory / f"{name}_receptor.maegz"
        reference = directory / f"{name}_reference_ligand.maegz"
        complex_file = directory / f"{name}_complex.maegz"
        grid = directory / f"{name}_grid.zip"
        metadata = directory / "complex_metadata.json"
        nheavy, residue, ligand_indices = max(candidates, key=lambda item: item[0])
        xyz = np.array([[st.atom[i].x, st.atom[i].y, st.atom[i].z] for i in ligand_indices], dtype=float)
        if force or not all(p.exists() for p in (receptor, reference, complex_file, metadata)):
            lig = st.extract(ligand_indices)
            lig.title = name
            lig.property["s_user_reference_complex"] = name
            with structure.StructureWriter(str(receptor)) as writer:
                writer.append(st.extract(protein))
            with structure.StructureWriter(str(reference)) as writer:
                writer.append(lig)
            with structure.StructureWriter(str(complex_file)) as writer:
                writer.append(st)
            metadata.write_text(json.dumps({
                "job_name": name, "source_archive": str(archive.resolve()), "source_CT": ct_number,
                "ligand_chain": residue[0], "ligand_resnum": residue[1], "ligand_resname": residue[2],
                "ligand_heavy_atoms": nheavy, "grid_center": xyz.mean(axis=0).tolist()
            }, indent=2))
        meta = json.loads(metadata.read_text())
        jobs.append({"name": name, "dir": directory, "receptor": receptor, "reference": reference,
                     "complex": complex_file, "grid": grid, "center": meta["grid_center"]})
    return jobs


@register_task("build_glide_grid", category="Glide", description="Build a Glide grid around a GlideEM reference pose")
def build_glide_grid(backend, ligand, config, **kwargs):
    job = kwargs["complex_job"]
    grid_cfg = config["grid"]
    if job["grid"].exists() and not config.get("force", False):
        return job["grid"]
    x, y, z = job["center"]
    inp = job["dir"] / f"{job['name']}_grid.in"
    inp.write_text(
        f"GRIDFILE {job['grid']}\nRECEP_FILE {job['receptor']}\n"
        f"GRID_CENTER {x:.4f},{y:.4f},{z:.4f}\n"
        f"INNERBOX {grid_cfg.get('innerbox', 10)},{grid_cfg.get('innerbox', 10)},{grid_cfg.get('innerbox', 10)}\n"
        f"OUTERBOX {grid_cfg.get('outerbox', 30)},{grid_cfg.get('outerbox', 30)},{grid_cfg.get('outerbox', 30)}\n"
    )
    backend.run_glide(inp, job["dir"], job["dir"] / f"{job['name']}_grid.log")
    if not job["grid"].exists():
        raise RuntimeError(f"Grid missing: {job['grid']}")
    return job["grid"]


@register_task("dock_mcss", category="Glide", description="Core-constrained MCSS docking using GlideEM pose")
def dock_mcss(backend, ligand, config, **kwargs):
    job = kwargs["complex_job"]
    ligand_file = Path(kwargs["ligand_file"])
    out_dir = Path(kwargs["mcss_root"]) / job["name"]
    out_dir.mkdir(parents=True, exist_ok=True)
    inp = out_dir / f"{job['name']}_mcss.in"
    pv = out_dir / f"{job['name']}_mcss_pv.maegz"
    mcss = config.get("mcss", {})
    inp.write_text(
        "CORECONS_FALLBACK True\nCORE_DEFINITION mcssmarts\nCORE_FILTER True\nCORE_RESTRAIN True\n"
        f"GRIDFILE {job['grid']}\nLIGANDFILE {ligand_file}\n"
        f"NREPORT {mcss.get('nreport', 10000)}\nPOSES_PER_LIG {mcss.get('poses_per_ligand', 1)}\n"
        f"POSTDOCK_NPOSE {mcss.get('postdock_npose', 1)}\nPRECISION {mcss.get('precision', 'SP')}\n"
        f"REF_LIGAND_FILE {job['reference']}\nUSE_REF_LIGAND True\n"
    )
    if not pv.exists() or config.get("force", False):
        backend.run_glide(inp, out_dir, out_dir / f"{job['name']}_mcss.log")
    if not pv.exists():
        raise RuntimeError(f"MCSS PoseViewer missing: {pv}")
    return pv


def analyse_correlations(complex_jobs, poseviewers, ligand_file, activity_column, output_dir, min_pairs=5):
    activities = {}
    for st in structure.StructureReader(str(ligand_file)):
        try:
            value = float(property_value(st, activity_column))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            activities.setdefault(compound_id(st), value)
    if not activities:
        raise RuntimeError(f"No numeric values found for activity property '{activity_column}'")

    by_name = {j["name"]: j for j in complex_jobs}
    details = []
    for name, pv in poseviewers.items():
        best = {}
        for st in structure.StructureReader(str(pv)):
            cid = compound_id(st)
            score = glide_score(st)
            if cid in activities and score is not None and (cid not in best or score < best[cid]):
                best[cid] = score
        for cid, score in best.items():
            details.append({"complex_pose": name, "compound_id": cid, "glide_score": score,
                            "activity": activities[cid], "activity_column": activity_column,
                            "poseviewer_path": str(Path(pv).resolve()),
                            "reference_ligand_path": str(by_name[name]["reference"].resolve()),
                            "grid_path": str(by_name[name]["grid"].resolve())})
    if not details:
        raise RuntimeError("No matched docking score/activity pairs")

    output_dir.mkdir(parents=True, exist_ok=True)

    detail_df = pd.DataFrame(details)
    detail_df.to_csv(output_dir / "matched_scores_and_activities.csv", index=False)
    summaries = []
    for name, group in detail_df.groupby("complex_pose"):
        if len(group) < min_pairs:
            logger.warning("%s has only %d matched pairs; excluded", name, len(group))
            continue
        x = group["glide_score"].to_numpy(float)
        y = group["activity"].to_numpy(float)
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
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(group["glide_score"], group["activity"])
        ax.set(title=name, xlabel="Glide docking score", ylabel=activity_column)
        fig.tight_layout(); fig.savefig(output_dir / f"{name}_score_vs_activity.png", dpi=200); plt.close(fig)
    best = summary.iloc[0]
    (output_dir / "best_binding_pose.txt").write_text(
        f"Best complex pose: {best.complex_pose}\nN: {best.n}\nPearson r: {best.pearson_r:.6f}\n"
        f"Spearman rho: {best.spearman_rho:.6f}\nKendall tau: {best.kendall_tau:.6f}\n"
    )
    return summary
