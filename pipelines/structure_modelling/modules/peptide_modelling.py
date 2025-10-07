from pipeline.task_registry import register_task

import io
import os
import math
import collections
import csv
import logging
import shutil

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeForceField, MMFFGetMoleculeProperties

from modules.utils.utils import calc_rmsd, is_rotor_bond
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

@register_task("build_peptide", category="Peptide modelling", description="Build peptides from sequence file")
def build_peptide(backend, config, **kwargs):
    logger.info("Building peptides...")
    results = backend.run()
    return results

@register_task("generate_peptide_conformers", 
               category="Peptide modelling", 
               description="Generate and minimise (MMFF) peptide conformers.")
def generate_conformers(backend, config, **kwargs):
    conf_dir = Path(config["output_dir"]) / "conformers"
    conf_dir.mkdir(parents=True, exist_ok=True)

    num_confs = config.get("generate_peptide_conformers", {}).get("number_of_conformers", 100)

    results = {}
    energy_map = {}

    for name, paths in backend.run().items():
        pdb_file = paths["capped_pdb"]
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            logger.warning(f"Failed to load molecule from {pdb_file}")
            continue

        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.maxAttempts = 1000
        params.pruneRmsThresh = -1.0
        params.numThreads = 10

        ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
        logger.info(f"[{name}] Generated {len(ids)} conformers.")

        conformer_paths = []
        sdf_writer = Chem.SDWriter(str(conf_dir / f"{name}.sdf"))

        for i in tqdm(ids, desc=f"[{name}] Optimising conformers", unit="conf"):
            try:
                mp = MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
                ff = MMFFGetMoleculeForceField(mol, mp, confId=i)
                ff.Minimize()
                energy = ff.CalcEnergy()

                filename = f"{name}_conf_{i}.xyz"
                path = conf_dir / filename
                Chem.MolToXYZFile(mol, str(path), confId=i)
                conformer_paths.append(str(path))

                mol.SetProp("_Name", f"{name}_conf_{i}")
                mol.SetProp("Energy", f"{energy:.4f}")
                sdf_writer.write(mol, confId=i)

                energy_map[filename] = energy
            except Exception as e:
                logger.warning(f"Skipping conformer {i} of {name}: {e}")
                continue

        sdf_writer.close()
        results[name] = conformer_paths

    energy_file = Path(config["output_dir"]) / "mmff_energies.txt"
    with open(energy_file, "w") as f:
        f.write("Conformer\tEnergy\n")
        for name, energy in energy_map.items():
            f.write(f"{name}\t{energy:.4f}\n")

    return energy_map


@register_task("analyse_peptide_flexibility", 
               category="Peptide modelling", 
               description="Compute peptide flexibility metrics.")
def analyse_flexibility(backend, config, **kwargs):

    logger = logging.getLogger(__name__)

    kT = 0.592  # kcal/mol at RT
    energy_map = kwargs.get("energy_map", {})
    if not energy_map:
        raise ValueError("Missing energy map for analysis.")

    output_dir = Path(config["output_dir"])
    grouped = collections.defaultdict(list)

    for conf_name, energy in energy_map.items():
        peptide_name = conf_name.split("_conf_")[0]
        grouped[peptide_name].append((conf_name, energy))

    def get_dihedral_indices(mol):
        """Identify rotatable dihedrals."""
        dihedrals = []
        for bond in mol.GetBonds():
            if not is_rotor_bond(bond):
                continue
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

            begin_atom = mol.GetAtomWithIdx(begin_idx)
            end_atom = mol.GetAtomWithIdx(end_idx)

            neighbors_begin = [nbr.GetIdx() for nbr in begin_atom.GetNeighbors() if nbr.GetIdx() != end_idx]
            neighbors_end = [nbr.GetIdx() for nbr in end_atom.GetNeighbors() if nbr.GetIdx() != begin_idx]

            if not neighbors_begin or not neighbors_end:
                continue

            i = neighbors_begin[0]
            j = begin_idx
            k = end_idx
            l = neighbors_end[0]

            dihedrals.append((i, j, k, l))
        return dihedrals

    summary = {}

    for peptide_name, items in grouped.items():
        conformer_names, energies = zip(*items)
        energies = np.array(energies)
        min_energy = np.min(energies)
        rel_energies = energies - min_energy
        boltzmann_weights = np.exp(-rel_energies / kT)
        boltzmann_weights /= boltzmann_weights.sum()

        entropy = -np.sum(boltzmann_weights * np.log(boltzmann_weights + 1e-12))

        logger.info(f"{peptide_name} flexibility (Boltzmann entropy): {entropy:.3f} kcal/mol")

        # Load conformers from SDF (annoying Rdkit workaround)
        sdf_path = output_dir / "conformers" / f"{peptide_name}.sdf"
        if not sdf_path.exists():
            logger.warning(f"SDF file not found for {peptide_name}")
            continue

        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            logger.warning(f"No valid molecules found in {sdf_path}")
            continue

        min_index = np.argmin(energies)
        ref_mol = mols[min_index]

        rmsds = []
        for mol in mols:
            try:
                rmsd = rdMolAlign.AlignMol(mol, ref_mol, maxIters=100)
            except Exception as e:
                logger.warning(f"AlignMol failed for peptide {peptide_name}: {e}")
                rmsd = calc_rmsd(mol, ref_mol)  # fallback to manual RMSD
            rmsds.append(rmsd)
        rmsds = np.array(rmsds)
        mean_rmsd = np.average(rmsds, weights=boltzmann_weights)

        coords = np.array([mol.GetConformer().GetPositions() for mol in mols])  # (n_confs, n_atoms, 3)
        flat_coords = coords.reshape((len(coords), -1))  # (n_confs, n_atoms * 3)
        flat_coords -= flat_coords.mean(axis=0)
        cov = np.cov(flat_coords, rowvar=False)
        eigvals = np.linalg.eigvalsh(cov)[::-1]
        total_var = np.sum(eigvals)
        pca_var_explained = np.sum(eigvals[:3]) / total_var if total_var > 0 else 0.0

        # Dihedral angle variance (circular variance weighted by Boltzmann weights)
        dihedral_indices = get_dihedral_indices(ref_mol)
        if dihedral_indices:
            dihedral_angles = []
            for mol in mols:
                conf = mol.GetConformer()
                angles = [rdMolTransforms.GetDihedralDeg(conf, *indices) for indices in dihedral_indices]
                dihedral_angles.append(angles)

            dihedral_angles = np.array(dihedral_angles)  
            angles_rad = np.deg2rad(dihedral_angles)

            sin_vals = np.sin(angles_rad)
            cos_vals = np.cos(angles_rad)

            mean_sin = np.average(sin_vals, axis=0, weights=boltzmann_weights)
            mean_cos = np.average(cos_vals, axis=0, weights=boltzmann_weights)

            R = np.sqrt(mean_sin**2 + mean_cos**2)
            circular_var = 1 - R

            dihedral_var = float(np.mean(circular_var))
        else:
            dihedral_var = float("nan")

        out_file = output_dir / f"{peptide_name}_boltzmann_weights.txt"
        np.savetxt(out_file,
                   np.column_stack((conformer_names, energies, boltzmann_weights)),
                   fmt="%s", header="Conformer\tEnergy\tBoltzmann_Weight")

        summary[peptide_name] = {
            "entropy": entropy,
            "normalised_entropy": 0.0,  
            "weights": dict(zip(conformer_names, boltzmann_weights.tolist())),
            "num_confs": len(conformer_names),
            "mean_rmsd": mean_rmsd,
            "pca_var_explained": pca_var_explained,
            "dihedral_var": dihedral_var,
        }

    entropies = [v["entropy"] for v in summary.values()]
    max_entropy = max(entropies) if entropies else 1.0

    csv_path = output_dir / "peptide_flexibility_summary.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Peptide", "Entropy (kcal/mol)", "Normalised Entropy", "Num Conformers",
            "Mean RMSD", "PCA Variance (Top 3)", "Dihedral Angle Variance"
        ])
        for pep, data in summary.items():
            norm_entropy = data["entropy"] / max_entropy if max_entropy else 0.0
            data["normalised_entropy"] = norm_entropy
            writer.writerow([
                pep,
                f"{data['entropy']:.4f}",
                f"{norm_entropy:.4f}",
                data["num_confs"],
                f"{data['mean_rmsd']:.4f}",
                f"{data['pca_var_explained']:.4f}",
                f"{data['dihedral_var']:.4f}"
            ])

    logger.info(f"Saved flexibility summary: {csv_path}")
    
    # cleanup conformers directory
    conf_dir = Path(config["output_dir"]) / "conformers"
    if config.get("cleanup", False) and conf_dir.exists():
        try:
            shutil.rmtree(conf_dir)
            logger.info(f"Cleaned up conformers directory: {conf_dir}")
        except Exception as e:
            logger.warning(f"Failed to clean up conformers directory: {e}")
        
    return summary

@register_task("plot_peptide_flexibility", category="Peptide modelling", description="Plot flexibility data")
def plot_flexibility(backend, config, **kwargs):
    output_dir = Path(config["output_dir"])
    txt_files = list(output_dir.glob("*_boltzmann_weights.txt"))

    if not txt_files:
        raise FileNotFoundError("No Boltzmann weight files found.")

    for txt_file in txt_files:
        peptide_name = txt_file.stem.replace("_boltzmann_weights", "")
        data = np.loadtxt(txt_file, skiprows=1, dtype=str)

        conformers = data[:, 0]
        energies = data[:, 1].astype(float)
        weights = data[:, 2].astype(float)

        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.bar(conformers, energies, color='tab:blue')
        ax1.set_ylabel("Energy (kcal/mol)")
        ax1.set_xlabel("Conformer")
        ax1.set_title(f"Flexibility for {peptide_name}")
        ax1.set_xticklabels(conformers, rotation=90)

        ax2 = ax1.twinx()
        ax2.plot(conformers, weights, 'r--', marker='o', label='Boltzmann weight')
        ax2.set_ylabel("Boltzmann Weight")

        fig.tight_layout()
        plot_path = output_dir / f"{peptide_name}_flexibility_plot.png"
        plt.savefig(plot_path, dpi=300)
        plt.close()

        logger.info(f"Saved flexibility plot for {peptide_name}: {plot_path}")

