from pipeline.task_registry import register_task

import os
import math

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Polypeptide import is_aa
from rdkit.Chem import rdDistGeom
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeForceField, MMFFGetMoleculeProperties

from modules.utils.utils import get_backend_from_config

import logging

logger = logging.getLogger(__name__)

@register_task("build_peptide", category="Peptide modelling", description="Build peptides from sequence file")
def build_peptide(backend, config, **kwargs):
    logger.info("Building peptides...")
    results = backend.run()
    return results

@register_task("generate_peptide_conformers", category="Peptide modelling", description="Generate peptide conformers and MMFF energies")
def generate_conformers(backend, config, **kwargs):

    conf_dir = Path(config["output_dir"]) / "conformers"
    conf_dir.mkdir(parents=True, exist_ok=True)

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

        ids = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)
        logger.info(f"[{name}] Generated {len(ids)} conformers.")

        conformer_paths = []

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

                energy_map[filename] = energy
            except Exception as e:
                logger.warning(f"Skipping conformer {i} of {name}: {e}")
                continue

        results[name] = conformer_paths

    # Save energies for debugging
    energy_file = Path(config["output_dir"]) / "mmff_energies.txt"
    with open(energy_file, "w") as f:
        f.write("Conformer\tEnergy\n")
        for name, energy in energy_map.items():
            f.write(f"{name}\t{energy:.4f}\n")

    return energy_map

# @register_task("optimise_peptide_conformers", category="Peptide modelling", description="Optimise conformers using xtb")
# def optimise_conformers(backend, config, **kwargs):
#     xtb_backend = get_backend_from_config("optimise_conformers", config)

#     conf_dir = Path(config["output_dir"]) / "conformers"
#     opt_dir = Path(config["output_dir"]) / "optimised"
#     opt_dir.mkdir(parents=True, exist_ok=True)

#     energy_map = {}

#     for file in os.listdir(conf_dir):
#         if not file.endswith(".xyz"):
#             continue

#         xyz_path = conf_dir / file
#         try:
#             optimised_xyz, energy = xtb_backend.optimise(str(xyz_path), config.get("xtb", {}))

#             base_name = file.replace(".xyz", "")
#             output_path = opt_dir / f"{base_name}_opt.xyz"
#             with open(output_path, "w") as f:
#                 f.write(optimised_xyz)

#             if energy is not None:
#                 energy_map[file] = energy
#             else:
#                 logger.warning(f"Could not extract energy for {file}")

#         except Exception as e:
#             logger.warning(f"Skipping {file}: {e}")
#             continue

#     return energy_map


@register_task("analyse_peptide_flexibility", category="Peptide modelling", description="Compute peptide flexibility metrics")
def analyze_flexibility(backend, config, **kwargs):
    import collections
    import csv

    kT = 0.592  # kcal/mol at room temp (~298K)
    energy_map = kwargs.get("energy_map", {})
    if not energy_map:
        raise ValueError("Missing energy map for analysis.")

    output_dir = Path(config["output_dir"])
    grouped = collections.defaultdict(list)

    # Group conformers by peptide
    for conf_name, energy in energy_map.items():
        peptide_name = conf_name.split("_conf_")[0]
        grouped[peptide_name].append((conf_name, energy))

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

        # Save individual Boltzmann data
        out_file = output_dir / f"{peptide_name}_boltzmann_weights.txt"
        np.savetxt(out_file,
                   np.column_stack((conformer_names, energies, boltzmann_weights)),
                   fmt="%s", header="Conformer\tEnergy\tBoltzmann_Weight")

        summary[peptide_name] = {
            "entropy": entropy,
            "weights": dict(zip(conformer_names, boltzmann_weights.tolist())),
            "num_confs": len(conformer_names)
        }

    # Compute normalized entropies
    entropies = [v["entropy"] for v in summary.values()]
    max_entropy = max(entropies) if entropies else 1.0

    csv_path = output_dir / "peptide_flexibility_summary.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Peptide", "Entropy (kcal/mol)", "Normalized Entropy", "Num Conformers"])
        for pep, data in summary.items():
            norm_entropy = data["entropy"] / max_entropy if max_entropy else 0.0
            writer.writerow([pep, f"{data['entropy']:.4f}", f"{norm_entropy:.4f}", data["num_confs"]])
            summary[pep]["normalized_entropy"] = norm_entropy

    logger.info(f"Saved flexibility summary: {csv_path}")
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

