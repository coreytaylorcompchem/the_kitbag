import subprocess
import tempfile
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina
import numpy as np
from tqdm import tqdm

class LigandPreparer:
    def __init__(self, smiles: str, name: str, xtb_path: str = "xtb"):
        self.smiles = smiles
        self.name = name
        self.xtb_path = xtb_path
        self.mol = None
        self.conformers = []

    def standardize(self):
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {self.smiles}")
        self.mol = Chem.AddHs(mol)

    def generate_conformers(self, n_confs: int = 250):
        if self.mol is None:
            raise RuntimeError("Molecule must be standardized before generating conformers.")

        params = AllChem.ETKDGv3()
        ids = AllChem.EmbedMultipleConfs(self.mol, numConfs=n_confs, params=params)

        results = AllChem.MMFFOptimizeMoleculeConfs(self.mol)
        # Store conformer IDs + energies
        self.conformer_energies = [
            (conf_id, result[1]) for conf_id, result in zip(ids, results)
        ]
        self.conformers = [x[0] for x in self.conformer_energies]

    def get_lowest_energy_conformer(self):
        if not hasattr(self, 'conformer_energies') or not self.conformer_energies:
            raise RuntimeError("No conformer energies available.")
        return min(self.conformer_energies, key=lambda x: x[1])[0]

    def cluster_and_select(self, final_n: int = 5, rmsd_threshold: float = 0.75, min_energy_gap: float = 0.5):
        """
        Clusters conformers using RMSD and selects top ones based on xTB energy.
        Applies a minimum energy gap between selected conformers.
        """
        if not self.conformers:
            raise RuntimeError("No conformers generated. Run generate_conformers() first.")

        # Calculate RMSD matrix
        rmslist = AllChem.GetConformerRMSMatrix(self.mol, prealigned=False)
        clusters = [[0]]  # Start first cluster with the first conformer

        for i in range(1, len(self.conformers)):
            added = False
            for cluster in clusters:
                rmsds = [rmslist[min(i, j) * (max(i, j) - 1) // 2] for j in cluster]
                if all(r < rmsd_threshold for r in rmsds):
                    cluster.append(i)
                    added = True
                    break
            if not added:
                clusters.append([i])

        print(f"[INFO] Found {len(clusters)} conformer clusters at RMSD threshold {rmsd_threshold}")

        # For each cluster, pick the lowest energy conformer (fake energy here for demo)
        selected = []
        energies = [AllChem.MMFFGetMoleculeForceField(self.mol, AllChem.MMFFGetMoleculeProperties(self.mol)).CalcEnergy() for _ in self.conformers]

        used = set()
        for cluster in clusters:
            sorted_cluster = sorted(cluster, key=lambda idx: energies[idx])
            for idx in sorted_cluster:
                if all(abs(energies[idx] - energies[u]) >= min_energy_gap for u in used):
                    selected.append(idx)
                    used.add(idx)
                    break
            if len(selected) >= final_n:
                break

        self.conformers = selected
        print(f"[INFO] Selected {len(selected)} conformers based on energy and RMSD")
        for i, conf_id in enumerate(self.conformers):
            print(f"  - Conf {conf_id:3d}: Energy = {energies[conf_id]:.4f}")

    def optimize_with_xtb(self, output_dir: Path):
        """
        Optimizes each conformer using xTB and stores energies.
        """
        self.conformer_energies = []

        for idx, conf_id in enumerate(tqdm(self.conformers, desc=f"[xTB] Optimizing {self.name}", unit="conf")):
            mol_block = Chem.MolToMolBlock(self.mol, confId=conf_id)
            temp_dir = tempfile.mkdtemp()
            mol_file = Path(temp_dir) / "input.mol"

            with open(mol_file, 'w') as f:
                f.write(mol_block)

            # Run xtb optimization
            cmd = [self.xtb_path, str(mol_file), "--opt", "--gfn", "1", "--chrg", "0", "--uhf", "0"]
            subprocess.run(cmd, cwd=temp_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Save optimized structure
            xtb_xyz = Path(temp_dir) / "xtbopt.xyz"
            xtb_log = Path(temp_dir) / "xtbopt.log"

            if xtb_xyz.exists():
                dest = output_dir / f"{self.name}_conf{idx}.xyz"
                shutil.copy(xtb_xyz, dest)

            # Parse xTB energy
            energy = None
            if xtb_log.exists():
                with open(xtb_log, 'r') as f:
                    for line in f:
                        if "TOTAL ENERGY" in line.upper():
                            try:
                                energy = float(line.strip().split()[-1])
                            except:
                                pass

            if energy is None:
                print(f"[WARNING] Energy not found for conf {conf_id}, skipping.")
            else:
                self.conformer_energies.append((conf_id, energy))

            shutil.rmtree(temp_dir)

        if not self.conformer_energies:
            raise RuntimeError("xTB optimization failed for all conformers.")

    def save_final_conformers(self, output_dir: Path) -> Path:
        if not self.conformers:
            raise RuntimeError("No conformers to save.")
        sdf_path = output_dir / f"{self.name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        for conf_id in self.conformers:
            writer.write(self.mol, confId=conf_id)
        writer.close()
        return sdf_path

    def convert_to_pdbqt(self, output_dir: Path, mode: str = "ensemble"):
        """
        Converts selected conformers to PDBQT files using Open Babel.
        
        Parameters:
            output_dir (Path): Directory where files will be written.
            mode (str): "ensemble" or "lowest_energy"
        
        Returns:
            List[Path]: Paths to generated PDBQT files.
        """
        pdbqt_paths = []

        if not self.conformers:
            raise RuntimeError("No conformers selected. Run conformer generation and selection first.")

        # Determine which conformers to convert
        if mode == "ensemble":
            conformers_to_process = self.conformers

        elif mode == "lowest_energy":
            conformers_to_process = [self.get_lowest_energy_conformer()]

        else:
            raise ValueError(f"Unknown docking mode: '{mode}' (expected 'ensemble' or 'lowest_energy')")

        for idx, conf_id in enumerate(conformers_to_process):
            sdf_path = output_dir / f"{self.name}_conf{idx}.sdf"
            pdbqt_path = output_dir / f"{self.name}_conf{idx}.pdbqt"

            # Save the conformer as an SDF file
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(self.mol, confId=conf_id)
            writer.close()

            # Convert to PDBQT using Open Babel
            cmd = [
                "obabel",
                str(sdf_path),
                "-O", str(pdbqt_path),
                "--partialcharge", "gasteiger"
            ]
            subprocess.run(cmd, check=True)

            pdbqt_paths.append(pdbqt_path)

        return pdbqt_paths


