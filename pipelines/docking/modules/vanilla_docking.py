import subprocess
import tempfile
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

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

    def cluster_and_select(self, final_n: int = 35):
        if not self.conformers:
            raise RuntimeError("No conformers generated. Run generate_conformers() first.")
        self.conformers = self.conformers[:final_n]

    def optimize_with_xtb(self, output_dir: Path):
        """
        Optimizes each conformer using xtb and writes out xyz files.
        """
        xtb_confs = []

        for idx, conf_id in enumerate(self.conformers):
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
            if xtb_xyz.exists():
                dest = output_dir / f"{self.name}_conf{idx}.xyz"
                shutil.copy(xtb_xyz, dest)
                xtb_confs.append(dest)

            shutil.rmtree(temp_dir)

        return xtb_confs

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
        Converts conformers to individual PDBQT files.
        """
        pdbqt_paths = []

        conformers_to_process = []

        if mode == "ensemble":
            conformers_to_process = self.conformers

        elif mode == "lowest_energy":
            lowest_conf = self.get_lowest_energy_conformer()
            conformers_to_process = [lowest_conf]

        else:
            raise ValueError(f"Unknown docking mode: {mode}")

        for idx, conf_id in enumerate(conformers_to_process):
            sdf_path = output_dir / f"{self.name}_conf{idx}.sdf"
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(self.mol, confId=conf_id)
            writer.close()

            pdbqt_path = output_dir / f"{self.name}_conf{idx}.pdbqt"
            cmd = [
                "obabel",
                str(sdf_path),
                "-O", str(pdbqt_path),
                "--partialcharge", "gasteiger"
            ]
            subprocess.run(cmd, check=True)

            pdbqt_paths.append(pdbqt_path)

        return pdbqt_paths

