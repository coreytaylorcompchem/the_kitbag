import subprocess
from pathlib import Path
from openbabel import pybel

class ProteinPreparer:
    def __init__(self, pdb_path: Path, work_dir: Path, pH: float = 7.4):
        self.pdb_path = pdb_path
        self.work_dir = work_dir
        self.pH = pH

        self.protonated_pdb = self.work_dir / "protonated.pdb"
        self.pdbqt_path = self.work_dir / "protein.pdbqt"
        self.reference_ligand_path = self.work_dir / "reference_ligand.pdb"

    def extract_crystal_ligand(self):
        with open(self.pdb_path, 'r') as infile, open(self.reference_ligand_path, 'w') as outfile:
            for line in infile:
                if line.startswith("HETATM"):
                    outfile.write(line)
        return self.reference_ligand_path

    def run_propka(self):
        try:
            # Attempt to run PROPKA
            cmd = ["propka", str(self.pdb_path), "--pH", str(self.pH)]
            subprocess.run(cmd, check=True)

            propka_out = self.pdb_path.with_suffix(".pka")
            if not propka_out.exists():
                raise RuntimeError("PROPKA output not found")

            # If successful, protonate using pdb2pqr
            protonated = self.work_dir / "protonated.pdb"
            cmd = [
                "pdb2pqr",
                "--ff", "PARSE",
                "--with-ph", str(self.pH),
                str(self.pdb_path),
                str(protonated)
            ]
            subprocess.run(cmd, check=True)
            self.protonated_pdb = protonated
            print("[INFO] Protonation with PROPKA + pdb2pqr succeeded.")
            return True
        except Exception as e:
            print(f"[WARNING] PROPKA protonation failed, falling back to Open Babel: {e}")
            return False

    def fallback_protonation_openbabel(self):
        """
        Fallback protonation using Open Babel.
        """
        molecule = list(pybel.readfile("pdb", str(self.pdb_path)))[0]
        molecule.OBMol.CorrectForPH(self.pH)
        molecule.addh()
        fallback_pdb = self.work_dir / "protonated_openbabel.pdb"
        molecule.write("pdb", str(fallback_pdb), overwrite=True)
        self.protonated_pdb = fallback_pdb
        return fallback_pdb

    def convert_to_pdbqt(self):
        molecule = list(pybel.readfile("pdb", str(self.protonated_pdb)))[0]
        molecule.OBMol.CorrectForPH(self.pH)
        molecule.addh()
        for atom in molecule.atoms:
            atom.OBAtom.GetPartialCharge()
        molecule.write("pdbqt", str(self.pdbqt_path), overwrite=True)
        return self.pdbqt_path

    def prepare(self):
        self.extract_crystal_ligand()

        # Try PROPKA + pdb2pqr
        success = self.run_propka()

        # Fallback to Open Babel if needed
        if not success:
            self.fallback_protonation_openbabel()

        self.convert_to_pdbqt()
        return self.pdbqt_path
