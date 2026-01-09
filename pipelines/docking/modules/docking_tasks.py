import subprocess
import multiprocessing
import tempfile
import shutil
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

def get_protein_preparer(backend, config):
    if "protein_preparer" not in backend.cache:
        protein_path = Path(config["protein"]["pdb_path"])
        backend.cache["protein_preparer"] = ProteinPreparer(pdb_path=protein_path)
    return backend.cache["protein_preparer"]


class ProteinPreparer:
    """
    Prepares receptor PDB files for docking / AmberTools:
      - Removes crystal ligands and water
      - Extracts ligand for binding site definition
      - Protonates protein using pdb4amber
      - Cleans PDB for tleap (removes CONECTs, renames disulfides)
    """

    def __init__(self, pdb_path: Path, name: str = "receptor", pH: float = 7.4):
        self.pdb_path = Path(pdb_path)
        self.name = name
        self.pH = pH
        self.cleaned_pdb = None
        self.extracted_ligand_path = None
        self.protonated_pdb = None
        self.pdbqt_path = None

    def prepare(self, output_dir: Path):
        output_dir.mkdir(parents=True, exist_ok=True)

        self.cleaned_pdb = output_dir / f"{self.name}_clean.pdb"
        self.extracted_ligand_path = output_dir / f"{self.name}_crystal_ligand.pdb"
        # Alias for downstream consumers (docking box, etc.)
        self.reference_ligand_path = self.extracted_ligand_path
        self.protonated_pdb = output_dir / f"{self.name}_protonated.pdb"
        self.pdbqt_path = output_dir / f"{self.name}.pdbqt"

        # Reuse if already exists
        if self.protonated_pdb.exists() and self.pdbqt_path.exists():
            logger.debug(f"[ProteinPreparer] Using existing prepared receptor files.")
            return self.pdbqt_path

        logger.debug(f"[ProteinPreparer] Cleaning protein: {self.pdb_path.name}")

        # Step 1: Separate protein vs. ligand atoms
        hetatm_blocks = {}
        with open(self.pdb_path, 'r') as infile, \
             open(self.cleaned_pdb, 'w') as protein_out, \
             open(self.extracted_ligand_path, 'w') as ligand_out:

            for line in infile:
                if line.startswith("ATOM"):
                    protein_out.write(line)
                elif line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    if res_name not in ("HOH", "WAT"):
                        hetatm_blocks.setdefault(res_name, []).append(line)
                elif line.startswith(("TER", "END")):
                    protein_out.write(line)

            if hetatm_blocks:
                largest_res = max(hetatm_blocks, key=lambda r: len(hetatm_blocks[r]))
                for l in hetatm_blocks[largest_res]:
                    ligand_out.write(l)
                logger.info(f"[ProteinPreparer] Extracted crystal ligand: {largest_res}")
            else:
                logger.warning("[ProteinPreparer] No non-water HETATM groups found.")

        # Step 2: Protonate protein with pdb4amber
        cmd_pdb4amber = [
            "pdb4amber",
            "-i", str(self.cleaned_pdb),
            "-o", str(self.protonated_pdb),
            # "-a", # keep amber compatible residues"
            "--add-missing-atoms",
            "--reduce",
            # "--pH", str(self.pH),
            "--most-populous",
            "--keep-altlocs",
            "--no-conect",
        ]

        result = subprocess.run(cmd_pdb4amber, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"pdb4amber failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
            raise RuntimeError("pdb4amber failed")

        logger.info(f"[ProteinPreparer] Protonated protein saved at: {self.protonated_pdb}")

        # Step 3: Strip CONECT records
        # self._remove_conect_records(self.protonated_pdb)

        # Step 4: Convert to PDBQT
        cmd_pdbqt = [
            "obabel",
            str(self.protonated_pdb),
            "-O", str(self.pdbqt_path),
            "--partialcharge", "gasteiger",
            "-xh"
        ]

        result = subprocess.run(cmd_pdbqt, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Open Babel PDB → PDBQT failed:\n{result.stderr}")
            raise RuntimeError("Protein PDBQT conversion failed.")

        logger.info(f"[ProteinPreparer] Receptor PDBQT saved at: {self.pdbqt_path}")
        return self.pdbqt_path

    # @staticmethod
    # def _remove_conect_records(pdb_path: Path):
    #     """Strip all CONECT lines from a PDB file in-place."""
    #     tmp_path = pdb_path.with_suffix(".tmp")
    #     with open(pdb_path) as f_in, open(tmp_path, 'w') as f_out:
    #         for line in f_in:
    #             if not line.startswith("CONECT"):
    #                 f_out.write(line)
    #     shutil.move(tmp_path, pdb_path)
    #     logger.debug(f"[ProteinPreparer] Removed CONECT records from {pdb_path.name}")


def get_ligand_preparer(backend, ligand):
    if "ligand_preparers" not in backend.cache:
        backend.cache["ligand_preparers"] = {}

    if ligand["name"] not in backend.cache["ligand_preparers"]:
        backend.cache["ligand_preparers"][ligand["name"]] = LigandPreparer(
            smiles=ligand["smiles"], name=ligand["name"]
        )
    return backend.cache["ligand_preparers"][ligand["name"]]

class LigandPreparer:
    def __init__(self, smiles: str, name: str, xtb_path: str = "xtb"):
        self.smiles = smiles
        self.name = name
        self.xtb_path = xtb_path
        self.mol = None
        self.conformers = []
        self.conformer_energies = []

    def standardise(self):
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {self.smiles}")
        self.mol = Chem.AddHs(mol)
    
    def protonate(self, pH: float = 7.4):
        if self.mol is None:
            raise RuntimeError("Ligand must be standardised before protonation.")

        mol = Chem.Mol(self.mol)

        # Step 1: Normalise functional groups
        normalizer = rdMolStandardize.Normalizer()
        mol = normalizer.normalize(mol)

        # Step 2: Use Rdkit reionizer (pH-aware where possible)
        reionizer = rdMolStandardize.Reionizer()
        mol = reionizer.reionize(mol)

        # Step 3: Explicit hydrogens & sanitisation
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)

        self.mol = mol

        logger.debug(f"[LigandPreparer] Protonated ligand {self.name} at pH {pH}")

    def generate_conformers(self, n_confs: int = 250):
        if self.mol is None:
            raise RuntimeError("Molecule must be standardised before generating conformers.")

        params = AllChem.ETKDGv3()
        ids = AllChem.EmbedMultipleConfs(self.mol, numConfs=n_confs, params=params)

        results = AllChem.MMFFOptimizeMoleculeConfs(self.mol)
        self.conformer_energies = [(conf_id, result[1]) for conf_id, result in zip(ids, results)]
        self.conformers = [x[0] for x in self.conformer_energies]

    def get_lowest_energy_conformer(self):
        if not self.conformer_energies:
            raise RuntimeError("No conformer energies available.")
        return min(self.conformer_energies, key=lambda x: x[1])[0]

    def cluster_and_select(self, final_n: int = 5, rmsd_threshold: float = 0.75, min_energy_gap: float = 0.5):
        if not self.conformers:
            raise RuntimeError("No conformers generated. Run generate_conformers() first.")

        rmslist = AllChem.GetConformerRMSMatrix(self.mol, prealigned=False)
        clusters = [[0]]

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

        logger.info(f"Found {len(clusters)} conformer clusters at RMSD threshold {rmsd_threshold}")

        # Use MMFF energies for clustering selection
        energies = [AllChem.MMFFGetMoleculeForceField(self.mol, AllChem.MMFFGetMoleculeProperties(self.mol)).CalcEnergy() for _ in self.conformers]

        selected = []
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
        logger.info(f"Selected {len(selected)} conformers based on energy and RMSD")
        for i, conf_id in enumerate(self.conformers):
            logger.info(f"  - Conf {conf_id:3d}: Energy = {energies[conf_id]:.4f}")

    def optimize_with_xtb(self, output_dir: Path):
        self.conformer_energies = []

        for idx, conf_id in enumerate(tqdm(self.conformers, desc=f"[xTB] Optimizing {self.name}", unit="conf")):
            mol_block = Chem.MolToMolBlock(self.mol, confId=conf_id)
            temp_dir = tempfile.mkdtemp()
            mol_file = Path(temp_dir) / "input.mol"

            with open(mol_file, 'w') as f:
                f.write(mol_block)

            cmd = [self.xtb_path, str(mol_file), "--opt", "--gfn", "1", "--chrg", "0", "--uhf", "0"]
            subprocess.run(cmd, cwd=temp_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            xtb_xyz = Path(temp_dir) / "xtbopt.xyz"
            xtb_log = Path(temp_dir) / "xtbopt.log"

            if xtb_xyz.exists():
                dest = output_dir / f"{self.name}_conf{idx}.xyz"
                shutil.copy(xtb_xyz, dest)

            energy = None
            if xtb_log.exists():
                with open(xtb_log, 'r') as f:
                    for line in f:
                        if "TOTAL ENERGY" in line.upper():
                            try:
                                energy = float(line.strip().split()[-1])
                            except Exception:
                                pass

            if energy is None:
                logger.warning(f"Could not read energy for conformer {conf_id}")

            self.conformer_energies.append((conf_id, energy))
            shutil.rmtree(temp_dir)

        # Sort conformers by energy after optimisation
        self.conformer_energies.sort(key=lambda x: (x[1] if x[1] is not None else float('inf')))
        self.conformers = [conf_id for conf_id, _ in self.conformer_energies]

    def convert_to_pdbqt(self, output_dir: Path, mode: str = "ensemble"):
        if self.mol is None:
            raise RuntimeError(f"[{self.name}] Molecule not initialised. Run standardise() before convert_to_pdbqt.")

        # Fallback: if no conformers, generate one
        if not self.conformers:
            logger.debug(f"No conformers available — generating a single default conformer for {self.name}.")
            conf_id = AllChem.EmbedMolecule(self.mol)
            if conf_id < 0:
                raise RuntimeError(f"Failed to embed 3D conformer for {self.name}.")
            self.conformers = [conf_id]

        pdbqt_paths = []

        if mode == "ensemble":
            conformers_to_convert = self.conformers
        elif mode == "lowest_energy":
            conformers_to_convert = [self.get_lowest_energy_conformer()]
        else:
            raise ValueError(f"Unknown docking mode: {mode}")

        for idx, conf_id in enumerate(conformers_to_convert):
            sdf_path = output_dir / f"{self.name}_conf{idx}.sdf"
            pdbqt_path = output_dir / f"{self.name}_conf{idx}.pdbqt"

            writer = Chem.SDWriter(str(sdf_path))
            writer.write(self.mol, confId=conf_id)
            writer.close()

            # Convert using Open Babel
            cmd = [
                "obabel", str(sdf_path), "-xh", "-p", "7.4", "-O", str(pdbqt_path)
            ]
            subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)

            pdbqt_paths.append(pdbqt_path)

        return pdbqt_paths

@register_task("prepare_receptor_pdbqt", 
               category="Receptor preparation",
               description="Prepare the receptor for docking. Format: pdbqt")
def prepare_receptor_pdbqt(backend, ligand, config, **kwargs):

    if ligand is not None:
        raise RuntimeError(
            "prepare_receptor_pdbqt must be run in global context (ligand=None)"
        )

    # Reuse preparer instance if it already exists
    if "protein_preparer" not in backend.cache:
        preparer = ProteinPreparer(
            pdb_path=Path(config["protein"]["pdb_path"]),
            name="receptor",
            pH=config["protein"].get("pH", 7.4)
        )
        backend.cache["protein_preparer"] = preparer
    else:
        preparer = backend.cache["protein_preparer"]

    if "receptor_pdbqt" in backend.cache:
        logger.debug("Receptor already prepared — reusing.")
        return backend.cache["receptor_pdbqt"]

    output_dir = Path(config["output_dir"])
    pdbqt_path = preparer.prepare(output_dir)

    backend.cache["receptor_pdbqt"] = pdbqt_path
    return pdbqt_path

@register_task("standardise_ligand", 
               category='Ligand preparation',
               description="Prepare and generate 3D coords for each ligand.")
def standardise_ligand(backend, ligand, config, **kwargs):
    preparer = get_ligand_preparer(backend, ligand)
    preparer.standardise()
    logger.debug(f"Standardised ligand {ligand['name']}.")

@register_task(
    "protonate_ligand",
    category="Ligand preparation",
    description="Assign ligand protonation state at target pH."
)
def protonate_ligand(backend, ligand, config, **kwargs):
    preparer = get_ligand_preparer(backend, ligand)
    pH = config.get("protein", {}).get("pH", 7.4)
    preparer.protonate(pH=pH)

@register_task("generate_conformers", 
               category='Ligand preparation',
               description="Generate multiple feasible conformers from ligands.")
def generate_conformers(backend, ligand, config, **kwargs):
    preparer = get_ligand_preparer(backend, ligand)
    n_confs = config.get("docking", {}).get("initial_n_conformers", 250)
    preparer.generate_conformers(n_confs=n_confs)
    logger.info(f"Generated {n_confs} conformers for ligand {ligand['name']}.")

@register_task("cluster_conformers", 
               category='Ligand preparation',
               description="Cluster conformers by specified energy and RMSD criteria.")
def cluster_conformers(backend, ligand, config, **kwargs):
    preparer = get_ligand_preparer(backend, ligand)
    final_n = config.get("docking", {}).get("final_n_conformers", 5)
    rmsd = config.get("docking", {}).get("rmsd_threshold", 0.75)
    energy_gap = config.get("docking", {}).get("min_energy_gap", 0.5)
    preparer.cluster_and_select(final_n=final_n, rmsd_threshold=rmsd, min_energy_gap=energy_gap)

@register_task("save_final_conformers", 
               category='Ligand preparation',
               description="Save conformers that meet energy and RMSD criteria.")
def save_final_conformers(backend, ligand, config, **kwargs):
    preparer = get_ligand_preparer(backend, ligand)
    out_dir = Path(config.get("output_dir", "output"))
    for idx, conf_id in enumerate(preparer.conformers):
        sdf_path = out_dir / f"{ligand['name']}_final_conf{idx}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(preparer.mol, confId=conf_id)
        writer.close()
    logger.info(f"Saved final conformers for ligand {ligand['name']}.")

@register_task("convert_to_pdbqt", 
               category='Ligand preparation',
               description="Convert ligands to pdbqt for docking.")
def convert_to_pdbqt(backend, ligand, config, **kwargs):
    preparer = get_ligand_preparer(backend, ligand)
    output_dir = Path(config.get("output_dir", "output"))
    output_dir.mkdir(parents=True, exist_ok=True)
    mode = config.get("docking", {}).get("mode", "ensemble")
    pdbqt_paths = preparer.convert_to_pdbqt(output_dir, mode=mode)
    ligand['pdbqt_paths'] = [str(p) for p in pdbqt_paths]
    logger.debug(f"Converted conformers to PDBQT for ligand: {ligand['name']}.")

@register_task("dock", category='Docking', description="Dock with backend specified.")
def dock(backend, ligand, config, **kwargs):
    pdbqt_paths = ligand.get("pdbqt_paths", [])
    if not pdbqt_paths:
        logger.warning(f"No PDBQT paths found for ligand {ligand['name']}, skipping docking.")
        return []

    docking_mode = config.get("docking", {}).get("docking_mode", "ensemble").lower()
    output_dir = Path(config['output_dir'])
    pocket_id = kwargs.get("pocket_id", None)

    docking_outputs = []

    if docking_mode == "lowest_energy":
        pdbqt_path = pdbqt_paths[0]
        ligand["pdbqt_path"] = pdbqt_path
        output_filename = f"{ligand['name']}_pocket{pocket_id}_docked.sdf" if pocket_id else f"{ligand['name']}_docked.sdf"
        output_path = output_dir / output_filename

        try:
            output_path = backend.dock(ligand, config, output_path=output_path)
            docking_outputs.append({
                "conformer": 0,
                "pdbqt": str(pdbqt_path),
                "docked_output": output_path,
            })
        except Exception as e:
            logger.error(f"❌ Docking failed for {ligand['name']}: {e}")

    elif docking_mode == "ensemble":
        logger.debug(f"Docking ligand {ligand['name']} with {len(pdbqt_paths)} conformers...")

        for idx, pdbqt_path in enumerate(pdbqt_paths):
            ligand["pdbqt_path"] = pdbqt_path
            output_filename = f"{ligand['name']}_conf{idx}_pocket{pocket_id}_docked.sdf" if pocket_id is not None else f"{ligand['name']}_conf{idx}_docked.sdf"
            output_path = output_dir / output_filename

            try:
                output_path = backend.dock(ligand, config, output_path=output_path)
                docking_outputs.append({
                    "conformer": idx,
                    "pdbqt": str(pdbqt_path),
                    "docked_output": output_path,
                })
            except Exception as e:
                logger.error(f"❌ Docking failed for ligand {ligand['name']} conformer {idx}: {e}")

    else:
        raise ValueError(f"Unknown docking_mode: {docking_mode}")

    return docking_outputs