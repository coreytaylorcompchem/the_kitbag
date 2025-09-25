import subprocess
import tempfile
import shutil
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem

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
    def __init__(self, pdb_path: Path, name: str = "receptor", method: str = "autodocktools"):
        self.pdb_path = Path(pdb_path)
        self.name = name
        self.method = method.lower()
        self.pdbqt_path = None

    def prepare(self, output_dir: Path):
        output_dir.mkdir(parents=True, exist_ok=True)
        self.pdbqt_path = output_dir / f"{self.name}.pdbqt"

        if self.pdbqt_path.exists():
            logger.info(f"Using existing receptor PDBQT at: {self.pdbqt_path}")
            return self.pdbqt_path

        logger.info(f"Preparing receptor using Open Babel")

        cmd = [
            "obabel",
            str(self.pdb_path),
            "-O", str(self.pdbqt_path),
            "-xr",  # Remove water
            "--partialcharge", "gasteiger"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Open Babel failed:\n{result.stderr}")
            raise RuntimeError("Protein preparation failed.")

        logger.info(f"[Protein] Receptor PDBQT saved at: {self.pdbqt_path}")
        return self.pdbqt_path

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

        # Sort conformers by energy after optimization
        self.conformer_energies.sort(key=lambda x: (x[1] if x[1] is not None else float('inf')))
        self.conformers = [conf_id for conf_id, _ in self.conformer_energies]

    def convert_to_pdbqt(self, output_dir: Path, mode: str = "ensemble"):
        if not self.conformers:
            raise RuntimeError("No conformers available to convert. Run generation and selection first.")

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

            # Convert using Open Babel with Gasteiger charges
            cmd = ["obabel", str(sdf_path), "-O", str(pdbqt_path), "--partialcharge", "gasteiger"]
            subprocess.run(cmd, check=True)

            pdbqt_paths.append(pdbqt_path)

        return pdbqt_paths

@register_task("prepare_receptor_pdbqt", 
               category='Receptor preparation',
               description="Prepare the receptor for docking. Format: pdbqt")
def prepare_receptor_pdbqt(backend, ligand, config, **kwargs):
    preparer = get_protein_preparer(backend, config)
    output_dir = Path(config.get("output_dir", "output"))
    pdbqt_path = preparer.prepare(output_dir)

    backend.cache["receptor_pdbqt"] = pdbqt_path

@register_task("standardise_ligand", 
               category='Ligand preparation',
               description="Prepare and generate 3D coords for each ligand.")
def standardise_ligand(backend, ligand, config, **kwargs):
    preparer = get_ligand_preparer(backend, ligand)
    preparer.standardise()
    logger.info(f"Standardised ligand {ligand['name']}.")

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
    logger.info(f"Converted conformers to PDBQT for ligand: {ligand['name']}.")

import pandas as pd
from rdkit import Chem

@register_task("dock", category='Docking', description="Dock with backend specified.")
def dock(backend, ligand, config, **kwargs):
    pdbqt_paths = ligand.get("pdbqt_paths", [])
    if not pdbqt_paths:
        logger.warning(f"No PDBQT paths found for ligand {ligand['name']}, skipping docking.")
        return

    docking_mode = config.get("docking", {}).get("docking_mode", "ensemble").lower()
    pocket_id = kwargs.get("pocket_id", None)
    output_dir = Path(config['output_dir'])

    docking_results = []

    if docking_mode == "lowest_energy":
        pdbqt_path = pdbqt_paths[0]
        ligand["pdbqt_path"] = pdbqt_path

        try:
            output_filename = f"{ligand['name']}_pocket{pocket_id}_docked.sdf" if pocket_id else f"{ligand['name']}_docked.sdf"
            output_path = output_dir / output_filename

            output_path = backend.dock(ligand, config, output_path=output_path)
            score = extract_gnina_score(output_path)

            docking_results.append({
                "ligand": ligand["name"],
                "conformer": 0,
                "pdbqt": str(pdbqt_path),
                "docked_sdf": str(output_path),
                "score": score
            })

            ligand["docking_results"] = docking_results
            logger.info(f"Docked lowest-energy conformer for ligand {ligand['name']} into pocket {pocket_id}.")

        except Exception as e:
            logger.error(f"❌ Docking failed for ligand {ligand['name']}: {e}")
            ligand["docking_results"] = []

    elif docking_mode == "ensemble":
        logger.info(f"Docking ligand {ligand['name']} with {len(pdbqt_paths)} conformers...")

        for idx, pdbqt_path in enumerate(pdbqt_paths):
            ligand["pdbqt_path"] = pdbqt_path

            try:
                output_filename = f"{ligand['name']}_conf{idx}_pocket{pocket_id}_docked.sdf" if pocket_id is not None else f"{ligand['name']}_conf{idx}_docked.sdf"
                output_path = output_dir / output_filename

                output_path = backend.dock(ligand, config, output_path=output_path)
                score = extract_gnina_score(output_path)

                docking_results.append({
                    "ligand": ligand["name"],
                    "conformer": idx,
                    "pdbqt": str(pdbqt_path),
                    "docked_sdf": str(output_path),
                    "score": score
                })

            except Exception as e:
                logger.error(f"❌ Docking failed for ligand {ligand['name']} conformer {idx}: {e}")

        ligand["docking_results"] = docking_results

    else:
        raise ValueError(f"❌❌ Unknown docking_mode: {docking_mode}")

    if docking_results:
        pocket_id = kwargs.get("pocket_id", 0)
        pocket_dir = Path(config["output_dir"]) / f"pocket_{pocket_id}"
        pocket_dir.mkdir(exist_ok=True, parents=True)

        score_csv_path = pocket_dir / "docking_scores.csv"
        score_rows = []

        for result in docking_results:
            sdf_path = result["docked_sdf"]
            suppl = Chem.SDMolSupplier(str(sdf_path))
            for pose_idx, mol in enumerate(suppl):
                if mol is None:
                    continue
                try:
                    cnn_score = float(mol.GetProp("CNNscore"))
                    affinity = float(mol.GetProp("minimizedAffinity"))
                except KeyError:
                    logger.warning(f"Missing scores in pose {pose_idx} for ligand {ligand['name']}")
                    continue

                score_rows.append({
                    "ligand": ligand["name"],
                    "conformer": result["conformer"],
                    "pose_idx": pose_idx,
                    "score": affinity,
                    "pose_rank": cnn_score,
                    "pdbqt": result["pdbqt"],
                    "docked_sdf": sdf_path,
                    "pocket": pocket_id
                })

        # Append if file exists
        if score_csv_path.exists():
            existing_df = pd.read_csv(score_csv_path)
            df = pd.concat([existing_df, pd.DataFrame(score_rows)], ignore_index=True)
        else:
            df = pd.DataFrame(score_rows)

        df.to_csv(score_csv_path, index=False)

def extract_gnina_score(sdf_path):
    """
    Extracts the GNINA docking score from the first molecule in an SDF file.
    """
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next((m for m in suppl if m is not None), None)

        if mol is None:
            logger.warning(f"Could not parse SDF: {sdf_path}")
            return None

        for field in ["CNNscore", "minimizedAffinity", "score"]:
            if mol.HasProp(field):
                return float(mol.GetProp(field))

        logger.warning(f"No recognized score field found in {sdf_path}")
        return None

    except Exception as e:
        logger.warning(f"Failed to extract score from {sdf_path}: {e}")
        return None



