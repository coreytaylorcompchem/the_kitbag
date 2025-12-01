import os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    'site_pka',
    description="Enumerate acidic sites, generate protonated and deprotonated geometries for pKa calculations.",
    modifies_geometry=False,
    category='Property'
)
def run(backend, csv_file, step_config, global_config=None):
    df = pd.read_csv(csv_file)

    if 'SMILES' not in df.columns or 'name' not in df.columns:
        raise ValueError("CSV must contain 'SMILES' and 'name' columns.")

    all_results = {}

    output_dir = step_config.get('output_directory', 'outputs/site_pka')
    os.makedirs(output_dir, exist_ok=True)

    for _, row in df.iterrows():
        smiles = row['SMILES']
        name = str(row['name'])

        logger.info(f"Processing molecule: {name} ({smiles})")

        try:
            mol = load_molecule_from_smiles(smiles)
        except Exception as e:
            logger.warning(f"Skipping {name}: failed to parse or embed molecule. Error: {e}")
            continue

        acidic_sites = find_acidic_sites(mol)
        if not acidic_sites:
            logger.warning(f"No acidic sites found for molecule: {name}")
            continue

        site_data = {}
        for i, atom_idx in enumerate(acidic_sites):
            try:
                ha_mol = protonate_site(mol, atom_idx)
                a_mol = deprotonate_site(mol, atom_idx)

                ha_xyz = mol_to_xyz_string(ha_mol)
                ha_charge = Chem.GetFormalCharge(ha_mol)
                ha_mult = 1   
                a_xyz = mol_to_xyz_string(a_mol)
                a_charge = Chem.GetFormalCharge(a_mol)
                a_mult = 1    

                site_data[i] = {
                    'atom_idx': atom_idx,
                    'acidic_atom_symbol': mol.GetAtomWithIdx(atom_idx).GetSymbol(),
                    'ha_xyz': ha_xyz,
                    'a_xyz': a_xyz,
                    'ha_mol': ha_mol,
                    'a_mol': a_mol,
                    'ha_charge': ha_charge,
                    'ha_mult': ha_mult,
                    'a_charge': a_charge,
                    'a_mult': a_mult,
                }
            except Exception as e:
                logger.warning(f"Failed to process site {atom_idx} in {name}: {e}")
                continue

        all_results[name] = site_data
        logger.info(f"Found {len(site_data)} acidic site(s) for {name}.")

        # Save geometries to disk
        save_site_geometries(site_data, name, output_dir)

    # Save all_results as pickle
    output_path = os.path.join(output_dir, "site_pka_data.pkl")
    with open(output_path, 'wb') as f:
        pickle.dump(all_results, f)

    logger.info(f"Saved all site pKa data to {output_path}")

    return all_results


def load_molecule_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    return mol

def find_acidic_sites(mol):
    """
    Identify acidic *functional groups*, regardless of protonation.
    Returns heavy atom indices that are protonatable/deprotonatable acidic sites.
    Does NOT require explicit hydrogens.
    Excludes C–H acidity.
    """
    mol = Chem.AddHs(mol)

    # Acidic fragments (stealing approach from EPIK)
    FRAGMENTS = [
        # Carboxylic acid / carboxylate
        ("carboxylic",       "[CX3](=O)[OX2H1]",        5.0),
        ("carboxylate",      "[CX3](=O)[O-]",           5.0),

        # Sulfonic acid / sulfonate
        ("sulfonic",         "S(=O)(=O)[OX2H1]",        5.0),
        ("sulfonate",        "S(=O)(=O)[O-]",           5.0),

        # Phosphoric acid series (mono/di/triprotic)
        ("phosphate_OH",     "P(=O)(O)(O)[OX2H1]",      4.5),
        ("phosphate_O-",     "P(=O)(O)(O)[O-]",         4.5),

        # Pyro/Meta/Poly-phosphate (any protonation)
        ("polyphosphate_OH", "P(=O)(O)O[OX2H1]",        4.5),
        ("polyphosphate_O-", "P(=O)(O)O[O-]",           4.5),

        # Boric acid / borate
        ("boric_OH",         "B(O)(O)[OX2H1]",          4.0),
        ("borate_O-",        "B(O)(O)[O-]",             4.0),

        # Arsenic acid / arsenate
        ("arsenic_OH",       "As(=O)(O)(O)[OX2H1]",     4.0),
        ("arsenate_O-",      "As(=O)(O)(O)[O-]",        4.0),

        # Chromic acid / chromate
        ("chromate_OH",      "Cr(=O)(=O)[OX2H1]",       3.5),
        ("chromate_O-",      "Cr(=O)(=O)[O-]",          3.5),

        # Phenols / Enols
        ("phenol",           "[OX2H1][cX3]",            2.0),
        ("phenolate",        "[O-][cX3]",               2.0),
        ("enol",             "[OX2H1][C]=C",            1.5),
        ("enolate",          "[O-][C]=C",               1.5),

        # Thiols / Thiolates
        ("thiol",            "[SX2H1]",                 2.0),
        ("thiolate",         "[S-]",                    2.0),

        # Amides / imides / heterocycles
        ("amide_NH",         "N([CX3](=O))[H]",         1.0),
        ("imide_NH",         "[CX3](=O)N([CX3](=O))[H]",1.2),
        ("heterocycle_nH",   "[nH]",                    1.0),
        ("protonated_amine", "[N+H3]",                  1.0),

        # Hydroxylamines / oximes
        ("oxime",            "N(=C)O",                  1.0),
    ]

    acidic_sites = {}

    for name, smarts, weight in FRAGMENTS:
        patt = Chem.MolFromSmarts(smarts)
        if not patt:
            continue

        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)

                # NO carbon acidity
                if atom.GetAtomicNum() == 6:
                    continue

                acidic_sites[idx] = acidic_sites.get(idx, 0.0) + weight

    # Additional generic rule:
    # Any O, S, N, P, B, As with a single-bond to electronegative center
    for atom in mol.GetAtoms():
        Z = atom.GetAtomicNum()
        if Z in (8, 7, 16, 15, 5, 33):  # O, N, S, P, B, As
            # ignore carbon-bound neutral oxygen unless matched above
            if Z == 8 and atom.GetDegree() == 1 and atom.GetNeighbors()[0].GetAtomicNum() == 6:
                continue
            acidic_sites.setdefault(atom.GetIdx(), 0.5)

    return sorted(acidic_sites.keys(), key=lambda i: -acidic_sites[i])



def protonate_site(mol, atom_idx):
    """
    Add a proton (H) to the acidic site atom, return new mol.
    If atom already has a bonded H, return copy of original mol.
    """
    mol_copy = Chem.Mol(mol)
    mol_copy = Chem.AddHs(mol_copy)
    atom = mol_copy.GetAtomWithIdx(atom_idx)

    # Check if atom already has H neighbor
    has_H = any(nbr.GetAtomicNum() == 1 for nbr in atom.GetNeighbors())
    if has_H:
        # Already protonated, just return as-is
        return mol_copy

    # Add new H bonded to this atom
    new_h = Chem.Atom(1)
    rw_mol = Chem.RWMol(mol_copy)
    new_h_idx = rw_mol.AddAtom(new_h)
    rw_mol.AddBond(atom_idx, new_h_idx, order=Chem.rdchem.BondType.SINGLE)

    mol_copy = rw_mol.GetMol()

    # Sanitise and optimise
    Chem.SanitizeMol(mol_copy)
    AllChem.EmbedMolecule(mol_copy)
    AllChem.UFFOptimizeMolecule(mol_copy)

    return mol_copy


def deprotonate_site(mol, atom_idx):
    """
    Remove acidic proton bonded to atom_idx, return new mol.
    """
    mol_copy = Chem.RWMol(mol)
    Chem.SanitizeMol(mol_copy)

    # Find one H neighbor to remove
    h_idx = None
    atom = mol_copy.GetAtomWithIdx(atom_idx)
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            h_idx = nbr.GetIdx()
            break

    if h_idx is not None:
        mol_copy.RemoveAtom(h_idx)
    else:
        logger.warning(f"No proton found on acidic site atom idx {atom_idx} for deprotonation.")

    # Adjust formal charge to -1 (TODO: make more betterer)
    atom = mol_copy.GetAtomWithIdx(atom_idx)
    atom.SetFormalCharge(atom.GetFormalCharge() - 1)

    mol_copy = mol_copy.GetMol()
    Chem.SanitizeMol(mol_copy)
    mol_copy = Chem.AddHs(mol_copy)
    AllChem.EmbedMolecule(mol_copy, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol_copy)

    return mol_copy


def mol_to_xyz_string(mol):
    """
    Convert RDKit mol to xyz format string.
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), "",]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
    return "\n".join(lines)


def save_site_geometries(site_data, name, output_dir):
    """
    Save the HA and A- XYZ geometries from site_pka data to disk.

    Args:
        site_data (dict): Result of the site_pka task.
        name (str): Molecule name (used in filenames).
        output_dir (str): Where to write files.
    """
    os.makedirs(output_dir, exist_ok=True)

    for site_idx, data in site_data.items():
        ha_path = os.path.join(output_dir, f"{name}_site_{site_idx}_HA.xyz")
        a_path = os.path.join(output_dir, f"{name}_site_{site_idx}_A.xyz")

        with open(ha_path, "w") as f:
            f.write(data['ha_xyz'])

        with open(a_path, "w") as f:
            f.write(data['a_xyz'])

        logger.info(f"Wrote geometry files for site {site_idx} of molecule {name} to {output_dir}")
