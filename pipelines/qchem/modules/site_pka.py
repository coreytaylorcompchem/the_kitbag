import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem.rdmolops import GetFormalCharge
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    'site_pka',
    description="Enumerate acidic sites, generate protonated and deprotonated geometries for pKa calculations.",
    modifies_geometry=False,
    category='Property'
)
def run(backend, xyz_file, step_config, global_config=None):
    """
    Given an input molecule in xyz_file, enumerate acidic sites,
    generate protonated (HA) and deprotonated (A-) forms as RDKit mols,
    convert to xyz strings, and return a dict of site info.

    Returns dict:
        {
          site_idx: {
            'atom_idx': int,
            'acidic_atom_symbol': str,
            'ha_mol': RDKit mol (protonated),
            'a_mol': RDKit mol (deprotonated),
            'ha_xyz': str,
            'a_xyz': str
          },
          ...
        }
    """
    mol = load_molecule_from_xyz(xyz_file)
    acidic_sites = find_acidic_sites(mol)

    if not acidic_sites:
        logger.warning("No acidic sites found.")
        return {}

    site_data = {}
    for i, atom_idx in enumerate(acidic_sites):
        ha_mol = protonate_site(mol, atom_idx)
        a_mol = deprotonate_site(mol, atom_idx)

        ha_xyz = mol_to_xyz_string(ha_mol)
        a_xyz = mol_to_xyz_string(a_mol)

        site_data[i] = {
            'atom_idx': atom_idx,
            'acidic_atom_symbol': mol.GetAtomWithIdx(atom_idx).GetSymbol(),
            'ha_mol': ha_mol,
            'a_mol': a_mol,
            'ha_xyz': ha_xyz,
            'a_xyz': a_xyz
        }

    logger.info(f"Found {len(site_data)} acidic sites for pKa calculation.")

    return site_data

def load_molecule_from_xyz(xyz_path):
    """
    Load an RDKit molecule from an XYZ file using the RDKit XYZ parser.
    """
    with open(xyz_path) as f:
        xyz_content = f.read()

    # RDKit doesn't natively parse XYZ directly to Mol with connectivity
    # Use workaround: embed mol from xyz coordinates
    mol = Chem.MolFromXYZBlock(xyz_content, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError("Failed to parse molecule from XYZ file.")

    # Add hydrogens and compute connectivity
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)

    return mol

def find_acidic_sites(mol):
    """
    Identify acidic protons' heavy atoms indices (sites for deprotonation).
    Uses SMARTS patterns for common acidic groups.
    Returns list of atom indices (heavy atoms bonded to acidic H).
    """
    acidic_smarts = {
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "phenol": "c[OX2H]",
        "alcohol": "[CX4][OX2H]",  
        "thiol": "[#16X2H]",  # sulfur with one H
        "primary_amine": "[NX3;H2]",  # -NH2
        "secondary_amine": "[NX3;H1]", # -NH
        # Add more patterns if needed
    }

    acidic_sites = set()

    for name, smarts in acidic_smarts.items():
        patt = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            # Find the acidic heavy atom index in the match
            # Usually the heavy atom bonded to acidic H
            heavy_atom_idx = None
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() > 1 and any(
                    nbr.GetAtomicNum() == 1 and nbr.GetTotalNumHs() > 0
                    for nbr in atom.GetNeighbors()
                ):
                    heavy_atom_idx = idx
                    break
            if heavy_atom_idx is not None:
                acidic_sites.add(heavy_atom_idx)

    return sorted(list(acidic_sites))

def protonate_site(mol, atom_idx):
    """
    Add a proton (H) to the acidic site atom, return new mol.
    """
    mol_copy = Chem.Mol(mol)
    mol_copy = Chem.AddHs(mol_copy)
    atom = mol_copy.GetAtomWithIdx(atom_idx)

    # Add new H bonded to this atom
    mol_copy.GetAtomWithIdx(atom_idx).SetFormalCharge(0)
    new_h = Chem.Atom(1)
    mol_copy.AddAtom(new_h)
    new_h_idx = mol_copy.GetNumAtoms() - 1
    mol_copy.AddBond(atom_idx, new_h_idx, order=Chem.rdchem.BondType.SINGLE)

    # Sanitize and optimize
    Chem.SanitizeMol(mol_copy)
    AllChem.EmbedMolecule(mol_copy)
    AllChem.UFFOptimizeMolecule(mol_copy)

    return mol_copy

def deprotonate_site(mol, atom_idx):
    """
    Remove acidic proton bonded to atom_idx, return new mol.
    """
    mol_copy = Chem.Mol(mol)
    mol_copy = Chem.AddHs(mol_copy)
    atom = mol_copy.GetAtomWithIdx(atom_idx)

    # Find an H neighbor to remove
    h_idx = None
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            h_idx = nbr.GetIdx()
            break
    if h_idx is not None:
        mol_copy = Chem.RWMol(mol_copy)
        mol_copy.RemoveAtom(h_idx)
        mol_copy = mol_copy.GetMol()
    else:
        # No H found, return original molecule
        pass

    # Adjust formal charge for deprotonation if needed (assume -1)
    atom = mol_copy.GetAtomWithIdx(atom_idx)
    atom.SetFormalCharge(atom.GetFormalCharge() - 1)

    Chem.SanitizeMol(mol_copy)
    AllChem.EmbedMolecule(mol_copy)
    AllChem.UFFOptimizeMolecule(mol_copy)

    return mol_copy

def mol_to_xyz_string(mol):
    """
    Convert RDKit mol to xyz format string.
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), "",]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
    return "\n".join(lines)
