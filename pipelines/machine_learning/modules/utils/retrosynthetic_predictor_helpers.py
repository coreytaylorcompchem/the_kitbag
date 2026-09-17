from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

from modules.utils.retrosynth_parameter_dicts import REACTION_HANDLE_SMARTS


def _process_building_block_chunk(records, smiles_col):
    """
    Process a chunk of building-block records in a worker process.

    Individual molecules that cause RDKit processing errors are marked
    invalid rather than terminating the entire chunk.
    """
    processed = []

    for row in records:
        smiles = row[smiles_col]

        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                processed.append({
                    **row,
                    "valid_smiles": False,
                })
                continue

            # Canonicalise and calculate descriptors
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            inchikey = Chem.MolToInchiKey(mol)

            elements = sorted({
                atom.GetSymbol()
                for atom in mol.GetAtoms()
            })

            processed.append({
                **row,
                "smiles": canonical_smiles,
                "canonical_smiles": canonical_smiles,
                "inchikey": inchikey,
                "valid_smiles": True,
                "elements": elements,
                "molecular_weight": Descriptors.MolWt(mol),
                "logp": Crippen.MolLogP(mol),
                "tpsa": Descriptors.TPSA(mol),
                "hbd": Lipinski.NumHDonors(mol),
                "hba": Lipinski.NumHAcceptors(mol),
                "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
                "ring_count": Lipinski.RingCount(mol),
                "heavy_atom_count": Lipinski.HeavyAtomCount(mol),
            })

        except Exception:
            # RDKit can fail on chemically problematic structures even
            # after MolFromSmiles() has returned a molecule.
            processed.append({
                **row,
                "valid_smiles": False,
            })

    return processed


REACTION_HANDLE_PATTERNS = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in REACTION_HANDLE_SMARTS.items()
}


def _annotate_building_block_chunk(records):
    """
    Annotate a chunk of building blocks with reaction handles.

    Individual RDKit failures are caught so that one problematic molecule
    cannot terminate the entire chunk.
    """

    annotated = []

    for row in records:

        smiles = row.get("canonical_smiles") or row.get("smiles")

        try:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                annotated.append({
                    **row,
                    "reaction_handles": "",
                    "n_reaction_handles": 0,
                    "has_reaction_handle": False,
                })
                continue

            handles = []

            for handle_name, pattern in REACTION_HANDLE_PATTERNS.items():

                if mol.HasSubstructMatch(pattern):
                    handles.append(handle_name)

            annotated.append({
                **row,
                "reaction_handles": ";".join(handles),
                "n_reaction_handles": len(handles),
                "has_reaction_handle": bool(handles),
            })

        except Exception:
            annotated.append({
                **row,
                "reaction_handles": "",
                "n_reaction_handles": 0,
                "has_reaction_handle": False,
            })

    return annotated