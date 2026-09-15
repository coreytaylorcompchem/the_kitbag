from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski


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

# ---------------------------------------------------------------------------
# Reaction-handle SMARTS
# ---------------------------------------------------------------------------

REACTION_HANDLE_SMARTS = {
    "amine": "[N;H1,H2;!$(N-C=O);!$(N-S=O)]",
    "alcohol": "[O;H1;!$(O-C=O)]",
    "phenol": "[O;H1]-[c]",
    "carboxylic_acid": "C(=O)[O;H1]",
    "acid_chloride": "C(=O)Cl",
    "aldehyde": "[CX3H1](=O)",
    "ketone": "[#6][CX3](=O)[#6]",
    "alkyl_halide": "[CX4][Cl,Br,I]",
    "aryl_vinyl_halide": "[c,C;X3,X2][Cl,Br,I]",
    "boronic_acid": "[B;$(B(O)O)]",
    "boronate_ester": "[B;X3]([O;X2])[O;X2]",
    "sulfonyl_chloride": "S(=O)(=O)Cl",
    "nitrile": "C#N",
}


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