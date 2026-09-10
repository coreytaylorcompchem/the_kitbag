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

            processed.append({
                **row,
                "smiles": canonical_smiles,
                "canonical_smiles": canonical_smiles,
                "inchikey": inchikey,
                "valid_smiles": True,
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