from rdkit import Chem

def sanitize_smiles_list(ligands):
    """Filter ligands by RDKit MolFromSmiles and full sanitization (SanitizeMol)."""
    valid_ligands = []
    invalid_ligands = []

    for lig in ligands:
        try:
            mol = Chem.MolFromSmiles(lig['smiles'])
            if mol is None:
                raise ValueError("MolFromSmiles failed")
            Chem.SanitizeMol(mol)  # Raises exception on failure
            valid_ligands.append(lig)
        except Exception as e:
            lig['sanitization_error'] = str(e)
            invalid_ligands.append(lig)

    return valid_ligands, invalid_ligands

