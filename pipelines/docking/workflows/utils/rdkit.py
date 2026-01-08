from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, QED

def compute_rdkit_descriptors(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}

    return {
        "MW": Descriptors.MolWt(mol),
        "QED": QED.qed(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "LogP": Crippen.MolLogP(mol),
        "RotB": Lipinski.NumRotatableBonds(mol),
        "Rings": Lipinski.RingCount(mol),
    }

