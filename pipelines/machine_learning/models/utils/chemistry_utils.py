from rdkit import Chem
from rdkit.Chem import QED, Descriptors
import numpy as np

def calculate_properties(smiles_list):
    qed_scores, mol_weights, logps, tpsas = [], [], [], []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            qed_scores.append(QED.qed(mol))
            mol_weights.append(Descriptors.MolWt(mol))
            logps.append(Descriptors.MolLogP(mol))
            tpsas.append(Descriptors.TPSA(mol))
        else:
            qed_scores += [0]
            mol_weights += [0]
            logps += [0]
            tpsas += [0]

    return {
        "qed_mean": np.mean(qed_scores),
        "mol_wt_mean": np.mean(mol_weights),
        "logp_mean": np.mean(logps),
        "tpsa_mean": np.mean(tpsas)
    }
