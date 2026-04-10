from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import defaultdict
import random

def generate_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)


def scaffold_split(graphs, smiles_list, val_fraction=0.2, seed=42):
    random.seed(seed)

    scaffold_to_indices = defaultdict(list)

    for i, smi in enumerate(smiles_list):
        scaffold = generate_scaffold(smi)
        scaffold_to_indices[scaffold].append(i)

    # sort scaffolds by size (largest first)
    scaffold_sets = sorted(
        scaffold_to_indices.values(),
        key=lambda x: len(x),
        reverse=True
    )

    train_idx, val_idx = [], []
    total = len(graphs)
    val_target = int(total * val_fraction)

    for scaffold in scaffold_sets:
        if len(val_idx) + len(scaffold) <= val_target:
            val_idx.extend(scaffold)
        else:
            train_idx.extend(scaffold)

    train_graphs = [graphs[i] for i in train_idx]
    val_graphs = [graphs[i] for i in val_idx]

    return train_graphs, val_graphs