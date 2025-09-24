from rdkit import Chem
from rdkit.Chem import BRICS
from collections import Counter

# Example SMILES — replace these with your dataset for better accuracy
smiles_list = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",      # Aspirin
    "CCN(CC)CCCC(C)NC1=NC2=CC=CC=C2C(=N1)N",  # Drug-like
    "CC(C)NCC(O)COc1ccc2nc(S(N)(=O)=O)sc2c1", # Drug-like
    "CC1=C(C(=O)NC2=CC=CC=C2)N=C(N1)N",      # Allopurinol-like
    "C1=CC=C(C=C1)C=O",                      # Benzaldehyde
]

# Convert to RDKit Mol objects
mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]

# Count BRICS fragments
fragment_counter = Counter()
for mol in mols:
    try:
        fragments = BRICS.BRICSDecompose(mol)
        fragment_counter.update(fragments)
    except Exception as e:
        print(f"BRICSDecompose failed: {e}")

# Helper to count attachment points (i.e., "*")
def count_attachment_points(frag_smiles):
    return frag_smiles.count("*")

# Split into scaffolds and R-groups
scaffolds = []
rgroups = []

for frag, count in fragment_counter.items():
    attachment_points = count_attachment_points(frag)
    if attachment_points >= 2:
        scaffolds.append((frag, count))
    elif attachment_points == 1:
        rgroups.append((frag, count))

# Output result
print("Top Scaffolds:")
for s in sorted(scaffolds, key=lambda x: -x[1])[:10]:
    print(f"  {s[0]} — {s[1]} uses")

print("\nTop R-groups:")
for r in sorted(rgroups, key=lambda x: -x[1])[:10]:
    print(f"  {r[0]} — {r[1]} uses")

