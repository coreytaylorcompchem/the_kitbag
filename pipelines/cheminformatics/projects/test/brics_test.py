from rdkit import Chem
from rdkit.Chem.BRICS import FindBRICSBonds

smi = "CCCCOCCOc1ccc(-c2ccc3c(c2)/C=C(/C(=O)Nc2ccc([S@@+]([O-])Cc4cncn4CCC)cc2)CCCN3CC(C)C)cc1"
mol = Chem.MolFromSmiles(smi)
brics_bonds = FindBRICSBonds(mol)
#print(len(brics_bonds))
for bond in brics_bonds:
    print(bond)
