from rdkit import Chem

# Relaxed SMARTS
core_smarts = "[*]c1cc(=O)oc2[nH]c(=S)[nH]c(=O)c12"

# Your lead molecule
lead_smiles = "CCCCc1cc(=O)oc2[nH]c(=S)[nH]c(=O)c12"

core = Chem.MolFromSmarts(core_smarts)
lead = Chem.MolFromSmiles(lead_smiles)

match = lead.GetSubstructMatch(core)
if match:
    print("✅ Match found!")
else:
    print("❌ Still no match.")

