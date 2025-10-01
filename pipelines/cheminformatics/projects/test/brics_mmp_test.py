from rdkit import Chem
from rdkit.Chem import BRICS

# Example molecule
smi = "Cc1ccnc(Nc2nc(-c3cn[nH]c3)c(C)s2)n1"
mol = Chem.MolFromSmiles(smi)

# BRICS Decomposition
frags = list(BRICS.BRICSDecompose(mol))
print("BRICS Fragments:")
for f in frags:
    print(f)

# Try reassembling using your function
def reassemble_from_brics(core_smi, sub_smi):
    core = Chem.MolFromSmiles(core_smi)
    sub = Chem.MolFromSmiles(sub_smi)
    if core is None or sub is None:
        print("Failed to parse core or sub.")
        return None

    core_dummies = [a for a in core.GetAtoms() if a.GetAtomicNum() == 0]
    sub_dummies = [a for a in sub.GetAtoms() if a.GetAtomicNum() == 0]
    print(f"Core dummies: {len(core_dummies)}, Sub dummies: {len(sub_dummies)}")

    if len(core_dummies) < 1 or len(sub_dummies) < 1:
        return None

    core_idx = core_dummies[0].GetIdx()
    sub_idx = sub_dummies[0].GetIdx() + core.GetNumAtoms()

    combo = Chem.CombineMols(core, sub)
    em = Chem.EditableMol(combo)
    em.AddBond(core_idx, sub_idx, Chem.rdchem.BondType.SINGLE)
    newmol = em.GetMol()

    dummy_idxs = [a.GetIdx() for a in newmol.GetAtoms() if a.GetAtomicNum() == 0]
    rw = Chem.RWMol(newmol)
    for idx in sorted(dummy_idxs, reverse=True):
        rw.RemoveAtom(idx)
    newmol = rw.GetMol()
    
    try:
        Chem.SanitizeMol(newmol)
    except Exception as e:
        print("Sanitization failed:", e)
        return None

    print("Reassembled:", Chem.MolToSmiles(newmol))
    return newmol

# Example: pick first fragment as core, second as substituent
if len(frags) >= 2:
    frag_list = sorted(frags, key=lambda x: len(x), reverse=True)
    reassemble_from_brics(frag_list[0], frag_list[1])

