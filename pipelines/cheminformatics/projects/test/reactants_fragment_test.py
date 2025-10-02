from rdkit import Chem

def reassemble_lead_substituent_fixed(lead_smi, core_smarts, sub_smarts):
    lead_mol = Chem.MolFromSmiles(lead_smi)
    core_mol = Chem.MolFromSmarts(core_smarts)
    sub_mol = Chem.MolFromSmarts(sub_smarts)

    if not (lead_mol and core_mol and sub_mol):
        print("Invalid molecules or SMARTS.")
        return None

    match = lead_mol.GetSubstructMatch(core_mol)
    if not match:
        print("No core match in lead")
        return None

    dummy_core_idx = [a.GetIdx() for a in core_mol.GetAtoms() if a.GetAtomicNum() == 0]
    if not dummy_core_idx:
        print("No dummy atom in core")
        return None
    dummy_core_idx = dummy_core_idx[0]

    lead_dummy_idx = match[dummy_core_idx]

    lead_rw = Chem.RWMol(lead_mol)
    lead_dummy_atom = lead_rw.GetAtomWithIdx(lead_dummy_idx)
    lead_neighbors = [nbr.GetIdx() for nbr in lead_dummy_atom.GetNeighbors()]
    if not lead_neighbors:
        print("Lead dummy atom has no neighbors")
        return None
    lead_neighbor_idx = lead_neighbors[0]

    # Remove dummy atom from lead first
    lead_rw.RemoveAtom(lead_dummy_idx)
    lead_no_dummy = lead_rw.GetMol()

    # Similarly for substituent
    sub_dummy_atoms = [a.GetIdx() for a in sub_mol.GetAtoms() if a.GetAtomicNum() == 0]
    if not sub_dummy_atoms:
        print("No dummy atom in substituent")
        return None
    sub_dummy_idx = sub_dummy_atoms[0]

    sub_rw = Chem.RWMol(sub_mol)
    sub_dummy_atom = sub_rw.GetAtomWithIdx(sub_dummy_idx)
    sub_neighbors = [nbr.GetIdx() for nbr in sub_dummy_atom.GetNeighbors()]
    if not sub_neighbors:
        print("Substituent dummy atom has no neighbors")
        return None
    sub_neighbor_idx = sub_neighbors[0]

    # Remove dummy atom from substituent
    sub_rw.RemoveAtom(sub_dummy_idx)
    sub_no_dummy = sub_rw.GetMol()

    # Combine lead and substituent without dummy atoms
    combo = Chem.CombineMols(lead_no_dummy, sub_no_dummy)
    em = Chem.EditableMol(combo)

    lead_num_atoms = lead_no_dummy.GetNumAtoms()

    # Add bond between lead_neighbor_idx and substituent neighbor shifted by lead atoms count
    em.AddBond(lead_neighbor_idx, lead_num_atoms + sub_neighbor_idx, Chem.rdchem.BondType.SINGLE)

    combined = em.GetMol()

    try:
        Chem.SanitizeMol(combined)
    except Exception as e:
        print(f"Sanitize failed: {e}")
        return None

    combined_smi = Chem.MolToSmiles(combined, canonical=True)
    print(f"Final SMILES: {combined_smi}")
    return combined_smi


if __name__ == "__main__":
    lead = "C[C@H]1CN(Cc2c(Cl)cccc2C(F)F)C[C@@]1(CC(=O)O)C(=O)NC1CCN(CC2CCCCC2)CC1"
    core = "[*:1]C1CNCC1"
    substituent = "[*:1]Clc1ccccc1"

    reassemble_lead_substituent_fixed(lead, core, substituent)

