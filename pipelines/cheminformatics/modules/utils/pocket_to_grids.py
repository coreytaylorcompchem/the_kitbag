import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, Superimposer

def extract_tm_ca_atoms(structure, pdb_to_uni, gpcrdb_segments):
    """
    Returns dict:
      uni_residue_number -> CA Atom
    Only for TM residues
    """
    atoms = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.get_id()
                if hetflag.strip():
                    continue

                key = (chain.id.strip(), resseq)
                if key not in pdb_to_uni:
                    continue

                _, uni_res = pdb_to_uni[key]
                seg = gpcrdb_segments.get(uni_res)

                if seg and seg.startswith("TM") and "CA" in residue:
                    atoms[uni_res] = residue["CA"]

    return atoms

def compute_alignment(reference_pdb, mobile_pdb, pdb_to_uni, gpcrdb_segments):
    """
    Aligns structures using common TM Cα atoms mapped by UniProt residue number.
    Returns Bio.PDB Superimposer
    """
    parser = PDBParser(QUIET=True)

    ref_struct = parser.get_structure("ref", reference_pdb)
    mob_struct = parser.get_structure("mob", mobile_pdb)

    ref_atoms_dict = extract_tm_ca_atoms(
        ref_struct, pdb_to_uni, gpcrdb_segments
    )
    mob_atoms_dict = extract_tm_ca_atoms(
        mob_struct, pdb_to_uni, gpcrdb_segments
    )

    # Intersection of UniProt residues
    common_residues = sorted(
        set(ref_atoms_dict.keys()) & set(mob_atoms_dict.keys())
    )

    if len(common_residues) < 20:
        raise ValueError(
            f"Not enough common TM residues for alignment "
            f"({len(common_residues)} found)"
        )

    ref_atoms = [ref_atoms_dict[r] for r in common_residues]
    mob_atoms = [mob_atoms_dict[r] for r in common_residues]

    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)

    return sup


def transform_pocket_atoms(pocket_pdb, superimposer):
    """
    Returns Nx3 numpy array of transformed atom coordinates
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", pocket_pdb)

    coords = []

    for atom in structure.get_atoms():
        if atom.element != "H":
            coord = atom.get_coord()
            new_coord = np.dot(coord, superimposer.rotran[0]) + superimposer.rotran[1]
            coords.append(new_coord)

    return np.array(coords)


def voxelize(coords, spacing=1.0):
    """
    Convert coordinates into a density grid
    """
    min_xyz = coords.min(axis=0) - spacing
    max_xyz = coords.max(axis=0) + spacing

    shape = np.ceil((max_xyz - min_xyz) / spacing).astype(int)
    grid = np.zeros(shape, dtype=np.float32)

    indices = ((coords - min_xyz) / spacing).astype(int)
    for x, y, z in indices:
        grid[x, y, z] += 1

    return grid, min_xyz, spacing

def write_dx(grid, origin, spacing, out_file):
    nx, ny, nz = grid.shape

    with open(out_file, "w") as f:
        f.write("object 1 class gridpositions counts ")
        f.write(f"{nx} {ny} {nz}\n")
        f.write(f"origin {origin[0]} {origin[1]} {origin[2]}\n")
        f.write(f"delta {spacing} 0 0\n")
        f.write(f"delta 0 {spacing} 0\n")
        f.write(f"delta 0 0 {spacing}\n")
        f.write(f"object 2 class gridconnections counts ")
        f.write(f"{nx} {ny} {nz}\n")
        f.write(f"object 3 class array type double rank 0 items {grid.size} data follows\n")

        flat = grid.flatten(order="C")
        for i in range(0, len(flat), 3):
            f.write(" ".join(f"{v:.3f}" for v in flat[i:i+3]) + "\n")

        f.write("object \"density\" class field\n")
        f.write("component \"positions\" value 1\n")
        f.write("component \"connections\" value 2\n")
        f.write("component \"data\" value 3\n")

def build_pocket_cluster_grid(
    cluster_rows,
    structure_dir,
    reference_pdb_id,
    pdb_to_uni,
    gpcrdb_segments,
    out_dir,
    spacing=1.0
):
    """
    cluster_rows: rows belonging to ONE pocket_cluster_id
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_pdb = structure_dir / f"{reference_pdb_id}.pdb"

    all_coords = []

    for row in cluster_rows:
        pdb_id = row["pdb_id"]
        pocket = row["pocket"]

        pocket_pdb = structure_dir / f"{pdb_id}_out" / f"{pocket}_env_atm.pdb"
        mob_pdb = structure_dir / f"{pdb_id}.pdb"

        sup = compute_alignment(ref_pdb, mob_pdb, pdb_to_uni, gpcrdb_segments)
        coords = transform_pocket_atoms(pocket_pdb, sup)
        all_coords.append(coords)

    all_coords = np.vstack(all_coords)
    grid, origin, spacing = voxelize(all_coords, spacing=spacing)

    out_file = out_dir / f"{cluster_rows[0]['pocket_cluster_id']}_density.dx"
    write_dx(grid, origin, spacing, out_file)

    return out_file