import numpy as np
from Bio.PDB import PDBParser, PDBIO

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

## NONE OF THIS SEEMS TO BE NEEDED ANY MORE SINCE WE'RE NO LONGER USING DSSP
# # Simple Kyte–Doolittle hydrophobicity scale - for orientation
# HYDROPHOBICITY = {
#     "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
#     "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
#     "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2,
#     "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5
# }

# def helix_hydrophobicity(helix_residues):
#     """Simple hydrophobicity score to find alpha helices
#     """
#     scores = [HYDROPHOBICITY.get(res.resname, 0.0) for res in helix_residues]
#     return np.mean(scores) if scores else 0.0

# def compute_helix_axis(residues):
#     """Compute axis of a helix via PCA of Cα coordinates."""
#     coords = np.array([atom.coord for res in residues for atom in res if atom.name == "CA"])
#     if coords.shape[0] < 2:
#         return np.array([0, 0, 1])
#     coords -= coords.mean(axis=0)
#     _, _, vh = np.linalg.svd(coords)
#     return vh[0]

# def compute_bundle_axis(helices):
#     axes = [compute_helix_axis(h) for h in helices if h]
#     if not axes:
#         return np.array([0, 0, 1])
#     avg_axis = np.mean(axes, axis=0)
#     return avg_axis / np.linalg.norm(avg_axis)

# def align_to_z(coords, axis):
#     """Rotate coordinates so that 'axis' aligns with Z = (0,0,1).
#     This is a necessary condition for OpenMMs addMembrane() module (coords will fall to nan otherwise).
#     """
#     z = np.array([0, 0, 1])
#     axis = axis / np.linalg.norm(axis)
#     v = np.cross(axis, z)
#     c = np.dot(axis, z)
#     if np.linalg.norm(v) < 1e-6:
#         return coords
#     vx = np.array([[0, -v[2], v[1]],
#                    [v[2], 0, -v[0]],
#                    [-v[1], v[0], 0]])
#     rot = np.eye(3) + vx + vx @ vx * (1 / (1 + c))
#     return np.dot(coords, rot.T)

# def detect_helices_fallback(structure, window_size=6, elongation_ratio=2.0):
#     """Fallback helix detection using sliding window PCA along backbone."""
#     helices = []
#     residues = [res for res in structure.get_residues()]
#     for i in range(len(residues) - window_size + 1):
#         window = residues[i:i+window_size]
#         ca_coords = np.array([atom.coord for res in window for atom in res if atom.name == "CA"])
#         if ca_coords.shape[0] < 2:
#             continue
#         coords_centered = ca_coords - ca_coords.mean(axis=0)
#         _, s, _ = np.linalg.svd(coords_centered)
#         if s[0] / max(s[1], 1e-6) > elongation_ratio:
#             helices.append(window)
#     return helices

def orient_gpcr_with_ligand(pdb_path, output_path, ligand_resnames=None, center_z=True):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)

    # CA atoms
    protein_atoms = [a for a in structure.get_atoms() if a.get_parent().get_id()[0] == " " and a.name=="CA"]
    coords_protein = np.array([a.coord for a in protein_atoms])

    # PCA : protein principal axis
    coords_centered = coords_protein - coords_protein.mean(axis=0)
    _, _, vh = np.linalg.svd(coords_centered)
    axis = vh[0]
    if coords_protein[:,2].mean() < 0:
        axis = -axis  # point "up"

    # Rotation matrix to align axis with Z
    z_axis = np.array([0,0,1])
    v = np.cross(axis, z_axis)
    s = np.linalg.norm(v)
    if s < 1e-6:
        R = np.eye(3)
    else:
        c = np.dot(axis, z_axis)
        vx = np.array([[0, -v[2], v[1]],
                       [v[2], 0, -v[0]],
                       [-v[1], v[0], 0]])
        R = np.eye(3) + vx + np.dot(vx, vx) * ((1-c)/(s**2))

    # COM of entire complex (protein + ligand)
    if ligand_resnames:
        complex_atoms = [a for a in structure.get_atoms() if (a.get_parent().get_id()[0] == " " or a.get_parent().resname in ligand_resnames)]
    else:
        complex_atoms = [a for a in structure.get_atoms() if a.get_parent().get_id()[0] == " "]
    com = np.mean([a.coord for a in complex_atoms], axis=0)

    # Rotate all atoms around COM
    for atom in complex_atoms:
        atom.coord = np.dot(R, atom.coord - com) + com

    if center_z:
        # Center TM bundle midplane at Z=0 using protein CA atoms
        ca_coords = np.array([a.coord for a in protein_atoms])
        z_mid = (ca_coords[:,2].max() + ca_coords[:,2].min()) / 2
        for atom in complex_atoms:
            atom.coord[2] -= z_mid

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)

    ligand_atoms_debug = [
        a for a in structure.get_atoms()
        if a.get_parent().resname in (ligand_resnames or [])
    ]
    
    logger.debug(f"Ligand atoms seen by orienter: {len(ligand_atoms_debug)}")

    return output_path
