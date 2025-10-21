# utils/graph_construction.py

import os
import warnings
import contextlib
import numpy as np
import torch
from torch_geometric.data import Data

import MDAnalysis as mda
from MDAnalysis.analysis import distances, rms
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
try:
    from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
except ImportError:
    SASAAnalysis = None
    warnings.warn("mdakit-sasa SASAAnalysis not available; ΔSASA metric will be disabled.")
from rdkit import Chem
from rdkit.Chem import AllChem

# --- Atomic feature helpers ---

def atom_features(atom):
    """Return a richer set of atomic descriptors."""
    feats = [
        atom.GetAtomicNum(),
        atom.GetTotalDegree(),
        atom.GetTotalValence(),
        atom.GetFormalCharge(),
        int(atom.GetIsAromatic()),
        atom.GetTotalNumHs(includeNeighbors=True),
        atom.GetImplicitValence(),
    ]
    hyb = atom.GetHybridization()
    hyb_onehot = [int(hyb == h) for h in (
        Chem.rdchem.HybridizationType.SP,
        Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3,
        Chem.rdchem.HybridizationType.SP3D,
        Chem.rdchem.HybridizationType.SP3D2,
    )]
    feats.extend(hyb_onehot)
    return feats


def atomic_number(atom):
    try:
        return Chem.GetPeriodicTable().GetAtomicNumber(atom.element)
    except Exception:
        return 0


# --- File discovery ---

def find_topology_and_trajectory(struct_dir):
    topo_exts = (".pdb", ".psf", ".gro", ".prmtop", ".top")
    traj_exts = (".dcd", ".xtc", ".trr", ".nc", ".mdcrd")

    topo_files, traj_files = [], []
    for root, _, files in os.walk(struct_dir):
        for f in files:
            path = os.path.join(root, f)
            if f.lower().endswith(topo_exts):
                topo_files.append(path)
            elif f.lower().endswith(traj_exts):
                traj_files.append(path)

    pick = lambda lst: max(lst, key=os.path.getmtime) if lst else None
    return pick(topo_files), pick(traj_files)


# --- Trajectory → Graphs conversion ---

def load_trajectory_as_graphs(pdb_path, traj_path, ligand_resname="UNK", cutoff=5.0):
    """
    Load MD trajectory as per-frame ligand graphs with enriched features.
    """
    try:
        u = mda.Universe(pdb_path, traj_path)
    except Exception as e:
        warnings.warn(f"Failed to load Universe from {pdb_path} + {traj_path}: {e}")
        return []

    ligand = u.select_atoms(f"resname {ligand_resname}")
    protein = u.select_atoms("protein")
    if len(ligand) == 0:
        warnings.warn(f"No ligand atoms found in {pdb_path} (resname {ligand_resname}); skipping.")
        return []

    # Load RDKit molecule for features
    rdkit_mol = None
    try:
        rdkit_mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        if rdkit_mol:
            AllChem.ComputeGasteigerCharges(rdkit_mol)
    except Exception as e:
        warnings.warn(f"RDKit feature extraction failed for {pdb_path}: {e}")
        rdkit_mol = None

    try:
        ligand.guess_bonds()
    except Exception as e:
        warnings.warn(f"Failed to guess bonds for {pdb_path}: {e}")
        return []

    atom_index_map = {atom.index: i for i, atom in enumerate(ligand.atoms)}
    graphs = []

    for ts in u.trajectory:
        coords = ligand.positions
        n_atoms = len(ligand)

        # Base atomic features
        if rdkit_mol and rdkit_mol.GetNumAtoms() == n_atoms:
            atom_feats = np.array([atom_features(a) for a in rdkit_mol.GetAtoms()], dtype=float)
        else:
            atom_feats = np.array([[atomic_number(a)] for a in ligand], dtype=float)

        # Protein context
        context = np.zeros((n_atoms, 3), dtype=float)
        if len(protein) > 0:
            try:
                dist = distances.distance_array(ligand.positions, protein.positions)
                close = dist < cutoff
                for i in range(n_atoms):
                    contacts = np.where(close[i])[0]
                    if len(contacts) == 0:
                        continue
                    mean_dist = np.mean(dist[i, contacts])
                    n_contacts = len(contacts)
                    resnames = protein[contacts].resnames
                    polar_count = sum(r in ("SER", "THR", "ASN", "GLN", "HIS") for r in resnames)
                    context[i] = [n_contacts, mean_dist, polar_count]
            except Exception as e:
                warnings.warn(f"Protein context failed for {pdb_path}: {e}")

        # Concatenate
        x = torch.tensor(np.concatenate([atom_feats, context], axis=1), dtype=torch.float)

        # Build edges
        edge_index = []
        for bond in ligand.bonds:
            try:
                i, j = [atom_index_map[a.index] for a in bond.atoms]
                edge_index.extend([[i, j], [j, i]])
            except KeyError:
                continue
        if not edge_index:
            continue

        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        pos = torch.tensor(coords, dtype=torch.float)

        graphs.append(Data(x=x, edge_index=edge_index, pos=pos))

    return graphs


# --- Pose metrics and helpers ---

@contextlib.contextmanager
def suppress_freesasa_output():
    old_stdout_fd = os.dup(1)
    old_stderr_fd = os.dup(2)
    try:
        with open(os.devnull, 'w') as devnull:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
    finally:
        os.dup2(old_stdout_fd, 1)
        os.dup2(old_stderr_fd, 2)
        os.close(old_stdout_fd)
        os.close(old_stderr_fd)

def interface_rmsd(current_ligand, ref_ligand, ref_protein, cutoff=5.0):
    """Compute interface RMSD using only ligand atoms near the protein in the reference."""
    ref_contacts = distances.distance_array(ref_ligand.positions, ref_protein.positions) < cutoff
    interface_atoms = np.any(ref_contacts, axis=1)
    if not np.any(interface_atoms):
        return np.nan
    return rms.rmsd(
        current_ligand.positions[interface_atoms],
        ref_ligand.positions[interface_atoms],
        superposition=True
    )

def geometric_hbond_analysis(u, donors_sel, acceptors_sel, d_a_cutoff=3.5, angle_cutoff=120.0):
    """
    Robust geometric H-bond detector for cases without explicit bonding info.
    Falls back to nearest hydrogen within 1.2 Å of donor.
    """
    donors = u.select_atoms(donors_sel)
    acceptors = u.select_atoms(acceptors_sel)
    hydrogens = u.select_atoms("name H* or element H")

    n_frames = len(u.trajectory)
    hb_counts = np.zeros(n_frames)
    hb_mean_dist = np.full(n_frames, np.nan)
    hb_mean_angle = np.full(n_frames, np.nan)

    for ts in u.trajectory:
        hbonds_in_frame = []
        dist_da = distances.distance_array(donors.positions, acceptors.positions)
        donor_idx, acceptor_idx = np.where(dist_da < d_a_cutoff)

        for d_i, a_i in zip(donor_idx, acceptor_idx):
            donor = donors[d_i]
            acceptor = acceptors[a_i]
            d_a_dist = dist_da[d_i, a_i]

            # Find the closest hydrogen near donor (<1.2 Å)
            dH_dists = distances.distance_array(donor.position[None, :], hydrogens.positions)[0]
            near_H_idx = np.where(dH_dists < 1.2)[0]

            if len(near_H_idx) > 0:
                Hs = hydrogens[near_H_idx]
                angles = []
                for H in Hs:
                    v1 = donor.position - H.position
                    v2 = acceptor.position - H.position
                    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    cos_angle = np.clip(cos_angle, -1.0, 1.0)
                    angles.append(np.degrees(np.arccos(cos_angle)))
                mean_angle = np.mean(angles)
            else:
                mean_angle = np.nan  # No H found nearby

            if np.isnan(mean_angle) or mean_angle > angle_cutoff:
                hbonds_in_frame.append((d_i, a_i, d_a_dist, mean_angle))

        if hbonds_in_frame:
            arr = np.array(hbonds_in_frame, dtype=object)
            hb_counts[ts.frame] = len(arr)
            hb_mean_dist[ts.frame] = np.nanmean(arr[:, 2].astype(float))
            hb_mean_angle[ts.frame] = np.nanmean(arr[:, 3].astype(float))
        else:
            hb_counts[ts.frame] = 0

    hb_mean_dist = np.nan_to_num(hb_mean_dist, nan=0.0)
    hb_mean_angle = np.nan_to_num(hb_mean_angle, nan=0.0)

    return hb_counts, hb_mean_dist, hb_mean_angle



def compute_pose_metrics(u, ligand_sel="resname UNK", cutoff=4.5, ref_frame=0):
    """Compute pose quality metrics with NaN handling and diagnostics."""

    ligand = u.select_atoms(ligand_sel)
    protein = u.select_atoms("protein")
    if len(ligand) == 0 or len(protein) == 0:
        warnings.warn("Could not find ligand or protein atoms")
        return np.empty((0, 9)), ["drmsd","irmsd","fnat","qscore","dsasa","energy_proxy","hbonds","hb_mean_dist","hb_mean_angle"]

    # Reference config
    u.trajectory[ref_frame]
    ref_coords = ligand.positions.copy()
    ref_dists = distances.distance_array(ref_coords, protein.positions)
    ref_contacts = ref_dists < cutoff

    ref_ligand = ligand.copy()
    ref_protein = protein.copy()

    # SASA
    try:
        with suppress_freesasa_output():
            sasa_analysis = SASAAnalysis(u, select=ligand_sel, probe_radius=1.4)
            sasa_analysis.run()
        ref_sasa_lig = sasa_analysis.results.total_area[ref_frame]
    except Exception as e:
        warnings.warn(f"SASA unavailable or failed: {e}")
        sasa_analysis = None
        ref_sasa_lig = 0.0

    ref_qdist = ref_dists.copy()
    
    # Hydrogen bonds (enhanced)
    n_frames = len(u.trajectory)
    hb_counts = np.zeros(n_frames)
    hb_mean_dist = np.zeros(n_frames)
    hb_mean_angle = np.zeros(n_frames)

    try:
        donors_sel = "protein and (name N* or name NE* or name ND* or name NZ or name OG* or name OH*)"
        acceptors_sel = f"{ligand_sel} and (element O or element N or name O* or name N*)"
        \
        # Check if charges exist
        has_charges = hasattr(u.atoms, "charges") and u.atoms.charges is not None

        if not has_charges:
            print("Universe missing charges — creating a new Universe with dummy zero charges for H-bond detection.")
            # Create new Universe from original topology and trajectory files
            topo_file = u.filename  # Usually the topology file path
            traj_file = None
            if hasattr(u.trajectory, "filename"):
                traj_file = u.trajectory.filename  # Trajectory file path if exists

            if traj_file:
                u_copy = mda.Universe(topo_file, traj_file)
            else:
                u_copy = mda.Universe(topo_file)

            u_copy.add_TopologyAttr("charges", np.zeros(len(u_copy.atoms)))
        else:
            u_copy = u

        h = HydrogenBondAnalysis(
            u_copy,
            donors_sel=donors_sel,
            acceptors_sel=acceptors_sel,
            d_h_cutoff=1.5,
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=120,
        )
        h.run()

        # We try HBonds from MDA first but, if that fails, we have a backup.

        if hasattr(h.results, "hbonds") and len(h.results.hbonds) > 0:
            hb_array = np.array(h.results.hbonds)
            frames, counts = np.unique(hb_array[:, 0].astype(int), return_counts=True)
            hb_counts[frames] = counts
            hb_mean_dist[:] = np.mean(hb_array[:, 7]) if hb_array.shape[1] > 7 else 0
            hb_mean_angle[:] = np.mean(hb_array[:, 8]) if hb_array.shape[1] > 8 else 0
        else:
            print("No H-bonds detected by HydrogenBondAnalysis — falling back to geometric method.")
            hb_counts, hb_mean_dist, hb_mean_angle = geometric_hbond_analysis(
                u_copy, donors_sel, acceptors_sel
            )

    except Exception as e:
        print(f"H-bond analysis failed: {e}")
        hb_counts, hb_mean_dist, hb_mean_angle = geometric_hbond_analysis(
            u_copy, donors_sel, acceptors_sel
        )

    metrics = []
    for ts in u.trajectory:
        try:
            drmsd = rms.rmsd(ligand.positions, ref_coords, superposition=True)
        except Exception:
            drmsd = np.nan

        # Interface RMSD
        try:
            ref_contacts_local = distances.distance_array(ref_ligand.positions, ref_protein.positions) < cutoff
            interface_atoms = np.any(ref_contacts_local, axis=1)
            if not np.any(interface_atoms):
                irmsd = 10.0  # no interface → large fallback
            else:
                irmsd = rms.rmsd(
                    ligand.positions[interface_atoms],
                    ref_ligand.positions[interface_atoms],
                    superposition=True
                )
        except Exception:
            irmsd = 10.0

        dist = distances.distance_array(ligand.positions, protein.positions)
        dist = np.where(dist < 1e-3, 1e-3, dist)  # avoid div by 0
        contacts = dist < cutoff

        fnat = (contacts & ref_contacts).sum() / (ref_contacts.sum() + 1e-8)
        qscore = np.mean(np.exp(-((dist - ref_qdist) ** 2) / (2 * 1.0 ** 2)))

        if sasa_analysis is not None:
            try:
                sasa_lig = sasa_analysis.results.total_area[ts.frame]
                delta_sasa = ref_sasa_lig - sasa_lig
            except Exception:
                delta_sasa = 0.0
        else:
            delta_sasa = 0.0

        try:
            inv_dist = np.clip(1.0 / dist, 0, 1e2)
            energy_proxy = np.sum(-inv_dist**6 + inv_dist**12)
        except Exception:
            energy_proxy = 0.0

        frame = ts.frame
        n_hbonds = hb_counts[frame] if frame < len(hb_counts) else 0.0
        mean_hb_dist = hb_mean_dist[frame] if frame < len(hb_mean_dist) else 0.0
        mean_hb_angle = hb_mean_angle[frame] if frame < len(hb_mean_angle) else 0.0

        frame_metrics = [drmsd, irmsd, fnat, qscore, delta_sasa, energy_proxy, n_hbonds, mean_hb_dist, mean_hb_angle]
        metrics.append(frame_metrics)

    metrics = np.array(metrics, dtype=float)
    colnames = ["drmsd", "irmsd", "fnat", "qscore", "dsasa", "energy_proxy", "hbonds", "hb_mean_dist", "hb_mean_angle"]

    # Clean and report
    nan_counts = np.isnan(metrics).sum(axis=0)
    if np.any(nan_counts):
        bad_info = {colnames[i]: int(nc) for i, nc in enumerate(nan_counts) if nc > 0}
        warnings.warn(f"NaNs detected in metrics: {bad_info}")
    metrics = np.nan_to_num(metrics, nan=0.0, posinf=1e6, neginf=-1e6)

    return metrics, colnames

def attach_metrics_to_graphs(pdb_path, traj_path, graphs, ligand_resname="UNK"):
    if not graphs:
        warnings.warn(f"No graphs to attach metrics for {pdb_path}")
        return []

    try:
        u = mda.Universe(pdb_path, traj_path)
    except Exception as e:
        warnings.warn(f"Universe creation failed for {pdb_path}: {e}")
        return []

    metrics, colnames = compute_pose_metrics(u, ligand_sel=f"resname {ligand_resname}")
    if metrics.size == 0:
        warnings.warn(f"No valid metrics for {pdb_path}, skipping.")
        return []

    n = min(len(graphs), len(metrics))
    graphs = graphs[:n]
    metrics = metrics[:n]

    valid_mask = np.isfinite(metrics).all(axis=1)
    valid_count = np.sum(valid_mask)
    if valid_count < n:
        warnings.warn(f"{pdb_path}: Dropping {n - valid_count} frames with invalid metrics.")
        metrics = metrics[valid_mask]
        graphs = [g for g, ok in zip(graphs, valid_mask) if ok]

    for g, m in zip(graphs, metrics):
        g.y = torch.tensor(m, dtype=torch.float32).unsqueeze(0)
        g.metric_names = colnames

    return graphs


def load_all_targets(base_dir, ligand_resname="UNK"):
    all_graphs = []
    skipped = []
    loaded = 0

    for target in os.listdir(base_dir):
        target_path = os.path.join(base_dir, target)
        if not os.path.isdir(target_path):
            continue

        for pdb_struct in os.listdir(target_path):
            struct_dir = os.path.join(target_path, pdb_struct)
            if not os.path.isdir(struct_dir):
                continue

            try:
                topo, traj = find_topology_and_trajectory(struct_dir)
                if not topo or not traj:
                    skipped.append((struct_dir, "missing topology/trajectory"))
                    continue
                if os.path.getsize(traj) == 0:
                    skipped.append((struct_dir, "empty trajectory"))
                    continue

                graphs = load_trajectory_as_graphs(topo, traj, ligand_resname)
                if not graphs:
                    skipped.append((struct_dir, "no graphs generated"))
                    continue

                graphs = attach_metrics_to_graphs(topo, traj, graphs, ligand_resname)
                if graphs:
                    all_graphs.extend(graphs)
                    loaded += 1
                else:
                    skipped.append((struct_dir, "metrics attachment failed"))

            except Exception as e:
                warnings.warn(f"Failed to process {struct_dir}: {e}")
                skipped.append((struct_dir, str(e)))
                continue

    print(f"\nLoaded {len(all_graphs)} graphs from {loaded} valid trajectories.")
    if skipped:
        print(f"Skipped {len(skipped)} trajectories. Examples:")
        for s in skipped[:10]:
            print(f"  - {s[0]} ({s[1]})")

    return all_graphs