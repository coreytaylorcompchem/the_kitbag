import requests
import csv

import numpy as np
from pathlib import Path
from collections import Counter
import pandas as pd
from Bio.PDB import PDBParser, Superimposer

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)
parser = PDBParser(QUIET=True)

##### convert pockets to grids

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


#### calculate persistence

def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def cluster_pockets_across_structures(pocket_rows, default_threshold=0.3):
    """
    pocket_rows: list of dicts, each with keys:
      - pdb_id
      - pocket
      - region
      - uni_residues (set of UniProt residue numbers)

    Returns:
      updated pocket_rows (with pocket_cluster_id assigned),
      clusters list
    """

    # Region-specific Jaccard thresholds
    REGION_THRESHOLDS = {
        "7TM": 0.30,
        "Loops": 0.20,
        "ECD": 0.25,
        "H8": 0.25,
        "Unknown": 0.35,
    }

    clusters = []
    cluster_counter = 0

    for row in pocket_rows:
        assigned = False
        res = row["uni_residues"]
        region = row.get("region", "Unknown")

        # Skip empty residue sets safely
        if not res:
            cluster_counter += 1
            cid = f"PC{cluster_counter:02d}"
            row["pocket_cluster_id"] = cid
            clusters.append({
                "id": cid,
                "region": region,
                "members": [row],
                "residues": set()
            })
            continue

        threshold = REGION_THRESHOLDS.get(region, default_threshold)

        for cluster in clusters: 
            if cluster["region"] != region:
                continue

            score = jaccard(res, cluster["residues"])
            if score >= threshold:
                cluster["members"].append(row)
                cluster["residues"] |= res
                row["pocket_cluster_id"] = cluster["id"]
                assigned = True
                break

        if not assigned:
            cluster_counter += 1
            cid = f"PC{cluster_counter:02d}"
            row["pocket_cluster_id"] = cid
            clusters.append({
                "id": cid,
                "region": region,
                "members": [row],
                "residues": set(res)
            })

    return pocket_rows, clusters


def extract_pocket_residues(pocket_atm_pdb):
    """
    Returns set of (chain_id, pdb_resseq) compatible with SIFTS mapping
    """
    structure = parser.get_structure("pocket", pocket_atm_pdb)
    residues = set()

    for res in structure.get_residues():
        hetflag, resseq, icode = res.get_id()

        # Skip hetero / water
        if hetflag.strip():
            continue

        chain = res.get_parent().id.strip()
        # if not chain:
        #     chain = "A"  # fpocket often uses blank chain

        residues.add((chain, resseq))

    return residues

def fetch_sifts_mapping(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def build_pdb_to_uniprot_map(sifts_data):
    """
    Map (chain_id, pdb_resseq) -> (uniprot_id, uniprot_resseq)
    using PDBe SIFTS residue_number (coordinate numbering).
    """
    mapping = {}

    for pdb_id, entry in sifts_data.items():
        for uniprot_id, udata in entry.get("UniProt", {}).items():
            for m in udata.get("mappings", []):
                chain = (m.get("chain_id") or "").strip()
                if not chain:
                    continue

                pdb_start = m["start"].get("residue_number")
                pdb_end = m["end"].get("residue_number")
                uni_start = m.get("unp_start")
                uni_end = m.get("unp_end")

                if None in (pdb_start, pdb_end, uni_start, uni_end):
                    continue

                for pdb_res, uni_res in zip(
                    range(int(pdb_start), int(pdb_end) + 1),
                    range(int(uni_start), int(uni_end) + 1),
                ):
                    mapping[(chain, pdb_res)] = (uniprot_id, uni_res)

    return mapping

def fetch_gpcrdb_protein_from_structure(pdb_id):
    """
    Returns GPCRdb protein entry name (e.g. pe2r4_human)
    """
    url = f"https://gpcrdb.org/services/structure/{pdb_id.upper()}/"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()["protein"]


def fetch_gpcrdb_segments(protein_name):
    """
    Returns mapping:
    UniProt residue number -> GPCRdb protein_segment
    """
    url = f"https://gpcrdb.org/services/residues/{protein_name}/"
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    segments = {}
    for entry in r.json():
        try:
            seqnum = int(entry["sequence_number"])
            seg = entry.get("protein_segment")
            if seg:
                segments[seqnum] = seg
        except (KeyError, ValueError, TypeError):
            continue

    return segments

def classify_pocket_topology(pocket_residues, pdb_to_uni, gpcrdb_segments):
    """
    Collapse GPCRdb residue annotations into biologically meaningful pocket regions.
    Only uses residues that successfully map to GPCRdb via UniProt.
    Returns: (region, confidence, raw_counts)
    """

    # Map each pocket residue to gpcrdb segment (if possible)
    raw_labels = []
    for chain_res in pocket_residues:
        if chain_res not in pdb_to_uni:
            continue

        uniprot_id, uni_res = pdb_to_uni[chain_res]
        # fetch segment for this residue number
        seg = gpcrdb_segments.get(uni_res)
        if seg:
            raw_labels.append(seg)

    if not raw_labels:
        return "Unknown", 0.0, {}

    raw_counts = Counter(raw_labels)

    collapsed = {
        "7TM": 0,
        "Loops": 0,
        "ECD": 0,
        "H8": 0,
    }

    for seg, count in raw_counts.items():
        if seg.startswith("TM"):
            collapsed["7TM"] += count
        elif "CL" in seg:              # ECL / ICL
            collapsed["Loops"] += count
        elif seg in {"N-term", "ECD"}:
            collapsed["ECD"] += count
        elif seg == "H8":
            collapsed["H8"] += count

    total = sum(collapsed.values())
    if total == 0:
        return "Unknown", 0.0, dict(raw_counts)

    region = max(collapsed, key=collapsed.get)
    confidence = round(collapsed[region] / total, 3)

    return region, confidence, dict(raw_counts)



def load_fpocket_descriptors(descriptor_file):
    """
    Load fpocket descriptor table.
    Handles whitespace-delimited fpocket output (-d flag).
    Returns DataFrame indexed by pocket number (string).
    """
    # IMPORTANT: fpocket -d output is whitespace-delimited, not CSV
    df = pd.read_csv(
        descriptor_file,
        delim_whitespace=True,
        comment="#"
    )

    # fpocket column IDs
    pocket_cols = [
        "cav_id",          # most common with -d parameter
        "pocket",
        "pocket_number",
        "pocket_id",
        "ID",
        "id"
    ]

    pocket_col = next((c for c in pocket_cols if c in df.columns), None)

    if pocket_col is None:
        raise ValueError(
            f"Could not find pocket ID column in {descriptor_file.name}. "
            f"Columns found: {list(df.columns)}"
        )

    df[pocket_col] = df[pocket_col].astype(str)
    df = df.set_index(pocket_col)

    return df

def assign_pocket_role(row, orthosteric_pocket):
    if row["pocket"] == orthosteric_pocket:
        return "orthosteric"

    if row["region"] == "Loops":
        return "known_allosteric"

    if row["region"] == "7TM":
        return "candidate_allosteric"

    return "orphan"

@register_task(
    "pocket_classification",
    category="Pocket detection",
    description="Classify GPCR fpocket pockets using GPCRdb topology"
)
def gpcr_pocket_classification(config: dict, data: dict = None) -> dict:

    fpocket_dir = Path(
        data["output_directory"] if data and "output_directory" in data
        else config["run_fpocket"]["output_directory"]
    )

    output_cfg = config.get("output", {})
    output_dir = Path(output_cfg.get("directory", "outputs/gpcr_pocket_classification"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / output_cfg.get(
        "filename", "classified_gpcr_pockets.csv"
    )

    records = []

    for struct_dir in fpocket_dir.glob("*_out"):
        pdb_id = struct_dir.name.replace("_out", "")
        logger.info(f"Classifying pockets for {pdb_id}")

        pocket_files = sorted(struct_dir.glob("pocket*_env_atm.pdb"))
        if not pocket_files:
            logger.warning(f"No env pockets found for {pdb_id}")
            continue

        logger.info(f"Found {len(pocket_files)} pockets for {pdb_id}")

        # Load fpocket descriptors (NOTE: whitespace-delimited)
        descriptor_file = struct_dir / f"{pdb_id}_pocket_descriptors.csv"
        fpocket_df = (
            load_fpocket_descriptors(descriptor_file)
            if descriptor_file.exists()
            else None
        )

        sifts = fetch_sifts_mapping(pdb_id)
        pdb_to_uni = build_pdb_to_uniprot_map(sifts)

        sifts_chains = {c for (c, _) in pdb_to_uni.keys()}

        if len(sifts_chains) == 1:
            canonical_chain = next(iter(sifts_chains))
        else:
            canonical_chain = None

        protein_name = fetch_gpcrdb_protein_from_structure(pdb_id)
        gpcrdb_segments = fetch_gpcrdb_segments(protein_name)

        for pocket_pdb in pocket_files:

            # Make sure fpocket indices map correctly: pocket10_env_atm.pdb → pocket10, idx=10
            stem = pocket_pdb.stem
            pocket_idx = stem.replace("pocket", "").replace("_env_atm", "")
            pocket_name = f"pocket{pocket_idx}"

            residues = extract_pocket_residues(pocket_pdb)

            normalized_residues = set()
            for chain, resseq in residues:
                if not chain and canonical_chain:
                    normalized_residues.add((canonical_chain, resseq))
                else:
                    normalized_residues.add((chain, resseq))
            
            # UniProt residue numbers for persistence clustering
            uni_residues = set()
            for chain_res in normalized_residues:
                if chain_res in pdb_to_uni:
                    _, uni_res = pdb_to_uni[chain_res]
                    uni_residues.add(uni_res)

            # DEBUG: how many pocket residues map to UniProt via SIFTS?
            mapped = sum(1 for r in residues if r in pdb_to_uni)
            logger.debug(
                f"{pdb_id} {pocket_name}: {mapped}/{len(normalized_residues)} residues mapped via SIFTS"
            )

            if mapped == 0:
                logger.debug(
                    f"Example pocket residues: {list(residues)[:5]}"
                )
                logger.debug(
                    f"Example SIFTS keys: {list(pdb_to_uni.keys())[:5]}"
                )

            region, confidence, raw = classify_pocket_topology(
                residues, pdb_to_uni, gpcrdb_segments
            )

            fpocket_metrics = {}
            if fpocket_df is not None and pocket_idx in fpocket_df.index:
                fpocket_metrics = fpocket_df.loc[pocket_idx]

            records.append({
                "pdb_id": pdb_id,
                "pocket": pocket_name,
                "region": region,
                "confidence": confidence,
                "n_residues": len(residues),
                "fpocket_drug_score": fpocket_metrics.get("drug_score"),
                "fpocket_volume": fpocket_metrics.get("volume"),
                "fpocket_hydrophobicity": fpocket_metrics.get("hydrophobicity_score"),
                "fpocket_polarity": fpocket_metrics.get("polarity_score"),
                "fpocket_inter_chain": fpocket_metrics.get("inter_chain"),
                "raw_segment_counts": raw,
                "uni_residues": uni_residues,
            })

    df = pd.DataFrame(records)

    # Assign regions in each structure
    final_rows = []

    for pdb_id, sdf in df.groupby("pdb_id"):
        
        # Identify orthosteric pocket: top 7TM by drug_score then volume
        seven_tm = sdf[sdf["region"] == "7TM"]

        if not seven_tm.empty:
            orthosteric_pocket = (
                seven_tm
                .sort_values(
                    ["fpocket_drug_score", "fpocket_volume"],
                    ascending=False
                )
                .iloc[0]["pocket"]
            )
        else:
            orthosteric_pocket = None

        for _, row in sdf.iterrows():
            row = row.copy()
            row["pocket_role"] = assign_pocket_role(row, orthosteric_pocket)
            final_rows.append(row)

    final_df = pd.DataFrame(final_rows)
    rows = final_df.to_dict(orient="records")

    clustered_rows, clusters = cluster_pockets_across_structures(rows)

    # Compute persistence stats
    total_structures = final_df["pdb_id"].nunique()
    cluster_stats = {}

    for cluster in clusters:
        pdbs = {r["pdb_id"] for r in cluster["members"]}
        cluster_stats[cluster["id"]] = {
            "persistence": round(len(pdbs) / total_structures, 3),
            "n_structures": len(pdbs)
        }

    # Annotate rows
    for row in clustered_rows:
        cid = row["pocket_cluster_id"]
        row["pocket_persistence"] = cluster_stats[cid]["persistence"]
        row["pocket_cluster_size"] = cluster_stats[cid]["n_structures"]
        row.pop("uni_residues", None)  # not for CSV

    final_df = pd.DataFrame(clustered_rows)

    # ---- Build density grids for persistent pocket clusters ----

    # Infer structure directory from fpocket input
    structure_dir = Path(
        config.get("run_fpocket", {}).get("output_directory", "")
    )

    if not structure_dir or not structure_dir.exists():
        raise KeyError(
            "Structure directory not found. Please define "
            "run_fpocket.structure_directory in your config."
        )
    reference_pdb_id = sorted(final_df["pdb_id"].unique())[0]

    grid_dir = output_dir / "pocket_density_grids"
    grid_dir.mkdir(exist_ok=True)

    for cluster in clusters:
        if len(cluster["members"]) < 2:
            continue  # skip non-persistent pockets

        try:
            dx_file = build_pocket_cluster_grid(
                cluster_rows=cluster["members"],
                structure_dir=structure_dir,
                reference_pdb_id=reference_pdb_id,
                pdb_to_uni=pdb_to_uni,
                gpcrdb_segments=gpcrdb_segments,
                out_dir=grid_dir,
                spacing=1.0
            )

            logger.info(f"Built density grid for {cluster['id']} → {dx_file}")

        except Exception as e:
            logger.warning(
                f"Failed grid build for {cluster['id']}: {e}"
            )

    final_df.to_csv(output_file, index=False, quoting=csv.QUOTE_ALL)

    logger.info(f"Wrote {len(final_df)} classified pockets → {output_file}")

    return {
        "classified_pockets_file": str(output_file),
        "df": final_df,
    }