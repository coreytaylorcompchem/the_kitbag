import pandas as pd
from collections import Counter

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
        sep=r"\s+",
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