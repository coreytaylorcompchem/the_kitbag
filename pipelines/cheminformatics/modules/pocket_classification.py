import requests
import csv
from pathlib import Path
from collections import Counter
import pandas as pd
from Bio.PDB import PDBParser

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)
parser = PDBParser(QUIET=True)

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
    # fpocket -d output is whitespace-delimited, not CSV
    df = pd.read_csv(
        descriptor_file,
        delim_whitespace=True,
        comment="#"
    )

    # Known pocket ID column names used by fpocket
    pocket_cols = [
        "cav_id",          # most common with -d
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

        # 🔹 Load fpocket descriptors (whitespace-delimited)
        descriptor_file = struct_dir / f"{pdb_id}_pocket_descriptors.csv"
        fpocket_df = (
            load_fpocket_descriptors(descriptor_file)
            if descriptor_file.exists()
            else None
        )

        # 🔹 Topology mapping
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
            # pocket10_env_atm.pdb → pocket10, idx=10
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

            # 🔍 DEBUG: how many pocket residues map to UniProt via SIFTS?
            mapped = sum(1 for r in residues if r in pdb_to_uni)
            logger.info(
                f"{pdb_id} {pocket_name}: {mapped}/{len(normalized_residues)} residues mapped via SIFTS"
            )

            if mapped == 0:
                logger.info(
                    f"Example pocket residues: {list(residues)[:5]}"
                )
                logger.info(
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

                # transparency
                "raw_segment_counts": raw,
            })

    df = pd.DataFrame(records)
    df.to_csv(output_file, index=False, quoting=csv.QUOTE_ALL)

    logger.info(f"Wrote {len(df)} classified pockets → {output_file}")

    return {
        "classified_pockets_file": str(output_file),
        "df": df,
    }