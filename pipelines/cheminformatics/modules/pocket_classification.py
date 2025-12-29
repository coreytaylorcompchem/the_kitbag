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
    Returns set of (chain_id, pdb_resseq)
    """
    structure = parser.get_structure("pocket", pocket_atm_pdb)
    residues = set()

    for atom in structure.get_atoms():
        res = atom.get_parent()
        chain = res.get_parent().id
        resseq = res.get_id()[1]
        residues.add((chain, resseq))

    return residues

def fetch_sifts_mapping(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def build_pdb_to_uniprot_map(sifts_data):
    mapping = {}

    for pdb_id, entry in sifts_data.items():
        for uniprot_id, udata in entry["UniProt"].items():
            for m in udata["mappings"]:
                chain = m["chain_id"]
                pdb_start = m["start"]["residue_number"]
                pdb_end = m["end"]["residue_number"]
                uni_start = m["unp_start"]
                uni_end = m["unp_end"]

                for pdb_res, uni_res in zip(
                    range(pdb_start, pdb_end + 1),
                    range(uni_start, uni_end + 1),
                ):
                    mapping[(chain, pdb_res)] = (uniprot_id, uni_res)

    return mapping

def fetch_gpcrdb_segments(uniprot_id):
    url = f"https://gpcrdb.org/services/residues/{uniprot_id}/"
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    return {
        r["sequence_number"]: r["protein_segment"]
        for r in r.json()
    }

def classify_pocket_topology(pocket_residues, pdb_to_uni, gpcrdb_segments):
    labels = []

    for key in pocket_residues:
        if key not in pdb_to_uni:
            continue

        _, uni_res = pdb_to_uni[key]
        segment = gpcrdb_segments.get(uni_res, "Unknown")
        labels.append(segment)

    if not labels:
        return "Unknown", 0.0, {}

    counts = Counter(labels)

    collapsed = {
        "7TM": sum(v for k, v in counts.items() if k.startswith("TM")),
        "ECL-ICL": sum(v for k, v in counts.items() if "CL" in k),
        "ECD": counts.get("N-term", 0),
        "H8": counts.get("H8", 0),
    }

    total = sum(collapsed.values())
    if total == 0:
        return "Unknown", 0.0, dict(counts)

    region = max(collapsed, key=collapsed.get)
    confidence = collapsed[region] / total

    return region, round(confidence, 3), dict(counts)

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

        uniprot_id = next(iter(sifts[pdb_id.lower()]["UniProt"]))
        gpcrdb_segments = fetch_gpcrdb_segments(uniprot_id)

        for pocket_pdb in pocket_files:
            # pocket10_env_atm.pdb → pocket10, idx=10
            stem = pocket_pdb.stem
            pocket_idx = stem.replace("pocket", "").replace("_env_atm", "")
            pocket_name = f"pocket{pocket_idx}"

            residues = extract_pocket_residues(pocket_pdb)

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

                # ✅ Correct fpocket columns
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