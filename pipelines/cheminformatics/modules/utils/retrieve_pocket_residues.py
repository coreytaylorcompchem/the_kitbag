import requests
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)

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