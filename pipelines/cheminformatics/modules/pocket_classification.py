import csv

from pathlib import Path
import pandas as pd
from tqdm import tqdm

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

from modules.utils.pocket_to_grids import build_pocket_cluster_grid
from modules.utils.retrieve_pocket_residues import ( 
                                                    fetch_sifts_mapping, 
                                                    fetch_gpcrdb_protein_from_structure, 
                                                    fetch_gpcrdb_segments,
                                                    build_pdb_to_uniprot_map,
                                                    extract_pocket_residues
                                                    )
from modules.utils.extract_metrics_and_classify_pockets import (
                                                                load_fpocket_descriptors, 
                                                                cluster_pockets_across_structures, 
                                                                classify_pocket_topology, 
                                                                assign_pocket_role
                                                                )

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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

    structure_dirs = sorted(fpocket_dir.glob("*_out"))

    for struct_dir in tqdm(
        structure_dirs,
        desc="Classifying structures",
        unit="structure"
    ):
        pdb_id = struct_dir.name.replace("_out", "")
        logger.debug(f"Classifying pockets for {pdb_id}")

        pocket_files = sorted(struct_dir.glob("pocket*_env_atm.pdb"))
        if not pocket_files:
            logger.warning(f"No env pockets found for {pdb_id}")
            continue

        logger.debug(f"Found {len(pocket_files)} pockets for {pdb_id}")

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

    persistent_clusters = [
        c for c in clusters if len(c["members"]) >= 2
    ]

    for cluster in tqdm(
        persistent_clusters,
        desc="Building pocket density grids",
        unit="cluster"
    ):
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