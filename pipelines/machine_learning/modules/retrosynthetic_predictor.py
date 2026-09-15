# import json
# import math

from pathlib import Path

import pandas as pd
import numpy as np

from tqdm import tqdm

import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

from modules.utils.retrosynth_rxn_rules import REACTION_RULES, COMPILED_REACTION_RULES, _canonicalise_smiles, _apply_reaction_rule

from pipeline.task_registry import register_task

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    "load_retrosynthesis_targets",
    category="SYNTHESIS",
    description="Load candidate molecules for retrosynthetic analysis."
)
def load_retrosynthesis_targets(config, context):

    input_file = Path(config["input_file"])

    if not input_file.exists():
        raise FileNotFoundError(
            f"Retrosynthesis target file not found: {input_file}"
        )

    file_type = input_file.suffix.lower()

    if file_type == ".csv":
        df = pd.read_csv(input_file)

    elif file_type == ".parquet":
        df = pd.read_parquet(input_file)

    else:
        raise ValueError(
            f"Unsupported target file format: {file_type}. "
            "Use CSV or Parquet."
        )

    smiles_col = config.get("smiles_col", "smiles")

    if smiles_col not in df.columns:
        raise ValueError(
            f"SMILES column '{smiles_col}' not found in "
            f"target dataset. Available columns: {list(df.columns)}"
        )

    df = df.copy()

    # Keep the original column but provide a standard name for
    # downstream retrosynthesis tasks.
    if smiles_col != "target_smiles":
        df["target_smiles"] = df[smiles_col]
    else:
        df["target_smiles"] = df[smiles_col]

    id_col = config.get("id_col")

    if id_col is not None:

        if id_col not in df.columns:
            raise ValueError(
                f"ID column '{id_col}' not found in target dataset."
            )

        df["target_id"] = df[id_col]

    else:
        df["target_id"] = range(len(df))

    df = df[
        df["target_smiles"].notna()
        & (df["target_smiles"].astype(str).str.strip() != "")
    ].copy()

    df.reset_index(drop=True, inplace=True)

    logger.info(
        f"Loaded {len(df):,} retrosynthesis targets "
        f"from {input_file}"
    )

    context["retrosynthesis_targets"] = df

    return {
        "retrosynthesis_targets": df,
    }

@register_task(
    "generate_reaction_disconnections",
    category="SYNTHESIS",
    description="Generate retrosynthetic disconnections using configured reaction rules."
)
def generate_reaction_disconnections(config, context):

    targets_df = context.get("retrosynthesis_targets")

    if targets_df is None:
        raise ValueError(
            "No retrosynthesis targets found in pipeline context. "
            "Run 'load_retrosynthesis_targets' before "
            "'generate_reaction_disconnections'."
        )

    if "target_smiles" not in targets_df.columns:
        raise ValueError(
            "Retrosynthesis target dataset must contain "
            "'target_smiles'."
        )

    enabled_rules = config.get(
        "rules",
        list(REACTION_RULES.keys()),
    )

    unknown_rules = [
        rule
        for rule in enabled_rules
        if rule not in COMPILED_REACTION_RULES
    ]

    if unknown_rules:
        raise ValueError(
            f"Unknown reaction rules: {unknown_rules}. "
            f"Available rules: {list(REACTION_RULES)}"
        )

    all_disconnections = []

    for _, target_row in targets_df.iterrows():

        target_smiles = target_row["target_smiles"]
        target_id = target_row["target_id"]

        canonical_target = _canonicalise_smiles(target_smiles)

        if canonical_target is None:
            logger.warning(
                f"Skipping invalid target {target_id}: "
                f"{target_smiles}"
            )
            continue

        for rule_name in enabled_rules:

            reaction = COMPILED_REACTION_RULES[rule_name]

            candidates = _apply_reaction_rule(
                canonical_target,
                rule_name,
                reaction,
            )

            for candidate in candidates:
                candidate["target_id"] = target_id

            all_disconnections.extend(candidates)

    # Assign stable IDs.
    for index, candidate in enumerate(
        all_disconnections,
        start=1,
    ):
        candidate["disconnection_id"] = index

    logger.info(
        f"Generated {len(all_disconnections):,} "
        f"retrosynthetic disconnections from "
        f"{len(targets_df):,} target(s)"
    )

    for rule_name in enabled_rules:

        count = sum(
            candidate["reaction_rule"] == rule_name
            for candidate in all_disconnections
        )

        logger.info(
            f"  {rule_name}: {count:,}"
        )

    context["reaction_disconnections"] = all_disconnections

    return {
        "reaction_disconnections": all_disconnections,
    }

@register_task(
    "search_commercial_building_blocks",
    category="SYNTHESIS",
    description="Search commercial building blocks using reaction-handle constrained molecular similarity."
)
def search_commercial_building_blocks(config, context):

    disconnections = context.get("reaction_disconnections")

    if disconnections is None:
        raise ValueError(
            "No reaction disconnections found in pipeline context. "
            "Run 'generate_reaction_disconnections' before "
            "'search_commercial_building_blocks'."
        )

    index_dir = Path(
        config.get(
            "index_dir",
            "outputs/commercial_bbs/search_index"
        )
    )

    if not index_dir.exists():
        raise FileNotFoundError(
            f"Building-block search index not found: {index_dir}"
        )

    max_candidates = int(
        config.get("max_candidates_per_fragment", 20)
    )

    similarity_threshold = float(
        config.get("similarity_threshold", 0.40)
    )

    fingerprint_radius = int(
        config.get("fingerprint_radius", 2)
    )

    fingerprint_bits = int(
        config.get("fingerprint_bits", 2048)
    )

    # ------------------------------------------------------------------
    # Reaction-rule → possible commercial reaction handles
    # ------------------------------------------------------------------

    rule_handles = {
        "amide_formation": {
            "amine": "amine",
            "carboxylic_acid": "carboxylic_acid",
            "acid_chloride": "acid_chloride",
        },

        "ester_formation": {
            "alcohol": "alcohol",
            "phenol": "phenol",
            "carboxylic_acid": "carboxylic_acid",
            "acid_chloride": "acid_chloride",
        },

        "suzuki_coupling": {
            "aryl_vinyl_halide": "aryl_vinyl_halide",
            "boronic_acid": "boronic_acid",
            "boronate_ester": "boronate_ester",
        },

        "reductive_amination": {
            "amine": "amine",
            "aldehyde": "aldehyde",
            "ketone": "ketone",
        },

        "ether_formation": {
            "alcohol": "alcohol",
            "phenol": "phenol",
            "alkyl_halide": "alkyl_halide",
        },
    }

    # ------------------------------------------------------------------
    # Cache commercial indexes and fingerprints
    # ------------------------------------------------------------------

    index_cache = {}
    fingerprint_cache = {}

    def load_handle_index(handle):

        if handle in index_cache:
            return index_cache[handle]

        index_file = index_dir / f"{handle}.parquet"

        if not index_file.exists():
            logger.warning(
                f"No commercial building-block index exists for "
                f"reaction handle '{handle}'"
            )

            index_cache[handle] = None
            fingerprint_cache[handle] = None

            return None

        handle_df = pd.read_parquet(index_file)

        if "canonical_smiles" not in handle_df.columns:
            raise ValueError(
                f"Commercial index '{index_file}' does not contain "
                "'canonical_smiles'."
            )

        index_cache[handle] = handle_df

        logger.info(
            f"Loaded commercial index '{handle}': "
            f"{len(handle_df):,} building blocks"
        )

        return handle_df

    def build_fingerprints(handle, handle_df):

        if handle in fingerprint_cache:
            return fingerprint_cache[handle]

        logger.info(
            f"Generating fingerprints for '{handle}' index"
        )

        fingerprints = []

        for smiles in tqdm(
            handle_df["canonical_smiles"],
            desc=f"Fingerprints: {handle}",
            leave=False,
        ):

            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                fingerprints.append(None)
                continue

            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol,
                radius=fingerprint_radius,
                nBits=fingerprint_bits,
            )

            fingerprints.append(fp)

        fingerprint_cache[handle] = fingerprints

        return fingerprints

    # ------------------------------------------------------------------
    # Similarity search
    # ------------------------------------------------------------------

    def search_index(
        precursor_smiles,
        handle,
    ):

        canonical_precursor = _canonicalise_smiles(
            precursor_smiles
        )

        if canonical_precursor is None:
            return []

        precursor_mol = Chem.MolFromSmiles(
            canonical_precursor
        )

        if precursor_mol is None:
            return []

        handle_df = load_handle_index(handle)

        if handle_df is None or handle_df.empty:
            return []

        # --------------------------------------------------------------
        # Exact match first
        # --------------------------------------------------------------

        exact_matches = handle_df[
            handle_df["canonical_smiles"] == canonical_precursor
        ].copy()

        exact_records = []

        for _, row in exact_matches.iterrows():

            record = row.to_dict()

            record["similarity"] = 1.0
            record["exact_match"] = True
            record["expected_reaction_handle"] = handle
            record["precursor_smiles"] = canonical_precursor

            exact_records.append(record)

        # --------------------------------------------------------------
        # Similarity search
        # --------------------------------------------------------------

        query_fp = AllChem.GetMorganFingerprintAsBitVect(
            precursor_mol,
            radius=fingerprint_radius,
            nBits=fingerprint_bits,
        )

        fingerprints = build_fingerprints(
            handle,
            handle_df,
        )

        similarities = []

        for index, fp in enumerate(fingerprints):

            if fp is None:
                continue

            similarity = DataStructs.TanimotoSimilarity(
                query_fp,
                fp,
            )

            if similarity >= similarity_threshold:

                similarities.append(
                    (
                        similarity,
                        index,
                    )
                )

        similarities.sort(
            key=lambda x: x[0],
            reverse=True,
        )

        # --------------------------------------------------------------
        # Convert similarity results into records
        # --------------------------------------------------------------

        similarity_records = []

        exact_inchikeys = {
            record["inchikey"]
            for record in exact_records
            if "inchikey" in record
        }

        for similarity, index in similarities:

            row = handle_df.iloc[index]

            # Don't return exact matches twice.
            if row["inchikey"] in exact_inchikeys:
                continue

            record = row.to_dict()

            record["similarity"] = float(similarity)
            record["exact_match"] = False
            record["expected_reaction_handle"] = handle
            record["precursor_smiles"] = canonical_precursor

            similarity_records.append(record)

            if len(similarity_records) >= max_candidates:
                break

        # Exact matches always appear first.
        results = (
            exact_records
            + similarity_records
        )

        return results[:max_candidates]

    # ------------------------------------------------------------------
    # Search all disconnections
    # ------------------------------------------------------------------

    searched_disconnections = []

    for disconnection in tqdm(
        disconnections,
        desc="Searching commercial building blocks",
    ):

        rule_name = disconnection["reaction_rule"]
        precursor_smiles = disconnection["precursor_smiles"]

        possible_handles = rule_handles.get(rule_name)

        if possible_handles is None:

            logger.warning(
                f"No commercial-search mapping defined for "
                f"reaction rule '{rule_name}'"
            )

            searched_disconnections.append({
                **disconnection,
                "commercial_candidates": [],
                "commercial_search_complete": False,
            })

            continue

        fragment_results = []

        for precursor in precursor_smiles:

            precursor_candidates = []

            # Search every chemically compatible handle index.
            for role, handle in possible_handles.items():

                matches = search_index(
                    precursor,
                    handle,
                )

                for match in matches:

                    match["precursor_role"] = role

                    precursor_candidates.append(
                        match
                    )

            # ----------------------------------------------------------
            # Deduplicate commercial molecules.
            #
            # A molecule can appear in multiple handle indexes, so
            # retain its best similarity result.
            # ----------------------------------------------------------

            unique_candidates = {}

            for candidate in precursor_candidates:

                key = candidate["inchikey"]

                if key not in unique_candidates:

                    unique_candidates[key] = candidate

                else:

                    existing = unique_candidates[key]

                    if candidate["similarity"] > existing["similarity"]:
                        unique_candidates[key] = candidate

            precursor_candidates = list(
                unique_candidates.values()
            )

            precursor_candidates.sort(
                key=lambda x: x["similarity"],
                reverse=True,
            )

            precursor_candidates = (
                precursor_candidates[:max_candidates]
            )

            fragment_results.append({
                "precursor_smiles": precursor,
                "commercial_candidates": precursor_candidates,
                "n_commercial_candidates": len(
                    precursor_candidates
                ),
                "commercial_match": bool(
                    precursor_candidates
                ),
            })

        searched_disconnections.append({
            **disconnection,
            "commercial_candidates": fragment_results,
            "commercial_search_complete": True,
        })

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------

    n_disconnections = len(
        searched_disconnections
    )

    n_with_commercial_matches = 0

    for disconnection in searched_disconnections:

        if any(
            fragment["commercial_match"]
            for fragment in disconnection["commercial_candidates"]
        ):
            n_with_commercial_matches += 1

    logger.info(
        "Commercial building-block search complete"
    )

    logger.info(
        f"  Disconnections searched: "
        f"{n_disconnections:,}"
    )

    logger.info(
        f"  Disconnections with at least one "
        f"commercial precursor match: "
        f"{n_with_commercial_matches:,}"
    )

    logger.info(
        f"  Similarity threshold: "
        f"{similarity_threshold:.2f}"
    )

    logger.info(
        f"  Maximum candidates per fragment: "
        f"{max_candidates}"
    )

    context["commercial_building_block_matches"] = (
        searched_disconnections
    )

    return {
        "commercial_building_block_matches": (
            searched_disconnections
        ),
    }

@register_task(
"assemble_candidate_routes",
category="SYNTHESIS",
description="Assemble candidate retrosynthetic routes from reaction disconnections and commercial building-block matches."
)
def assemble_candidate_routes(config, context):

    commercial_matches = context.get(
        "commercial_building_block_matches"
    )

    if commercial_matches is None:
        raise ValueError(
            "No commercial building-block matches found in "
            "pipeline context. Run "
            "'search_commercial_building_blocks' before "
            "'assemble_candidate_routes'."
        )

    max_routes_per_disconnection = int(
        config.get("max_routes_per_disconnection", 20)
    )

    candidate_routes = []

    route_id = 1

    for disconnection in tqdm(
        commercial_matches,
        desc="Assembling candidate routes",
    ):

        target_id = disconnection["target_id"]
        target_smiles = disconnection["target_smiles"]
        reaction_rule = disconnection["reaction_rule"]
        disconnection_id = disconnection["disconnection_id"]

        fragment_results = disconnection.get(
            "commercial_candidates",
            []
        )

        if not fragment_results:
            continue

        # Each precursor fragment may have zero or more
        # commercial candidates.
        #
        # We deliberately retain unmatched fragments rather
        # than discarding the entire disconnection. This allows
        # later route scoring to distinguish:
        #
        #   - fully commercial routes
        #   - partially commercial routes
        #   - completely non-commercial routes
        candidate_lists = []

        for fragment in fragment_results:

            candidates = fragment.get(
                "commercial_candidates",
                []
            )

            if candidates:
                candidate_lists.append(
                    candidates
                )
            else:
                # No commercial candidate for this precursor.
                #
                # None is used as a placeholder so that the
                # disconnection can still become a candidate
                # route.
                candidate_lists.append(
                    [None]
                )

        from itertools import product

        route_combinations = product(
            *candidate_lists
        )

        n_routes = 0

        for combination in route_combinations:

            if n_routes >= max_routes_per_disconnection:
                break

            precursors = []

            n_commercial_precursors = 0
            precursor_similarities = []

            for fragment, candidate in zip(
                fragment_results,
                combination,
            ):

                precursor = {
                    "precursor_smiles": (
                        fragment["precursor_smiles"]
                    ),
                    "precursor_role": None,
                    "required_reaction_handle": None,
                    "commercial_smiles": None,
                    "inchikey": None,
                    "mcule_id": None,
                    "similarity": None,
                    "exact_match": False,
                    "commercially_available": False,
                }

                if candidate is not None:

                    precursor["precursor_role"] = (
                        candidate.get(
                            "precursor_role"
                        )
                    )

                    precursor[
                        "required_reaction_handle"
                    ] = candidate.get(
                        "expected_reaction_handle"
                    )

                    precursor["commercial_smiles"] = (
                        candidate.get(
                            "canonical_smiles"
                        )
                    )

                    precursor["inchikey"] = (
                        candidate.get(
                            "inchikey"
                        )
                    )

                    precursor["mcule_id"] = (
                        candidate.get(
                            "mcule_id"
                        )
                    )

                    precursor["similarity"] = float(
                        candidate.get(
                            "similarity",
                            0.0
                        )
                    )

                    precursor["exact_match"] = bool(
                        candidate.get(
                            "exact_match",
                            False
                        )
                    )

                    precursor[
                        "commercially_available"
                    ] = True

                    n_commercial_precursors += 1

                    precursor_similarities.append(
                        precursor["similarity"]
                    )

                precursors.append(precursor)

            n_precursors = len(precursors)

            all_precursors_commercial = (
                n_commercial_precursors
                == n_precursors
            )

            any_precursors_commercial = (
                n_commercial_precursors > 0
            )

            all_precursors_exact = (
                all_precursors_commercial
                and all(
                    precursor["exact_match"]
                    for precursor in precursors
                )
            )

            if precursor_similarities:
                min_precursor_similarity = float(
                    min(precursor_similarities)
                )

                mean_precursor_similarity = float(
                    np.mean(
                        precursor_similarities
                    )
                )
            else:
                min_precursor_similarity = None
                mean_precursor_similarity = None

            route = {
                "route_id": route_id,

                "target_id": target_id,
                "target_smiles": target_smiles,

                "disconnection_id": (
                    disconnection_id
                ),

                "reaction_rule": reaction_rule,

                "n_steps": 1,

                "precursors": precursors,

                "n_precursors": n_precursors,

                "n_commercial_precursors": (
                    n_commercial_precursors
                ),

                "all_precursors_commercial": (
                    all_precursors_commercial
                ),

                "any_precursors_commercial": (
                    any_precursors_commercial
                ),

                "all_precursors_exact": (
                    all_precursors_exact
                ),

                "min_precursor_similarity": (
                    min_precursor_similarity
                ),

                "mean_precursor_similarity": (
                    mean_precursor_similarity
                ),
            }

            candidate_routes.append(route)

            route_id += 1
            n_routes += 1

    logger.info(
        "Candidate route assembly complete"
    )

    logger.info(
        f"  Disconnections considered: "
        f"{len(commercial_matches):,}"
    )

    logger.info(
        f"  Candidate routes generated: "
        f"{len(candidate_routes):,}"
    )

    n_fully_commercial = sum(
        route["all_precursors_commercial"]
        for route in candidate_routes
    )

    n_partially_commercial = sum(
        route["any_precursors_commercial"]
        and not route["all_precursors_commercial"]
        for route in candidate_routes
    )

    n_non_commercial = sum(
        not route["any_precursors_commercial"]
        for route in candidate_routes
    )

    logger.info(
        f"  Fully commercial routes: "
        f"{n_fully_commercial:,}"
    )

    logger.info(
        f"  Partially commercial routes: "
        f"{n_partially_commercial:,}"
    )

    logger.info(
        f"  Non-commercial routes: "
        f"{n_non_commercial:,}"
    )

    context["candidate_routes"] = candidate_routes

    return {
        "candidate_routes": candidate_routes,
    }

@register_task(
"save_candidate_routes",
category="SYNTHESIS",
description="Save assembled retrosynthetic routes and precursor details to Parquet."
)
def save_candidate_routes(config, context):

    candidate_routes = context.get("candidate_routes")

    if candidate_routes is None:
        raise ValueError(
            "No candidate routes found in pipeline context. "
            "Run 'assemble_candidate_routes' before "
            "'save_candidate_routes'."
        )

    output_dir = Path(
        config.get(
            "output_dir",
            "outputs/retrosynthesis/routes"
        )
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    route_output_file = output_dir / config.get(
        "route_output_file",
        "candidate_routes.parquet"
    )

    precursor_output_file = output_dir / config.get(
        "precursor_output_file",
        "candidate_route_precursors.parquet"
    )

    # ------------------------------------------------------------------
    # Route-level table
    # ------------------------------------------------------------------

    route_records = []

    for route in candidate_routes:

        route_records.append({
            "route_id": route["route_id"],
            "target_id": route["target_id"],
            "target_smiles": route["target_smiles"],
            "disconnection_id": route["disconnection_id"],
            "reaction_rule": route["reaction_rule"],
            "n_steps": route["n_steps"],
            "n_precursors": route["n_precursors"],
            "n_commercial_precursors": (
                route["n_commercial_precursors"]
            ),
            "all_precursors_commercial": (
                route["all_precursors_commercial"]
            ),
            "any_precursors_commercial": (
                route["any_precursors_commercial"]
            ),
            "all_precursors_exact": (
                route["all_precursors_exact"]
            ),
            "min_precursor_similarity": (
                route["min_precursor_similarity"]
            ),
            "mean_precursor_similarity": (
                route["mean_precursor_similarity"]
            ),
        })

    routes_df = pd.DataFrame(
        route_records
    )

    routes_df.to_parquet(
        route_output_file,
        index=False,
    )

    # ------------------------------------------------------------------
    # Precursor-level table
    # ------------------------------------------------------------------

    precursor_records = []

    for route in candidate_routes:

        for precursor_index, precursor in enumerate(
            route["precursors"],
            start=1,
        ):

            precursor_records.append({
                "route_id": route["route_id"],
                "target_id": route["target_id"],
                "target_smiles": route["target_smiles"],
                "disconnection_id": (
                    route["disconnection_id"]
                ),
                "reaction_rule": (
                    route["reaction_rule"]
                ),
                "precursor_index": precursor_index,
                "precursor_smiles": (
                    precursor["precursor_smiles"]
                ),
                "precursor_role": (
                    precursor["precursor_role"]
                ),
                "required_reaction_handle": (
                    precursor[
                        "required_reaction_handle"
                    ]
                ),
                "commercial_smiles": (
                    precursor["commercial_smiles"]
                ),
                "inchikey": precursor["inchikey"],
                "mcule_id": precursor["mcule_id"],
                "similarity": precursor["similarity"],
                "exact_match": precursor["exact_match"],
                "commercially_available": (
                    precursor[
                        "commercially_available"
                    ]
                ),
            })

    precursors_df = pd.DataFrame(
        precursor_records
    )

    precursors_df.to_parquet(
        precursor_output_file,
        index=False,
    )

    logger.info(
        "Candidate route files saved"
    )

    logger.info(
        f"  Routes: {route_output_file}"
    )

    logger.info(
        f"  Route records: {len(routes_df):,}"
    )

    logger.info(
        f"  Precursors: {precursor_output_file}"
    )

    logger.info(
        f"  Precursor records: {len(precursors_df):,}"
    )

    context["candidate_routes_df"] = routes_df
    context["candidate_route_precursors_df"] = precursors_df

    return {
        "candidate_routes_df": routes_df,
        "candidate_route_precursors_df": precursors_df,
    }

@register_task(
"validate_candidate_routes",
category="SYNTHESIS",
description="Validate and visualise assembled retrosynthetic candidate routes."
)
def validate_candidate_routes(config, context):

    route_file = Path(
        config.get(
            "route_file",
            "outputs/retrosynthesis/routes/candidate_routes.parquet"
        )
    )

    precursor_file = Path(
        config.get(
            "precursor_file",
            "outputs/retrosynthesis/routes/candidate_route_precursors.parquet"
        )
    )

    output_dir = Path(
        config.get(
            "output_dir",
            "outputs/retrosynthesis/routes/validation"
        )
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not route_file.exists():
        raise FileNotFoundError(
            f"Candidate route file not found: {route_file}"
        )

    if not precursor_file.exists():
        raise FileNotFoundError(
            f"Candidate precursor file not found: {precursor_file}"
        )

    routes_df = pd.read_parquet(
        route_file
    )

    precursors_df = pd.read_parquet(
        precursor_file
    )

    if routes_df.empty:
        logger.warning(
            "Candidate route file is empty. "
            "Nothing to validate."
        )
        return {}

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------

    summary_records = [
        {
            "metric": "n_routes",
            "value": len(routes_df),
        },
        {
            "metric": "n_targets",
            "value": routes_df["target_id"].nunique(),
        },
        {
            "metric": "n_disconnections",
            "value": routes_df["disconnection_id"].nunique(),
        },
        {
            "metric": "n_fully_commercial_routes",
            "value": int(
                routes_df[
                    "all_precursors_commercial"
                ].sum()
            ),
        },
        {
            "metric": "n_partially_commercial_routes",
            "value": int(
                (
                    routes_df[
                        "any_precursors_commercial"
                    ]
                    & ~routes_df[
                        "all_precursors_commercial"
                    ]
                ).sum()
            ),
        },
        {
            "metric": "n_non_commercial_routes",
            "value": int(
                (
                    ~routes_df[
                        "any_precursors_commercial"
                    ]
                ).sum()
            ),
        },
        {
            "metric": "n_precursor_records",
            "value": len(precursors_df),
        },
        {
            "metric": "n_commercial_precursor_records",
            "value": int(
                precursors_df[
                    "commercially_available"
                ].sum()
            ),
        },
        {
            "metric": "mean_precursor_similarity",
            "value": float(
                precursors_df[
                    "similarity"
                ].dropna().mean()
            ),
        },
        {
            "metric": "median_precursor_similarity",
            "value": float(
                precursors_df[
                    "similarity"
                ].dropna().median()
            ),
        },
    ]

    summary_df = pd.DataFrame(
        summary_records
    )

    summary_df.to_csv(
        output_dir / "route_summary.csv",
        index=False,
    )

    with open(
        output_dir / "route_summary.txt",
        "w",
    ) as handle:

        handle.write(
            "RETROSYNTHESIS ROUTE VALIDATION\n"
        )
        handle.write(
            "================================\n\n"
        )

        for _, row in summary_df.iterrows():

            handle.write(
                f"{row['metric']}: "
                f"{row['value']}\n"
            )

    # ------------------------------------------------------------------
    # Plot 1: commercial coverage
    # ------------------------------------------------------------------

    commercial_counts = pd.Series({
        "Fully commercial": int(
            routes_df[
                "all_precursors_commercial"
            ].sum()
        ),
        "Partially commercial": int(
            (
                routes_df[
                    "any_precursors_commercial"
                ]
                & ~routes_df[
                    "all_precursors_commercial"
                ]
            ).sum()
        ),
        "Non-commercial": int(
            (
                ~routes_df[
                    "any_precursors_commercial"
                ]
            ).sum()
        ),
    })

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    commercial_counts.plot.bar(
        ax=ax
    )

    ax.set_ylabel(
        "Number of routes"
    )

    ax.set_title(
        "Commercial availability of candidate routes"
    )

    ax.tick_params(
        axis="x",
        rotation=0,
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / "commercial_coverage.svg"
    )

    plt.close(fig)

    # ------------------------------------------------------------------
    # Plot 2: routes by reaction rule
    # ------------------------------------------------------------------

    rule_counts = (
        routes_df[
            "reaction_rule"
        ]
        .value_counts()
        .sort_values(
            ascending=False
        )
    )

    fig, ax = plt.subplots(
        figsize=(9, 5)
    )

    rule_counts.plot.bar(
        ax=ax
    )

    ax.set_ylabel(
        "Number of routes"
    )

    ax.set_title(
        "Candidate routes by reaction rule"
    )

    ax.tick_params(
        axis="x",
        rotation=45,
    )

    fig.tight_layout()

    fig.savefig(
        output_dir / "routes_by_reaction.svg"
    )

    plt.close(fig)

    # ------------------------------------------------------------------
    # Plot 3: precursor similarity distribution
    # ------------------------------------------------------------------

    similarity_values = (
        precursors_df[
            "similarity"
        ]
        .dropna()
    )

    if not similarity_values.empty:

        fig, ax = plt.subplots(
            figsize=(8, 5)
        )

        ax.hist(
            similarity_values,
            bins=20,
        )

        ax.set_xlabel(
            "Morgan/Tanimoto similarity"
        )

        ax.set_ylabel(
            "Number of precursor matches"
        )

        ax.set_title(
            "Commercial precursor similarity"
        )

        fig.tight_layout()

        fig.savefig(
            output_dir / "precursor_similarity.svg"
        )

        plt.close(fig)

    # ------------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------------

    logger.info(
        "Candidate route validation complete"
    )

    logger.info(
        f"  Routes: {len(routes_df):,}"
    )

    logger.info(
        f"  Targets: "
        f"{routes_df['target_id'].nunique():,}"
    )

    logger.info(
        f"  Fully commercial: "
        f"{int(routes_df['all_precursors_commercial'].sum()):,}"
    )

    logger.info(
        f"  Partially commercial: "
        f"{int((routes_df['any_precursors_commercial'] & ~routes_df['all_precursors_commercial']).sum()):,}"
    )

    logger.info(
        f"  Validation output: {output_dir}"
    )

    return {
        "candidate_route_validation": summary_df,
    }