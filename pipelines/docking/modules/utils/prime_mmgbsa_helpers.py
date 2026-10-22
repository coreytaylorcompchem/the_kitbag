import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path

from schrodinger import structure

from modules.utils.glideem_helpers import glide_score, compound_id

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

MMGBSA_POSE_KEY_PROPERTY = (
    "s_user_kitbag_mmgbsa_pose_key"
)

MMGBSA_COMPOUND_ID_PROPERTY = (
    "s_user_kitbag_compound_id"
)

MMGBSA_POSE_RANK_PROPERTY = (
    "i_user_kitbag_mcss_pose_rank"
)

def _first_numeric_property(
    st,
    candidates: tuple[str, ...],
) -> float | None:
    """
    Return the first usable numeric structure property.

    Prime property names can vary slightly between Suite
    versions, so use an explicit candidate list rather than
    relying on one property spelling.
    """

    for name in candidates:
        value = st.property.get(
            name
        )

        if value is None:
            continue

        try:
            value = float(
                value
            )

        except (
            TypeError,
            ValueError,
        ):
            continue

        if np.isfinite(
            value
        ):
            return value

    return None


def _find_mmgbsa_dg_bind(
    st,
) -> float | None:
    """
    Locate Prime's MM-GBSA dG Bind property.

    First try known property names, then make a controlled
    fallback search over numeric MM-GBSA binding properties.
    """

    known_names = (
        "r_psp_MMGBSA_dG_Bind",
        "r_psp_MMGBSA_dG_bind",
        "r_psp_MMGBSA_dG_Bind_Coulomb",
    )

    # Only the total dG-bind value is acceptable here.
    value = _first_numeric_property(
        st,
        known_names[:2],
    )

    if value is not None:
        return value

    matches = []

    for name, raw_value in (
        st.property.items()
    ):
        normalized = (
            name.lower()
            .replace("_", "")
            .replace(" ", "")
        )

        if "mmgbsa" not in normalized:
            continue

        if "dgbind" not in normalized:
            continue

        # Exclude component and no-strain values.
        excluded_tokens = (
            "coulomb",
            "covalent",
            "vdw",
            "lipo",
            "solv",
            "hbond",
            "packing",
            "selfcont",
            "strain",
        )

        if any(
            token in normalized
            for token in excluded_tokens
        ):
            continue

        try:
            numeric_value = float(
                raw_value
            )

        except (
            TypeError,
            ValueError,
        ):
            continue

        if np.isfinite(
            numeric_value
        ):
            matches.append(
                (
                    name,
                    numeric_value,
                )
            )

    if len(matches) == 1:
        return matches[0][1]

    if len(matches) > 1:
        raise RuntimeError(
            "Multiple possible MM-GBSA dG Bind "
            "properties were found: "
            + ", ".join(
                name
                for name, _ in matches
            )
        )

    return None


def _write_mmgbsa_input(
    poseviewer: Path,
    output_file: Path,
) -> pd.DataFrame:
    """
    Write a dedicated Prime MM-GBSA input PoseViewer.

    Each ligand pose receives a unique title and persistent
    custom properties. This prevents MM-GBSA results from
    being joined back to docking data by row order.
    """

    receptor = None
    ligand_records = []
    pose_count_by_compound = defaultdict(
        int
    )

    for st in structure.StructureReader(
        str(poseviewer)
    ):
        if st.property.get(
            "b_glide_receptor",
            False,
        ):
            if receptor is not None:
                raise RuntimeError(
                    "Multiple receptor CTs found in "
                    f"{poseviewer}"
                )

            receptor = st
            continue

        score = glide_score(
            st
        )

        if score is None:
            continue

        cid = compound_id(
            st
        )

        score = glide_score(
            st
        )

        if cid:
            cid = str(
                cid
            )

        # Must match _write_mmgbsa_input():
        # only scored ligand CTs participate
        # in pose ranking.
        if (
            cid
            and score is not None
        ):
            pose_count_by_compound[
                cid
            ] = (
                pose_count_by_compound.get(
                    cid,
                    0,
                )
                + 1
            )

            mcss_pose_rank = (
                pose_count_by_compound[
                    cid
                ]
            )

        else:
            mcss_pose_rank = None

        pose_count_by_compound[
            cid
        ] += 1

        pose_rank = (
            pose_count_by_compound[
                cid
            ]
        )

        pose_key = (
            f"{cid}__pose{pose_rank:04d}"
        )

        st = st.copy()

        st.title = pose_key

        st.property[
            MMGBSA_POSE_KEY_PROPERTY
        ] = pose_key

        st.property[
            MMGBSA_COMPOUND_ID_PROPERTY
        ] = cid

        st.property[
            MMGBSA_POSE_RANK_PROPERTY
        ] = pose_rank

        ligand_records.append(
            {
                "mmgbsa_pose_key": (
                    pose_key
                ),
                "compound_id": cid,
                "mcss_pose_rank": (
                    pose_rank
                ),
                "docking_score": float(
                    score
                ),
                "structure": st,
            }
        )

    if receptor is None:
        raise RuntimeError(
            "No receptor CT found in MCSS "
            f"PoseViewer: {poseviewer}"
        )

    if not ligand_records:
        raise RuntimeError(
            "No scored ligand poses found in MCSS "
            f"PoseViewer: {poseviewer}"
        )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with structure.StructureWriter(
        str(output_file)
    ) as writer:
        writer.append(
            receptor
        )

        for record in ligand_records:
            writer.append(
                record["structure"]
            )

    return pd.DataFrame(
        [
            {
                "mmgbsa_pose_key": (
                    record[
                        "mmgbsa_pose_key"
                    ]
                ),
                "compound_id": (
                    record[
                        "compound_id"
                    ]
                ),
                "mcss_pose_rank": (
                    record[
                        "mcss_pose_rank"
                    ]
                ),
                "docking_score": (
                    record[
                        "docking_score"
                    ]
                ),
            }
            for record in ligand_records
        ]
    )


def _parse_mmgbsa_output(
    output_file: Path,
    input_manifest: pd.DataFrame,
) -> pd.DataFrame:
    """
    Parse ligand results from Prime MM-GBSA output.

    Returns one row per input MCSS pose.
    """

    rows = []

    for st in structure.StructureReader(
        str(output_file)
    ):
        pose_key = st.property.get(
            MMGBSA_POSE_KEY_PROPERTY
        )

        if not pose_key:
            title = (
                st.title or ""
            ).strip()

            if title:
                pose_key = title

        if not pose_key:
            continue

        dg_bind = _find_mmgbsa_dg_bind(
            st
        )

        if dg_bind is None:
            logger.warning(
                "No MM-GBSA dG Bind property found "
                "for output structure '%s'",
                pose_key,
            )
            continue

        rows.append(
            {
                "mmgbsa_pose_key": str(
                    pose_key
                ),
                "mmgbsa_dg_bind": float(
                    dg_bind
                ),
            }
        )

    results = pd.DataFrame(
        rows
    )

    if results.empty:
        raise RuntimeError(
            "Prime MM-GBSA output contained no "
            "parseable dG Bind results. Inspect the "
            "output properties and Prime log."
        )

    duplicates = results[
        "mmgbsa_pose_key"
    ].duplicated(
        keep=False
    )

    if duplicates.any():
        duplicate_keys = sorted(
            results.loc[
                duplicates,
                "mmgbsa_pose_key",
            ].unique()
        )

        raise RuntimeError(
            "Duplicate MM-GBSA pose keys found: "
            + ", ".join(
                duplicate_keys
            )
        )

    merged = input_manifest.merge(
        results,
        on="mmgbsa_pose_key",
        how="left",
        validate="one_to_one",
    )

    missing = merged[
        "mmgbsa_dg_bind"
    ].isna()

    if missing.any():
        logger.warning(
            "MM-GBSA returned no usable dG Bind "
            "for %d of %d MCSS poses",
            int(
                missing.sum()
            ),
            len(
                merged
            ),
        )

    return merged

MMGBSA_POSE_KEY_PROPERTY = (
    "s_user_kitbag_mmgbsa_pose_key"
)

MMGBSA_COMPOUND_ID_PROPERTY = (
    "s_user_kitbag_compound_id"
)

MMGBSA_POSE_RANK_PROPERTY = (
    "i_user_kitbag_mcss_pose_rank"
)


def _write_mmgbsa_input(
    poseviewer: Path,
    output_file: Path,
) -> pd.DataFrame:
    """
    Create a Prime MM-GBSA input PoseViewer from an MCSS
    Glide PoseViewer.

    Each docked ligand CT is assigned:

        compound_id
        mcss_pose_rank
        mmgbsa_pose_key

    The pose key survives into the Prime output, allowing
    MM-GBSA results to be joined back to the exact MCSS pose
    rather than relying on output row order.
    """

    poseviewer = (
        Path(poseviewer)
        .expanduser()
        .resolve()
    )

    output_file = (
        Path(output_file)
        .expanduser()
        .resolve()
    )

    if not poseviewer.is_file():
        raise FileNotFoundError(
            f"MCSS PoseViewer not found: "
            f"{poseviewer}"
        )

    receptor = None
    ligand_records = []

    pose_count_by_compound = defaultdict(
        int
    )

    for st in structure.StructureReader(
        str(poseviewer)
    ):
        # ---------------------------------------------
        # Receptor CT
        # ---------------------------------------------

        if st.property.get(
            "b_glide_receptor",
            False,
        ):
            if receptor is not None:
                raise RuntimeError(
                    "Multiple receptor CTs found in "
                    f"{poseviewer}"
                )

            receptor = st.copy()
            continue

        # ---------------------------------------------
        # Docked ligand CT
        # ---------------------------------------------

        score = glide_score(
            st
        )

        if score is None:
            continue

        cid = compound_id(
            st
        )

        if not cid:
            logger.warning(
                "Skipping scored MCSS pose without "
                "compound ID in %s",
                poseviewer,
            )
            continue

        cid = str(
            cid
        )

        pose_count_by_compound[
            cid
        ] += 1

        pose_rank = (
            pose_count_by_compound[
                cid
            ]
        )

        pose_key = (
            f"{cid}__pose{pose_rank:04d}"
        )

        ligand_st = st.copy()

        # Give Prime a human-readable unique title as
        # a backup identifier.
        ligand_st.title = pose_key

        ligand_st.property[
            MMGBSA_POSE_KEY_PROPERTY
        ] = pose_key

        ligand_st.property[
            MMGBSA_COMPOUND_ID_PROPERTY
        ] = cid

        ligand_st.property[
            MMGBSA_POSE_RANK_PROPERTY
        ] = pose_rank

        ligand_records.append(
            {
                "mmgbsa_pose_key": pose_key,
                "compound_id": cid,
                "mcss_pose_rank": pose_rank,
                "docking_score": float(
                    score
                ),
                "structure": ligand_st,
            }
        )

    if receptor is None:
        raise RuntimeError(
            "No receptor CT found in MCSS "
            f"PoseViewer: {poseviewer}"
        )

    if not ligand_records:
        raise RuntimeError(
            "No scored docked ligand poses found in "
            f"MCSS PoseViewer: {poseviewer}"
        )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with structure.StructureWriter(
        str(output_file)
    ) as writer:

        writer.append(
            receptor
        )

        for record in ligand_records:
            writer.append(
                record[
                    "structure"
                ]
            )

    manifest = pd.DataFrame(
        [
            {
                "mmgbsa_pose_key": (
                    record[
                        "mmgbsa_pose_key"
                    ]
                ),
                "compound_id": (
                    record[
                        "compound_id"
                    ]
                ),
                "mcss_pose_rank": (
                    record[
                        "mcss_pose_rank"
                    ]
                ),
                "docking_score": (
                    record[
                        "docking_score"
                    ]
                ),
            }
            for record in ligand_records
        ]
    )

    logger.info(
        "Prepared Prime MM-GBSA input with "
        "%d docked poses: %s",
        len(manifest),
        output_file,
    )

    return manifest


def _mmgbsa_dg_bind(
    st,
) -> float | None:
    """
    Extract the total Prime MM-GBSA dG Bind value.

    The first entries cover commonly used Maestro property
    names. A conservative fallback searches MMGBSA/dG Bind
    properties while explicitly excluding decomposed energy
    terms and no-strain values.
    """

    known_properties = (
        "r_psp_MMGBSA_dG_Bind",
        "r_psp_MMGBSA_dG_bind",
    )

    for property_name in (
        known_properties
    ):
        value = st.property.get(
            property_name
        )

        if value is None:
            continue

        try:
            value = float(
                value
            )

        except (
            TypeError,
            ValueError,
        ):
            continue

        if np.isfinite(
            value
        ):
            return value

    # ---------------------------------------------
    # Conservative fallback
    # ---------------------------------------------

    candidates = []

    excluded_tokens = (
        "coulomb",
        "covalent",
        "vdw",
        "lipo",
        "solv",
        "hbond",
        "packing",
        "selfcont",
        "strain",
        "nostrain",
    )

    for property_name, raw_value in (
        st.property.items()
    ):
        normalized = (
            str(property_name)
            .lower()
            .replace("_", "")
            .replace(" ", "")
        )

        if "mmgbsa" not in normalized:
            continue

        if "dgbind" not in normalized:
            continue

        if any(
            token in normalized
            for token in excluded_tokens
        ):
            continue

        try:
            value = float(
                raw_value
            )

        except (
            TypeError,
            ValueError,
        ):
            continue

        if np.isfinite(
            value
        ):
            candidates.append(
                (
                    property_name,
                    value,
                )
            )

    if len(candidates) == 1:
        return candidates[0][1]

    if len(candidates) > 1:
        raise RuntimeError(
            "Multiple candidate total MM-GBSA "
            "dG Bind properties found: "
            + ", ".join(
                property_name
                for property_name, _ in candidates
            )
        )

    return None


def _parse_mmgbsa_output(
    output_file: Path,
    input_manifest: pd.DataFrame,
) -> pd.DataFrame:
    """
    Extract MM-GBSA dG Bind values from Prime output and
    join them to the input MCSS pose manifest.

    The returned table contains one row for every original
    MCSS input pose, including NaN where Prime did not
    produce a usable dG Bind value.
    """

    output_file = (
        Path(output_file)
        .expanduser()
        .resolve()
    )

    if not output_file.is_file():
        raise FileNotFoundError(
            f"Prime MM-GBSA output not found: "
            f"{output_file}"
        )

    rows = []

    for st in structure.StructureReader(
        str(output_file)
    ):
        # Receptor-only CTs do not represent a ligand result.
        if st.property.get(
            "b_glide_receptor",
            False,
        ):
            continue

        pose_key = st.property.get(
            MMGBSA_POSE_KEY_PROPERTY
        )

        # Backup in case Prime preserved the title but not
        # our custom Maestro structure property.
        if not pose_key:
            title = (
                st.title or ""
            ).strip()

            if title:
                pose_key = title

        if not pose_key:
            continue

        dg_bind = _mmgbsa_dg_bind(
            st
        )

        if dg_bind is None:
            logger.warning(
                "No total MM-GBSA dG Bind value "
                "found for pose %s",
                pose_key,
            )
            continue

        rows.append(
            {
                "mmgbsa_pose_key": str(
                    pose_key
                ),
                "mmgbsa_dg_bind": float(
                    dg_bind
                ),
            }
        )

    if not rows:
        raise RuntimeError(
            "Prime MM-GBSA output contained no "
            "parseable ligand dG Bind values: "
            f"{output_file}"
        )

    result_df = pd.DataFrame(
        rows
    )

    duplicated = result_df[
        "mmgbsa_pose_key"
    ].duplicated(
        keep=False
    )

    if duplicated.any():
        duplicate_keys = sorted(
            result_df.loc[
                duplicated,
                "mmgbsa_pose_key",
            ].unique()
        )

        raise RuntimeError(
            "Duplicate MM-GBSA output pose keys "
            "were found: "
            + ", ".join(
                duplicate_keys
            )
        )

    merged = input_manifest.merge(
        result_df,
        on="mmgbsa_pose_key",
        how="left",
        validate="one_to_one",
    )

    missing_count = int(
        merged[
            "mmgbsa_dg_bind"
        ].isna().sum()
    )

    if missing_count:
        logger.warning(
            "Prime MM-GBSA produced usable dG Bind "
            "values for %d/%d input poses",
            (
                len(merged)
                - missing_count
            ),
            len(merged),
        )

    return merged