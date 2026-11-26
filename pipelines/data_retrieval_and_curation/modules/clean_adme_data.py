import os
import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from tqdm.auto import tqdm
from functools import reduce
import seaborn as sns
from collections import Counter

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from sklearn.manifold import TSNE
from PIL import Image, ImageDraw
import networkx as nx

from pipeline.task_registry import register_task

from modules.utils.convert_adme import (
    convert_permeability,
    convert_unscaled_clint,
    convert_ppb,
    convert_solubility,
    convert_logp_logd,
    convert_cyp_activity,
    convert_vd,
    # scale_clint_to_liver_normalised,
    # scale_liver_clint_to_l_h_kg,
    # scale_clint_to_l_h_kg,
    # convert_metstab_half_life,
    # convert_metstab_percent_remaining,
    # convert_extraction_ratio,
    # convert_absolute_clearance,
)

from modules.utils.detect_adme import (
    detect_papp_direction,
    detect_metstab_system,
    extract_species,
    classify_permeability_system,
    classify_bioavailability_record,
    flatten_pivot_columns,
    classify_bsep_record,
    extract_inhibition_concentration,
    normalise_conc,
    extract_cyp3a4_substrate,
    train_cyp_classifier,
    extract_cyp3a4_substrate_hybrid,
    # classify_metstab_record,
    classify_metstab_endpoint_family,
    infer_metstab_matrix_with_context,
    extract_species_with_context,
    convert_bsep_activity,
    classify_oatp_record,
    convert_oatp_activity,
    classify_pxr_record,
    convert_pxr_activity,
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

# Assay categories
permeability_assays = ["caco", "mdck", "pampa", "p-gp", "bcrp", "mrp"]
permeability_assays += ["bbb", "blood brain barrier", "brain penetration", "brain/plasma", "logbb", "permeability"]
solubility_assays = ["solubility", "logs", "logS"]
metstab_assays = ["microsomal", "microsomes", "hepato", "hepatocyte", "hepatic"]
ppb_assays = ["plasma protein binding", "ppb", "protein binding"]
logp_assays = ["logp", "log p", "partition coefficient"]
logd_assays = ["logd", "distribution coefficient"]
cyp_assays = ["cyp1a2","cyp2c9","cyp2c19","cyp2d6","cyp3a4","cyp3a5","cyp2b6","cyp2e1", "cyp2c8"]
bioavailability_assays = ["bioavailability", "f%", "fraction absorbed", "oral bioavailability"]
vd_assays = ["vd", "vdss", "volume of distribution"]
mate_assays = ["mate", "multidrug and toxin extrusion"]
oatp_assays = ["oatp", "slco", "organic anion transporting polypeptide"]
oct_assays = ["oct", "slc22"]
bsep_assays = ["bsep", "bile salt export pump", "abcd11"]

@register_task(
    "harmonise_units",
    category="ADME",
    description="Harmonise units for ADME and type-B tox endpoints with before/after plots"
)
def harmonise_units(config, enriched=None):
    output_dir = config.get("output", {}).get("directory", "outputs/adme")
    plots_dir = os.path.join(output_dir, "_plots")
    os.makedirs(plots_dir, exist_ok=True)

    cleaned = {}

    # Tox assay names from config
    config_targets = config.get("targets", {})
    tox_assays = {k.strip().lower() for k in config_targets.keys()}

    for assay_name, records in tqdm(
            enriched.items(),
            desc="Harmonising assays",
            unit="assay"
        ):
        logger.debug(f"Processing assay: {assay_name}")
        records = [r for r in records if isinstance(r, dict)]
        if not records:
            continue

        lname = assay_name.lower()
        logger.debug(f"Processing assay: '{assay_name}' -> lname: '{lname}'")
        logger.debug(f"Permeability assays list: {permeability_assays}")

        units_before = [r.get("standard_units") for r in records if r.get("standard_units")]
        new_records = []

        # -----------------------------
        # ADME assays
        # -----------------------------
        if any(x in lname for x in permeability_assays):

            logger.debug(f"Matched permeability assay: {assay_name}")
            logger.debug(f"Number of raw records: {len(records)}")

            for r in records:

                system, directionality = classify_permeability_system(r)

                # Skip unclassified permeability assays
                if system is None:
                    continue

                assay_family = assay_name.lower()

                if "caco" in assay_family and system != "caco2":
                    continue

                if "mdck" in assay_family and system != "mdck":
                    continue

                if "pampa" in assay_family and system != "pampa":
                    continue

                r["permeability_system"] = system

                # Only directional systems get AB/BA labels
                if directionality == "directional":
                    r["papp_direction"] = detect_papp_direction(
                        r.get("assay_description")
                    )
                else:
                    r["papp_direction"] = "NA"

                val, unit = r.get("standard_value"), r.get("standard_units")

                new_val, new_unit = convert_permeability(val, unit)

                if new_val is None:
                    continue

                r["standard_value"] = new_val
                r["standard_units"] = new_unit

                new_records.append(r)

            logger.debug(
                f"{assay_name}: kept {len(new_records)} permeability records"
            )

        elif any(x in lname for x in solubility_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                smiles = r.get("smiles")
                new_val = convert_solubility(val, unit, smiles)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, "log10(mol/L)"
                    new_records.append(r)

        elif any(x in lname for x in logp_assays + logd_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val = convert_logp_logd(val, unit)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, "log_unitless"
                    new_records.append(r)

        elif any(x in lname for x in metstab_assays):

            # =====================================================
            # Reported unscaled intrinsic clearance only
            #
            # Retained canonical representations:
            #
            #   microsomes:
            #       uL/min/mg protein
            #
            #   hepatocytes:
            #       uL/min/10^6 cells
            #
            # No physiological scaling or derived records are
            # created in this baseline workflow.
            # =====================================================

            drop_reasons = Counter()
            unit_counter = Counter()
            stype_counter = Counter()
            species_counter = Counter()
            matrix_counter = Counter()
            converted_unit_counter = Counter()
            relation_counter = Counter()

            rejected_examples = []

            stage_counts = Counter()

            # Records rejected before Clint conversion.
            rejected_non_clint_unit_counter = Counter()
            rejected_non_clint_type_counter = Counter()

            # Records classified as Clint but rejected by the
            # strict reported-unscaled unit converter.
            rejected_clint_unit_counter = Counter()
            rejected_clint_type_counter = Counter()

            # Raw units successfully harmonised as reported
            # unscaled microsomal or hepatocyte Clint.
            accepted_clint_unit_counter = Counter()

            clint_preconversion_records = []

            def record_rejection(
                reason,
                record,
                species=None,
                matrix=None,
                endpoint_family=None,
                converted_unit=None,
            ):
                """
                Count every rejection and retain up to 1000 example rows.

                Summary tables are written only after all records in the
                retrieval bucket have been processed.
                """
                drop_reasons[reason] += 1

                if len(rejected_examples) >= 1000:
                    return

                rejected_examples.append({
                    "reason": reason,
                    "retrieval_assay_bucket": assay_name,
                    "activity_id": record.get("activity_id"),
                    "activity_chembl_id":
                        record.get("activity_chembl_id"),
                    "assay_chembl_id":
                        record.get("assay_chembl_id"),
                    "document_chembl_id":
                        record.get("document_chembl_id"),
                    "smiles": record.get("smiles"),
                    "standard_type":
                        record.get("standard_type"),
                    "standard_value":
                        record.get("standard_value"),
                    "standard_units":
                        record.get("standard_units"),
                    "relation":
                        record.get("relation"),
                    "assay_description":
                        record.get("assay_description"),
                    "target_organism":
                        record.get("target_organism"),
                    "species": species,
                    "matrix": matrix,
                    "endpoint_family": endpoint_family,
                    "converted_unit": converted_unit,
                })

            for r in records:

                stage_counts["raw_records"] += 1

                raw_value = r.get("standard_value")
                raw_unit = r.get("standard_units")
                raw_type = r.get("standard_type")

                unit_counter[str(raw_unit)] += 1
                stype_counter[str(raw_type)] += 1
                relation_counter[str(r.get("relation"))] += 1

                # ---------------------------------------------
                # 1. Keep intrinsic-clearance records only.
                # ---------------------------------------------

                endpoint_family = classify_metstab_endpoint_family(
                    r,
                    assay_name=assay_name,
                )

                stage_counts[
                    f"classified_as_{endpoint_family}"
                ] += 1                            

                if endpoint_family != "clint":

                    stage_counts[
                        "rejected_not_clint"
                    ] += 1

                    rejected_non_clint_unit_counter[
                        str(raw_unit)
                    ] += 1

                    rejected_non_clint_type_counter[
                        str(raw_type)
                    ] += 1

                    record_rejection(
                        reason=(
                            "not_intrinsic_clearance_"
                            f"{endpoint_family}"
                        ),
                        record=r,
                        endpoint_family=endpoint_family,
                    )
                    continue

                stage_counts[
                    "clint_before_unit_filter"
                ] += 1                

                clint_preconversion_records.append({
                    "retrieval_assay_bucket": assay_name,
                    "activity_id": r.get("activity_id"),
                    "activity_chembl_id":
                        r.get("activity_chembl_id"),
                    "assay_chembl_id":
                        r.get("assay_chembl_id"),
                    "document_chembl_id":
                        r.get("document_chembl_id"),
                    "smiles": r.get("smiles"),
                    "standard_type": raw_type,
                    "standard_value": raw_value,
                    "standard_units": raw_unit,
                    "relation": r.get("relation"),
                    "assay_description":
                        r.get("assay_description"),
                    "target_organism":
                        r.get("target_organism"),
                    "species_from_description":
                        extract_species(
                            r.get("assay_description")
                        ),
                    "species_from_target":
                        extract_species(
                            r.get("target_organism")
                        ),
                })

                # ---------------------------------------------
                # 2. Harmonise explicitly reported unscaled
                #    intrinsic clearance.
                #
                # The converter returns:
                #   - harmonised value
                #   - canonical unit
                #   - matrix implied by the denominator
                # ---------------------------------------------

                new_value, new_unit, matrix = (
                    convert_unscaled_clint(
                        raw_value,
                        raw_unit,
                    )
                )

                if new_value is None:

                    stage_counts[
                        "failed_unscaled_unit_conversion"
                    ] += 1

                    rejected_clint_unit_counter[
                        str(raw_unit)
                    ] += 1

                    rejected_clint_type_counter[
                        str(raw_type)
                    ] += 1

                    record_rejection(
                        reason=(
                            "not_explicitly_unit_resolved_"
                            "unscaled_clint"
                        ),
                        record=r,
                        endpoint_family=endpoint_family,
                    )
                    continue

                # Defensive validation. Only these two biological
                # representations may leave convert_unscaled_clint().
                expected_unit_by_matrix = {
                    "microsome": "uL/min/mg protein",
                    "hepatocyte": "uL/min/10^6 cells",
                }

                expected_unit = expected_unit_by_matrix.get(
                    matrix
                )

                if expected_unit is None:

                    stage_counts[
                        "failed_matrix_unit_validation"
                    ] += 1

                    rejected_clint_unit_counter[
                        str(raw_unit)
                    ] += 1

                    rejected_clint_type_counter[
                        str(raw_type)
                    ] += 1

                    record_rejection(
                        reason=(
                            "unexpected_unscaled_clint_matrix_"
                            f"{matrix}"
                        ),
                        record=r,
                        matrix=matrix,
                        endpoint_family=endpoint_family,
                        converted_unit=new_unit,
                    )
                    continue

                if new_unit != expected_unit:

                    stage_counts[
                        "failed_matrix_unit_validation"
                    ] += 1

                    rejected_clint_unit_counter[
                        str(raw_unit)
                    ] += 1

                    rejected_clint_type_counter[
                        str(raw_type)
                    ] += 1

                    record_rejection(
                        reason=(
                            "matrix_unit_mismatch_"
                            f"{matrix}"
                        ),
                        record=r,
                        matrix=matrix,
                        endpoint_family=endpoint_family,
                        converted_unit=new_unit,
                    )
                    continue

                # ---------------------------------------------
                # 5. Determine species independently.
                # ---------------------------------------------

                species = extract_species_with_context(r)

                if species is None:
                    species = "unknown"

                species_counter[str(species)] += 1
                matrix_counter[str(matrix)] += 1
                converted_unit_counter[str(new_unit)] += 1

                # ---------------------------------------------
                # 6. Construct exactly one retained record for
                #    exactly one reported ChEMBL activity.
                # ---------------------------------------------

                rec = r.copy()

                rec["standard_value"] = new_value
                rec["standard_units"] = new_unit

                rec["species"] = species
                rec["species_from_assay_description"] = extract_species(
                    r.get("assay_description")
                )

                rec["species_from_target_organism"] = extract_species(
                    r.get("target_organism")
                )                
                rec["metstab_matrix"] = matrix
                rec["metstab_endpoint_family"] = "clint"
                rec["clint_representation"] = (
                    "unscaled_in_vitro"
                )

                rec["metstab_system"] = "in_vitro"
                rec["metstab_value_origin"] = "reported"

                rec["retrieval_assay_bucket"] = assay_name

                # Preserve raw source values for auditing.
                rec["source_standard_value"] = raw_value
                rec["source_standard_units"] = raw_unit
                rec["source_standard_type"] = raw_type

                # Stable identifier used later to remove duplicate
                # activities retrieved through more than one broad
                # keyword bucket.
                activity_id = r.get("activity_id")
                activity_chembl_id = r.get(
                    "activity_chembl_id"
                )

                if activity_id is not None:
                    rec["metstab_source_id"] = (
                        f"activity_id:{activity_id}"
                    )

                elif activity_chembl_id is not None:
                    rec["metstab_source_id"] = (
                        "activity_chembl_id:"
                        f"{activity_chembl_id}"
                    )

                else:
                    rec["metstab_source_id"] = "|".join([
                        str(r.get("assay_chembl_id", "")),
                        str(r.get("document_chembl_id", "")),
                        str(r.get("smiles", "")),
                        str(raw_type),
                        str(raw_value),
                        str(raw_unit),
                    ])

                stage_counts["retained_rows"] += 1
                accepted_clint_unit_counter[str(raw_unit)] += 1

                new_records.append(rec)

            # -----------------------------------------------------
            # Validate stage accounting before writing outputs.
            # -----------------------------------------------------

            classified_clint = stage_counts[
                "classified_as_clint"
            ]

            converter_input = stage_counts[
                "clint_before_unit_filter"
            ]

            failed_conversion = stage_counts[
                "failed_unscaled_unit_conversion"
            ]

            failed_validation = stage_counts[
                "failed_matrix_unit_validation"
            ]

            retained = stage_counts[
                "retained_rows"
            ]

            if classified_clint != converter_input:
                raise RuntimeError(
                    f"{assay_name}: Clint classifier/converter "
                    "count mismatch: "
                    f"classified_as_clint={classified_clint:,}, "
                    f"clint_before_unit_filter={converter_input:,}"
                )

            accounted_clint = (
                failed_conversion
                + failed_validation
                + retained
            )

            if converter_input != accounted_clint:
                raise RuntimeError(
                    f"{assay_name}: incomplete Clint accounting: "
                    f"converter_input={converter_input:,}, "
                    f"failed_conversion={failed_conversion:,}, "
                    f"failed_validation={failed_validation:,}, "
                    f"retained={retained:,}, "
                    f"accounted={accounted_clint:,}"
                )

            if len(clint_preconversion_records) != converter_input:
                raise RuntimeError(
                    f"{assay_name}: pre-conversion audit length "
                    "does not match converter input: "
                    f"audit_rows="
                    f"{len(clint_preconversion_records):,}, "
                    f"converter_input={converter_input:,}"
                )

            if len(new_records) != retained:
                raise RuntimeError(
                    f"{assay_name}: retained-record count mismatch: "
                    f"new_records={len(new_records):,}, "
                    f"retained_counter={retained:,}"
                )

            if clint_preconversion_records:

                clint_preconversion_df = pd.DataFrame(
                    clint_preconversion_records
                )

                clint_unit_summary = (
                    clint_preconversion_df
                    .groupby(
                        [
                            "standard_type",
                            "standard_units",
                        ],
                        dropna=False,
                    )
                    .agg(
                        rows=("standard_value", "size"),
                        unique_smiles=("smiles", "nunique"),
                        unique_assays=(
                            "assay_chembl_id",
                            "nunique",
                        ),
                    )
                    .reset_index()
                    .sort_values(
                        "rows",
                        ascending=False,
                    )
                )

                clint_unit_summary.to_csv(
                    os.path.join(
                        output_dir,
                        (
                            f"{assay_name.replace(' ', '_')}"
                            "_clint_preconversion_unit_summary.csv"
                        ),
                    ),
                    index=False,
                )                

                clint_preconversion_path = os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_clint_preconversion.csv"
                    ),
                )

                clint_preconversion_df.to_csv(
                    clint_preconversion_path,
                    index=False,
                )

                logger.info(
                    f"{assay_name}: wrote "
                    f"{len(clint_preconversion_df):,} "
                    "pre-conversion Clint candidates to "
                    f"{clint_preconversion_path}"
                )                

            # -----------------------------------------------------
            # Save the harmonised long-form rows for this retrieval
            # bucket before aggregation by SMILES.
            # -----------------------------------------------------

            if new_records:

                metstab_long_df = pd.DataFrame(new_records)

                long_path = os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_unscaled_clint_long.csv"
                    ),
                )

                metstab_long_df.to_csv(
                    long_path,
                    index=False,
                )

                coverage_summary = (
                    metstab_long_df
                    .groupby(
                        [
                            "species",
                            "metstab_matrix",
                            "standard_units",
                        ],
                        dropna=False,
                    )
                    .agg(
                        rows=("standard_value", "size"),
                        unique_activities=(
                            "metstab_source_id",
                            "nunique",
                        ),
                        unique_smiles=(
                            "smiles",
                            "nunique",
                        ),
                    )
                    .reset_index()
                    .sort_values(
                        "rows",
                        ascending=False,
                    )
                )

                summary_path = os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_unscaled_clint_summary.csv"
                    ),
                )

                coverage_summary.to_csv(
                    summary_path,
                    index=False,
                )

                logger.info(
                    f"{assay_name}: wrote reported unscaled "
                    f"Clint rows to {long_path}"
                )

                logger.info(
                    f"{assay_name}: wrote unscaled Clint "
                    f"summary to {summary_path}"
                )

            logger.debug(
                f"{assay_name}: raw records = {len(records)}"
            )
            logger.debug(
                f"{assay_name}: retained reported unscaled "
                f"Clint records = {len(new_records)}"
            )
            logger.debug(
                f"{assay_name}: drop reasons = "
                f"{dict(drop_reasons)}"
            )
            logger.debug(
                f"{assay_name}: raw unit counts = "
                f"{unit_counter.most_common(30)}"
            )
            logger.debug(
                f"{assay_name}: standard-type counts = "
                f"{stype_counter.most_common(30)}"
            )
            logger.debug(
                f"{assay_name}: species counts = "
                f"{dict(species_counter)}"
            )
            logger.debug(
                f"{assay_name}: matrix counts = "
                f"{dict(matrix_counter)}"
            )
            logger.debug(
                f"{assay_name}: canonical unit counts = "
                f"{dict(converted_unit_counter)}"
            )
            logger.debug(
                f"{assay_name}: relation counts = "
                f"{relation_counter.most_common(20)}"
            )

            if rejected_examples:

                reject_path = os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_unscaled_clint_rejected.csv"
                    ),
                )

                pd.DataFrame(
                    rejected_examples
                ).to_csv(
                    reject_path,
                    index=False,
                )

                logger.info(
                    f"{assay_name}: wrote rejected metabolic-"
                    f"stability rows to {reject_path}"
                )        

            # -----------------------------------------------------
            # Complete summaries. Unlike rejected_examples, these
            # are not capped at 1000 rows because they are counters.
            # -----------------------------------------------------

            accepted_clint_unit_summary = pd.DataFrame(
                list(accepted_clint_unit_counter.items()),
                columns=["raw_unit", "count"],
            )

            if not accepted_clint_unit_summary.empty:
                accepted_clint_unit_summary = (
                    accepted_clint_unit_summary
                    .sort_values(
                        "count",
                        ascending=False,
                    )
                )

            accepted_clint_unit_summary.to_csv(
                os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_accepted_clint_units.csv"
                    ),
                ),
                index=False,
            )

            rejected_clint_unit_summary = pd.DataFrame(
                list(rejected_clint_unit_counter.items()),
                columns=["raw_unit", "count"],
            )

            if not rejected_clint_unit_summary.empty:
                rejected_clint_unit_summary = (
                    rejected_clint_unit_summary
                    .sort_values(
                        "count",
                        ascending=False,
                    )
                )

            rejected_clint_unit_summary.to_csv(
                os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_rejected_clint_units.csv"
                    ),
                ),
                index=False,
            )

            rejected_clint_type_summary = pd.DataFrame(
                list(rejected_clint_type_counter.items()),
                columns=["standard_type", "count"],
            )

            if not rejected_clint_type_summary.empty:
                rejected_clint_type_summary = (
                    rejected_clint_type_summary
                    .sort_values(
                        "count",
                        ascending=False,
                    )
                )

            rejected_clint_type_summary.to_csv(
                os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_rejected_clint_standard_types.csv"
                    ),
                ),
                index=False,
            )

            non_clint_unit_summary = pd.DataFrame(
                list(rejected_non_clint_unit_counter.items()),
                columns=["raw_unit", "count"],
            )

            if not non_clint_unit_summary.empty:
                non_clint_unit_summary = (
                    non_clint_unit_summary
                    .sort_values(
                        "count",
                        ascending=False,
                    )
                )

            non_clint_unit_summary.to_csv(
                os.path.join(
                    output_dir,
                    (
                        f"{assay_name.replace(' ', '_')}"
                        "_non_clint_units.csv"
                    ),
                ),
                index=False,
            )

            logger.info(
                f"{assay_name}: stage counts = "
                f"{dict(stage_counts)}"
            )

            logger.info(
                f"{assay_name}: top accepted Clint raw units = "
                f"{accepted_clint_unit_counter.most_common(30)}"
            )

            logger.info(
                f"{assay_name}: top rejected Clint raw units = "
                f"{rejected_clint_unit_counter.most_common(50)}"
            )

            logger.info(
                f"{assay_name}: top rejected Clint standard types = "
                f"{rejected_clint_type_counter.most_common(30)}"
            )

            logger.debug(
                f"{assay_name}: top non-Clint raw units = "
                f"{rejected_non_clint_unit_counter.most_common(30)}"
            )

            logger.debug(
                f"{assay_name}: top non-Clint standard types = "
                f"{rejected_non_clint_type_counter.most_common(30)}"
            )        

        elif any(x in lname for x in ppb_assays):
            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val = convert_ppb(val, unit)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, "fraction_unbound"
                    r["species"] = extract_species(r.get("target_organism") or r.get("assay_description"))
                    new_records.append(r)

        elif any(x in lname for x in cyp_assays):  # CYPs are ADME (type A)
            
            sample = [
                r.get("assay_description")
                for r in records
            ]

            # Existing crude sanity checks
            logger.debug(f"midazolam raw mentions: {sum('midazolam' in str(t).lower() for t in sample)}")
            logger.debug(f"testosterone raw mentions: {sum('testosterone' in str(t).lower() for t in sample)}")

            # NEW: test substrate classifier directly on raw records
            # Build weak labels using rules
            texts = []
            labels = []

            for r in records:
                desc = r.get("assay_description")
                if not desc:
                    continue

                label = extract_cyp3a4_substrate(desc)

                if label in ["midazolam_like", "testosterone_like"]:
                    texts.append(desc)
                    labels.append(label)

            # Train ML model if enough data
            if len(texts) > 50:
                train_cyp_classifier(texts, labels)
                logger.debug(f"Trained CYP classifier on {len(texts)} samples")

            logger.debug(f"substrate classification counts: {Counter(labels)}")

            for r in records:
                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, endpoint, new_unit = convert_cyp_activity(val, unit)

                if new_val is None:
                    continue

                r["standard_value"], r["standard_units"] = new_val, new_unit
                r["cyp_endpoint"] = endpoint

                def map_substrate_label(label):
                    return {
                        "midazolam_like": "mdz",
                        "testosterone_like": "tes"
                    }.get(label, label)

                # Only split CYP3A4
                if "cyp3a4" in lname:
                    raw_label = extract_cyp3a4_substrate_hybrid(
                    r.get("assay_description")
                )
                    
                    r["cyp_substrate"] = map_substrate_label(raw_label)
                else:
                    r["cyp_substrate"] = "na"   # important: keeps grouping stable

                if endpoint == "IC50":
                    r["IC50"] = new_val
                    r["IC50_unit"] = new_unit
                elif endpoint == "inhibition":
                    r["inhibition"] = new_val
                    r["inhibition_unit"] = new_unit

                r["tox_endpoint"] = endpoint
                new_records.append(r)
            logger.debug(f"{assay_name}: kept after conversion = {len(new_records)}")
        # elif any(x in lname for x in bbb_assays):
        #     for r in records:
        #         val, unit = r.get("standard_value"), r.get("standard_units")
        #         try:
        #             val = float(val)
        #         except (TypeError, ValueError):
        #             continue
                
        #         # Convert units to 10^-6 cm/s if possible
        #         if unit is not None:
        #             val, unit = convert_permeability(val, unit)
                
        #         if val is not None:
        #             r["standard_value"], r["standard_units"] = val, unit or "10^-6 cm/s"
        #             r["harmonised_endpoint"] = "BBB_permeability"
        #             new_records.append(r)

        elif any(x in lname for x in bioavailability_assays):
            for r in records:
                val = r.get("standard_value")
                try:
                    val = float(val)
                except (TypeError, ValueError):
                    continue
                if not classify_bioavailability_record(r):
                    continue
                # Assume %; store species if available
                r["standard_value"], r["standard_units"] = val, "%"
                # r["species"] = extract_species(r.get("target_organism") or r.get("assay_description"))
                r["species"] = extract_species(
                    r.get("target_organism") or r.get("assay_description"),
                    context="pk"
                )
                r["harmonised_endpoint"] = "F_percent"
                new_records.append(r)

        elif any(x in lname for x in vd_assays):
            for r in records:
                val = r.get("standard_value")
                unit = r.get("standard_units")
                try:
                    val = float(val)
                except (TypeError, ValueError):
                    continue
                r["species"] = extract_species(
                    r.get("target_organism") or r.get("assay_description"),
                    context="pk"
                )
                val, unit = convert_vd(val, unit)
                r["standard_value"], r["standard_units"] = val, unit
                r["harmonised_endpoint"] = "Vd"
                new_records.append(r)
        
        elif any(x in lname for x in bsep_assays):

            for r in records:

                endpoint = classify_bsep_record(r)

                if endpoint is None:
                    continue

                val = r.get("standard_value")
                unit = r.get("standard_units")

                new_val, endpoint_type, new_unit = convert_bsep_activity(
                    val,
                    unit,
                    endpoint
                )

                if new_val is None:
                    continue

                r["standard_value"] = new_val
                r["standard_units"] = new_unit

                r["bsep_endpoint"] = endpoint_type
                r["species"] = extract_species(
                    r.get("target_organism") or r.get("assay_description")
                )

                # Optional but strongly recommended
                desc = str(r.get("assay_description", "")).lower()

                if "taurocholate" in desc:
                    r["bsep_substrate"] = "taurocholate"
                else:
                    r["bsep_substrate"] = "unknown"

                new_records.append(r)


        elif any(x in lname for x in mate_assays):
            for r in records:
                val = r.get("standard_value")
                unit = r.get("standard_units")
                try:
                    val = float(val)
                except (TypeError, ValueError):
                    continue
                # Optionally: convert units to nM / ng/mL if needed
                r["standard_value"], r["standard_units"] = val, unit
                r["harmonised_endpoint"] = "MATE_exposure"
                r["species"] = extract_species(r.get("target_organism") or r.get("assay_description"))
                new_records.append(r)

        elif any(x in lname for x in oatp_assays):

            for r in records:

                endpoint = classify_oatp_record(r)

                if endpoint is None:
                    continue

                val = r.get("standard_value")
                unit = r.get("standard_units")

                new_val, endpoint_type, new_unit = convert_oatp_activity(
                    val,
                    unit,
                    endpoint
                )

                if new_val is None:
                    continue

                r["standard_value"] = new_val
                r["standard_units"] = new_unit

                r["oatp_endpoint"] = endpoint_type

                # ---------------------------------
                # Parse inhibition concentration
                # ---------------------------------

                if endpoint_type == "inhibition":

                    conc = extract_inhibition_concentration(
                        r.get("assay_description")
                    )

                    if conc is not None:
                        conc = normalise_conc(conc)

                    else:
                        # transporter assays are often 10 uM
                        # but we should only use that as a fallback
                        conc = 10

                    r["inhibition_conc_uM"] = conc

                r["species"] = extract_species(
                    r.get("target_organism") or r.get("assay_description")
                )

                # # Optional transporter subtype split
                # desc = str(r.get("assay_description", "")).lower()

                # if "oatp1b1" in desc or "slco1b1" in desc:
                #     r["oatp_subtype"] = "1B1"

                # elif "oatp1b3" in desc or "slco1b3" in desc:
                #     r["oatp_subtype"] = "1B3"

                # elif "oatp2b1" in desc or "slco2b1" in desc:
                #     r["oatp_subtype"] = "2B1"

                # else:
                #     r["oatp_subtype"] = "unknown"

                new_records.append(r)

        elif any(x in lname for x in oct_assays):
            for r in records:
                val = r.get("standard_value")
                try:
                    val = float(val)
                except (TypeError, ValueError):
                    val = None  # Keep the entry but numeric value may be missing
                r["standard_value"] = val
                r["standard_units"] = r.get("standard_units")
                r["harmonised_endpoint"] = "OCT_transport"
                r["species"] = extract_species(r.get("target_organism") or r.get("assay_description"))
                new_records.append(r)                
        # -----------------------------
        # Type-B tox endpoints (hERG, Nav1.5, etc.) from config
        # -----------------------------
        elif any(lname == t for t in tox_assays):

            # -----------------------------
            # PXR special handling
            # -----------------------------

            if lname == "pxr":

                new_records = []

                for r in records:

                    endpoint = classify_pxr_record(r)

                    if endpoint is None:
                        continue

                    val = r.get("standard_value")
                    unit = r.get("standard_units")

                    new_val, endpoint_type, new_unit = convert_pxr_activity(
                        val,
                        unit,
                        endpoint
                    )

                    if new_val is None:
                        continue

                    r["standard_value"] = new_val
                    r["standard_units"] = new_unit
                    r["tox_type"] = endpoint_type

                    # ---------------------------------
                    # Parse concentration from assay text
                    # ---------------------------------

                    conc = extract_inhibition_concentration(
                        r.get("assay_description")
                    )

                    if conc is not None:
                        conc = normalise_conc(conc)
                        r["inhibition_conc_uM"] = conc

                    # activation assays are experimental
                    r["is_estimated"] = False

                    new_records.append(r)

                cleaned[assay_name] = new_records
                continue

            for r in records:
                val = r.get("standard_value")
                unit = r.get("standard_units")
                stype = str(r.get("standard_type", "")).strip().upper()

                if val is None:
                    continue

                try:
                    val = float(val)
                except:
                    continue

                endpoint_type = None

                # concentration endpoints
                if stype in {"IC50", "EC50", "KI", "KD"}:
                    endpoint_type = stype

                # percent inhibition
                elif "%" in str(unit).lower() or "INHIBITION" in stype:
                    endpoint_type = "inhibition"

                    conc = extract_inhibition_concentration(r.get("assay_description"))
                    if conc is not None:
                        conc = normalise_conc(conc)
                        r["inhibition_conc_uM"] = conc
                    else:
                        continue

                    r["is_estimated"] = False

                if endpoint_type is None:
                    continue

                # unit normalisation
                if endpoint_type in {"IC50", "EC50", "KI", "KD"}:
                    u = str(unit).lower()

                    factor = 1
                    if "um" in u or "μm" in u:
                        factor = 1e3
                    elif "mm" in u:
                        factor = 1e6

                    val = val * factor
                    unit = "nM"

                r["standard_value"] = val
                r["standard_units"] = unit
                r["tox_type"] = endpoint_type
                r["endpoint"] = assay_name

                # debug log for verification
                # logger.debug({
                #     "smiles": r.get("smiles"),
                #     "assay_name": assay_name,
                #     "tox_type": r["tox_type"],
                #     "standard_value": r["standard_value"],
                #     "standard_units": r["standard_units"]
                # })

                # -----------------------------
                # IC50 → pseudo % inhibition
                # -----------------------------
                if endpoint_type in {"IC50", "EC50", "KI", "KD"}:

                    # choose reference concentration (uM)
                    ref_conc_uM = 10.0

                    # convert IC50 (currently nM) → uM
                    ic50_uM = val / 1000.0 if val is not None else None

                    if ic50_uM and ic50_uM > 0:
                        inh = 100.0 * ref_conc_uM / (ref_conc_uM + ic50_uM)

                        # create record derived from DR data
                        inh_record = r.copy()
                        inh_record["standard_value"] = inh
                        inh_record["standard_units"] = "%"
                        inh_record["tox_type"] = "inhibition"
                        inh_record["inhibition_conc_uM"] = ref_conc_uM

                        inh_record["is_estimated"] = True

                        new_records.append(inh_record)
                
                new_records.append(r)

        # -----------------------------
        # Fallback
        # -----------------------------
        else:
            for r in records:
                val = r.get("standard_value")
                if val is None:
                    continue
                r["endpoint"] = assay_name
                new_records.append(r)
    
        # -----------------------------
        # Plot before/after proportions
        # -----------------------------

        # logger.debug(f"{assay_name} endpoint types: {Counter(r['tox_type'] for r in new_records)}")

        units_after = [r.get("standard_units") for r in new_records if r.get("standard_units")]
        if units_before or units_after:
            all_units = sorted(set(units_before) | set(units_after))
            before_counts = Counter(units_before)
            after_counts = Counter(units_after)
            n_before = sum(before_counts.values()) or 1
            n_after = sum(after_counts.values()) or 1
            proportions_before = [before_counts.get(u,0)/n_before for u in all_units]
            proportions_after = [after_counts.get(u,0)/n_after for u in all_units]

            x = range(len(all_units))
            width = 0.35

            plt.figure(figsize=(8,4))
            plt.bar([i - width/2 for i in x], proportions_before, width=width, color='skyblue', label='Before')
            plt.bar([i + width/2 for i in x], proportions_after, width=width, color='lightgreen', label='After')
            plt.xticks(x, all_units, rotation=45, fontsize=6)
            plt.ylabel("Proportion of measurements")
            plt.title(f"{assay_name} unit harmonisation (before/after)")
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"{assay_name}_units_before_after.png"), dpi=300)
            plt.close()
        
        # debug checks of assay data

        # if lname in tox_assays:
        #     logger.debug(f"{assay_name} example records after harmonisation:")
        #     for r in new_records[:5]:
        #         logger.debug({
        #             "standard_type": r.get("standard_type"),
        #             "endpoint": r.get("endpoint"),
        #             "standard_units": r.get("standard_units"),
        #             "standard_value": r.get("standard_value")
        #         })

        cleaned[assay_name] = new_records

    return cleaned

@register_task(
    "build_multitask_dataset",
    category="ADME",
    description="Construct multi-task dataset including ADME (type-A) and type-B tox endpoints"
)
def build_multitask_dataset(config, cleaned=None):
    output_cfg = config.get("output", {})
    out_dir = output_cfg.get("directory", "outputs/adme")
    filename = output_cfg.get("filename", "chembl_mtl_dataset.csv")
    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, filename)

    dfs = []
    unscaled_clint_frames = []
    # allowed_species = {"Human", "Mouse", "Rat", "Monkey", "unknown"}

    # Helper to recover tox values (IC50, inhibition, etc.)
    def recover_tox_value(row, endpoint_col):
        val = row["standard_value"]
        unit = str(row.get("standard_units", "")).lower() if row.get("standard_units") else ""
        ep = row[endpoint_col]
        if unit in ["", "none", "no unit"] and ep in {"IC50", "EC50", "Ki", "Kd"}:
            val = val * 1000
            unit = "nM"
        elif "um" in unit:
            val = val * 1000
            unit = "nM"
        elif unit == "m":
            val = val * 1e9
            unit = "nM"
        return val, ep, unit

    for assay_name, records in tqdm(
            cleaned.items(),
            desc="Building multitask dataset",
            unit="assay"
        ):
        records = [r for r in records if isinstance(r, dict)]
        if not records:
            continue

        df = pd.DataFrame(records)

        # if (
        #     "metstab_endpoint_family" in df.columns
        #     and "clint_scale" in df.columns
        # ):
        #     clint_mask = df["metstab_endpoint_family"].eq("clint")

        #     invalid_scale_mask = (
        #         clint_mask
        #         & ~df["clint_scale"].isin({
        #             "unscaled_in_vitro",
        #             "liver_normalised",
        #             "body_normalised",
        #             "absolute",
        #         })
        #     )

        #     if invalid_scale_mask.any():
        #         bad_values = (
        #             df.loc[invalid_scale_mask, "clint_scale"]
        #             .value_counts(dropna=False)
        #             .to_dict()
        #         )

        #         raise ValueError(
        #             f"{assay_name}: non-canonical Clint scale labels "
        #             f"reached dataset construction: {bad_values}"
        #         )

        # debug checks of assay data

        logger.debug(f"{assay_name} DataFrame columns after harmonisation: {df.columns.tolist()}")

        # logger.debug(f"Processing assay: {assay_name}")
        # logger.debug(f"df columns: {df.columns.tolist()}")
        # if len(records) > 0:
        #     logger.debug(f"Sample record: {records[0]}")

        # if assay_name.lower() in {"herg", "nav1.5"}:
        #     logger.debug(f"{assay_name} columns in DataFrame: {df.columns.tolist()}")
        #     logger.debug(df[["standard_type","endpoint","standard_units"]].head())
    
        df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
        df = df.dropna(subset=["standard_value"])
        if df.empty:
            continue

        # ---------------------------------------------------------
        # Collect reported unscaled Clint records across the broad
        # retrieval buckets. They will be deduplicated and pivoted
        # once after the assay loop.
        # ---------------------------------------------------------

        if {
            "metstab_endpoint_family",
            "clint_representation",
            "metstab_matrix",
            "standard_units",
            "species",
        }.issubset(df.columns):

            is_unscaled_clint = (
                df["metstab_endpoint_family"].eq("clint")
                & df["clint_representation"].eq(
                    "unscaled_in_vitro"
                )
            )

            df_unscaled_clint = df.loc[
                is_unscaled_clint
            ].copy()

            if not df_unscaled_clint.empty:
                unscaled_clint_frames.append(
                    df_unscaled_clint
                )

            # This assay bucket contains only the metabolic-
            # stability rows created by the replacement branch.
            # Do not also pass them through the generic ADME pivot.
            continue

        # -----------------------------
        # Type-B assays (tox endpoints: hERG, Nav, etc.)
        # -----------------------------
        if "tox_type" in df.columns:

            logger.debug(f"Building pivot for tox assay: {assay_name}")

            # --- SPLIT DATA ---

            # percent inhibition style endpoints
            df_inh = df[
                df["tox_type"].isin([
                    "inhibition",
                    "inhibition_percent"
                ])
            ].copy()

            # activation / induction style endpoints
            df_act = df[
                df["tox_type"].isin([
                    "activation",
                    "activation_percent",
                    "activation_fold"
                ])
            ].copy()

            # concentration potency endpoints
            df_aff = df[
                df["tox_type"].isin([
                    "IC50",
                    "EC50",
                    "KI",
                    "KD",
                    "activation_EC50",
                    "inhibition_EC50",
                ])
            ].copy()

            pivots = []

            # -----------------------------
            # 1. Inhibition block
            # -----------------------------
            if not df_inh.empty:

                if "inhibition_conc_uM" not in df_inh.columns:
                    df_inh["inhibition_conc_uM"] = -1
                else:
                    df_inh["inhibition_conc_uM"] = (
                        df_inh["inhibition_conc_uM"]
                        .fillna(-1)
                    )

                if "is_estimated" not in df_inh.columns:
                    df_inh["is_estimated"] = False

                df_grouped_inh = (
                    df_inh.groupby(
                        ["smiles", "inhibition_conc_uM", "is_estimated"]
                    )["standard_value"]
                    .agg(["mean", "std", "count"])
                    .reset_index()
                )

                def make_inh_label(row):

                    label = (
                        f"inhibition_{int(row['inhibition_conc_uM'])}uM"
                        if row["inhibition_conc_uM"] > 0
                        else "inhibition"
                    )

                    if row["is_estimated"]:
                        label += "_est"
                    else:
                        label += "_exp"

                    return label

                df_grouped_inh["tox_label"] = (
                    df_grouped_inh.apply(make_inh_label, axis=1)
                )

                df_pivot_inh = df_grouped_inh.pivot(
                    index="smiles",
                    columns="tox_label"
                )

                df_pivot_inh = flatten_pivot_columns(
                    df_pivot_inh,
                    assay_name
                )

                pivots.append(df_pivot_inh)

            # -----------------------------
            # 1b. Activation block
            # -----------------------------
            if not df_act.empty:

                df_act["inhibition_conc_uM"] = (
                    df_act["inhibition_conc_uM"]
                    if "inhibition_conc_uM" in df_act.columns
                    else -1
                )

                df_act["inhibition_conc_uM"] = (
                    pd.to_numeric(
                        df_act["inhibition_conc_uM"],
                        errors="coerce"
                    )
                    .fillna(-1)
                )

                if "is_estimated" not in df_act.columns:
                    df_act["is_estimated"] = False

                df_grouped_act = (
                    df_act.groupby(
                        ["smiles", "inhibition_conc_uM", "is_estimated"]
                    )["standard_value"]
                    .agg(["mean", "std", "count"])
                    .reset_index()
                )

                def make_act_label(row):

                    label = (
                        f"activation_{int(row['inhibition_conc_uM'])}uM"
                        if row["inhibition_conc_uM"] > 0
                        else "activation"
                    )

                    if row["is_estimated"]:
                        label += "_est"
                    else:
                        label += "_exp"

                    return label

                df_grouped_act["tox_label"] = (
                    df_grouped_act.apply(make_act_label, axis=1)
                )

                df_pivot_act = df_grouped_act.pivot(
                    index="smiles",
                    columns="tox_label"
                )

                df_pivot_act = flatten_pivot_columns(
                    df_pivot_act,
                    assay_name
                )

                pivots.append(df_pivot_act)

            # -----------------------------
            # 2. Affinity block
            # -----------------------------
            if not df_aff.empty:

                df_grouped_aff = (
                    df_aff.groupby(["smiles", "tox_type"])["standard_value"]
                    .agg(["mean", "std", "count"])
                    .reset_index()
                )

                df_pivot_aff = df_grouped_aff.pivot(
                    index="smiles",
                    columns="tox_type"
                )

                df_pivot_aff = flatten_pivot_columns(
                    df_pivot_aff,
                    assay_name
                )

                pivots.append(df_pivot_aff)

            # -----------------------------
            # 3. Merge all tox blocks
            # -----------------------------
            if pivots:

                df_pivot = reduce(
                    lambda l, r: pd.merge(
                        l,
                        r,
                        on="smiles",
                        how="outer"
                    ),
                    pivots
                )

                dfs.append(df_pivot.reset_index())

            continue

        # -----------------------------
        # Type-A assays (ADME, CYPs)
        # -----------------------------
        group_cols = []

        def add_group_col(col):
            if col in df.columns and col not in group_cols:
                group_cols.append(col)

        add_group_col("species")
        add_group_col("permeability_system")

        if (
            "papp_direction" in df.columns
            and df["papp_direction"].nunique(dropna=True) > 1
        ):
            add_group_col("papp_direction")

        add_group_col("cyp_endpoint")

        add_group_col("bsep_endpoint")
        add_group_col("oatp_endpoint")

        if (
            "inhibition_conc_uM" in df.columns
            and df["inhibition_conc_uM"].notna().any()
        ):
            add_group_col("inhibition_conc_uM")

        if "cyp_substrate" in df.columns:
            if df["cyp_substrate"].nunique(dropna=True) > 1:
                add_group_col("cyp_substrate")

        if "bsep_substrate" in df.columns:
            if df["bsep_substrate"].nunique(dropna=True) > 1:
                add_group_col("bsep_substrate")
        # if "oatp_subtype" in df.columns:
        #     if df["oatp_subtype"].nunique() > 1:
        #         group_cols.append("oatp_subtype")

        if group_cols:

            group_cols = list(dict.fromkeys(group_cols))

            for col in group_cols:
                if df[col].dtype == object:
                    df[col] = df[col].fillna("unknown")

            logger.debug(f"{assay_name} raw rows entering grouping: {len(df)}")
            logger.debug(f"{assay_name} unique smiles before grouping: {df['smiles'].nunique()}")

            if "cyp_substrate" in df.columns:
                logger.debug(
                    f"{assay_name} substrate counts before grouping: "
                    f"{df['cyp_substrate'].value_counts(dropna=False).to_dict()}"
                )

            logger.debug(f"group_cols = {group_cols}")

            df_grouped = df.groupby(
                ["smiles"] + group_cols
            )["standard_value"].agg(
                ["mean","std","count"]
            ).reset_index()

            logger.debug(f"{assay_name} rows after grouping: {len(df_grouped)}")

            df_pivot = df_grouped.pivot(index="smiles", columns=group_cols)

            logger.debug(
                f"{assay_name} unique smiles after pivot: {len(df_pivot)}"
            )

            # Properly join MultiIndex columns without splitting letters
            # df_pivot.columns = [
            #     f"{assay_name}_{'_'.join(map(str, col[1:]))}_{col[0]}"
            #     for col in df_pivot.columns
            # ]
            df_pivot = flatten_pivot_columns(df_pivot, assay_name)
            dfs.append(df_pivot.reset_index())
            substr_cols = [
                c for c in df_pivot.columns
                if "midazolam" in c or "testosterone" in c
            ]

            logger.debug(
                f"Compounds with any substrate data: "
                f"{df_pivot[substr_cols].notna().any(axis=1).sum()}"
            )
            logger.debug(
                "Non-null counts per CYP3A4 column:\n"
                f"{df_pivot.notna().sum().sort_values(ascending=False)}"
            )
            continue

        # -----------------------------
        # Fallback for remaining assays
        # -----------------------------
        df_subset = df.groupby("smiles")["standard_value"].mean().reset_index()
        df_subset.rename(columns={"standard_value": assay_name}, inplace=True)
        dfs.append(df_subset)

    # =============================================================
    # Build the reported unscaled Clint endpoints once across all
    # retrieval buckets.
    # =============================================================

    if unscaled_clint_frames:

        unscaled_clint_all = pd.concat(
            unscaled_clint_frames,
            ignore_index=True,
        )

        # ---------------------------------------------------------
        # Remove the same ChEMBL activity if broad retrieval terms
        # placed it in both Microsomal Stability and Hepatocyte
        # Stability buckets.
        # ---------------------------------------------------------

        if "metstab_source_id" in unscaled_clint_all.columns:

            before_dedup = len(unscaled_clint_all)

            unscaled_clint_all = (
                unscaled_clint_all
                .drop_duplicates(
                    subset=[
                        "metstab_source_id",
                        "smiles",
                        "species",
                        "metstab_matrix",
                        "standard_units",
                    ],
                    keep="first",
                )
            )

            removed = (
                before_dedup - len(unscaled_clint_all)
            )

            if removed:
                logger.info(
                    "Removed "
                    f"{removed} duplicate unscaled Clint "
                    "activities retrieved through multiple "
                    "assay keyword buckets"
                )

        # ---------------------------------------------------------
        # Save the authoritative combined long-form table before
        # aggregation by SMILES.
        # ---------------------------------------------------------

        combined_long_path = os.path.join(
            out_dir,
            "reported_unscaled_clint_long.csv",
        )

        unscaled_clint_all.to_csv(
            combined_long_path,
            index=False,
        )

        logger.info(
            "Saved combined reported unscaled Clint rows to "
            f"{combined_long_path}"
        )

        # ---------------------------------------------------------
        # Convert canonical units to filenames/column-safe labels.
        #
        # Unit information remains explicit in the final task name.
        # ---------------------------------------------------------

        unit_label_map = {
            "uL/min/mg protein":
                "uL_min_mg_protein",
            "uL/min/10^6 cells":
                "uL_min_1e6_cells",
        }

        unscaled_clint_all["clint_unit_label"] = (
            unscaled_clint_all["standard_units"]
            .map(unit_label_map)
        )

        if unscaled_clint_all[
            "clint_unit_label"
        ].isna().any():

            bad_units = (
                unscaled_clint_all.loc[
                    unscaled_clint_all[
                        "clint_unit_label"
                    ].isna(),
                    "standard_units",
                ]
                .value_counts(dropna=False)
                .to_dict()
            )

            raise ValueError(
                "Unexpected unit reached the reported unscaled "
                f"Clint pivot: {bad_units}"
            )

        # ---------------------------------------------------------
        # Aggregate repeated reported activities for each compound.
        #
        # Species, matrix, and canonical unit all define the task.
        # ---------------------------------------------------------

        clint_group_cols = [
            "species",
            "metstab_matrix",
            "clint_unit_label",
        ]

        for col in clint_group_cols:
            unscaled_clint_all[col] = (
                unscaled_clint_all[col]
                .fillna("unknown")
            )

        unscaled_clint_grouped = (
            unscaled_clint_all
            .groupby(
                ["smiles"] + clint_group_cols,
                dropna=False,
            )["standard_value"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )

        unscaled_clint_pivot = (
            unscaled_clint_grouped
            .pivot(
                index="smiles",
                columns=clint_group_cols,
            )
        )

        unscaled_clint_pivot = flatten_pivot_columns(
            unscaled_clint_pivot,
            "Clint_unscaled_in_vitro",
        )

        dfs.append(
            unscaled_clint_pivot.reset_index()
        )

    if not dfs:
        return pd.DataFrame()

    # Merge all assay tables
    
    mtl_df = reduce(lambda l, r: pd.merge(l, r, on="smiles", how="outer"), dfs)
    
    duplicate_columns = (
        pd.Index(mtl_df.columns)
        [pd.Index(mtl_df.columns).duplicated()]
        .unique()
        .tolist()
    )

    if duplicate_columns:
        raise ValueError(
            "Duplicate multitask output columns were generated after "
            f"flattening/sanitisation: {duplicate_columns}"
        )
    
    mtl_df.to_csv(output_path, index=False)
    logger.debug(f"Saved multi-task dataset CSV to: {output_path}")

    return mtl_df

@register_task(
    "generate_diagnostics",
    category="ADME",
    description="Generate diagnostic plots and tables for ADME dataset"
)
def generate_diagnostics(config, mtl_df=None):
    """
    Generates:
    1. Counts of data points per assay
    2. Unit distributions per assay
    3. Pairwise assay overlap table and heatmap
    4. Network plot to show task connectivity
    5. t-SNE of Morgan fingerprints
    6. Scaffold frequency histogram and top scaffold grid
    
    Saves outputs in <output_dir>_plots
    """
    
    # Define endpoint groups using exhaustive YAML keywords
    
    # Helper to get columns containing any keyword (case-insensitive)
    assigned_cols = set()

    def get_group_cols(mtl_df, keywords):
        cols = []

        for c in mtl_df.columns:

            if c in assigned_cols:
                continue

            if (
                c.endswith("_std")
                or c.endswith("_count")
                or "inhibition_conc_uM" in c
                or "tox_type" in c
                or "is_estimated" in c
            ):
                continue

            for kw in keywords:
                if kw.lower() in c.lower():
                    cols.append(c)
                    assigned_cols.add(c)
                    break

        return cols
    
    endpoint_groups = {
        "PhysChem": get_group_cols(mtl_df, ["logp", "logd", "logS", "solubility"]),
        "Permeability": get_group_cols(mtl_df, [
            "caco", "mdck", "pampa", "p-gp", "bcrp", "mrp", "bbb",
            "blood brain barrier", "brain penetration", "logbb", "brain/plasma"
        ]),
        "Metabolism": get_group_cols(
            mtl_df,
            [
                "metstab",
                "microsomal",
                "microsome",
                "hepatocyte",
                "s9",
                "clint",
                "clearance",
                "scaled_L_h_kg",
                "liver_normalised",
                "t_half",
                "percent_remaining",
                "extraction_ratio",
                "cyp",
            ]
        ),
        "Tox": get_group_cols(mtl_df, ["herg", "nav1.5", "cav1.2", "5-ht1a", "5-ht1b", "5-ht2a", "5-ht2b", "5-ht2c", "pxr"]),
        "Transporters": get_group_cols(mtl_df, [
            "p-gp", "pgp",
            "bcrp",
            "mrp",
            "bsep",
            "oatp", "slco",
            "oct", "slc22",
            "mate"
        ]),
        "PK": get_group_cols(mtl_df, ["vd", "bioavailability", "f%", "fraction_unbound", "ppb"])
    }

    if mtl_df is None or mtl_df.empty:
        logger.warning("mtl_df is None or empty - skipping diagnostics")
        return {"plots_dir": None}

    output_dir = config.get("output", {}).get("directory", "outputs/adme")
    plots_dir = os.path.join(output_dir, "_plots")
    os.makedirs(plots_dir, exist_ok=True)

    task_cols = [
        c for c in mtl_df.columns
        if (
            not (c.endswith("_std") or c.endswith("_count"))
            and "inhibition_conc_uM" not in c
            and np.issubdtype(mtl_df[c].dtype, np.number)
            and not pd.api.types.is_bool_dtype(mtl_df[c])
            and "is_estimated" not in c
        )
    ]
    numeric_cols = task_cols.copy()
    mask = mtl_df[task_cols].notna()

    # ------------------------
    # Number of data points per assay (panelled by endpoint group)
    # ------------------------
    n_groups = len(endpoint_groups)

    ncols = 2
    nrows = math.ceil(n_groups / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(24, 6 * nrows)
    )

    axes = np.array(axes).flatten()

    for idx, (group_name, cols) in enumerate(endpoint_groups.items()):

        ax = axes[idx]

        if not cols:
            ax.text(
                0.5,
                0.5,
                "No endpoints found",
                ha='center',
                va='center',
                fontsize=12
            )
            ax.set_title(group_name)
            ax.axis("off")
            continue

        counts = (
            mtl_df[cols]
            .notna()
            .sum()
            .sort_values(ascending=False)
        )

        labels = (
            counts.index
            .str.replace("_mean", "", regex=False)
            .str.replace("_std", "", regex=False)
            .str.replace("_count", "", regex=False)
        )

        sns.barplot(
            x=labels,
            y=counts.values,
            ax=ax,
            palette="Blues_d",
            errorbar=None
        )

        # ---------------------------------
        # Dynamic label/font scaling
        # ---------------------------------

        n_labels = len(labels)

        if n_labels > 40:
            fontsize = 5
            rotation = 90

        elif n_labels > 25:
            fontsize = 6
            rotation = 75

        elif n_labels > 15:
            fontsize = 7
            rotation = 60

        else:
            fontsize = 9
            rotation = 45

        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=rotation,
            ha="right",
            fontsize=fontsize
        )

        ax.set_ylabel("Number of Data Points")
        ax.set_title(f"{group_name} ({n_labels} assays)")

    # hide unused panels
    for j in range(idx + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    plt.savefig(
        os.path.join(plots_dir, "assay_data_counts_grid.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    # ------------------------
    # Pairwise overlap
    # ------------------------

    overlap_counts = pd.DataFrame(index=task_cols, columns=task_cols, dtype=int)
    overlap_frac = pd.DataFrame(index=task_cols, columns=task_cols, dtype=float)

    for t1 in task_cols:
        for t2 in task_cols:
            both = (mask[t1] & mask[t2]).sum()
            overlap_counts.loc[t1, t2] = both
            min_count = min(mask[t1].sum(), mask[t2].sum())
            overlap_frac.loc[t1, t2] = both / min_count if min_count > 0 else np.nan

    # Save raw counts and fraction tables
    overlap_counts.to_csv(os.path.join(plots_dir, "pairwise_overlap_counts.csv"))
    overlap_frac.to_csv(os.path.join(plots_dir, "pairwise_overlap_fraction.csv"))

    # ------------------------
    # Clean labels for plotting
    # ------------------------
    plot_labels = [c.replace("_mean", "") for c in task_cols]
    overlap_plot = overlap_frac.copy()
    overlap_plot.index = plot_labels
    overlap_plot.columns = plot_labels

    # Make NaN-safe for clustering
    overlap_plot_finite = overlap_plot.dropna(axis=0, how="all").dropna(axis=1, how="all").fillna(0)
    n_tasks = len(overlap_plot_finite)

    if n_tasks == 0:
        logger.warning("No valid tasks for overlap plotting after removing all-NaN rows/columns")
    else:
        if n_tasks < 10:
            figsize = max(6, n_tasks * 0.8)
            annot = True
            annot_fontsize = 12
            g = sns.heatmap(
                overlap_plot_finite.astype(float),
                cmap="viridis",
                annot=annot,
                fmt=".2f",
                cbar_kws={"label": "Overlap fraction"},
                linewidths=0.5
            )
            plt.xticks(rotation=45, ha="right", fontsize=annot_fontsize)
            plt.yticks(rotation=0, fontsize=annot_fontsize)
            plt.title("Pairwise Assay Overlap Fraction", fontsize=14)
            plt.tight_layout()
        elif n_tasks <= 35:
            cell_size = 0.45
            figsize = max(8, n_tasks * cell_size + 3)
            annot = True
            annot_fontsize = max(6, 12 - n_tasks * 0.15)
            g = sns.clustermap(
                overlap_plot_finite.astype(float),
                cmap="viridis",
                annot=annot,
                fmt=".2f",
                figsize=(figsize, figsize),
                annot_kws={"size": annot_fontsize},
                cbar_kws={"label": "Overlap fraction"},
                dendrogram_ratio=(0.15, 0.15)
            )
            plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
            plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.suptitle("Pairwise Assay Overlap Fraction (clustered)", y=1.02)
        else:
            cell_size = 0.45
            figsize = max(12, n_tasks * cell_size + 3)
            g = sns.clustermap(
                overlap_plot_finite.astype(float),
                cmap="viridis",
                annot=False,
                figsize=(figsize, figsize),
                cbar_kws={"label": "Overlap fraction"},
                dendrogram_ratio=(0.1, 0.1)
            )
            plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
            plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.suptitle("Pairwise Assay Overlap Fraction (clustered)", y=1.02)

        plt.savefig(os.path.join(plots_dir, "pairwise_overlap_heatmap.png"), dpi=300, bbox_inches="tight")
        plt.close()

    # ------------------------
    # Task–Task correlation heatmap
    # ------------------------
    # Keep only numeric columns, exclude *_std and *_count
    numeric_cols = [
        c for c in mtl_df.columns
        if (
            not (c.endswith("_std") or c.endswith("_count"))
            and "inhibition_conc_uM" not in c
            and np.issubdtype(mtl_df[c].dtype, np.number)
            and not pd.api.types.is_bool_dtype(mtl_df[c])   # ← ADD
            and "is_estimated" not in c                     # ← ADD
        )
    ]
    corr_df = mtl_df[numeric_cols].corr()

    corr_df_plot = corr_df.copy()
    plot_labels = [c.replace("_mean", "") for c in numeric_cols]
    corr_df_plot.index = plot_labels
    corr_df_plot.columns = plot_labels

    # Make NaN-safe
    corr_df_plot_finite = corr_df_plot.fillna(0)

    plt.figure(figsize=(12,10))
    sns.heatmap(corr_df_plot_finite, cmap="PiYG", center=0, annot=False)
    plt.title("Task–Task Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "task_correlation_heatmap.png"), dpi=300)
    plt.close()

    # ------------------------
    # Correlation vs Overlap Scatter
    # ------------------------
    max_overlap = overlap_frac.where(~np.eye(len(overlap_frac), dtype=bool)).max(axis=0)

    cor_vs_overlap = pd.DataFrame({
        "task": numeric_cols,
        "max_overlap": [max_overlap.get(c, np.nan) for c in numeric_cols],
        "mean_correlation": [corr_df[c].drop(c).mean() for c in numeric_cols]
    })

    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=cor_vs_overlap,
        x="max_overlap",
        y="mean_correlation"
    )
    plt.xlabel("Maximum Fraction Overlap with Any Task")
    plt.ylabel("Mean Correlation with Other Tasks")
    plt.title("Correlation vs Overlap per Task")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "correlation_vs_overlap.png"), dpi=300)
    plt.close()

    # ------------------------
    # Task connectivity graph
    # ------------------------
    G = nx.Graph()

    for task in task_cols:
        G.add_node(task)

    # Only add edges with finite overlap
    for t1 in task_cols:
        for t2 in task_cols:
            if t1 == t2:
                continue
            weight = overlap_frac.loc[t1, t2]
            if pd.notna(weight) and np.isfinite(weight) and weight > 0.05:
                G.add_edge(t1, t2, weight=weight)

    plt.figure(figsize=(10,8))
    pos = nx.spring_layout(G, seed=42)
    edges = G.edges(data=True)
    weights = [d["weight"]*5 for (_,_,d) in edges]

    nx.draw_networkx_nodes(G, pos, node_size=1500, node_color="lightblue")
    nx.draw_networkx_labels(G, pos, font_size=9)
    nx.draw_networkx_edges(G, pos, width=weights, alpha=0.6)

    plt.title("Task Connectivity Graph (shared compound overlap)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "task_connectivity_graph.png"), dpi=300)
    plt.close()

    # TSNE of all fingerprints
    mols = [
        Chem.MolFromSmiles(sm)
        for sm in tqdm(mtl_df.smiles, desc="Parsing SMILES", unit="mol")
    ]
    fps = [
        rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius=2, nBits=1024)
        if m is not None else None
        for m in tqdm(mols, desc="Generating fingerprints", unit="mol")
    ]

    def fp_to_array(fp):
        arr = np.zeros((fp.GetNumBits(),), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    fp_arrays = np.array([fp_to_array(fp) for fp in fps if fp is not None])
    tsne = TSNE(n_components=2, random_state=42, perplexity=5)
    fp_tsne = tsne.fit_transform(fp_arrays)

    plt.figure(figsize=(8,6))
    sns.scatterplot(x=fp_tsne[:,0], y=fp_tsne[:,1])
    plt.title("t-SNE of Molecular Fingerprints")
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "tsne_fingerprints.png"), dpi=300)
    plt.close()

    # t-SNE grid with highlighted endpoints
    # Debug: print detected columns for sanity check (trim suffixes)
    for group_name, cols in endpoint_groups.items():
        clean_cols = [c.replace("_mean","").replace("_std","").replace("_count","") for c in cols]
        logger.debug(f"{group_name}: {clean_cols}")

    # Create 6-panel grid
    fig, axes = plt.subplots(3, 2, figsize=(18,18))
    axes = axes.flatten()

    for idx, (group_name, cols) in enumerate(endpoint_groups.items()):
        ax = axes[idx]

        if not cols:
            ax.scatter(fp_tsne[:,0], fp_tsne[:,1], color="lightgrey", s=30, alpha=0.5)
            ax.set_title(f"t-SNE colored by {group_name} (no endpoints found)")
            ax.set_xlabel("t-SNE 1")
            ax.set_ylabel("t-SNE 2")
            continue

        # Determine which compounds have at least one value present in this group
        present_mask = mtl_df[cols].notna().any(axis=1)

        # Plot highlighted points
        ax.scatter(
            fp_tsne[~present_mask,0], 
            fp_tsne[~present_mask,1], 
            color="lightgrey", s=15, alpha=0.1, label="Absent"
        )
        ax.scatter(
            fp_tsne[present_mask,0], 
            fp_tsne[present_mask,1], 
            color="dodgerblue", s=20, alpha=0.35, label="Present"
        )

        ax.set_title(f"t-SNE colored by {group_name}")
        ax.set_xlabel("t-SNE 1")
        ax.set_ylabel("t-SNE 2")
        ax.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "tsne_grid_highlighted.png"), dpi=300)
    plt.close()

    # Histogram of scaffold frequency

    scaffolds = [
        MurckoScaffold.GetScaffoldForMol(m)
        for m in tqdm(mols, desc="Extracting scaffolds", unit="mol")
        if m is not None
    ]
    scaffold_smiles = [Chem.MolToSmiles(s) for s in scaffolds if s is not None]
    scaffold_counts = Counter(scaffold_smiles)

    plt.figure(figsize=(8,5))
    plt.hist(list(scaffold_counts.values()), bins=30, log=True)
    plt.xlabel("Scaffold Frequency")
    plt.ylabel("Number of Scaffolds (log scale)")
    plt.title("Distribution of Scaffold Frequencies")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "scaffold_frequency.png"), dpi=300)
    plt.close()

    # Top 10 scaffolds grid
    
    top_n = 10
    top_scaffolds = scaffold_counts.most_common(top_n)
    scaffold_mols = [Chem.MolFromSmiles(smi) for smi, _ in top_scaffolds]
    scaffold_labels = [f"{count} molecules" for _, count in top_scaffolds]

    mol_size = (200, 200)
    images = [Draw.MolToImage(mol, size=mol_size) for mol in scaffold_mols]

    def add_label(img, label, size=(200, 220)):
        new_img = Image.new("RGB", size, "white")
        new_img.paste(img, (0, 0))
        draw = ImageDraw.Draw(new_img)
        draw.text((10, 200), label, fill="black")
        return new_img

    labeled_images = [add_label(img, label) for img, label in zip(images, scaffold_labels)]
    cols = 5
    rows = (len(labeled_images) + cols - 1) // cols
    grid_width = mol_size[0] * cols
    grid_height = 220 * rows
    grid_img = Image.new("RGB", (grid_width, grid_height), "white")

    for idx, img in enumerate(labeled_images):
        x = (idx % cols) * mol_size[0]
        y = (idx // cols) * 220
        grid_img.paste(img, (x, y))

    grid_img.save(os.path.join(plots_dir, "top_10_scaffolds.png"))

    logger.info(f"Plots saved to: {plots_dir}")

    return {"df": mtl_df, "plots_dir": plots_dir}