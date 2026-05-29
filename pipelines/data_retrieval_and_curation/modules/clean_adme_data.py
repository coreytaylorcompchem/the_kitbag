import os

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
    convert_met_stab,
    convert_ppb,
    convert_solubility,
    convert_logp_logd,
    convert_cyp_activity,
    convert_vd,
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
    classify_metstab_record,
    convert_bsep_activity,
    classify_oatp_record,
    convert_oatp_activity,
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

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
            for r in records:
                
                if not classify_metstab_record(r): # new met stab classifier
                    continue

                val, unit = r.get("standard_value"), r.get("standard_units")
                new_val, new_unit = convert_met_stab(val, unit)
                if new_val is not None:
                    r["standard_value"], r["standard_units"] = new_val, new_unit
                    r["species"] = extract_species(r.get("target_organism") or r.get("assay_description"))
                    system = detect_metstab_system(r.get("assay_description"))
                    if system == "in_vitro":
                        r["metstab_system"] = system
                        new_records.append(r)

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
    allowed_species = {"Human", "Mouse", "Rat", "Monkey", "unknown"}

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

        # -----------------------------
        # Type-A assays (ADME, CYPs)
        # -----------------------------
        group_cols = []
        if "species" in df.columns:
            group_cols.append("species")
        if "permeability_system" in df.columns:
            group_cols.append("permeability_system")

        if (
            "papp_direction" in df.columns
            and df["papp_direction"].nunique() > 1
        ):
            group_cols.append("papp_direction")
        if "cyp_endpoint" in df.columns:
            group_cols.append("cyp_endpoint")
        if "bsep_endpoint" in df.columns:
            group_cols.append("bsep_endpoint")
        if "oatp_endpoint" in df.columns:
            group_cols.append("oatp_endpoint")
        if (
            "inhibition_conc_uM" in df.columns
            and df["inhibition_conc_uM"].notna().any()
        ):
            group_cols.append("inhibition_conc_uM")
        if "cyp_substrate" in df.columns:
            if df["cyp_substrate"].nunique() > 1:
                group_cols.append("cyp_substrate")
        if "bsep_substrate" in df.columns:
            if df["bsep_substrate"].nunique() > 1:
                group_cols.append("bsep_substrate")
        # if "oatp_subtype" in df.columns:
        #     if df["oatp_subtype"].nunique() > 1:
        #         group_cols.append("oatp_subtype")

        if group_cols:

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
        # Type-B assays (tox endpoints: hERG, Nav, etc.)
        # -----------------------------
        if "tox_type" in df.columns:

            logger.debug(f"Building pivot for tox assay: {assay_name}")

            # --- SPLIT DATA ---
            df_inh = df[df["tox_type"] == "inhibition"].copy()
            df_aff = df[df["tox_type"].isin(["IC50", "EC50", "KI", "KD"])].copy()

            pivots = []

            # -----------------------------
            # 1. Inhibition block (measured + estimated)
            # -----------------------------
            if not df_inh.empty:

                df_inh["inhibition_conc_uM"] = df_inh.get("inhibition_conc_uM", -1).fillna(-1)

                if "is_estimated" not in df_inh.columns:
                    df_inh["is_estimated"] = False

                df_grouped_inh = (
                    df_inh.groupby(["smiles", "inhibition_conc_uM", "is_estimated"])["standard_value"]
                    .agg(["mean", "std", "count"])
                    .reset_index()
                )

                def make_inh_label(row):
                    label = f"inhibition_{int(row['inhibition_conc_uM'])}uM" if row["inhibition_conc_uM"] > 0 else "inhibition"

                    if row["is_estimated"]:
                        label += "_est"
                    else:
                        label += "_exp"

                    return label

                df_grouped_inh["tox_label"] = df_grouped_inh.apply(make_inh_label, axis=1)

                df_pivot_inh = df_grouped_inh.pivot(index="smiles", columns="tox_label")
                # df_pivot_inh.columns = [
                #     f"{assay_name}_{tox}_{stat}"
                #     for stat, tox in df_pivot_inh.columns
                # ]
                df_pivot_inh = flatten_pivot_columns(df_pivot_inh, assay_name)

                pivots.append(df_pivot_inh)

            # -----------------------------
            # 2. Affinity block (IC50, EC50, etc.)
            # -----------------------------
            if not df_aff.empty:

                df_grouped_aff = (
                    df_aff.groupby(["smiles", "tox_type"])["standard_value"]
                    .agg(["mean", "std", "count"])
                    .reset_index()
                )

                df_pivot_aff = df_grouped_aff.pivot(index="smiles", columns="tox_type")

                # df_pivot_aff.columns = [
                #     f"{assay_name}_{tox}_{stat}"
                #     for stat, tox in df_pivot_aff.columns
                # ]
                df_pivot_aff = flatten_pivot_columns(df_pivot_aff, assay_name)

                pivots.append(df_pivot_aff)

            # -----------------------------
            # 3. Merge both (if both exist)
            # -----------------------------
            if pivots:
                df_pivot = reduce(lambda l, r: pd.merge(l, r, on="smiles", how="outer"), pivots)
                dfs.append(df_pivot.reset_index())

            continue

        # -----------------------------
        # Fallback for remaining assays
        # -----------------------------
        df_subset = df.groupby("smiles")["standard_value"].mean().reset_index()
        df_subset.rename(columns={"standard_value": assay_name}, inplace=True)
        dfs.append(df_subset)

    if not dfs:
        return pd.DataFrame()

    # Merge all assay tables
    
    mtl_df = reduce(lambda l, r: pd.merge(l, r, on="smiles", how="outer"), dfs)
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
        "Metabolism": get_group_cols(mtl_df, ["metstab", "microsomal", "hepatocyte", "cyp"]),
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
            and not pd.api.types.is_bool_dtype(mtl_df[c])   # ← ADD
            and "is_estimated" not in c                     # ← ADD
        )
    ]
    numeric_cols = task_cols.copy()
    mask = mtl_df[task_cols].notna()

    # ------------------------
    # Number of data points per assay (panelled by endpoint group)
    # ------------------------
    fig, axes = plt.subplots(3, 2, figsize=(18, 18))
    axes = axes.flatten()

    for idx, (group_name, cols) in enumerate(endpoint_groups.items()):
        ax = axes[idx]

        if not cols:
            ax.text(0.5, 0.5, "No endpoints found", ha='center', va='center', fontsize=12)
            ax.set_title(f"{group_name}")
            ax.axis("off")
            continue

        # Count non-NaN values per assay
        counts = mtl_df[cols].notna().sum().sort_values(ascending=False)

        sns.barplot(
            x=counts.index.str.replace("_mean","").str.replace("_std","").str.replace("_count",""),
            y=counts.values,
            ax=ax,
            palette="Blues_d",
            errorbar=None
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
        ax.set_ylabel("Number of Data Points")
        ax.set_title(f"{group_name}")

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "assay_data_counts_grid.png"), dpi=300)
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