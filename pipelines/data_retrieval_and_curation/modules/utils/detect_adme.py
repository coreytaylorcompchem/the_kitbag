import re
import numpy as np

def build_dual_endpoint_dataset(df, assay_name, endpoint_col, value_cols):
    """
    Generic handler for assays with two endpoint types (e.g. IC50 and % inhibition).

    Parameters
    ----------
    df : dataframe
    assay_name : str
    endpoint_col : column indicating endpoint type (e.g. IC50 / inhibition)
    value_cols : dict mapping endpoint → value column

    Returns
    -------
    pivoted dataframe ready to merge into multitask dataset
    """

    df = df.dropna(subset=[endpoint_col])

    # choose correct value column based on endpoint
    df["endpoint_value"] = np.nan

    for endpoint, col in value_cols.items():
        df.loc[df[endpoint_col] == endpoint, "endpoint_value"] = df[col]

    df = df.dropna(subset=["endpoint_value"])

    grouped = (
        df.groupby(["smiles", endpoint_col])["endpoint_value"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )

    pivot = grouped.pivot(index="smiles", columns=endpoint_col)

    pivot.columns = [
        f"{assay_name}_{endpoint}_{stat}"
        for stat, endpoint in pivot.columns
    ]

    pivot = pivot.reset_index()

    return pivot

def detect_metstab_system(description):
    if not description:
        return "unknown"

    desc = description.lower()

    in_vitro_patterns = [
        "microsome",
        "microsomal",
        "hepatocyte",
        "hepatocytes",
        "liver s9",
        "s9 fraction",
        "intrinsic clearance"
    ]

    in_vivo_patterns = [
        "plasma clearance",
        "blood clearance",
        "pk",
        "pharmacokinetic",
        "half life",
        "t1/2",
        "in vivo",
        "iv clearance"
    ]

    for p in in_vitro_patterns:
        if p in desc:
            return "in_vitro"

    for p in in_vivo_patterns:
        if p in desc:
            return "in_vivo"

    return "unknown"

def detect_papp_direction(description):
        if not description:
            return None
        desc = description.lower().replace("→", "->").replace("–", "-").replace("—", "-")
        ab_patterns = [r"a\s*->\s*b", r"apical\s*(?:to|->|-)\s*basolateral", r"\bab\b"]
        ba_patterns = [r"b\s*->\s*a", r"basolateral\s*(?:to|->|-)\s*apical", r"\bba\b"]
        for p in ab_patterns:
            if re.search(p, desc):
                return "AB"
        for p in ba_patterns:
            if re.search(p, desc):
                return "BA"
        return None

def extract_species(text, context="generic"):

    if not text:
        return "unknown"

    t = str(text).lower()

    species_patterns = {
        "Human": [
            "human",
            "homo sapiens",
        ],
        "Mouse": [
            "mouse",
            "mice",
            "mus musculus",
            "cd-1",
            "c57bl/6",
        ],
        "Rat": [
            "rat",
            "rats",
            "rattus",
            "sprague-dawley",
            "wistar",
        ],
        "Monkey": [
            "monkey",
            "macaca",
            "cyno",
            "cynomolgus",
            "rhesus",
            "macaque",
        ],
        "Dog": [
            "dog",
            "dogs",
            "beagle",
        ]
    }

    # -----------------------------
    # PK-specific exclusions
    # -----------------------------
    if context == "pk":

        exclude_terms = [
            "microsome",
            "microsomal",
            "hepatocyte",
            "hepatocytes",
            "s9",
            "cytosol",
            "recombinant",
            "transfected",
            "cell line",
            "well stirred",
            "well-stirred",
            "predicted",
            "prediction",
            "simulated",
            "simulation",
            "model",
            "in vitro",
            "intrinsic clearance",
            "clint",
        ]

        if any(term in t for term in exclude_terms):
            return "unknown"

    # -----------------------------
    # Species matching
    # -----------------------------
    found = []

    for species, patterns in species_patterns.items():

        for pattern in patterns:

            if pattern in t:
                found.append(species)
                break

    found = list(set(found))

    if len(found) == 1:
        return found[0]

    return "unknown"

def classify_permeability_system(record):
    """
    Classify permeability assay system.

    Returns
    -------
    tuple:
        (system, directionality)

    system:
        - "caco2"
        - "mdck"
        - "pampa"
        - None

    directionality:
        - "directional"
        - "non_directional"
    """

    text = " ".join([
        str(record.get("assay_description", "")),
        str(record.get("assay_type", "")),
        str(record.get("standard_type", "")),
        str(record.get("document_chembl_id", "")),
    ]).lower()

    # -------------------------
    # PAMPA
    # -------------------------
    if any(x in text for x in [
        "pampa",
        "parallel artificial membrane",
        "artificial membrane permeability",
    ]):
        return "pampa", "non_directional"

    # -------------------------
    # Caco-2
    # -------------------------
    if any(x in text for x in [
        "caco",
        "caco-2",
        "caco2",
    ]):
        return "caco2", "directional"

    # -------------------------
    # MDCK
    # -------------------------
    if any(x in text for x in [
        "mdck",
        "mdcki",
        "mdr1-mdck",
        "mdckii",
    ]):
        return "mdck", "directional"

    return None, None

def flatten_pivot_columns(df_pivot, assay_name):
    """
    Flatten MultiIndex columns from grouped assay pivots.

    Expected MultiIndex structure:
        (stat, group1, group2, ...)

    Produces:
        assay_group1_group2_stat

    Example:
        ("mean", "caco2", "AB")
    ->
        Caco-2_caco2_AB_mean
    """

    new_cols = []

    for col in df_pivot.columns:

        # Non-MultiIndex fallback
        if not isinstance(col, tuple):
            new_cols.append(f"{assay_name}_{col}")
            continue

        stat = col[0]
        groups = [str(x) for x in col[1:] if x not in [None, "", "nan"]]

        label = "_".join([assay_name] + groups + [stat])

        # cleanup accidental duplicates
        label = label.replace("__", "_")

        new_cols.append(label)

    df_pivot.columns = new_cols

    return df_pivot

def classify_bioavailability_record(record):

    text = " ".join([
        str(record.get("assay_description", "")),
        str(record.get("standard_type", "")),
        str(record.get("assay_type", "")),
    ]).lower()

    exclude_terms = [
        "microsome",
        "microsomal",
        "hepatocyte",
        "well stirred",
        "well-stirred",
        "predicted",
        "prediction",
        "simulated",
        "calculated",
        "model",
        "in vitro",
        "intrinsic clearance",
        "clint",
        "hepatic extraction",
    ]

    if any(t in text for t in exclude_terms):
        return False

    include_terms = [
        "oral bioavailability",
        "absolute bioavailability",
        "bioavailability",
        "f%",
        "fraction absorbed",
    ]

    if not any(t in text for t in include_terms):
        return False

    return True