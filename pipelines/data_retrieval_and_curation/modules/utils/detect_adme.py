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

def extract_species(text):
    if not text:
        return "unknown"
    t = str(text).lower()
    species_map = {
        "human": "Human", "homo sapiens": "Human",
        "mouse": "Mouse", "mus musculus": "Mouse",
        "rat": "Rat", "rattus": "Rat",
        "monkey": "Monkey", "macaca": "Monkey", "cyno": "Monkey"
    }
    for k, v in species_map.items():
        if k in t:
            return v
    return "unknown"