import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors

def convert_cyp_activity(value, unit):
    """
    Standardize CYP inhibition endpoints:
      - IC50 → always return in nM
      - Percent inhibition (%) retained as-is
    Returns (numeric_value, endpoint_type, standard_units)
    """
    if value is None or unit is None:
        return None, None, None

    try:
        val = float(value)
    except:
        return None, None, None

    u = str(unit).lower().replace(" ", "").replace("µ", "u").strip()

    # Percent inhibition
    if "%" in u:
        return val, "inhibition", "%"

    # IC50 numeric
    factor = None
    if "nm" in u:
        factor = 1.0
    elif "um" in u or "μm" in u:
        factor = 1e3  # convert μM -> nM
    elif "m" == u:
        factor = 1e9  # M -> nM

    if factor is not None:
        return val * factor, "IC50", "nM"

    # Unhandled unit
    return None, None, None

def convert_vd(value, unit):
    """
    Convert volume of distribution to canonical L/kg.
    Returns (float_value, canonical_unit) or (None, None) if unconvertible.
    """
    if value is None:
        return None, None

    try:
        val = float(value)
    except (TypeError, ValueError):
        return None, None

    if not unit:
        return val, "L/kg"  # assume default if missing

    u = str(unit).lower().replace(" ", "").strip()

    # Conversion map
    conversions = {
        "l/kg": 1.0,
        "ml/kg": 1e-3,   # mL/kg → L/kg
        "l": 1.0,        # assume per kg if no further info
        "ml": 1e-3,
    }

    for k, factor in conversions.items():
        if k == u:
            return val * factor, "L/kg"

    # Unknown unit
    return val, u

def convert_logp_logd(value, unit):
    """
    Standardise LogP / LogD values (dimensionless).
    Returns float or None if invalid.
    """

    if value is None:
        return None

    try:
        val = float(value)
    except:
        return None

    if unit is None:
        return val

    u = str(unit).lower().strip()

    # Allowed annotations
    allowed = [
        "",
        "none",
        "unitless",
        "ph7.4",
        "ph 7.4"
    ]

    if u in allowed:
        return val

    # Reject incompatible units
    reject_patterns = [
        "%",
        "nm",
        "mm",
        "ug",
        "mol",
        "min",
        "hr",
        "c.p.m",
        "sigma",
        "ml",
        "l-1"
    ]

    if any(p in u for p in reject_patterns):
        return None

    # fallback: allow but treat as unitless
    return val

def convert_permeability(value, unit):
    """
    Convert various permeability units to canonical 10^-6 cm/s.
    Returns (float_value, canonical_unit) or (None, None) if unconvertible.
    """
    if value is None or unit is None:
        return None, None

    u = unit.lower().replace(" ", "").replace("'", "").replace("*", "").replace("^", "-").strip()

    # Canonical unit
    canonical = "10^-6 cm/s"

    # Map known variants directly
    direct_map = {
        "10-6cm/s": canonical,
        "10^-6cm/s": canonical,
        "10e-6cm/s": canonical,
        "10-6cms": canonical,
        "10^-6cms": canonical,
        "10e-6cms": canonical,
        "ucm/s": canonical,
        "papp10e6cm/s": canonical,
        "10^6cm/s": canonical,
        "10-7cm/s": canonical,   # convert to 0.1 * 10^-6 cm/s
        "10^-7cm/s": canonical,  # same here
        "10-8cm/s": canonical,   # 0.01 * 10^-6 cm/s
        "10^-8cm/s": canonical,
    }

    if u in direct_map:
        factor = 1.0
        if u in ["10-7cm/s", "10^-7cm/s"]:
            factor = 0.1
        elif u in ["10-8cm/s", "10^-8cm/s"]:
            factor = 0.01
        return float(value) * factor, canonical

    # Convert cm/s, nm/s, um/s
    if "cms" in u:
        return float(value) * 1e6, canonical
    if "nm/s" in u:
        return float(value) * 1e-1, canonical  # 1 nm/s = 0.1 * 10^-6 cm/s
    if "um/s" in u or "μm/s" in u:
        return float(value) * 1e-2, canonical
    if "cmhr-1" in u:
        return float(value) * 2.778e-6, canonical  # cm/hr → cm/s × 10^-6

    # Reject obviously incompatible units
    reject_patterns = ["%", "pmol", "nmol", "min", "sec", "umol", "ug", "a.u.", "uL/hr/cm2"]
    if any(r in u for r in reject_patterns):
        return None, None

    # Fallback: unhandled
    return None, None

def convert_solubility(value, unit, smiles):
    if value is None or unit is None:
        return None

    try:
        value = float(value)
    except:
        return None

    unit = unit.lower().strip()

    # nM -> mol/L
    if unit == "nm":
        mol_l = value * 1e-9

    # ug/mL -> mol/L (requires MW)
    elif unit in ["ug/ml", "ug.mL-1", "µg/ml"]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mw = Descriptors.MolWt(mol)
        mol_l = (value / mw) * 1e-3

    else:
        return None

    # Convert to logS
    if mol_l <= 0:
        return None

    return np.log10(mol_l)

def normalise_unit_string(unit):
    """
    Normalize common CHEMBL unit spellings.
    """
    if unit is None:
        return ""

    u = str(unit).lower().strip()

    u = (
        u
        .replace("µ", "u")
        .replace("μ", "u")
        .replace("micro", "u")
        .replace(" ", "")
        .replace("*", "")
        .replace("^", "-")
        .replace("(", "")
        .replace(")", "")
    )

    return u


def infer_metstab_matrix(record):
    """
    Infer whether a metabolic stability record is microsomal or hepatocyte.

    Returns:
      - "microsome"
      - "hepatocyte"
      - None
    """
    text = " ".join(
        str(record.get(k, "") or "")
        for k in [
            "assay_description",
            "standard_units",
            "standard_type",
            "assay_type",
            "target_organism",
        ]
    ).lower()

    if any(x in text for x in [
        "hepatocyte",
        "hepatocytes",
        "hepatic cell",
        "primary hepatocyte",
    ]):
        return "hepatocyte"

    if any(x in text for x in [
        "microsome",
        "microsomes",
        "microsomal",
        "hlm",
        "rlm",
        "mlm",
        "liver microsome",
        "liver microsomes",
    ]):
        return "microsome"

    # Unit-based fallback
    u = normalise_unit_string(record.get("standard_units"))

    if "10-6cells" in u or "1e6cells" in u or "106cells" in u or "millioncells" in u:
        return "hepatocyte"

    if "mgprotein" in u or "mgmicrosomalprotein" in u or "mg/ml" in u:
        return "microsome"

    return None


def convert_met_stab(value, unit):
    """
    Convert metabolic stability / clearance values.

    Returns
    -------
    tuple:
        new_value, canonical_unit, scale_state

    scale_state:
        - "unscaled"                    for uL/min/mg protein or uL/min/10^6 cells
        - "liver_normalised"            for mL/min/g liver
        - "scaled"                      for L/h/kg or converted mL/min/kg
        - None                          if unconvertible
    """

    if value is None or unit is None:
        return None, None, None

    try:
        val = float(value)
    except (TypeError, ValueError):
        return None, None, None

    u = normalise_unit_string(unit)

    # -----------------------------
    # Already body-weight-scaled Clint / clearance
    # -----------------------------
    scaled_l_h_kg_patterns = [
        "l/h/kg",
        "l/hr/kg",
        "lh-1kg-1",
        "lhr-1kg-1",
        "l.kg-1.h-1",
        "lkg-1h-1",
        "l/min/kg",
        "lmin-1kg-1",
    ]

    if any(p in u for p in scaled_l_h_kg_patterns):
        # If unit is L/min/kg, convert to L/h/kg.
        if "min" in u:
            return val * 60.0, "L/h/kg", "scaled"

        return val, "L/h/kg", "scaled"

    # mL/min/kg -> L/h/kg
    ml_min_kg_patterns = [
        "ml/min/kg",
        "mlmin-1kg-1",
        "ml.min-1.kg-1",
        "mlkg-1min-1",
        "ml/kg/min",
    ]

    if any(p in u for p in ml_min_kg_patterns):
        return val * 60.0 / 1000.0, "L/h/kg", "scaled"

    # -----------------------------
    # Liver-normalised clearance
    # mL/min/g liver
    # Needs liver weight to become L/h/kg.
    # -----------------------------
    ml_min_g_liver_patterns = [
        "ml/min/g",
        "ml/min/gliver",
        "mlmin-1g-1",
        "ml.min-1.g-1",
        "mlg-1min-1",
        "ml/g/min",
    ]

    if any(p in u for p in ml_min_g_liver_patterns):
        return val, "mL/min/g liver", "liver_normalised"

    # -----------------------------
    # Microsomal intrinsic clearance
    # uL/min/mg protein
    # -----------------------------
    microsome_ul_patterns = [
        "ul/min/mg",
        "ul/min/mgprotein",
        "ul/min/mgmicrosomalprotein",
        "ul/min/mgmicrosome",
        "ulmin-1mg-1",
        "ul.min-1.mg-1",
        "ulmin-1mgprotein-1",
        "ulmin-1mg-1protein",
        "ulmin-1mgmicrosomalprotein-1",
        "ulmin-1mgmicrosome-1",
    ]

    if any(p in u for p in microsome_ul_patterns):
        return val, "uL/min/mg protein", "unscaled"

    # microL/min/mg gets normalized to uL/min/mg above,
    # but keep this explicit for safety.
    if "ul" in u and "min" in u and "mg" in u:
        return val, "uL/min/mg protein", "unscaled"

    # l/min/mg -> uL/min/mg
    # Rare and maybe suspicious, but convert literally.
    l_min_mg_patterns = [
        "l/min/mg",
        "lmin-1mg-1",
        "l.min-1.mg-1",
    ]

    if any(p in u for p in l_min_mg_patterns):
        return val * 1_000_000.0, "uL/min/mg protein", "unscaled"

    # mL/min/mg -> uL/min/mg
    ml_min_mg_patterns = [
        "ml/min/mg",
        "mlmin-1mg-1",
        "ml.min-1.mg-1",
    ]

    if any(p in u for p in ml_min_mg_patterns):
        return val * 1000.0, "uL/min/mg protein", "unscaled"

    # -----------------------------
    # Hepatocyte intrinsic clearance
    # uL/min/10^6 cells
    # -----------------------------
    hepatocyte_patterns = [
        "ul/min/10-6cells",
        "ul/min/106cells",
        "ul/min/1e6cells",
        "ul/min/10e6cells",
        "ulmin-110-6cells-1",
        "ul.min-1.10-6cells-1",
        "ul.min-1.106cells-1",
        "ul.min-1.1e6cells-1",
        "ul/min/millioncells",
        "ulmin-1millioncells-1",
        "ul/106cells/min",
        "ul/1e6cells/min",
    ]

    if any(p in u for p in hepatocyte_patterns):
        return val, "uL/min/10^6 cells", "unscaled"

    if "ul" in u and "min" in u and "cell" in u:
        return val, "uL/min/10^6 cells", "unscaled"

    # -----------------------------
    # Ambiguous uL/min
    # Do not handle here because matrix context is needed.
    # The harmonisation block can resolve this using microsome/hepatocyte text.
    # -----------------------------
    ambiguous_ul_min_patterns = [
        "ul/min",
        "ulmin-1",
        "ul.min-1",
    ]

    if u in ambiguous_ul_min_patterns:
        return None, None, None

    # -----------------------------
    # Half-life / stability / percent
    # -----------------------------
    # hr or min are usually t1/2, not Clint. Need incubation conditions.
    if u in {"hr", "h", "min"}:
        return None, None, None

    if "%" in u:
        return None, None, None

    if u in {"nm", "um", "mm", "m"}:
        return None, None, None

    return None, None, None

def convert_scaled_clint_to_l_h_kg(value, unit):
    """
    Convert already-scaled Clint to L/h/kg where possible.
    """
    if value is None or unit is None:
        return None, None

    try:
        val = float(value)
    except (TypeError, ValueError):
        return None, None

    u = normalise_unit_string(unit)

    if any(p in u for p in [
        "l/h/kg",
        "l/hr/kg",
        "lhr-1kg-1",
        "lh-1kg-1",
        "l.kg-1.h-1",
        "lkg-1h-1",
    ]):
        return val, "L/h/kg"

    if any(p in u for p in [
        "ml/min/kg",
        "mlmin-1kg-1",
        "ml.kg-1.min-1",
        "mlkg-1min-1",
    ]):
        return val * 60.0 / 1000.0, "L/h/kg"

    return None, None

def scale_clint_to_l_h_kg(
    value,
    canonical_unit,
    matrix,
    species,
    scaling_cfg,
):
    """
    Scale metabolic stability / clearance to L/h/kg.

    Handles:
      microsomes:
        uL/min/mg protein

      hepatocytes:
        uL/min/10^6 cells

      liver-normalised clearance:
        mL/min/g liver
    """

    if value is None:
        return None, None

    try:
        val = float(value)
    except (TypeError, ValueError):
        return None, None

    if not scaling_cfg or not scaling_cfg.get("enabled", False):
        return None, None

    species_params = scaling_cfg.get("species_parameters", {})
    default_species = scaling_cfg.get("default_species")

    use_species = species

    if use_species not in species_params:
        if default_species in species_params:
            use_species = default_species
        else:
            return None, None

    params = species_params.get(use_species, {})
    liver_weight = params.get("liver_weight_g_per_kg")

    if liver_weight is None:
        return None, None

    # -----------------------------
    # Microsomal Clint
    # uL/min/mg protein -> L/h/kg
    # -----------------------------
    if canonical_unit == "uL/min/mg protein":
        mppgl = params.get("mppgl_mg_per_g_liver")

        if mppgl is None:
            return None, None

        scaled = val * mppgl * liver_weight * 60.0 / 1_000_000.0
        return scaled, "L/h/kg"

    # -----------------------------
    # Hepatocyte Clint
    # uL/min/10^6 cells -> L/h/kg
    # -----------------------------
    if canonical_unit == "uL/min/10^6 cells":
        hpgl = params.get("hpgl_million_cells_per_g_liver")

        if hpgl is None:
            return None, None

        scaled = val * hpgl * liver_weight * 60.0 / 1_000_000.0
        return scaled, "L/h/kg"

    # -----------------------------
    # Liver-normalised clearance
    # mL/min/g liver -> L/h/kg
    # -----------------------------
    if canonical_unit == "mL/min/g liver":
        scaled = val * liver_weight * 60.0 / 1000.0
        return scaled, "L/h/kg"

    # Already scaled
    if canonical_unit == "L/h/kg":
        return val, "L/h/kg"

    return None, None

def convert_ppb(value, unit):
    """
    Convert PPB to fraction unbound (0-1)
    """
    if value is None or unit is None:
        return None

    try:
        val = float(value)
    except:
        return None

    u = str(unit).lower().replace(" ", "").strip()

    if u in ["%", "percent", "%bound"]:
        return val / 100.0
    elif u in ["fu","fractionunbound","fraction_unbound"]:
        return val
    else:
        # units like ug/mL, pmol, etc are ignored
        return None