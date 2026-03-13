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

def convert_met_stab(value, unit):
    """
    Convert microsomal/hepatocyte clearance units to canonical mL/min/g liver or uL/min/10^6 cells.
    Returns (float_value, canonical_unit) or (None, None) if unconvertible.
    """
    if value is None or unit is None:
        return None, None

    u = str(unit).lower().replace(" ", "").replace("*","").replace("^","-").replace("(","").replace(")","").strip()

    # Canonical units
    if "cell" in u:
        canonical = "uL/min/10^6 cells"
    else:
        canonical = "mL/min/g liver"

    try:
        val = float(value)
    except:
        return None, None

    # Mapping known variants
    direct_map = {
        "ul/min": 1.0 if canonical=="uL/min/10^6 cells" else None,
        "ul.min-1.(10^6cells)-1":1.0,
        "uL/min/1e6cells":1.0,
        "ml.min-1.kg-1": 1.0,   # assume mL/min/g conversion will be handled later
        "ml.min-1.g-1": 1.0,
        "hr": None,  # needs half-life conversion, skip for now
        "min": 1.0,
    }

    for k,v in direct_map.items():
        if k in u and v is not None:
            return val * v, canonical

    # Reject obviously incompatible units
    reject_patterns = ["%", "nmol", "nM", "ug", "mM", "pmol"]
    if any(r in u for r in reject_patterns):
        return None, None

    return val, canonical

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