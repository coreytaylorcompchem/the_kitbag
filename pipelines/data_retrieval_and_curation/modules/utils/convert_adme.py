import re

from collections import Counter
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression

_CYP_MODEL = None
_CYP_VECTORIZER = None

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

    # return val, canonical
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
    
def extract_inhibition_concentration(description):
    """
    Extract concentration from assay description.
    Returns concentration in µM (float) or None.
    """

    if not description:
        return None

    desc = description.lower().replace("μ", "u")

    # patterns like:
    # "at 30 uM", "at 100uM", "at 1 um"
    patterns = [
        r'at\s+([\d\.]+)\s*u?m',
        r'([\d\.]+)\s*u?m',
        # r'([\d\.e-]+)\s*m',  # convert M → µM might add later
    ]

    for pattern in patterns:
        match = re.search(pattern, desc)
        if match:
            try:
                return float(match.group(1))
            except:
                continue

    return None

def normalise_conc(conc):
    allowed = [0.1, 1, 3, 10, 30, 100]
    return min(allowed, key=lambda x: abs(x - conc))


def normalise_text(t):
    if not t:
        return ""
    t = str(t).lower()
    t = t.replace("μ", "u")
    t = t.replace("β", "beta")
    t = re.sub(r"[’']", "", t)
    t = re.sub(r"[-_/]", " ", t)
    t = re.sub(r"\s+", " ", t)
    return t

# ----------------------------
# CYP3A4 substrate/probe groups
# ----------------------------

def score_cyp_keywords(t):
    """
    Weak keyword scoring system.
    Returns (label, score_dict)
    """

    scores = {
        "midazolam_like": 0,
        "testosterone_like": 0,
    }

    # MIDAZOLAM SIGNALS
    if any(k in t for k in [
        "midazolam", "mdz", "1 oh midazolam", "4 oh midazolam",
        "hydroxymidazolam", "triazolam", "alprazolam"
    ]):
        scores["midazolam_like"] += 3

    if any(k in t for k in [
        "benzodiazepine", "azolam", "cyp3a4 metabolism mdz"
    ]):
        scores["midazolam_like"] += 1

    # TESTOSTERONE SIGNALS
    if any(k in t for k in [
        "testosterone", "6beta", "6 beta", "6-b", "6β",
        "hydroxytestosterone"
    ]):
        scores["testosterone_like"] += 3

    if any(k in t for k in [
        "steroid", "androgen", "cortisol", "progesterone"
    ]):
        scores["testosterone_like"] += 1

    # DECISION
    best = max(scores, key=scores.get)
    if scores[best] >= 2:
        return best, scores

    return None, scores

MIDAZ_LIKE = [

    # canonical probe drugs
    r"\bmidazolam\b",
    r"\bmdz\b",
    r"\btriazolam\b",
    r"\bnifedipine\b",
    r"\berythromycin\b",

    # metabolite / hydroxylation forms
    r"hydroxy.*midazol",
    r"midazol.*hydroxy",
    r"midazol.*hydroxyl",
    r"hydroxylation.*midazol",

    r"hydroxy.*mdz",
    r"mdz.*hydroxy",

    # positional forms
    r"1.? ?hydroxy.*midazol",
    r"1.? ?oh.*midazol",
    r"4.? ?hydroxy.*midazol",
    r"4.? ?oh.*midazol",

    r"1.?hydroxymidazolam",
    r"midazolam.*hydroxylation",
]

TESTO_LIKE = [

    r"\btestosterone\b",
    r"\btst\b",
    r"\btestoster",

    r"hydroxy.*testosterone",
    r"testosterone.*hydroxy",

    # CYP3A4 signals
    r"6.?beta.*hydroxyl",
    r"6beta.*hydroxylation",
    r"6beta.*hydroxylase",
    r"formation.*6beta.*hydroxytestosterone",

    r"6.?b.*hydroxy.*testosterone",

    r"6.?beta.?hydroxytestosterone",
    r"6.?β.?hydroxytestosterone",
    r"testosterone.*6.?beta",
    r"6.?beta.*testosterone",
    r"6.?beta.*oh",
    r"6β",
    r"6 beta",
    r"6-b",
    r"hydroxytestosterone",
    r"steroid.*hydroxylation",
    r"androgen.*metabolism",
]


FLUOROGENIC = [

    r"\bdibenzylfluorescein\b",
    r"\bdbf\b",

    r"\bbfc\b",
    r"benzyloxy.*trifluoromethyl.*coumarin",
    r"trifluoromethyl.*coumarin",

    r"benzyloxyresorufin",
    r"\bbrod\b",

    r"luciferin",
]


REPORTER = [

    r"qhts assay.*cyp3a4",
    r"activators? of cytochrome p450 3a4",
    r"inhibitors? of cytochrome p450 3a4",
    r"pubchem_bioassay.*cyp3a4",
]

OTHER_PROBES = [

    # Fluorogenic probes
    r"\bdibenzylfluorescein\b",
    r"\bdbf\b",
    r"\bbfc\b",
    r"7-benzyloxy-4-\(trifluoromethyl\)-coumarin",
    r"7-benzyloxy.*coumarin",
    r"benzyloxyresorufin",
    r"luciferin[- ]?ipa",
    r"luciferin",
    r"p450-glo",

    # Major clinical CYP3A4 probes
    r"\bcyclosporin\b",
    r"\bcyclosporine\b",
    r"\btacrolimus\b",
    r"\bterfenadine\b",
    r"\bfelodipine\b",
    r"\bquinidine\b",
    r"\balprazolam\b",
    r"\bbuspirone\b",

    # Steroid-like probes related to testosterone
    r"\bcortisol\b",
    r"\bcortisone\b",
    r"\bprogesterone\b",

    # Generic metabolism readouts
    r"6beta hydroxylation",
    r"1.?hydroxymidazolam",
]

GENERIC_FUNCTIONAL = [

    # PubChem style
    r"pubchem_bioassay.*cyp3a4",
    r"qhts assay.*cyp3a4",

    # Functional assay language
    r"activators? of cytochrome p450 3a4",
    r"inhibitors? of cytochrome p450 3a4",
    r"cyp3a4 inhibition assay",
    r"cyp3a4 activity assay",
    r"recombinant cyp3a4",
    r"human liver microsome.*cyp3a4",
    r"time dependent inhibition.*3a4",
    r"mechanism based inhibition.*3a4",

    # Generic "substrate metabolism" wording
    r"cyp3a4 mediated metabolism",
    r"cyp3a4 catalyzed",
    r"substrate depletion.*3a4",
]

def extract_cyp3a4_substrate(text):

    if not text:
        return "unknown"

    t = normalise_text(text)

    # STRONG ANCHORS
    if "midazolam" in t or "mdz" in t:
        return "midazolam_like"

    if "testosterone" in t:
        return "testosterone_like"

    # WEAK KEYWORD SCORING
    weak_label, scores = score_cyp_keywords(t)
    if weak_label:
        return weak_label

    # RULES
    for p in REPORTER:
        if re.search(p, t):
            return "substrate_agnostic"

    for p in FLUOROGENIC:
        if re.search(p, t):
            return "substrate_agnostic"

    for p in OTHER_PROBES:
        if re.search(p, t):
            return "probe_substrate_other"

    for p in GENERIC_FUNCTIONAL:
        if re.search(p, t):
            return "substrate_agnostic"

    if "cyp3a4" in t:
        return "substrate_agnostic"

    return "unknown"
    
def train_cyp_classifier(texts, labels):
    global _CYP_MODEL, _CYP_VECTORIZER

    # Balance dataset
    counts = Counter(labels)
    min_class = min(counts.values())

    balanced_texts = []
    balanced_labels = []

    for cls in counts:
        cls_samples = [(t, l) for t, l in zip(texts, labels) if l == cls]
        cls_samples = cls_samples[:min_class]

        for t, l in cls_samples:
            balanced_texts.append(t)
            balanced_labels.append(l)

    vec = TfidfVectorizer(
        ngram_range=(1, 3),   # was (1,2) - broader match
        min_df=2,  
        max_features=10000
    )

    X = vec.fit_transform(balanced_texts)

    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    clf.fit(X, balanced_labels)

    _CYP_MODEL = clf
    _CYP_VECTORIZER = vec

def predict_substrate_ml(text):
    global _CYP_MODEL, _CYP_VECTORIZER

    if _CYP_MODEL is None:
        return None, 0.0

    t = normalise_text(text)
    X = _CYP_VECTORIZER.transform([t])

    probs = _CYP_MODEL.predict_proba(X)[0]
    idx = probs.argmax()

    return _CYP_MODEL.classes_[idx], probs[idx]

def extract_cyp3a4_substrate_hybrid(text, threshold=0.5):

    t = normalise_text(text)

    # 1. Rules first
    rule = extract_cyp3a4_substrate(t)

    if rule in ["midazolam_like", "testosterone_like"]:
        return rule

    # 2. Weak keyword override
    weak_label, scores = score_cyp_keywords(t)
    if weak_label:
        return weak_label

    # 3. ML fallback, which is more permissive
    pred, conf = predict_substrate_ml(t)

    if pred is not None:
        X = _CYP_VECTORIZER.transform([t])
        probs = _CYP_MODEL.predict_proba(X)[0]

        # loose acceptance - TODO: should tweak this
        if max(probs) > 0.45:
            return pred

    return rule

def classify_metstab_record(record):
    """
    Determine whether a record is genuinely a metabolic stability assay
    rather than a CYP inhibition/metabolism assay performed in microsomes.
    """

    desc = normalise_text(record.get("assay_description"))
    stype = normalise_text(record.get("standard_type"))
    unit = normalise_text(record.get("standard_units"))

    combined = f"{desc} {stype} {unit}"

    # -------------------------
    # Positive signals
    # -------------------------
    positive_patterns = [

        # intrinsic clearance
        "intrinsic clearance",
        "clint",

        # stability language
        "metabolic stability",
        "microsomal stability",
        "hepatocyte stability",
        "substrate depletion",

        # half-life
        "half life",
        "t1/2",

        # common units
        "ul/min/mg",
        "ml/min/g",
        "ul/min/10^6",

        # depletion style
        "% remaining",
        "percent remaining",
    ]

    # -------------------------
    # Negative signals
    # -------------------------
    negative_patterns = [

        # CYP inhibition
        "inhibition",
        "ic50",
        "ec50",
        "ki",
        "kd",

        # enzyme phenotyping
        "metabolite formation",
        "enzyme activity",
        "probe substrate",

        # common CYP wording
        "cyp1a2",
        "cyp2c9",
        "cyp2c19",
        "cyp2d6",
        "cyp3a4",
        "cyp3a5",

        # classic probe substrates
        "midazolam",
        "testosterone",
        "diclofenac",
        "dextromethorphan",
    ]

    positive = any(p in combined for p in positive_patterns)
    negative = any(n in combined for n in negative_patterns)

    return positive and not negative