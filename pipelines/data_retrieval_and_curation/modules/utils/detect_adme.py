import re
import numpy as np

from collections import Counter

from rdkit import Chem
from rdkit.Chem import Descriptors

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression

_CYP_MODEL = None
_CYP_VECTORIZER = None

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

import re


def extract_inhibition_concentration(description):
    """
    Extract inhibitor test concentration from assay text.

    Returns
    -------
    float
        concentration in µM
    """

    if not description:
        return None

    desc = str(description).lower()

    # unicode normalization
    desc = (
        desc.replace("μ", "u")
            .replace("µ", "u")
    )

    # ----------------------------------------
    # Supported forms:
    #
    # 10 uM
    # 10uM
    # 0.3 mM
    # 100 nM
    # 10 micromolar
    # 10 umol/l
    # at 30 uM
    # tested at 1 mM
    # ----------------------------------------

    pattern = re.compile(
        r'([\d\.]+)\s*'
        r'(nm|um|mm|micromolar|nanomolar|millimolar|umol/l|nmol/l|mmol/l)',
        re.IGNORECASE
    )

    matches = pattern.findall(desc)

    if not matches:
        return None

    # usually first concentration is the assay concentration
    value, unit = matches[0]

    try:
        value = float(value)
    except ValueError:
        return None

    unit = unit.lower()

    # ----------------------------------------
    # Convert everything -> µM
    # ----------------------------------------

    if unit in ["nm", "nanomolar", "nmol/l"]:
        value = value / 1000.0

    elif unit in ["um", "micromolar", "umol/l"]:
        value = value

    elif unit in ["mm", "millimolar", "mmol/l"]:
        value = value * 1000.0

    return value

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

def convert_bsep_activity(value, unit, endpoint_class):
    """
    Standardise BSEP endpoints.

    endpoint_class should come from classify_bsep_record():
        - "inhibition"
        - "exposure_ratio"
        - "transport"

    Returns:
        (numeric_value, endpoint_type, canonical_unit)

    Examples:
        IC50 uM -> nM
        inhibition % -> %
        ratio -> ratio
    """

    if value is None:
        return None, None, None

    try:
        val = float(value)
    except (TypeError, ValueError):
        return None, None, None

    u = (
        str(unit)
        .lower()
        .replace(" ", "")
        .replace("µ", "u")
        .strip()
        if unit is not None
        else ""
    )

    # =========================================================
    # EXPOSURE RATIO
    # =========================================================

    if endpoint_class == "exposure_ratio":

        # Ratios are unitless
        # Examples:
        #   Css/IC50
        #   Cmax/IC50
        #   exposure margin

        if val < 0:
            return None, None, None

        return val, "exposure_ratio", "ratio"

    # =========================================================
    # INHIBITION ASSAYS
    # =========================================================

    if endpoint_class == "inhibition":

        # -------------------------
        # Percent inhibition
        # -------------------------

        if "%" in u:
            return val, "inhibition", "%"

        # -------------------------
        # IC50/Ki concentration endpoints
        # Canonical unit = nM
        # -------------------------

        factor = None

        if u in ["nm", "nmol", "nanom", "nanomolar"]:
            factor = 1.0

        elif (
            "um" in u
            or "umol" in u
            or "micromolar" in u
        ):
            factor = 1e3

        elif u == "m":
            factor = 1e9

        elif (
            "mm" in u
            or "millimolar" in u
        ):
            factor = 1e6

        # Some ChEMBL records have missing units
        # but standard_type=IC50
        elif u in ["", "none", "null"]:
            # Too risky to guess generally
            return None, None, None

        if factor is not None:

            val_nM = val * factor

            # Remove pathological values
            if val_nM <= 0:
                return None, None, None

            return val_nM, "IC50", "nM"

        return None, None, None

    # =========================================================
    # TRANSPORT / SUBSTRATE ASSAYS
    # =========================================================

    if endpoint_class == "transport":

        # Keep transport-like readouts mostly raw for now.
        # You can later standardise:
        # uptake ratio
        # efflux ratio
        # % transport
        # ATPase activity

        if "%" in u:
            return val, "transport_percent", "%"

        if "ratio" in u or u == "":
            return val, "transport_ratio", "ratio"

        return val, "transport", u

    return None, None, None

def convert_oatp_activity(value, unit, endpoint_class):
    """
    Standardise OATP assay endpoints.

    Returns:
        (value, endpoint_type, canonical_unit)
    """

    if value is None:
        return None, None, None

    try:
        val = float(value)
    except (TypeError, ValueError):
        return None, None, None

    u = (
        str(unit)
        .lower()
        .replace(" ", "")
        .replace("µ", "u")
        .strip()
        if unit is not None
        else ""
    )

    # ---------------------------------
    # IC50/Ki concentration assays
    # ---------------------------------

    if endpoint_class == "IC50":

        factor = None

        if u in ["nm", "nanomolar"]:
            factor = 1.0

        elif "um" in u or "micromolar" in u:
            factor = 1e3

        elif "mm" in u:
            factor = 1e6

        elif u == "m":
            factor = 1e9

        if factor is None:
            return None, None, None

        val_nM = val * factor

        if val_nM <= 0:
            return None, None, None

        return val_nM, "IC50", "nM"

    # ---------------------------------
    # Percent inhibition
    # ---------------------------------

    if endpoint_class == "inhibition":

        return val, "inhibition", "%"

    # ---------------------------------
    # Transport / uptake
    # ---------------------------------

    if endpoint_class == "transport":

        if "%" in u:
            return val, "transport_percent", "%"

        return val, "transport", u or "raw"

    return None, None, None

def classify_bsep_record(record):

    desc = normalise_text(record.get("assay_description"))
    stype = normalise_text(record.get("standard_type"))
    unit = normalise_text(record.get("standard_units"))

    combined = f"{desc} {stype} {unit}"

    # -------------------------
    # Exposure/risk ratio
    # -------------------------
    ratio_patterns = [
        "ratio of drug concentration",
        "steady state",
        "css",
        "cmax",
        "auc",
        "to ic50",
        "exposure margin",
        "safety margin",
    ]

    if (
        stype == "ratio"
        or any(p in combined for p in ratio_patterns)
    ):
        return "exposure_ratio"

    # -------------------------
    # True inhibition assays
    # -------------------------
    inhibition_patterns = [
        "inhibition",
        "taurocholate transport",
        "atp dependent",
        "vesicle",
        "transport into vesicles",
        "bsep expressed",
    ]

    if (
        stype in ["ic50", "ki"]
        or any(p in combined for p in inhibition_patterns)
    ):
        return "inhibition"

    # -------------------------
    # Substrate/transport
    # -------------------------
    transport_patterns = [
        "substrate",
        "uptake",
        "transport assay",
    ]

    if any(p in combined for p in transport_patterns):
        return "transport"

    return None

def classify_oatp_record(record):
    """
    Classify OATP assay endpoint type.

    Returns:
        - "IC50"
        - "inhibition"
        - "transport"
        - None
    """

    desc = normalise_text(record.get("assay_description"))
    stype = normalise_text(record.get("standard_type"))
    unit = normalise_text(record.get("standard_units"))

    combined = f"{desc} {stype} {unit}"

    # ---------------------------------
    # IC50 / Ki style inhibition assays
    # ---------------------------------

    if stype in ["ic50", "ki", "kd", "ec50"]:
        return "IC50"

    # ---------------------------------
    # Percent inhibition assays
    # ---------------------------------

    if (
        "%" in unit
        or "inhibition" in stype
    ):
        return "inhibition"

    # ---------------------------------
    # Transport / uptake assays
    # ---------------------------------

    transport_patterns = [
        "uptake",
        "transport",
        "substrate",
        "efflux",
        "estrone",
        "estradiol",
        "pitavastatin",
        "taurocholate",
    ]

    if any(p in combined for p in transport_patterns):
        return "transport"

    return None