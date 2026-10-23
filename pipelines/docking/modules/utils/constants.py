import re

ELEMENT_COLOURS = {
    "N": "blue",
    "O": "red",
    "S": "yellow",
    "P": "orange",
    "H": "white",
    "F": "palecyan",
    "Cl": "green",
    "Br": "brown",
    "I": "purple",
}

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","HID","HIE","HIP","ASH","GLH","LYN","CYX"}
IGNORE = {"HOH","WAT","SOL","DOD","NA","K","CL","CA","MG","ZN","MN","FE","SO4","PO4","GOL","EDO","PEG"}
ID_FIELDS = ("s_user_compound_id", "Compound ID", "Compound_ID", "ID", "Name", "s_m_title")

ELEMENT_COLOURS = {
    "N": "blue",
    "O": "red",
    "S": "yellow",
    "P": "orange",
    "H": "white",
    "F": "palecyan",
    "Cl": "green",
    "Br": "brown",
    "I": "purple",
}


def element_colour_commands(
    object_name: str,
    carbon_colour: str,
) -> list[str]:
    """
    Return PyMOL commands for consistent element colouring.

    Carbon uses the supplied object-specific colour.
    Other supported elements use a conventional element colour.
    """

    commands = [
        (
            f"color {carbon_colour}, "
            f"({object_name} and elem C)"
        )
    ]

    for element, colour in ELEMENT_COLOURS.items():
        commands.append(
            f"color {colour}, "
            f"({object_name} and elem {element})"
        )

    return commands

THREE_TO_ONE = {
    # Standard residues
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",

    # Common prepared/protonation variants
    "ASH": "D",
    "GLH": "E",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "LYN": "K",
    "CYM": "C",
    "CYX": "C",

    # Common modified residue interpreted as parent residue
    "MSE": "M",
}

PYMOL_COLOURS = [
    "marine",
    "orange",
    "green",
    "magenta",
    "cyan",
    "salmon",
    "yellow",
    "violet",
    "teal",
    "wheat",
]

def _safe_object_name(
    value: str,
) -> str:
    value = re.sub(
        r"[^A-Za-z0-9_]+",
        "_",
        value,
    )

    return value.strip(
        "_"
    )
