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
