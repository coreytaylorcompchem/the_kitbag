# ---------------------------------------------------------------------------
# Reaction-handle SMARTS
# ---------------------------------------------------------------------------

REACTION_HANDLE_SMARTS = {
    "amine": "[N;H1,H2;!$(N-C=O);!$(N-S=O)]",
    "alcohol": "[O;H1;!$(O-C=O)]",
    "phenol": "[O;H1]-[c]",
    "carboxylic_acid": "C(=O)[O;H1]",
    "acid_chloride": "C(=O)Cl",
    "aldehyde": "[CX3H1](=O)",
    "ketone": "[#6][CX3](=O)[#6]",
    "alkyl_halide": "[CX4][Cl,Br,I]",
    "aryl_vinyl_halide": "[c,C;X3,X2][Cl,Br,I]",
    "boronic_acid": "[B;$(B(O)O)]",
    "boronate_ester": "[B;X3]([O;X2])[O;X2]",
    "sulfonyl_chloride": "S(=O)(=O)Cl",
    "nitrile": "C#N",
}