from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdChemReactions

# ---------------------------------------------------------------------------
# Reaction rules
# ---------------------------------------------------------------------------

REACTION_RULES = {
    "amide_formation": {
        "description": "Amide bond formation from an amine and carboxylic acid derivative.",
        "reaction_smarts": (
            "[C:1](=[O:2])-[N:3]"
            ">>"
            "[C:1](=[O:2])[O:4]."
            "[N:3]"
        ),
    },

    "ester_formation": {
        "description": "Ester bond formation from an alcohol and carboxylic acid derivative.",
        "reaction_smarts": (
            "[C:1](=[O:2])-[O:3]"
            ">>"
            "[C:1](=[O:2])[O:4]."
            "[O:3]"
        ),
    },

    "suzuki_coupling": {
        "description": "Suzuki coupling between an aryl/vinyl halide and a boronic acid/boronate.",
        "reaction_smarts": (
            "[c,C:1]-[Cl,Br,I:2]"
            ">>"
            "[c,C:1]-[Cl,Br,I:2]."
            "[B:3]"
        ),
    },

    "reductive_amination": {
        "description": "Reductive amination from an aldehyde/ketone and amine.",
        "reaction_smarts": (
            "[C:1](=[O:2])-[N:3]"
            ">>"
            "[C:1]=[O:2]."
            "[N:3]"
        ),
    },

    "ether_formation": {
        "description": "Ether formation from an alcohol/phenol and alkyl halide.",
        "reaction_smarts": (
            "[O:1]-[C:2]"
            ">>"
            "[O:1]."
            "[C:2]-[Cl,Br,I]"
        ),
    },
}


# Compile once when the module is imported.
COMPILED_REACTION_RULES = {}

for rule_name, rule_config in REACTION_RULES.items():
    reaction = rdChemReactions.ReactionFromSmarts(
        rule_config["reaction_smarts"]
    )

    if reaction is None:
        raise ValueError(
            f"Could not compile reaction rule: {rule_name}"
        )

    COMPILED_REACTION_RULES[rule_name] = reaction

def _canonicalise_smiles(smiles):
    """Return canonical SMILES or None if invalid."""
    try:
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            return None

        return Chem.MolToSmiles(mol, canonical=True)

    except Exception:
        return None


def _deduplicate_fragments(fragments):
    """Canonicalise and deduplicate a list of fragment SMILES."""

    unique = set()

    for smiles in fragments:
        canonical = _canonicalise_smiles(smiles)

        if canonical:
            unique.add(canonical)

    return sorted(unique)


def _apply_reaction_rule(
    target_smiles,
    rule_name,
    reaction,
):
    """
    Apply one retrosynthetic reaction rule to a target molecule.

    Returns a list of candidate disconnections.
    """

    target_mol = Chem.MolFromSmiles(target_smiles)

    if target_mol is None:
        return []

    candidates = []

    try:
        products = reaction.RunReactants((target_mol,))
    except Exception:
        return []

    for product_set in products:

        fragments = []

        for product in product_set:

            try:
                smiles = Chem.MolToSmiles(
                    product,
                    canonical=True,
                )
            except Exception:
                continue

            if smiles:
                fragments.append(smiles)

        fragments = _deduplicate_fragments(fragments)

        if not fragments:
            continue

        candidates.append({
            "reaction_rule": rule_name,
            "target_smiles": target_smiles,
            "precursor_smiles": fragments,
            "n_precursors": len(fragments),
        })

    # Deduplicate identical precursor sets.
    unique_candidates = []
    seen = set()

    for candidate in candidates:

        key = tuple(candidate["precursor_smiles"])

        if key in seen:
            continue

        seen.add(key)
        unique_candidates.append(candidate)

    return unique_candidates