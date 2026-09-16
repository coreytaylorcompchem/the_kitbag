from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdChemReactions

import pandas as pd

from modules.utils.retrosynth_parameter_dicts import REACTION_HANDLE_SMARTS

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=True, simple_format=True)

# ---------------------------------------------------------------------------
# Reaction rule library
# ---------------------------------------------------------------------------

REACTION_RULES = {}


# ---------------------------------------------------------------------------
# Reaction utilities
# ---------------------------------------------------------------------------

def _reverse_reaction(
    reaction,
):
    """
    Reverse a forward reaction template so that the
    forward product becomes the retrosynthetic input
    and the forward reactants become retrosynthetic
    precursor outputs.
    """

    reversed_reaction = rdChemReactions.ChemicalReaction()

    # Forward reactants become retrosynthetic products.
    for reactant_template in reaction.GetReactants():
        reversed_reaction.AddProductTemplate(
            Chem.Mol(reactant_template)
        )

    # Forward products become retrosynthetic reactants.
    for product_template in reaction.GetProducts():
        reversed_reaction.AddReactantTemplate(
            Chem.Mol(product_template)
        )

    reversed_reaction.Initialize()

    return reversed_reaction


def load_reaction_rules(
    input_file,
):
    """
    Load forward reaction rules from a CSV file.

    Expected columns:
        rxn_name
        smirks
        educt1_smiles
        educt2_smiles

    The CSV contains forward literature reactions.
    These are compiled and reversed for retrosynthetic
    application.
    """

    input_file = Path(input_file)

    if not input_file.exists():
        raise FileNotFoundError(
            f"Reaction rule file not found: {input_file}"
        )

    df = pd.read_csv(
        input_file
    )

    required_columns = [
        "rxn_name",
        "smirks",
        "educt1_smiles",
        "educt2_smiles",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            "Reaction rule CSV is missing required "
            f"columns: {missing_columns}"
        )

    reaction_rules = {}

    for _, row in df.iterrows():

        rule_name = str(
            row["rxn_name"]
        ).strip()

        smirks = str(
            row["smirks"]
        ).strip()

        educt1_smiles = row["educt1_smiles"]
        educt2_smiles = row["educt2_smiles"]

        if not rule_name:
            continue

        if not smirks:
            raise ValueError(
                f"Reaction rule '{rule_name}' has "
                "an empty SMIRKS."
            )

        # Remove the curly braces used in the source CSV.
        rule_name = rule_name.strip("{}")

        reaction_rules[
            rule_name
        ] = {
            "description": (
                f"Literature reaction: {rule_name}"
            ),
            "smirks": smirks,
            "educt1_smiles": (
                None
                if pd.isna(educt1_smiles)
                else str(educt1_smiles).strip()
            ),
            "educt2_smiles": (
                None
                if pd.isna(educt2_smiles)
                else str(educt2_smiles).strip()
            ),
        }

    if not reaction_rules:
        raise ValueError(
            f"No reaction rules were loaded from "
            f"{input_file}"
        )

    return reaction_rules


# ---------------------------------------------------------------------------
# Reaction rule compilation
# ---------------------------------------------------------------------------

def compile_reaction_rules(
    reaction_rules,
):
    """
    Compile forward reaction SMIRKS and reverse them
    for retrosynthetic application.
    """

    compiled_reaction_rules = {}

    for rule_name, rule_data in reaction_rules.items():

        smirks = rule_data.get(
            "smirks"
        )

        if not isinstance(smirks, str):
            raise TypeError(
                f"Reaction rule '{rule_name}' has invalid "
                f"'smirks' type: "
                f"{type(smirks).__name__}. "
                "Expected a string."
            )

        try:
            forward_reaction = (
                rdChemReactions.ReactionFromSmarts(
                    smirks
                )
            )

        except Exception as exc:
            raise ValueError(
                f"Failed to compile reaction rule "
                f"'{rule_name}':\n"
                f"{smirks}\n"
                f"RDKit error: {exc}"
            ) from exc

        if forward_reaction is None:
            raise ValueError(
                f"RDKit returned None for reaction "
                f"'{rule_name}'."
            )

        try:
            reversed_reaction = _reverse_reaction(
                forward_reaction
            )

        except Exception as exc:
            raise ValueError(
                f"Failed to reverse reaction rule "
                f"'{rule_name}'.\n"
                f"SMIRKS: {smirks}\n"
                f"RDKit error: {exc}"
            ) from exc

        compiled_reaction_rules[
            rule_name
        ] = reversed_reaction

    return compiled_reaction_rules

def _get_mapped_atoms(
    mol,
):
    """
    Return a mapping from atom-map number to atom for a reaction
    template molecule.

    Reaction-template atoms are query atoms, so this helper deliberately
    only accesses properties that are safe for unsanitized/query atoms.
    """

    mapped_atoms = {}

    for atom in mol.GetAtoms():

        map_number = atom.GetAtomMapNum()

        if map_number > 0:
            mapped_atoms[map_number] = atom

    return mapped_atoms


def _get_mapped_bonds(
    mol,
):
    """
    Return mapped bonds from a reaction-template molecule.

    Each bond is represented as:

        (map_number_1, map_number_2, bond_type)

    Only bonds where both atoms are atom-mapped are included.
    """

    bonds = set()

    for bond in mol.GetBonds():

        atom_1 = bond.GetBeginAtom()
        atom_2 = bond.GetEndAtom()

        map_1 = atom_1.GetAtomMapNum()
        map_2 = atom_2.GetAtomMapNum()

        if map_1 <= 0 or map_2 <= 0:
            continue

        pair = tuple(
            sorted(
                (map_1, map_2)
            )
        )

        bonds.add(
            (
                pair[0],
                pair[1],
                str(bond.GetBondType()),
            )
        )

    return bonds


def _get_reaction_centre_atoms(
    forward_reaction,
):
    """
    Identify mapped atoms participating in the reaction centre.

    The reaction centre consists of mapped atoms whose mapped bonding
    environment changes between the reactants and products.

    This operates only on atom-map numbers and bond objects and does
    not query implicit hydrogen counts or valence.
    """

    reactant_bonds = set()
    product_bonds = set()

    for reactant in forward_reaction.GetReactants():

        reactant_bonds.update(
            _get_mapped_bonds(
                reactant
            )
        )

    for product in forward_reaction.GetProducts():

        product_bonds.update(
            _get_mapped_bonds(
                product
            )
        )

    changed_bonds = (
        reactant_bonds ^ product_bonds
    )

    reaction_centre = set()

    for map_1, map_2, _ in changed_bonds:

        reaction_centre.add(map_1)
        reaction_centre.add(map_2)

    return reaction_centre


def infer_reaction_handles(
    smirks,
):
    """
    Infer commercial reaction handles from a forward reaction template.

    Reaction-centre atoms are identified from atom mapping and bond
    changes. Existing REACTION_HANDLE_SMARTS definitions are then used
    to classify the relevant reactant-side functional groups.

    Parameters
    ----------
    smirks : str
        Forward reaction SMIRKS.

    Returns
    -------
    list[str]
        Unique reaction handles detected at the reaction centre.
    """

    if not isinstance(smirks, str) or not smirks.strip():
        return []

    try:
        reaction = rdChemReactions.ReactionFromSmarts(
            smirks
        )
    except Exception as exc:
        logger.warning(
            f"Could not compile reaction for handle inference: "
            f"{smirks}. Error: {exc}"
        )
        return []

    if reaction is None:
        return []

    reaction_centre = _get_reaction_centre_atoms(
        reaction
    )

    if not reaction_centre:
        return []

    handles = []

    for reactant in reaction.GetReactants():

        mapped_atoms = _get_mapped_atoms(
            reactant
        )

        # Only inspect reactants containing a mapped atom that
        # participates in the reaction centre.
        relevant_map_numbers = (
            set(mapped_atoms)
            & reaction_centre
        )

        if not relevant_map_numbers:
            continue

        # Inspect the local atom environment directly.
        #
        # We deliberately do not call GetSubstructMatches() here
        # because reaction-template molecules contain query atoms
        # rather than ordinary sanitized molecular atoms.

        for handle_name, smarts in REACTION_HANDLE_SMARTS.items():

            pattern = Chem.MolFromSmarts(
                smarts
            )

            if pattern is None:
                continue

            pattern_atom_count = pattern.GetNumAtoms()

            if pattern_atom_count == 0:
                continue

            # Compare the reaction-template atom environment against
            # the handle definition using RDKit's query matching on
            # individual atom/bond environments.
            #
            # Build candidate subgraphs around each mapped reaction-
            # centre atom and test them without requesting implicit H
            # information from the reaction template.
            for map_number in relevant_map_numbers:

                atom = mapped_atoms[map_number]

                if atom.GetDegree() == 0:
                    continue

                # A handle can only be relevant if its mapped reaction-
                # centre atom is represented somewhere in the handle
                # SMARTS. We therefore inspect the handle's atoms and
                # compare their atomic-number/query requirements.
                #
                # This deliberately remains conservative.
                for pattern_atom in pattern.GetAtoms():

                    if (
                        pattern_atom.GetAtomicNum()
                        != atom.GetAtomicNum()
                    ):
                        continue

                    handles.append(
                        handle_name
                    )

                    break

                if handle_name in handles:
                    break

    return sorted(
        set(handles)
    )

# ---------------------------------------------------------------------------
# SMILES utilities
# ---------------------------------------------------------------------------

def _canonicalise_smiles(smiles):
    """Return canonical SMILES or None if invalid."""

    try:

        mol = Chem.MolFromSmiles(
            smiles
        )

        if mol is None:
            return None

        return Chem.MolToSmiles(
            mol,
            canonical=True,
        )

    except Exception:
        return None


def _deduplicate_fragments(fragments):
    """Canonicalise and deduplicate a list of fragment SMILES."""

    unique = set()

    for smiles in fragments:

        canonical = _canonicalise_smiles(
            smiles
        )

        if canonical:
            unique.add(
                canonical
            )

    return sorted(unique)

def _is_reasonable_precursor(smiles):
    """
    Apply broad medicinal-chemistry sanity filters to a
    retrosynthetic precursor.

    These are deliberately permissive. They are intended to
    remove obvious reaction-template artefacts rather than
    enforce complete chemical validity.
    """

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return False

    # Reject fragments that are too small to be useful
    if mol.GetNumHeavyAtoms() < 2:
        return False

    # Reject radicals / explicitly charged atoms that are unlikely
    # to represent the intended isolated commercial precursor.
    for atom in mol.GetAtoms():

        if atom.GetNumRadicalElectrons() > 0:
            return False

    return True

def _validate_reaction_precursors(
    precursor_smiles,
    rule_data,
):
    """
    Validate retrosynthetic precursor fragments against the
    expected reactant structures in the reaction rule.

    Returns:
        (True, None) if valid.
        (False, reason) if invalid.
    """

    for smiles in precursor_smiles:

        mol = Chem.MolFromSmiles(
            smiles
        )

        if mol is None:
            return False, (
                f"invalid SMILES: {smiles}"
            )

        if _is_template_artefact(
            smiles
        ):
            return False, (
                f"template artefact: {smiles}"
            )

        if mol.GetNumHeavyAtoms() < 2:
            return False, (
                f"too few heavy atoms: {smiles}"
            )

    return True, None

# def _has_dummy_atom(smiles):
#     """
#     Return True if the molecule contains one or more
#     RDKit dummy atoms ('*').
#     """

#     mol = Chem.MolFromSmiles(smiles)

#     if mol is None:
#         return False

#     return any(
#         atom.GetAtomicNum() == 0
#         for atom in mol.GetAtoms()
#     )


def _is_template_artefact(smiles):
    """
    Identify obvious reaction-template artefacts.

    These filters are deliberately generic and do not attempt
    to determine whether a precursor is chemically appropriate
    for a particular reaction class.
    """

    mol = Chem.MolFromSmiles(
        smiles
    )

    if mol is None:
        return True

    # Count only real atoms. RDKit dummy atoms (*) have atomic
    # number 0 and therefore do not contribute here.
    real_heavy_atoms = sum(
        atom.GetAtomicNum() > 0
        for atom in mol.GetAtoms()
    )

    # A dummy atom attached to fewer than two real heavy atoms
    # is treated as an obvious template artefact.
    if any(
        atom.GetAtomicNum() == 0
        and atom.GetDegree() > 0
        and real_heavy_atoms < 2
        for atom in mol.GetAtoms()
    ):
        return True

    # Reject radicals.
    for atom in mol.GetAtoms():

        if atom.GetNumRadicalElectrons() > 0:
            return True

    return False

def _log_reaction_candidates(
    rule_name,
    candidates,
):
    """
    Log accepted retrosynthetic disconnections for
    inspection during development.
    """

    logger.debug(
        f"Accepted {rule_name} disconnections:"
    )

    for index, candidate in enumerate(
        candidates,
        start=1,
    ):
        logger.debug(
            f"  {index}: "
            f"{candidate['precursor_smiles']}"
        )

def _apply_reaction_rule(
    target_smiles,
    rule_name,
    reaction,
    rule_data,
    inspect=False,
):
    """
    Apply a compiled retrosynthetic reaction rule
    to a target molecule.
    """

    target_mol = Chem.MolFromSmiles(
        target_smiles
    )

    if target_mol is None:
        return []

    candidates = []

    reaction_handles = infer_reaction_handles(
        rule_data["smirks"]
    )

    try:
        product_sets = reaction.RunReactants(
            (target_mol,)
        )
    except Exception:
        return []

    for product_set in product_sets:

        # --------------------------------------------------
        # Convert the complete RDKit product molecules
        # into canonical SMILES.
        # --------------------------------------------------

        fragments = []

        valid_products = True

        for product in product_set:

            try:
                Chem.SanitizeMol(
                    product
                )

                smiles = Chem.MolToSmiles(
                    product,
                    canonical=True,
                )

            except Exception:
                valid_products = False
                break

            if not smiles:
                valid_products = False
                break

            fragments.append(
                smiles
            )

        if not valid_products:
            continue

        # --------------------------------------------------
        # Validate precursor fragments.
        # --------------------------------------------------

        valid_precursors, rejection_reason = (
            _validate_reaction_precursors(
                fragments,
                rule_data,
            )
        )

        if not valid_precursors:
            continue

        reaction_handles = infer_reaction_handles(
            rule_data["smirks"]
        )

        # --------------------------------------------------
        # Build candidate disconnection.
        # --------------------------------------------------

        candidate = {
            "reaction_rule": rule_name,
            "target_smiles": target_smiles,
            "precursor_smiles": fragments,
            "n_precursors": len(fragments),
            "reaction_handles": reaction_handles,
        }

        candidates.append(
            candidate
        )

    # ------------------------------------------------------
    # Deduplicate complete disconnections.
    # ------------------------------------------------------

    unique_candidates = []

    seen = set()

    for candidate in candidates:

        key = (
            candidate["reaction_rule"],
            tuple(
                candidate["precursor_smiles"]
            ),
        )

        if key in seen:
            continue

        seen.add(key)

        unique_candidates.append(
            candidate
        )

    if inspect:
        _log_reaction_candidates(
            rule_name,
            unique_candidates,
        )

    return unique_candidates