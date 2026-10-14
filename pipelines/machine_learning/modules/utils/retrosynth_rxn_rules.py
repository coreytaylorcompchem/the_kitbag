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

def _handle_match_is_reaction_centre_related(
    mol,
    handle_match,
    reaction_centre_atom_indices,
):
    """
    Determine whether a matched commercial reaction handle is
    associated with the reaction centre.

    A handle is considered relevant when:

        1. one of its atoms is itself a reaction-centre atom, or
        2. one of its atoms is directly bonded to a reaction-centre atom.

    The second condition is important for handles such as boronic
    acid in Suzuki coupling, where the reaction-centre atom is the
    carbon attached to boron rather than the boron atom itself.
    """

    handle_atom_indices = set(
        handle_match
    )

    reaction_centre_atom_indices = set(
        reaction_centre_atom_indices
    )

    # Direct overlap.
    if (
        handle_atom_indices
        & reaction_centre_atom_indices
    ):
        return True

    # Directly attached to the reaction centre.
    for atom_index in handle_atom_indices:

        atom = mol.GetAtomWithIdx(
            atom_index
        )

        for neighbour in atom.GetNeighbors():

            if (
                neighbour.GetIdx()
                in reaction_centre_atom_indices
            ):
                return True

    return False


def _get_reaction_centre_atom_indices(
    reactant_template,
    template_match,
    reaction_centre,
):
    """
    Map reaction-centre atom-map numbers from a reaction template
    onto atom indices in an actual educt molecule.

    Parameters
    ----------
    reactant_template
        RDKit reaction-template molecule.

    template_match
        Atom-index mapping returned by matching the reaction template
        against an actual educt molecule.

    reaction_centre
        Set of reaction-centre atom-map numbers.

    Returns
    -------
    set[int]
        Atom indices in the actual educt molecule corresponding to
        reaction-centre atoms.
    """

    reaction_centre_atom_indices = set()

    for template_atom_index, actual_atom_index in enumerate(
        template_match
    ):

        template_atom = (
            reactant_template.GetAtomWithIdx(
                template_atom_index
            )
        )

        map_number = (
            template_atom.GetAtomMapNum()
        )

        if map_number in reaction_centre:

            reaction_centre_atom_indices.add(
                actual_atom_index
            )

    return reaction_centre_atom_indices

def _is_sensible_suzuki_boron_candidate(mol):
    """
    Check whether a commercial boron-containing molecule is a
    chemically sensible Suzuki coupling partner.
    """

    if mol is None:
        return False

    boron_query = Chem.MolFromSmarts(
        "[B]"
    )

    boron_atoms = mol.GetSubstructMatches(
        boron_query
    )

    if len(boron_atoms) != 1:
        return False

    boron_idx = boron_atoms[0][0]

    boron = mol.GetAtomWithIdx(
        boron_idx
    )

    carbon_neighbours = [
        neighbour
        for neighbour in boron.GetNeighbors()
        if neighbour.GetAtomicNum() == 6
    ]

    if len(carbon_neighbours) != 1:
        return False

    carbon = carbon_neighbours[0]

    if carbon.GetHybridization() != Chem.HybridizationType.SP2:
        return False

    if any(
        bond.GetBondType() == Chem.BondType.TRIPLE
        for bond in carbon.GetBonds()
    ):
        return False

    oxygen_neighbours = [
        neighbour
        for neighbour in boron.GetNeighbors()
        if neighbour.GetAtomicNum() == 8
    ]

    if len(oxygen_neighbours) != 2:
        return False

    return True


def infer_reaction_handles(
    rule_data,
):
    """
    Infer the required commercial reaction handles for each
    reactant/precursor position in a reaction rule.

    The returned list corresponds positionally to the forward
    reaction reactants and therefore to the retrosynthetic
    precursor fragments generated by the reversed reaction.

    For example, a Suzuki reaction may return:

        [
            ["boronic_acid"],
            ["aryl_vinyl_halide"],
        ]

    meaning:

        precursor 0 requires boronic_acid
        precursor 1 requires aryl_vinyl_halide

    Parameters
    ----------
    rule_data : dict
        Reaction-rule dictionary containing:

            smirks
            educt1_smiles
            educt2_smiles

    Returns
    -------
    list[list[str]]
        Reaction handles for each reactant position.
    """

    if not isinstance(
        rule_data,
        dict,
    ):
        return []

    smirks = rule_data.get(
        "smirks"
    )

    if (
        not isinstance(smirks, str)
        or not smirks.strip()
    ):
        return []

    try:

        reaction = (
            rdChemReactions.ReactionFromSmarts(
                smirks
            )
        )

    except Exception as exc:

        logger.warning(
            "Could not compile reaction for handle "
            f"inference: {smirks}. Error: {exc}"
        )

        return []

    if reaction is None:
        return []

    reaction_centre = (
        _get_reaction_centre_atoms(
            reaction
        )
    )

    if not reaction_centre:

        logger.debug(
            "No mapped reaction centre found for "
            f"reaction rule: {smirks}"
        )

        return [
            []
            for _ in reaction.GetReactants()
        ]

    # --------------------------------------------------------------
    # Get the actual example educts from the reaction-rule CSV.
    #
    # These are ordinary molecules, unlike the query molecules
    # contained in the reaction template.
    # --------------------------------------------------------------

    educt_smiles = []

    for column in (
        "educt1_smiles",
        "educt2_smiles",
    ):

        value = rule_data.get(
            column
        )

        if (
            value is None
            or pd.isna(value)
        ):
            educt_smiles.append(
                None
            )
        else:
            value = str(
                value
            ).strip()

            educt_smiles.append(
                value
                if value
                else None
            )

    reactant_templates = list(
        reaction.GetReactants()
    )

    # --------------------------------------------------------------
    # Infer handles independently for each reactant.
    # --------------------------------------------------------------

    reactant_handles = []

    for reactant_index, reactant_template in enumerate(
        reactant_templates
    ):

        handles = []

        if reactant_index >= len(
            educt_smiles
        ):

            logger.debug(
                "No example educt available for reactant "
                f"position {reactant_index}."
            )

            reactant_handles.append(
                handles
            )

            continue

        actual_educt_smiles = (
            educt_smiles[
                reactant_index
            ]
        )

        if not actual_educt_smiles:

            reactant_handles.append(
                handles
            )

            continue

        actual_educt = Chem.MolFromSmiles(
            actual_educt_smiles
        )

        if actual_educt is None:

            logger.warning(
                "Could not parse example educt for "
                f"reactant position {reactant_index}: "
                f"{actual_educt_smiles}"
            )

            reactant_handles.append(
                handles
            )

            continue

        # ----------------------------------------------------------
        # Match the reaction-template query against the actual
        # example educt.
        #
        # IMPORTANT:
        #
        # actual_educt.GetSubstructMatches(
        #     reactant_template
        # )
        #
        # is deliberate. The actual molecule is the target and
        # the reaction template is the query.
        # ----------------------------------------------------------

        try:

            template_matches = (
                actual_educt.GetSubstructMatches(
                    reactant_template
                )
            )

        except Exception as exc:

            logger.debug(
                "Could not map reaction template onto "
                f"example educt for reactant position "
                f"{reactant_index}: {exc}"
            )

            reactant_handles.append(
                handles
            )

            continue

        if not template_matches:

            logger.debug(
                "Reaction template did not match example "
                f"educt for reactant position "
                f"{reactant_index}: "
                f"{actual_educt_smiles}"
            )

            reactant_handles.append(
                handles
            )

            continue

        # ----------------------------------------------------------
        # Compile the commercial handle definitions once for this
        # reactant.
        # ----------------------------------------------------------

        handle_patterns = {}

        for handle_name, handle_smarts in (
            REACTION_HANDLE_SMARTS.items()
        ):

            pattern = Chem.MolFromSmarts(
                handle_smarts
            )

            if pattern is not None:
                handle_patterns[
                    handle_name
                ] = pattern

        # ----------------------------------------------------------
        # Test every valid template mapping.
        #
        # There may be multiple matches of a generic reaction
        # template against the example educt.
        # ----------------------------------------------------------

        for template_match in template_matches:

            reaction_centre_atom_indices = (
                _get_reaction_centre_atom_indices(
                    reactant_template,
                    template_match,
                    reaction_centre,
                )
            )

            if not reaction_centre_atom_indices:
                continue

            # ------------------------------------------------------
            # Search the actual molecule for each defined handle.
            # ------------------------------------------------------

            for handle_name, pattern in (
                handle_patterns.items()
            ):

                try:

                    handle_matches = (
                        actual_educt.GetSubstructMatches(
                            pattern
                        )
                    )

                except Exception as exc:

                    logger.debug(
                        "Could not match reaction handle "
                        f"'{handle_name}' against example "
                        f"educt: {exc}"
                    )

                    continue

                for handle_match in handle_matches:

                    if _handle_match_is_reaction_centre_related(
                        actual_educt,
                        handle_match,
                        reaction_centre_atom_indices,
                    ):

                        handles.append(
                            handle_name
                        )

                        break

        handles = sorted(
            set(handles)
        )

        reactant_handles.append(
            handles
        )

    return reactant_handles

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

    precursor_mols = []

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

        precursor_mols.append(
            mol
        )

    # ----------------------------------------------------------
    # Rule-specific chemical sanity checks.
    # ----------------------------------------------------------

    reaction_rule = rule_data.get(
        "description",
        ""
    )

    if reaction_rule.startswith(
        "Literature reaction: "
    ):
        reaction_rule = (
            reaction_rule.replace(
                "Literature reaction: ",
                "",
                1,
            )
        )

    if not validate_generated_reaction_precursors(
        reaction_rule,
        precursor_mols,
    ):
        return False, (
            f"failed rule-specific precursor validation: "
            f"{reaction_rule}"
        )

    return True, None

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

    # reaction_handles = infer_reaction_handles(
    #     rule_data["smirks"]
    # )

    reaction_handles = infer_reaction_handles(
        rule_data
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

            logger.debug(
                f"  {rule_name}: rejected precursor set: "
                f"{rejection_reason}"
            )

            continue

        # reaction_handles = infer_reaction_handles(
        #     rule_data["smirks"]
        # )

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

def validate_generated_reaction_precursors(
    reaction_rule,
    precursor_mols,
):
    """
    Apply lightweight chemical sanity checks to generated
    retrosynthetic precursors.

    Returns:
        True if the precursor set passes the rule-specific checks.
        False otherwise.

    These checks are intentionally conservative and are not intended
    to replace reaction-feasibility prediction.
    """

    if not precursor_mols:
        return False

    if reaction_rule == "Suzuki":

        # A Suzuki disconnection should generate exactly two
        # coupling partners.
        if len(precursor_mols) != 2:
            return False

        boron_precursors = []

        boronic_acid_query = Chem.MolFromSmarts(
            "[B;X3]([O;H1])[O;H1]"
        )

        boronate_ester_query = Chem.MolFromSmarts(
            "[B;X3]([O;X2][#6])[O;X2][#6]"
        )

        for precursor in precursor_mols:

            if precursor is None:
                continue

            if (
                precursor.HasSubstructMatch(
                    boronic_acid_query
                )
                or precursor.HasSubstructMatch(
                    boronate_ester_query
                )
            ):
                boron_precursors.append(
                    precursor
                )

        # Suzuki should have one boron-containing partner.
        if len(boron_precursors) != 1:
            return False

        if not _is_sensible_suzuki_boron_candidate(
            boron_precursors[0]
        ):
            logger.debug(
                "Rejected Suzuki precursor set because "
                "the boron-containing precursor does not "
                "contain a sensible C(sp2)-B coupling motif."
            )
            return False

    return True