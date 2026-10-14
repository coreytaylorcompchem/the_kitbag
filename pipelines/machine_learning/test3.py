from rdkit import Chem
import pandas as pd

from modules.utils.retrosynth_rxn_rules import (
    COMPILED_REACTION_RULES,
)


df = pd.read_csv(
    "projects/retrosynthetic_predictor/input/retrosynthesis_candidates.csv"
)


for _, row in df.iterrows():

    target_id = row["molecule_id"]
    smiles = row["SMILES"]

    mol = Chem.MolFromSmiles(smiles)

    print()
    print("=" * 80)
    print(f"TARGET: {target_id}")
    print(smiles)
    print("=" * 80)

    if mol is None:
        print("INVALID SMILES")
        continue

    for rule_name, reaction in COMPILED_REACTION_RULES.items():

        try:
            matches = reaction.RunReactants(
                (mol,)
            )
        except Exception as exc:
            print(
                f"{rule_name}: ERROR: {exc}"
            )
            continue

        print(
            f"{rule_name}: "
            f"{len(matches)} raw match(es)"
        )

        for product_set in matches[:3]:

            print(
                "   ",
                [
                    Chem.MolToSmiles(product)
                    for product in product_set
                ]
            )
