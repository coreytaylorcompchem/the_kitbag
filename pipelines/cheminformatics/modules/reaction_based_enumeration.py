import csv
# import re
# import gc
import os
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
from pipeline.parallel_runner import ParallelWorkflowRunner


logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def apply_reactions_to_scaffold(scaffold_row, reaction_info_list, reactant_mols, max_products=500):
    """
    reaction_info_list: list of dicts with keys 'smarts' and 'name'
    """
    products = []
    scaffold_smi = scaffold_row["smiles"]
    scaffold_mol = Chem.MolFromSmiles(scaffold_smi)
    if scaffold_mol is None:
        logger.warning(f"Invalid scaffold SMILES: {scaffold_smi}")
        return pd.DataFrame()

    logger.debug(f"Processing scaffold: {scaffold_smi}")

    for rxn_idx, rxn_info in enumerate(reaction_info_list):
        smarts = rxn_info['smarts']
        rxn_name = rxn_info['name']
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None:
                logger.warning(f"[{rxn_name}] Could not parse SMIRKS: {smarts}")
                continue
        except Exception as e:
            logger.warning(f"[{rxn_name}] Failed to parse SMIRKS: {smarts} | Error: {e}")
            continue

        product_count = 0
        for reactant_mol in reactant_mols:
            try:
                product_sets = rxn.RunReactants((scaffold_mol, reactant_mol))
                if not product_sets:
                    logger.debug(f"[{rxn_name}] No products generated for reactant.")

                for product_set in product_sets:
                    for product in product_set:
                        try:
                            Chem.SanitizeMol(product)
                            smi = Chem.MolToSmiles(product, canonical=True)
                            reactant_smi = Chem.MolToSmiles(reactant_mol, canonical=True)
                            products.append({
                                "smiles": smi,
                                "scaffold": scaffold_smi,
                                "reactant": reactant_smi,
                                "rxn_name": rxn_name
                            })
                            product_count += 1
                            if product_count >= max_products:
                                break
                        except Exception as e:
                            logger.debug(f"[{rxn_name}] Sanitization failed: {e}")
                    if product_count >= max_products:
                        break
            except Exception as e:
                logger.debug(f"[{rxn_name}] Reaction failed for scaffold/reactant: {e}")
                continue

        logger.debug(f"[{rxn_name}] Total products for scaffold: {product_count}")

    logger.info(f"[{scaffold_smi}] Total products generated: {len(products)}")
    return pd.DataFrame(products)


class ScaffoldReactantRunner:
    def __init__(self, reaction_info_list, max_products):
        self.reaction_info_list = reaction_info_list
        self.max_products = max_products

    def __call__(self, task):
        scaffold_row = task["scaffold"]
        reactant_chunk_file = task["reactant_chunk_file"]

        reactant_chunk_df = pd.read_csv(reactant_chunk_file)
        reactant_mols = []
        for smi in reactant_chunk_df["smiles"].dropna():
            mol = Chem.MolFromSmiles(smi)
            if mol:
                reactant_mols.append(mol)
            else:
                logger.warning(f"Invalid reactant SMILES in chunk: {smi}")

        if not reactant_mols:
            logger.warning(f"No valid reactants in chunk {reactant_chunk_file}. Skipping.")
            return pd.DataFrame()

        return apply_reactions_to_scaffold(
            scaffold_row=scaffold_row,
            reaction_info_list=self.reaction_info_list,
            reactant_mols=reactant_mols,
            max_products=self.max_products
        )


@register_task("reaction_based_enumeration", category="Library generation", description="Enumerate reaction-based product library from scaffold and fragment inputs.")
def reaction_based_enumeration(config: dict, data: dict = None) -> dict:
    params = config.get("reaction_based_enumeration", {})

    rxn_file = Path(params.get("input_file_rxns"))
    scaffold_file = Path(params.get("input_file_scaffolds"))
    reactant_file = Path(params.get("input_file_reactants"))
    max_products = params.get("max_products_per_scaffold", 500)

    logger.info(f"Loading input files...")

    if not rxn_file.exists() or not scaffold_file.exists() or not reactant_file.exists():
        raise FileNotFoundError("Missing one or more input files.")


    output_cfg = config.get("output", {})
    output_dir = Path(output_cfg.get("directory", "outputs/reaction_based_enumeration"))
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_filename = output_cfg.get("filename", "reaction_library.csv")

    rxn_df = pd.read_csv(rxn_file)
    scaffold_df = pd.read_csv(scaffold_file)
    reactant_df = pd.read_csv(reactant_file)

    logger.info(f"Loaded {len(rxn_df)} reactions, {len(scaffold_df)} scaffolds, {len(reactant_df)} reactants.")

    # Get reaction info list: list of dicts {'smarts': ..., 'name': ...}
    if "smirks" in rxn_df.columns:
        smarts_col = "smirks"
    elif "smarts" in rxn_df.columns:
        smarts_col = "smarts"
    else:
        raise ValueError("Reaction file must have either 'smirks' or 'smarts' column.")

    if "rxn_name" not in rxn_df.columns:
        raise ValueError("Reaction file must have 'rxn_name' column.")

    reaction_info_list = []
    for _, row in rxn_df.iterrows():
        smarts = row[smarts_col]
        name = row["rxn_name"]
        if pd.notna(smarts) and pd.notna(name):
            reaction_info_list.append({"smarts": smarts, "name": name})

    if "smiles" not in scaffold_df.columns:
        raise ValueError("Scaffold file must have 'smiles' column.")
    if "smiles" not in reactant_df.columns:
        raise ValueError("Reactant file must have 'smiles' column.")

    #  Chunk reactants for speed (not necessary usually, but good to have in case we have tons of reactants)
    chunk_size = params.get("chunk_size", 1000)
    chunk_dir = output_dir / "reactant_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_files = []
    for i in range(0, len(reactant_df), chunk_size):
        chunk = reactant_df.iloc[i:i + chunk_size].copy()
        chunk_file = chunk_dir / f"reactants_chunk_{i // chunk_size}.csv"
        chunk.to_csv(chunk_file, index=False, quoting=csv.QUOTE_ALL)
        chunk_files.append(str(chunk_file))
    logger.info(f"Split reactants into {len(chunk_files)} chunks.")

    scaffold_inputs = scaffold_df.to_dict(orient="records")

    task_list = []
    for chunk_file in chunk_files:
        for scaffold_row in scaffold_inputs:
            task_list.append({
                "scaffold": scaffold_row,
                "reactant_chunk_file": chunk_file
            })

    scaffold_runner = ScaffoldReactantRunner(
        reaction_info_list=reaction_info_list,
        max_products=max_products
    )

    def safe_filename(smi):
        import re
        return re.sub(r'[^A-Za-z0-9_.-]+', '_', smi)[:40]

    for i, task in enumerate(task_list):
        smi = task["scaffold"]["smiles"]
        task["safe_smi"] = safe_filename(smi)
        task["chunk_idx"] = i

    runner = ParallelWorkflowRunner(
        workflow_func=scaffold_runner,
        config={"manual": task_list, **config},
        input_key="manual",
        output_key="scaffold_product",
        filename_pattern="products_{safe_smi}_{chunk_idx}.csv",
        combined_filename=combined_filename,
        output_dir=str(output_dir),
        use_multiprocessing=True,
        reserve_cpus=4
    )

    logger.info(f"Starting parallel reaction-based enumeration over {len(task_list)} scaffold/reactant chunk pairs...")
    combined_df = runner.run()
    logger.info(f"Finished. Combined product count: {len(combined_df)}")

    # --- Cleanup chunk files ---
    cleanup_flag = config.get("cleanup", False)
    if cleanup_flag:
        logger.info(f"Cleanup enabled. Removing {len(chunk_files)} chunk files...")
        for chunk_file in chunk_files:
            try:
                os.remove(chunk_file)
            except Exception as e:
                logger.warning(f"Failed to remove chunk file {chunk_file}: {e}")

        # Try to remove the chunk directory if it's empty
        try:
            chunk_dir.rmdir()
            logger.info(f"Removed empty chunk directory: {chunk_dir}")
        except OSError:
            logger.warning(f"Chunk directory {chunk_dir} not empty or could not be removed.")

    return {"df": combined_df}
