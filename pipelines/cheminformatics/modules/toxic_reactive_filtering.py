from rdkit import Chem
from rdkit.Chem import rdFMCS
from pipeline.task_registry import register_task
from pathlib import Path
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def match_smarts(mol, smarts_list):
    toxic_matches = []
    reactive_matches = []
    
    for category, smarts in smarts_list.items():
        patt = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(patt):
            if category == 'PAINS':
                toxic_matches.append(smarts)
            elif category == 'REACTIVE':
                reactive_matches.append(smarts)
    
    return toxic_matches, reactive_matches

@register_task("toxic_reactive_filtering", category="Filtering", description="Toxic and reactive SMARTS matching for molecules")
def toxic_reactive_filtering(config, data=None):
    input_file = config.get("input_file")
    toxic_reactive_file = config.get("toxic_reactive_file")
    
    if input_file is None or toxic_reactive_file is None:
        raise ValueError("Both 'input_file' and 'toxic_reactive_file' must be specified in config")

    # Load input file and toxic/reactive SMARTS
    input_path = Path(input_file)
    toxic_reactive_path = Path(toxic_reactive_file)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")
    if not toxic_reactive_path.exists():
        raise FileNotFoundError(f"Toxic/reactive SMARTS file does not exist: {toxic_reactive_file}")
    
    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")
    
    # Read the toxic/reactive SMARTS
    toxic_reactive_df = pd.read_csv(toxic_reactive_file)
    smarts_list = {
        row['toxic_category']: row['smarts'] for _, row in toxic_reactive_df.iterrows()
    }

    logger.debug(f"Loaded {len(smarts_list)} toxic/reactive SMARTS patterns")

    # Prepare columns for toxic/reactive flags and matching SMARTS
    toxic_flags = []
    reactive_flags = []
    toxic_matches_list = []
    reactive_matches_list = []

    for smi in df["smiles"].dropna():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            toxic_flags.append("N")
            reactive_flags.append("N")
            toxic_matches_list.append("")
            reactive_matches_list.append("")
            continue

        toxic_matches, reactive_matches = match_smarts(mol, smarts_list)
        
        toxic_flags.append("Y" if toxic_matches else "N")
        reactive_flags.append("Y" if reactive_matches else "N")
        toxic_matches_list.append(",".join(toxic_matches))
        reactive_matches_list.append(",".join(reactive_matches))

    # Add new columns to dataframe
    df["toxic_flag"] = toxic_flags
    df["reactive_flag"] = reactive_flags
    df["toxic_matches"] = toxic_matches_list
    df["reactive_matches"] = reactive_matches_list

    # Return dataframe with additional columns
    output_dir = Path(config.get("output", {}).get("directory", "outputs/toxic_reactive"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename = output_dir / config.get("output", {}).get("filename", "filtered_toxic_reactive.csv")
    df.to_csv(output_filename, index=False)

    logger.debug(f"Saved filtered dataframe to {output_filename}")

    return {
        "mols": df['smiles'].tolist(),
        "df": df
    }
