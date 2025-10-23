import os
import re
# import itertools

import matplotlib.pyplot as plt
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

from sklearn.model_selection import train_test_split

from collections import Counter
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from pipeline.task_registry import register_task

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def is_valid_smiles(smi):
    try:
        return Chem.MolFromSmiles(smi) is not None
    except:
        return False

def standardise_smiles(smi):

    fragment_remover = rdMolStandardize.FragmentRemover()
    uncharger = rdMolStandardize.Uncharger()

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    # Remove salts
    mol = fragment_remover.remove(mol)
    # Remove charged moieties
    mol = uncharger.uncharge(mol)
    return Chem.MolToSmiles(mol)

def tokenize_smiles(smiles):
    # Basic regex for tokenizing SMILES
    pattern =  "(\[[^\[\]]{1,6}\])" + \
               "|Br|Cl" + \
               "|Si|Se|Na|Li|Ca|Fe|Zn|Cu" + \
               "|[B-Zb-z]" + \
               "|\d+" + \
               "|=|#|\(|\)|\.|:|\/|\\|\+|\-|\%|\@|\?"  # extended syntax
    regex = re.compile(pattern)
    tokens = regex.findall(smiles)
    return tokens

@register_task("load_preprocess_standardise_data", category="Molecular generation", description="Load CSV, prepare/standardise and generate data loaders for training.")
def load_data_and_standardise(config: dict, context: dict):
    """
    aaa
    """
    input_file = config.get("input_file")
    if input_file is None:
        raise ValueError("No input_file specified in config")

    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = pd.read_csv(input_file)
    if "smiles" not in df.columns or "standard_value" not in df.columns:
        raise ValueError("Input CSV must contain 'smiles' and 'standard_value' columns")
    
    df = df[df['smiles'].apply(is_valid_smiles)].reset_index(drop=True)
    
    logger.info(f"Loaded {len(df)} valid SMILES.")

    df['standardised_smiles'] = df['smiles'].apply(standardise_smiles)
    df = df[df['smiles'].notnull()].reset_index(drop=True)
    
    logger.info(f"After standardization: {len(df)} SMILES")

    token_counts = Counter()
    for smi in tqdm(df['standardised_smiles']):
        tokens = tokenize_smiles(smi)
        token_counts.update(tokens)

    # Add special tokens
    special_tokens = ['<pad>', '<bos>', '<eos>', '<unk>']
    vocab = special_tokens + sorted(token_counts.keys())
    token_to_idx = {token: idx for idx, token in enumerate(vocab)}
    idx_to_token = {idx: token for token, idx in token_to_idx.items()}
    vocab_size = len(vocab)

    logger.info(f"Vocab size: {vocab_size}")

    context.update({
        "token_to_idx": token_to_idx,
        "idx_to_token": idx_to_token,
        "vocab_size": vocab_size,
    })

    return {
        "token_to_idx": token_to_idx,
        "idx_to_token": idx_to_token,
        "vocab_size": vocab_size,
    }

# @register_task("prepare_dataloaders", category="Pose ranker", description="Normalise labels and create train/val DataLoaders")
# def prepare_dataloaders(config: dict, context: dict):

#     class SmilesDataset(Dataset):
#         def __init__(self, smiles_list, token_to_idx, max_len=128):
#             self.smiles_list = smiles_list
#             self.token_to_idx = token_to_idx
#             self.max_len = max_len

#         def __len__(self):
#             return len(self.smiles_list)

#         def __getitem__(self, idx):
#             smiles = self.smiles_list[idx]
#             tokens = tokenize_smiles(smiles)
#             tokens = ['<bos>'] + tokens + ['<eos>']
#             token_ids = [self.token_to_idx.get(tok, self.token_to_idx['<unk>']) for tok in tokens]

#             if len(token_ids) > self.max_len:
#                 token_ids = token_ids[:self.max_len]
#             else:
#                 token_ids += [self.token_to_idx['<pad>']] * (self.max_len - len(token_ids))

#             input_ids = torch.tensor(token_ids[:-1], dtype=torch.long)
#             target_ids = torch.tensor(token_ids[1:], dtype=torch.long)

#             return input_ids, target_ids

#     graphs = context.get("graphs")
#     train_smiles, val_smiles = train_test_split(df['standardised_smiles'].tolist(), test_size=0.1, random_state=42)

#     train_dataset = SmilesDataset(train_smiles, token_to_idx, max_len=128)
#     val_dataset = SmilesDataset(val_smiles, token_to_idx, max_len=128)

#     train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
#     val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False)

#     context.update({
#         "train_loader": train_loader,
#         "val_loader": val_loader,
#     })

#     return {
#         "train_loader": train_loader,
#         "val_loader": val_loader,
#     }
    