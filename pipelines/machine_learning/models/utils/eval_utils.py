import torch
from rdkit import Chem
import numpy as np
from tqdm import tqdm

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def evaluate_generated_smiles(samples, train_smiles, max_samples=1000, max_smi_len=200):
    if not samples:
        logger.warning("No samples provided to evaluate_generated_smiles.")
        return 0.0, 0.0, 0.0

    # Defensive limit on number of samples processed
    if len(samples) > max_samples:
        logger.warning(f"Too many samples ({len(samples)}), truncating to {max_samples}.")
        samples = samples[:max_samples]

    valid = []
    for i, smi in enumerate(samples):
        if not isinstance(smi, str):
            logger.warning(f"Sample {i} is not a string: {smi}. Skipping.")
            valid.append(False)
            continue

        if len(smi) > max_smi_len:
            logger.warning(f"Sample {i} is too long (> {max_smi_len} chars). Skipping.")
            valid.append(False)
            continue

        try:
            mol = Chem.MolFromSmiles(smi)
        except Exception as e:
            logger.warning(f"Exception on MolFromSmiles for sample {i}: {e}")
            valid.append(False)
            continue

        valid.append(mol is not None)

    validity = np.mean(valid) if valid else 0.0

    unique_samples = set(samples)
    uniqueness = len(unique_samples) / len(samples) if samples else 0.0
    novelty = len([s for s in unique_samples if s not in train_smiles]) / len(unique_samples) if unique_samples else 0.0
    return validity, uniqueness, novelty

def eval_epoch_enc_dec(model, val_loader, criterion, device):
    model.eval()
    total_loss = 0
    val_iter = tqdm(val_loader, desc="Validation", leave=False)

    with torch.no_grad():
        for src, tgt in val_iter:
            src, tgt = src.to(device), tgt.to(device)
            decoder_input = tgt[:, :-1]
            target_output = tgt[:, 1:]

            logits = model(src, decoder_input)
            loss = criterion(logits.reshape(-1, logits.size(-1)), target_output.reshape(-1))

            total_loss += loss.item()
            val_iter.set_postfix(val_loss=loss.item())

    return total_loss / len(val_loader)