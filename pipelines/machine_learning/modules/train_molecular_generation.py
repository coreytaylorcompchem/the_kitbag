import os
import re
import json

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

from models.molecular_generation import SmilesTransformerEncDec
from models.utils.train_utils import train_epoch_enc_dec
from models.utils.generation_utils import generate_smiles_beam
from models.utils.eval_utils import evaluate_generated_smiles, eval_epoch_enc_dec
from models.utils.chemistry_utils import calculate_properties
from models.utils.losses_and_schedulers import LabelSmoothingLoss, WarmupInverseSqrtScheduler

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

@register_task("load_preprocess_standardise_data", category="Molecular generation")
def load_data_and_standardise(config: dict, context: dict):
    input_file = config.get("input_file")
    max_len = config.get("max_len", 128)
    save_vocab = config.get("save_vocab", True)

    # --- UPDATED: get output config from top-level context if available ---
    # If 'output' is in context, use that; else fallback to config or default
    output_cfg = context.get("output") or config.get("output", {})
    output_dir = Path(output_cfg.get("directory", "outputs/molecular_generation"))
    overwrite = output_cfg.get("overwrite", True)

    output_dir.mkdir(parents=True, exist_ok=True)

    processed_path = output_dir / "processed_smiles.csv"
    vocab_path = output_dir / "vocab.json"

    if processed_path.exists() and not overwrite:
        logger.warning(f"{processed_path} exists and overwrite=False — skipping processing.")
        context.update({
            "processed_data_path": str(processed_path),
            "vocab_path": str(vocab_path),
            "output_dir": str(output_dir),
        })
        return context

    # --- main data load + preprocessing ---
    df = pd.read_csv(input_file)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")

    df = df[df['smiles'].apply(is_valid_smiles)].reset_index(drop=True)
    df['standardised_smiles'] = df['smiles'].apply(standardise_smiles)
    df = df[df['standardised_smiles'].notnull()].reset_index(drop=True)

    logger.info(f"After standardization: {len(df)} SMILES")

    token_counts = Counter()
    for smi in tqdm(df['standardised_smiles']):
        tokens = tokenize_smiles(smi)
        token_counts.update(tokens)

    special_tokens = ['<pad>', '<bos>', '<eos>', '<unk>']
    vocab = special_tokens + sorted(token_counts.keys())
    token_to_idx = {tok: i for i, tok in enumerate(vocab)}
    idx_to_token = {idx: token for token, idx in token_to_idx.items()}

    df.to_csv(processed_path, index=False)

    if save_vocab:
        import json
        with open(vocab_path, "w") as f:
            json.dump(token_to_idx, f, indent=2)

    logger.info(f"Saved processed data to: {processed_path}")
    logger.info(f"Saved vocab to: {vocab_path}")

    context.update({
        "processed_data_path": str(processed_path),
        "vocab_path": str(vocab_path),
        "vocab_size": len(vocab),
        "output_dir": str(output_dir),
        "token_to_idx": token_to_idx,
        "idx_to_token": idx_to_token
    })

    return context

class SmilesDataset(Dataset):
    def __init__(self, smiles_list, token_to_idx, max_len=128):
        self.smiles_list = smiles_list
        self.token_to_idx = token_to_idx
        self.max_len = max_len

    def __len__(self):
        return len(self.smiles_list)

    def __getitem__(self, idx):
        smiles = self.smiles_list[idx]
        tokens = tokenize_smiles(smiles)
        tokens = ['<bos>'] + tokens + ['<eos>']
        token_ids = [self.token_to_idx.get(tok, self.token_to_idx['<unk>']) for tok in tokens]
        if len(token_ids) > self.max_len:
            token_ids = token_ids[:self.max_len]
        else:
            token_ids += [self.token_to_idx['<pad>']] * (self.max_len - len(token_ids))
        return (
            torch.tensor(token_ids[:-1], dtype=torch.long),
            torch.tensor(token_ids[1:], dtype=torch.long)
        )

@register_task(
    "split_and_create_dataloaders",
    category="Molecular generation",
    description="Create train/val dataloaders for SMILES transformer training."
)
def split_and_create_dataloaders(config: dict, context: dict):
    """
    Split the standardized SMILES into train/val sets and create PyTorch DataLoaders.

    Args:
        config (dict): task configuration from YAML
        context (dict): shared context between workflow tasks
    """

    # === Retrieve context from previous task ===
    token_to_idx = context.get("token_to_idx")
    # idx_to_token = context.get("idx_to_token")

    output_cfg = context.get("output", {})
    output_dir = Path(output_cfg.get("directory", "outputs/molecular_generation"))
    overwrite = output_cfg.get("overwrite", True)

    processed_data_path = output_dir / "processed_smiles.csv"
    if not processed_data_path.exists():
        raise FileNotFoundError(f"Processed data not found: {processed_data_path}")

    df = pd.read_csv(processed_data_path)
    if "standardised_smiles" not in df.columns:
        raise ValueError("Processed CSV must contain 'standardised_smiles' column")

    smiles_list = df["standardised_smiles"].tolist()

    # === Config parameters ===
    test_size = config.get("test_size", 0.1)
    random_state = config.get("random_state", 42)
    batch_size = config.get("batch_size", 64)
    max_len = config.get("max_len", 128)

    # === Split data ===
    train_smiles, val_smiles = train_test_split(smiles_list, test_size=test_size, random_state=random_state)

    # === Save splits for reproducibility ===
    split_dir = output_dir / "splits"
    split_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"smiles": train_smiles}).to_csv(split_dir / "train_smiles.csv", index=False)
    pd.DataFrame({"smiles": val_smiles}).to_csv(split_dir / "val_smiles.csv", index=False)
    logger.info(f"Saved train/val splits to: {split_dir}")

    # === Create datasets/loaders ===
    train_dataset = SmilesDataset(train_smiles, token_to_idx, max_len=max_len)
    val_dataset = SmilesDataset(val_smiles, token_to_idx, max_len=max_len)

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    logger.info(f"Created dataloaders — Train: {len(train_loader.dataset)}, Val: {len(val_loader.dataset)}")

    # === Update context for downstream tasks ===
    context.update({
        "train_loader": train_loader,
        "val_loader": val_loader,
        "train_smiles": train_smiles,
        "val_smiles": val_smiles,
        # "token_to_idx": token_to_idx,
        # "idx_to_token": idx_to_token
    })

    return {
        "train_loader": train_loader,
        "val_loader": val_loader,
        "train_smiles": train_smiles,
        "val_smiles": val_smiles,
    }

@register_task("build_model", category="Molecular generation")
def build_model(config: dict, context: dict):
    vocab_size = context.get("vocab_size")
    if vocab_size is None:
        raise ValueError("vocab_size not found in context — did you run load_preprocess_standardise_data?")

    token_to_idx = context.get("token_to_idx")
    if token_to_idx is None:
        raise ValueError("token_to_idx not found in context — did you pass it from previous task?")

    device = torch.device(config.get("device", "cuda" if torch.cuda.is_available() else "cpu"))

    pad_idx = token_to_idx.get("<pad>", 0)

    model = SmilesTransformerEncDec(
        vocab_size=vocab_size,
        pad_idx=pad_idx,
        embed_dim=config.get("embed_dim", 512),
        num_heads=config.get("num_heads", 8),
        hidden_dim=config.get("hidden_dim", 2048),
        num_layers=config.get("num_layers", 6),
        max_len=config.get("max_len", 128),
        dropout=config.get("dropout", 0.1)
    ).to(device)

    context.update({
        "model": model,
        "device": str(device),
    })

    return context


@register_task(
    "train_molecular_generation",
    category="Molecular generation",
    description="Train SMILES Transformer model using configuration-driven hyperparameters."
)
def train_model(config: dict, context: dict):
    """
    Train a SMILES Transformer model based on YAML configuration.

    Args:
        config: Dict of hyperparameters from the YAML task block.
        context: Shared dictionary containing objects passed from previous tasks (model, dataloaders, vocab, etc.)
    """
    # def smiles_to_tensor(smiles, token_to_idx, max_len=128, device='cpu'):
    #     tokens = tokenize_smiles(smiles)  # your SMILES tokenizer
    #     tokens = ['<bos>'] + tokens + ['<eos>']
    #     token_ids = [token_to_idx.get(tok, token_to_idx['<unk>']) for tok in tokens]
        
    #     if len(token_ids) > max_len:
    #         token_ids = token_ids[:max_len]
    #     else:
    #         token_ids += [token_to_idx['<pad>']] * (max_len - len(token_ids))
        
    #     return torch.tensor(token_ids, dtype=torch.long).unsqueeze(0).to(device)  # shape: [1, max_len]

    logger.info("Starting model training...")

    # === Retrieve required context objects ===
    model = context.get("model")
    train_loader = context.get("train_loader")
    val_loader = context.get("val_loader")
    token_to_idx = context.get("token_to_idx")
    idx_to_token = context.get("idx_to_token")
    train_smiles = context.get("train_smiles")

    missing = [k for k, v in {
        "model": model,
        "train_loader": train_loader,
        "val_loader": val_loader,
        "token_to_idx": token_to_idx,
        "idx_to_token": idx_to_token,
        "train_smiles": train_smiles,
    }.items() if v is None]
    
    if missing:
        raise ValueError(f"Missing one or more objects in context ({', '.join(missing)}).")

    # === Task-specific config parameters from YAML ===
    num_epochs = config.get("num_epochs", 5)
    batch_size = config.get("batch_size", 64)
    learning_rate = float(config.get("learning_rate", 1e-3))
    warmup_steps = config.get("warmup_steps", 4000)
    loss_fn_name = config.get("loss_fn", "LabelSmoothingLoss")
    optimizer_name = config.get("optimizer", "Adam")
    scheduler_name = config.get("scheduler", "WarmupInverseSqrtScheduler")
    save_dir = config.get("save_dir", "output/models/")
    save_frequency = config.get("save_frequency", 1)

    # === Read global output directory from YAML ===
    output_cfg = context.get("output", {})
    base_output_dir = output_cfg.get("directory", "outputs/molecular_generation")
    overwrite = output_cfg.get("overwrite", True)

    # Final model save directory
    save_dir = os.path.join(base_output_dir, save_dir)
    os.makedirs(save_dir, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    vocab_size = len(token_to_idx)

    logger.info(f"Training on device: {device}")
    logger.info(f"Model save directory: {save_dir}")

    # === Initialize optimizer, loss, scheduler ===
    if optimizer_name.lower() == "adam":
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    else:
        raise ValueError(f"Unsupported optimizer: {optimizer_name}")

    if loss_fn_name == "LabelSmoothingLoss":
        criterion = LabelSmoothingLoss(classes=vocab_size, smoothing=0.1, ignore_index=token_to_idx["<pad>"])
    else:
        raise ValueError(f"Unsupported loss function: {loss_fn_name}")

    if scheduler_name == "WarmupInverseSqrtScheduler":
        scheduler = WarmupInverseSqrtScheduler(optimizer, warmup_steps=warmup_steps)
    else:
        scheduler = None

    # === Initialize metric history ===
    history = {
        "train_loss": [], "val_loss": [],
        "train_perplexity": [], "val_perplexity": [],
        "validity": [], "uniqueness": [], "novelty": [],
        "avg_length": [],
        "qed_mean": [], "mol_wt_mean": [], "logp_mean": [], "tpsa_mean": []
    }

    # === Training loop ===
    for epoch in range(num_epochs):
        logger.info(f">>>>>>>>>> Epoch {epoch + 1}/{num_epochs} <<<<<<<<<<")

        train_loss = train_epoch_enc_dec(model, train_loader, optimizer, scheduler, criterion, device)
        val_loss = eval_epoch_enc_dec(model, val_loader, criterion, device)

        train_perplexity = np.exp(train_loss)
        val_perplexity = np.exp(val_loss)

        logger.info(f"Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")
        logger.info(f"Perplexity — Train: {train_perplexity:.2f} | Val: {val_perplexity:.2f}")

        # === Generate and evaluate molecules ===
        dummy_src = torch.tensor([[token_to_idx["<bos>"]]], device=device)

        samples = generate_smiles_beam(
            model, dummy_src,
            token_to_idx=token_to_idx, idx_to_token=idx_to_token,
            n_samples=200, beam_width=5, max_len=100, temperature=1.0, device=device
        )

        # import random
        # seed_smiles = random.sample(train_smiles, 20)
        # generated_samples = []

        ### NEED TO DEBUG FEEDING IN REAL DATA, SEEMS TO BE SILENTLY FAILING SOMETIMES.
        
        # for smi in seed_smiles:
        #     src_tensor = smiles_to_tensor(smi, token_to_idx, max_len=128, device=device)
        #     samples = generate_smiles_beam(
        #         model, src_tensor,
        #         token_to_idx=token_to_idx,
        #         idx_to_token=idx_to_token,
        #         n_samples=5,       # number of beams per seed molecule
        #         beam_width=5,
        #         max_len=100,
        #         temperature=1.0,
        #         device=device
        #     )
        #     generated_samples.extend(samples)

        validity, uniqueness, novelty = evaluate_generated_smiles(samples, train_smiles)
        avg_length = np.mean([len(s) for s in samples]) if samples else 0
        properties = calculate_properties(samples)

        logger.info(
            f"Validity: {validity:.2%}, Uniqueness: {uniqueness:.2%}, Novelty: {novelty:.2%}, "
            f"QED={properties['qed_mean']:.3f}, MW={properties['mol_wt_mean']:.1f}, "
            f"LogP={properties['logp_mean']:.2f}, TPSA={properties['tpsa_mean']:.2f}"
        )

        # === Update history ===
        for k, v in [
            ("train_loss", train_loss), ("val_loss", val_loss),
            ("train_perplexity", train_perplexity), ("val_perplexity", val_perplexity),
            ("validity", validity), ("uniqueness", uniqueness), ("novelty", novelty),
            ("avg_length", avg_length),
            ("qed_mean", properties["qed_mean"]), ("mol_wt_mean", properties["mol_wt_mean"]),
            ("logp_mean", properties["logp_mean"]), ("tpsa_mean", properties["tpsa_mean"])
        ]:
            history[k].append(v)

        # === Save checkpoint ===
        if (epoch + 1) % save_frequency == 0:
            ckpt_path = os.path.join(save_dir, f"model_epoch_{epoch + 1}.pt")
            torch.save(model.state_dict(), ckpt_path)
            logger.info(f"Saved checkpoint: {ckpt_path}")

    # Save history
    history_path = os.path.join(save_dir, "training_history.npy")
    np.save(history_path, history)
    logger.info(f"Training history saved: {history_path}")

    # Update context
    context.update({
        "model": model,
        "training_history": history
    })

    return {"model": model, "training_history": history}

@register_task(
    "evaluate_model",
    category="Molecular generation",
    description="Plot training history and evaluation metrics."
)
def evaluate_model(config: dict, context: dict):
    """
    Plot training and evaluation metrics from history.

    Args:
        config: Dict of evaluation config parameters
        context: Dict containing training history and global output directory
    """

    history = context.get("training_history")
    if history is None:
        raise ValueError("training_history not found in context. Run train_model first.")

    # Get global output directory from context, default to current dir if missing
    global_output_dir = context.get("output_directory", ".")

    # Get local save_dir from task config, default to "evaluation"
    local_save_dir = config.get("save_dir", "evaluation")

    # Combine paths
    save_dir = os.path.join(global_output_dir, local_save_dir)
    os.makedirs(save_dir, exist_ok=True)

    epochs = range(1, len(history["train_loss"]) + 1)

    # Plot 1: Loss, Perplexity, Validity, Uniqueness, Novelty, Avg Length
    plt.figure(figsize=(12, 8))

    plt.subplot(2, 2, 1)
    plt.plot(epochs, history["train_loss"], label="Train Loss")
    plt.plot(epochs, history["val_loss"], label="Val Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(epochs, history["train_perplexity"], label="Train Perplexity")
    plt.plot(epochs, history["val_perplexity"], label="Val Perplexity")
    plt.xlabel("Epoch")
    plt.ylabel("Perplexity")
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(epochs, history["validity"], label="Validity")
    plt.plot(epochs, history["uniqueness"], label="Uniqueness")
    plt.plot(epochs, history["novelty"], label="Novelty")
    plt.xlabel("Epoch")
    plt.ylabel("Fraction")
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.plot(epochs, history["avg_length"], label="Avg SMILES Length")
    plt.xlabel("Epoch")
    plt.ylabel("Length")
    plt.legend()

    plt.tight_layout()
    plot1_path = os.path.join(save_dir, "training_and_quality_metrics.png")
    plt.savefig(plot1_path)
    plt.close()
    print(f"Saved plot to {plot1_path}")

    # Plot 2: Drug-like properties
    plt.figure(figsize=(12, 8))

    plt.subplot(2, 2, 1)
    plt.plot(epochs, history["qed_mean"], label="QED Mean", color="green")
    plt.xlabel("Epoch")
    plt.ylabel("QED")
    plt.ylim(0, 1)
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(epochs, history["mol_wt_mean"], label="Molecular Weight", color="purple")
    plt.xlabel("Epoch")
    plt.ylabel("Molecular Weight")
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(epochs, history["logp_mean"], label="LogP Mean", color="orange")
    plt.xlabel("Epoch")
    plt.ylabel("LogP")
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.plot(epochs, history["tpsa_mean"], label="TPSA Mean", color="red")
    plt.xlabel("Epoch")
    plt.ylabel("TPSA")
    plt.legend()

    plt.tight_layout()
    plot2_path = os.path.join(save_dir, "drug_like_properties.png")
    plt.savefig(plot2_path)
    plt.close()
    print(f"Saved plot to {plot2_path}")

    # Update context with actual save_dir and plot paths
    context.update({
        "evaluation_plots_dir": save_dir,
        "plot_training_metrics": plot1_path,
        "plot_drug_properties": plot2_path,
    })

    return {"evaluation_plots_dir": save_dir}
