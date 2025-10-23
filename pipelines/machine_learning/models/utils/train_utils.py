import torch

from tqdm import tqdm
import numpy as np

class EarlyStopping:
    def __init__(self, patience=10, min_delta=1e-4):
        self.patience = patience
        self.min_delta = min_delta
        self.best_loss = float("inf")
        self.counter = 0
        self.should_stop = False

    def step(self, val_loss):
        if val_loss < self.best_loss - self.min_delta:
            self.best_loss = val_loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.should_stop = True

def train_epoch(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        pred = model(batch)
        loss = criterion(pred, batch.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)

def eval_epoch(model, loader, criterion, device):
    model.eval()
    total_loss = 0
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            pred = model(batch)
            loss = criterion(pred, batch.y)
            total_loss += loss.item()
    return total_loss / len(loader)

def train_epoch_enc_dec(model, train_loader, optimizer, scheduler, criterion, device):
    model.train()
    total_loss = 0
    train_iter = tqdm(train_loader, desc="Training", leave=False)

    for src, tgt in train_iter:
        src, tgt = src.to(device), tgt.to(device)
        optimizer.zero_grad()

        decoder_input = tgt[:, :-1]
        target_output = tgt[:, 1:]

        logits = model(src, decoder_input)
        loss = criterion(logits.reshape(-1, logits.size(-1)), target_output.reshape(-1))
        loss.backward()

        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        scheduler.step()

        total_loss += loss.item()
        train_iter.set_postfix(loss=loss.item(), lr=optimizer.param_groups[0]['lr'])

    return total_loss / len(train_loader)

