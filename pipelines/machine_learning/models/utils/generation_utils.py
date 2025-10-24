import heapq
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

def beam_search_decoding(model, src_input, token_to_idx, idx_to_token,
                         beam_width=5, max_len=100, temperature=1.0, device='cpu'):
    model.eval()
    src_input = src_input.to(device)
    assert src_input.size(0) == 1, "Beam search currently supports batch size 1."

    with torch.no_grad():
        src_len = src_input.size(1)
        src_positions = torch.arange(src_len, device=device).unsqueeze(0)
        src_emb = model.token_embedding(src_input) + model.position_embedding(src_positions)
        src_emb = src_emb.permute(1, 0, 2)
        src_key_padding_mask = (src_input == token_to_idx['<pad>'])
        memory = model.encoder(src_emb, src_key_padding_mask=src_key_padding_mask)

        beams = [(0.0, [token_to_idx['<bos>']])]
        for _ in range(max_len):
            candidates = []
            for neg_log_prob, seq in beams:
                tgt_input = torch.tensor(seq, device=device).unsqueeze(0)
                tgt_len = tgt_input.size(1)
                tgt_positions = torch.arange(tgt_len, device=device).unsqueeze(0)
                tgt_emb = model.token_embedding(tgt_input) + model.position_embedding(tgt_positions)
                tgt_emb = tgt_emb.permute(1, 0, 2)
                tgt_mask = nn.Transformer.generate_square_subsequent_mask(tgt_len).to(device)

                output = model.decoder(tgt_emb, memory, tgt_mask=tgt_mask)
                output = output.permute(1, 0, 2)
                logits = model.fc_out(output)
                next_token_logits = logits[0, -1, :] / temperature
                probs = torch.softmax(next_token_logits, dim=-1)
                top_probs, top_indices = torch.topk(probs, beam_width)

                for prob, idx in zip(top_probs.tolist(), top_indices.tolist()):
                    new_seq = seq + [idx]
                    new_neg_log_prob = neg_log_prob - np.log(prob + 1e-9)
                    candidates.append((new_neg_log_prob, new_seq))

            beams = heapq.nsmallest(beam_width, candidates, key=lambda x: x[0])
            if any(seq[-1] == token_to_idx['<eos>'] for _, seq in beams):
                break

        best_seq = next((seq for _, seq in beams if seq[-1] == token_to_idx['<eos>']), beams[0][1])
        tokens = [idx_to_token[idx] for idx in best_seq[1:] if idx != token_to_idx['<eos>']]
        return "".join(tokens)


def generate_smiles_beam(model, src_input, token_to_idx, idx_to_token,
                         n_samples=5, beam_width=5, max_len=100, temperature=1.0, device='cpu'):
    return [
        beam_search_decoding(model, src_input, token_to_idx, idx_to_token, beam_width, max_len, temperature, device)
        for _ in range(n_samples)
    ]

def sample_smiles_topk(
    model,
    src_input,
    token_to_idx,
    idx_to_token,
    max_len=100,
    temperature=0.8,
    top_k=50,
    top_p=0.9,
    device='cpu'
):
    """
    Generate a SMILES string using top-k / top-p sampling from a trained encoder-decoder model.

    Args:
        model: Trained Transformer model (encoder-decoder)
        src_input: Tensor([1, seq_len]) source input tokens (SMILES)
        token_to_idx: Dict[str, int]
        idx_to_token: Dict[int, str]
        max_len: Maximum decoding length
        temperature: Softmax temperature (lower = more deterministic)
        top_k: Keep only top-k tokens per step
        top_p: Keep smallest set of tokens with cumulative prob >= top_p
        device: torch.device
    """
    model.eval()
    src_input = src_input.to(device)
    assert src_input.size(0) == 1, "Batch size must be 1 for sampling."

    with torch.no_grad():
        src_len = src_input.size(1)
        src_positions = torch.arange(src_len, device=device).unsqueeze(0)
        src_emb = model.token_embedding(src_input) + model.position_embedding(src_positions)
        src_emb = src_emb.permute(1, 0, 2)
        src_key_padding_mask = (src_input == token_to_idx['<pad>'])
        memory = model.encoder(src_emb, src_key_padding_mask=src_key_padding_mask)

        generated = [token_to_idx['<bos>']]

        for _ in range(max_len):
            tgt_input = torch.tensor(generated, device=device).unsqueeze(0)
            tgt_len = tgt_input.size(1)
            tgt_positions = torch.arange(tgt_len, device=device).unsqueeze(0)
            tgt_emb = model.token_embedding(tgt_input) + model.position_embedding(tgt_positions)
            tgt_emb = tgt_emb.permute(1, 0, 2)
            tgt_mask = torch.nn.Transformer.generate_square_subsequent_mask(tgt_len).to(device)

            output = model.decoder(tgt_emb, memory, tgt_mask=tgt_mask)
            logits = model.fc_out(output.permute(1, 0, 2))[:, -1, :] / temperature
            probs = F.softmax(logits, dim=-1).squeeze(0)

            # === Top-k filtering ===
            if top_k > 0:
                top_probs, top_idx = torch.topk(probs, top_k)
                probs = torch.zeros_like(probs).scatter_(0, top_idx, top_probs)
                probs = probs / probs.sum()

            # === Top-p (nucleus) filtering ===
            if 0 < top_p < 1.0:
                sorted_probs, sorted_idx = torch.sort(probs, descending=True)
                cumulative = torch.cumsum(sorted_probs, dim=0)
                mask = cumulative <= top_p
                if torch.any(mask):
                    cutoff = torch.where(mask)[0][-1]
                    keep_idx = sorted_idx[:cutoff + 1]
                    probs = torch.zeros_like(probs).scatter_(0, keep_idx, probs[keep_idx])
                    probs = probs / probs.sum()

            # Sample next token
            next_token = torch.multinomial(probs, 1).item()
            generated.append(next_token)

            # Stop if EOS
            if next_token == token_to_idx['<eos>']:
                break

        tokens = [idx_to_token[idx] for idx in generated[1:] if idx != token_to_idx['<eos>']]
        return "".join(tokens)
