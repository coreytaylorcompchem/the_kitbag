import torch
import torch.nn as nn

class SmilesTransformerEncDec(nn.Module):
    def __init__(self, vocab_size, pad_idx, embed_dim=512, num_heads=8, hidden_dim=2048,
                 num_layers=6, max_len=128, dropout=0.1):
        super().__init__()
        
        self.pad_idx = pad_idx
        
        self.token_embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=pad_idx)
        self.position_embedding = nn.Embedding(max_len, embed_dim)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embed_dim, nhead=num_heads,
            dim_feedforward=hidden_dim, dropout=dropout
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        decoder_layer = nn.TransformerDecoderLayer(
            d_model=embed_dim, nhead=num_heads,
            dim_feedforward=hidden_dim, dropout=dropout
        )
        self.decoder = nn.TransformerDecoder(decoder_layer, num_layers=num_layers)

        self.fc_out = nn.Linear(embed_dim, vocab_size)
        self.dropout = nn.Dropout(dropout)
        self.max_len = max_len

    def forward(self, src, tgt, pad_idx=None):
        if pad_idx is None:
            pad_idx = self.pad_idx
        batch_size, src_len = src.size()
        batch_size, tgt_len = tgt.size()

        # --- Encoder embeddings + positions ---
        src_positions = torch.arange(src_len, device=src.device).unsqueeze(0).expand(batch_size, -1)
        src_emb = self.token_embedding(src) + self.position_embedding(src_positions)
        src_emb = self.dropout(src_emb).permute(1, 0, 2)
        src_key_padding_mask = (src == pad_idx)

        memory = self.encoder(src_emb, src_key_padding_mask=src_key_padding_mask)

        # --- Decoder embeddings + positions ---
        tgt_positions = torch.arange(tgt_len, device=tgt.device).unsqueeze(0).expand(batch_size, -1)
        tgt_emb = self.token_embedding(tgt) + self.position_embedding(tgt_positions)
        tgt_emb = self.dropout(tgt_emb).permute(1, 0, 2)
        tgt_key_padding_mask = (tgt == pad_idx)
        tgt_mask = nn.Transformer.generate_square_subsequent_mask(tgt_len).to(tgt.device)

        output = self.decoder(
            tgt_emb, memory,
            tgt_mask=tgt_mask,
            tgt_key_padding_mask=tgt_key_padding_mask,
            memory_key_padding_mask=src_key_padding_mask
        )

        output = output.permute(1, 0, 2)
        logits = self.fc_out(output)
        return logits
