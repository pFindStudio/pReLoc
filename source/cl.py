import torch
import torch.nn as nn
import torch.utils.data as data
import pandas as pd
import numpy as np
import os
import random
import gc
from tqdm import tqdm
import logging
import esm
import argparse
from typing import List
import sys

def resource_path(relative_path):
    try:
        # PyInstaller stores the bundled resource directory in _MEIPASS.
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.dirname(os.path.abspath(__file__))
    
    return os.path.join(base_path, relative_path)

# ====================== Fix PyTorch 2.6 Weight Safety ======================
torch.serialization.add_safe_globals([argparse.Namespace])

# ====================== Prediction Config ======================
PREDICT_CONFIG = {
    "model_path": resource_path("best_model.pth"),
    "device": "auto",
    "esm_path": resource_path("esm2_t33_650M_UR50D.pt"),
    "seed": 42,
    "embedding_dim": 128,
    "temperature": 0.1,
    "batch_size": 1,
    "num_workers": 0,
}

# -------------------------- Constants --------------------------
AMINO_ACIDS = [
    "<pad>",
    "G", "A", "S", "P", "V", "T", "C", "L", "I", "N", "D", "Q",
    "K", "E", "M", "H", "F", "R", "Y", "W",
]
AA_TO_IDX = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}
NUM_AMINO_ACIDS = len(AMINO_ACIDS)
CONTEXT_WINDOW = 3

# ====================== Load ESM2 650M Locally ======================
esm_model = None
esm_alphabet = None
esm_batch_converter = None
_esm_path = None


def load_esm(path, device):
    global esm_model, esm_alphabet, esm_batch_converter, _esm_path
    path = os.path.abspath(path)
    if esm_model is None or _esm_path != path:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"CL ESM checkpoint not found: {path}")
        esm_model, esm_alphabet = esm.pretrained.load_model_and_alphabet_local(path)
        esm_model.eval()
        esm_batch_converter = esm_alphabet.get_batch_converter()
        _esm_path = path
    esm_model = esm_model.to(device)


def get_aa_context(pep_seq: str, target_pos: float) -> List[str]:
    target_pos = int(target_pos)
    pep_len = len(pep_seq)
    target_idx = target_pos - 1
    context = []

    for i in range(CONTEXT_WINDOW, 0, -1):
        pos = target_idx - i
        context.append(pep_seq[pos] if pos >= 0 else "<pad>")

    context.append(pep_seq[target_idx] if 0 <= target_idx < pep_len else "<pad>")

    for i in range(1, CONTEXT_WINDOW + 1):
        pos = target_idx + i
        context.append(pep_seq[pos] if pos < pep_len else "<pad>")

    return context


def get_esm_embedding(aa_seq: List[str]) -> np.ndarray:
    ESM_EMBED_DIM = esm_model.embed_dim
    seq_str = "".join(aa_seq)

    valid_chars = set(esm_alphabet.tok_to_idx.keys())
    filtered_seq = [c for c in seq_str if c in valid_chars]
    if not filtered_seq:
        return np.zeros(ESM_EMBED_DIM, dtype=np.float32)
    filtered_str = "".join(filtered_seq)

    data = [("context_seq", filtered_str)]
    batch_labels, batch_strs, batch_tokens = esm_batch_converter(data)
    model_device = next(esm_model.parameters()).device
    batch_tokens = batch_tokens.to(model_device)

    with torch.no_grad():
        results = esm_model(batch_tokens, repr_layers=[33])

    token_reprs = results["representations"][33]
    seq_emb = token_reprs[0, 1:-1, :].mean(dim=0).cpu().numpy()

    return seq_emb.astype(np.float32)


# -------------------------- Data Preprocessing --------------------------
def load_raw_data(df):
    if "pep" in df.columns and "peptide_sequence" not in df.columns:
        df["peptide_sequence"] = df["pep"].copy()

    required_cols = [
        "peptide_sequence", "bm25_score", "ion_ratio", "intensity_ratio",
        "can_bm25_score", "can_ion_ratio", "can_intensity_ratio",
        "aa_position", "aa_type", "pep_length",
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")

    return df


def process_peptide_features(df):
    df["aa_context"] = df.apply(
        lambda x: get_aa_context(x["peptide_sequence"], x["aa_position"]), axis=1
    )

    # FIXED: Removed conflicting progress_apply that caused pandas crash
    print("Extracting ESM embeddings...")
    embedding_list = []
    for ctx in tqdm(df["aa_context"], desc="Extracting ESM embeddings"):
        embedding_list.append(get_esm_embedding(ctx))
    df["esm_embedding"] = embedding_list

    df["aa_pos_norm"] = df["aa_position"] / pd.to_numeric(df["pep_length"])

    df = df[
        (df["pep_length"] > 0)
        & (df["aa_position"] >= 1)
        & (df["aa_position"] <= df["pep_length"])
        & (df["aa_context"].apply(len) == 2 * CONTEXT_WINDOW + 1)
    ].reset_index(drop=True)

    return df


# -------------------------- Dataset Class --------------------------
class ContrastiveDataset(data.Dataset):
    def __init__(self, df, has_label=False):
        self.df = df
        self.esm_embed_dim = esm_model.embed_dim
        self.has_label = has_label

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]

        pep_features = np.array(
            [row["bm25_score"], row["ion_ratio"], row["intensity_ratio"]],
            dtype=np.float32,
        )

        can_pep_features = np.array(
            [row["can_bm25_score"], row["can_ion_ratio"], row["can_intensity_ratio"]],
            dtype=np.float32,
        )

        extra_features = np.array(
            [row["aa_position"], row["aa_pos_norm"], row["pep_length"]],
            dtype=np.float32,
        )

        esm_feat = row["esm_embedding"].astype(np.float32)

        result = {
            "pep_features": torch.FloatTensor(pep_features),
            "can_pep_features": torch.FloatTensor(can_pep_features),
            "extra_features": torch.FloatTensor(extra_features),
            "esm_feat": torch.FloatTensor(esm_feat),
            "index": idx,
        }

        if self.has_label and "label" in self.df.columns:
            result["label"] = torch.FloatTensor([row["label"]])

        return result


# -------------------------- Model Structure --------------------------
class FeatureEncoder(nn.Module):
    def __init__(self, embedding_dim=128, esm_embed_dim=1280):
        super().__init__()
        self.main_feat_fc = nn.Sequential(
            nn.Linear(3, 128),
            nn.BatchNorm1d(128),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(128, 64),
            nn.BatchNorm1d(64),
            nn.GELU(),
            nn.Dropout(0.1),
        )

        self.esm_extra_fc = nn.Sequential(
            nn.Linear(3 + esm_embed_dim, 256),
            nn.BatchNorm1d(256),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.GELU(),
            nn.Dropout(0.1),
        )

        self.fusion = nn.Sequential(
            nn.Linear(64 + 128, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(256, embedding_dim),
            nn.LayerNorm(embedding_dim),
        )

    def forward(self, main_features, extra_features, esm_feat):
        x_main = self.main_feat_fc(main_features)
        esm_extra_combined = torch.cat([extra_features, esm_feat], dim=1)
        x_esm_extra = self.esm_extra_fc(esm_extra_combined)
        x = torch.cat([x_main, x_esm_extra], dim=1)
        x = self.fusion(x)
        return nn.functional.normalize(x, p=2, dim=1)


class ContrastiveModel(nn.Module):
    def __init__(self, feature_encoder, embedding_dim=128, temperature=0.1, num_classes=2):
        super().__init__()
        self.feature_encoder = feature_encoder
        self.temperature = temperature
        self.classifier = nn.Sequential(
            nn.Linear(embedding_dim, 64),
            nn.LayerNorm(64),
            nn.GELU(),
            nn.Linear(64, num_classes),
        )

    def forward(self, main_features, extra_features, esm_feat):
        emb = self.feature_encoder(main_features, extra_features, esm_feat)
        logits = self.classifier(emb)
        return emb, logits

    def info_nce_loss(self, pos_emb, neg_emb):
        batch_size = pos_emb.size(0)
        pos_sim = torch.bmm(pos_emb.unsqueeze(1), pos_emb.unsqueeze(2)).squeeze(-1)
        neg_sim = torch.bmm(pos_emb.unsqueeze(1), neg_emb.unsqueeze(2)).squeeze(-1)
        logits = torch.cat([pos_sim, neg_sim], dim=1) / self.temperature
        labels = torch.zeros(batch_size, dtype=torch.long, device=pos_emb.device)
        return nn.CrossEntropyLoss()(logits, labels)

    def supervised_loss(self, logits, labels):
        return nn.CrossEntropyLoss()(logits, labels.squeeze(1).long())


# -------------------------- Prediction --------------------------
def predict(model, dataloader, device):
    model.eval()
    all_probs = []
    all_indices = []

    with torch.no_grad():
        pbar = tqdm(dataloader, desc="Predicting, please wait...")
        for batch in pbar:
            pep_feat = batch["pep_features"].to(device)
            extra_feat = batch["extra_features"].to(device)
            esm_feat = batch["esm_feat"].to(device)
            indices = batch["index"]

            _, pos_logits = model(pep_feat, extra_feat, esm_feat)
            pos_probs = torch.softmax(pos_logits, dim=1)[:, 1].cpu().numpy()

            all_probs.extend(pos_probs)
            all_indices.extend(indices.numpy())

            del pep_feat, extra_feat, esm_feat, pos_logits
            torch.cuda.empty_cache() if device.type == "cuda" else None

    return all_probs, all_indices


# -------------------------- Main Prediction API --------------------------
def main_predict(test_df, config=PREDICT_CONFIG):
    config = {**PREDICT_CONFIG, **config}
    requested_device = config["device"].lower()
    if requested_device == "auto":
        requested_device = "cuda:0" if torch.cuda.is_available() else "cpu"
    device = torch.device(requested_device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise ValueError("CL requested CUDA, but CUDA is unavailable; use clDevice=auto or cpu")
    if not os.path.isfile(config["model_path"]):
        raise FileNotFoundError(f"CL checkpoint not found: {config['model_path']}")
    if config["batch_size"] < 1:
        raise ValueError("CL batch size must be at least 1")

    test_df = test_df.copy().reset_index(drop=True)
    if test_df.empty:
        test_df["pred_prob"] = pd.Series(dtype=float)
        return test_df
    invalid = (
        (test_df["pep_length"] <= 0)
        | (test_df["aa_position"] < 1)
        | (test_df["aa_position"] > test_df["pep_length"])
    )
    if invalid.any():
        raise ValueError("CL features contain invalid peptide lengths or site positions")

    random.seed(config["seed"])
    np.random.seed(config["seed"])
    torch.manual_seed(config["seed"])
    if torch.cuda.is_available():
        torch.cuda.manual_seed(config["seed"])

    load_esm(config["esm_path"], device)

    test_raw_df = load_raw_data(test_df)
    test_df = process_peptide_features(test_raw_df)
    test_df = test_df.reset_index(drop=True)

    test_dataset = ContrastiveDataset(test_df, has_label="label" in test_df.columns)
    test_loader = data.DataLoader(
        test_dataset,
        batch_size=config["batch_size"],
        shuffle=False,
        num_workers=config["num_workers"],
        pin_memory=True if device.type == "cuda" else False,
    )

    esm_embed_dim = esm_model.embed_dim
    feature_encoder = FeatureEncoder(
        embedding_dim=config["embedding_dim"], esm_embed_dim=esm_embed_dim
    )
    model = ContrastiveModel(
        feature_encoder=feature_encoder,
        embedding_dim=config["embedding_dim"],
        temperature=config["temperature"],
    ).to(device)

    checkpoint = torch.load(config["model_path"], map_location=device, weights_only=False)
    if "model_state_dict" in checkpoint:
        model.load_state_dict(checkpoint["model_state_dict"])
    else:
        model.load_state_dict(checkpoint)

    all_probs, all_indices = predict(model, test_loader, device)

    results_df = test_raw_df.copy()
    results_df["pred_prob"] = np.nan
    results_df.loc[all_indices, "pred_prob"] = all_probs

    print("\nSuccessfully predicted site scores.")

    return results_df
