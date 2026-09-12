"""Sequence-keyed ESM embeddings (mean-pooled)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from mpatlas.paths import PROCESSED

OUT = PROCESSED / "esm2_650m.npz"
MODEL_ID = "facebook/esm2_t33_650M_UR50D"
MAX_LEN = 1022


def seq_key(seq: str) -> str:
    seq = "".join(aa for aa in seq.upper() if aa.isalpha())[:MAX_LEN]
    return hashlib.sha1(seq.encode("ascii", errors="ignore")).hexdigest()


class EmbedIndex:
    def __init__(self, path: Path | None = None):
        path = path or OUT
        data = np.load(path, allow_pickle=True)
        keys = [str(k) for k in data["keys"]]
        self._ix = {k: i for i, k in enumerate(keys)}
        self.X = np.asarray(data["X"], dtype=np.float32)
        self.model = str(data["model"][0]) if "model" in data.files else MODEL_ID

    def has(self, seq: str) -> bool:
        return seq_key(seq) in self._ix

    def get(self, seq: str) -> np.ndarray | None:
        i = self._ix.get(seq_key(seq))
        if i is None:
            return None
        return self.X[i]

    def coverage(self, sequences: list[str]) -> float:
        if not sequences:
            return 0.0
        return sum(self.has(s) for s in sequences) / len(sequences)

    def matrix(self, sequences: list[str]) -> np.ndarray | None:
        rows = []
        for s in sequences:
            v = self.get(s)
            if v is None:
                return None
            rows.append(v)
        return np.vstack(rows)

    def __len__(self) -> int:
        return len(self._ix)


def load(path: Path | None = None) -> EmbedIndex | None:
    path = path or OUT
    if not path.exists():
        return None
    return EmbedIndex(path)
