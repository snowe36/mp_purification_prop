from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def random_split(
    n: int, y: np.ndarray, test_size: float = 0.3, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    pos = idx[y == 1]
    neg = idx[y == 0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    n_te_pos = max(1, int(round(len(pos) * test_size)))
    n_te_neg = max(1, int(round(len(neg) * test_size)))
    te = np.concatenate([pos[:n_te_pos], neg[:n_te_neg]])
    tr = np.concatenate([pos[n_te_pos:], neg[n_te_neg:]])
    rng.shuffle(tr)
    rng.shuffle(te)
    return tr, te


def hamming_clusters(sequences: Sequence[str], max_frac: float = 0.12) -> np.ndarray:
    """Greedy clusters: assign to first representative within max_frac Hamming."""
    if not sequences:
        return np.array([], dtype=int)
    length = len(sequences[0])
    if any(len(s) != length for s in sequences):
        labels = np.empty(len(sequences), dtype=int)
        reps: list[str] = []
        cutoff = max_frac * length
        for i, seq in enumerate(sequences):
            assigned = False
            for c, rep in enumerate(reps):
                if len(seq) != len(rep):
                    continue
                dist = sum(a != b for a, b in zip(seq, rep, strict=True))
                if dist <= cutoff:
                    labels[i] = c
                    assigned = True
                    break
            if not assigned:
                labels[i] = len(reps)
                reps.append(seq)
        return labels
    X = np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8).reshape(
        len(sequences), length
    )
    cutoff = max_frac * length
    labels = np.empty(len(sequences), dtype=int)
    reps: list[int] = []
    for i in range(len(sequences)):
        assigned = False
        for c, ri in enumerate(reps):
            if int((X[i] != X[ri]).sum()) <= cutoff:
                labels[i] = c
                assigned = True
                break
        if not assigned:
            labels[i] = len(reps)
            reps.append(i)
    return labels


def cluster_split(
    cluster_ids: np.ndarray, test_size: float = 0.3, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    uniq = np.unique(cluster_ids)
    rng.shuffle(uniq)
    n_te = max(1, int(round(len(uniq) * test_size)))
    te_c = set(uniq[:n_te].tolist())
    idx = np.arange(len(cluster_ids))
    te = idx[np.isin(cluster_ids, list(te_c))]
    tr = idx[~np.isin(cluster_ids, list(te_c))]
    return tr, te


def group_split(
    groups: Sequence[str], test_groups: Sequence[str]
) -> tuple[np.ndarray, np.ndarray]:
    te_set = set(test_groups)
    idx = np.arange(len(groups))
    mask = np.array([g in te_set for g in groups])
    return idx[~mask], idx[mask]


def split_ok(y: np.ndarray, tr: np.ndarray, te: np.ndarray, min_n: int = 30) -> bool:
    if len(tr) == 0 or len(te) == 0:
        return False
    ytr, yte = y[tr], y[te]
    return (
        int((ytr == 1).sum()) >= min_n
        and int((ytr == 0).sum()) >= min_n
        and int((yte == 1).sum()) >= min_n
        and int((yte == 0).sum()) >= min_n
    )
