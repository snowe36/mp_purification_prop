from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from mpatlas.features import matrix
from mpatlas.paths import PROCESSED, REPORTS, ensure_dirs
from mpatlas.splits import cluster_split, group_split, hamming_clusters, random_split, split_ok


def _stack(sequences: list[str], embeds, mode: str) -> np.ndarray:
    if mode == "esm":
        X = embeds.matrix(sequences)
        if X is None:
            raise ValueError("missing embeddings")
        return X
    Xc = matrix(sequences)
    if mode != "both" or embeds is None:
        return Xc
    Xe = embeds.matrix(sequences)
    if Xe is None:
        raise ValueError("missing embeddings")
    return np.hstack([Xc, Xe])


def _models(seed: int = 0) -> dict:
    return {
        "logreg": Pipeline(
            [
                ("sc", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(max_iter=400, class_weight="balanced", random_state=seed),
                ),
            ]
        ),
        "rf": RandomForestClassifier(
            n_estimators=200,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=seed,
            n_jobs=1,
        ),
    }


def evaluate_split(
    X: np.ndarray, y: np.ndarray, tr: np.ndarray, te: np.ndarray, split: str
) -> list[dict]:
    rows = []
    for name, model in _models().items():
        model.fit(X[tr], y[tr])
        if hasattr(model, "predict_proba"):
            p = model.predict_proba(X[te])[:, 1]
        else:
            p = model.decision_function(X[te])
        yte = y[te]
        auc = float(roc_auc_score(yte, p)) if len(np.unique(yte)) > 1 else float("nan")
        ap = float(average_precision_score(yte, p)) if len(np.unique(yte)) > 1 else float("nan")
        rows.append(
            {
                "split": split,
                "model": name,
                "n_train": int(len(tr)),
                "n_test": int(len(te)),
                "pos_rate_train": float(y[tr].mean()),
                "pos_rate_test": float(yte.mean()),
                "roc_auc": auc,
                "pr_auc": ap,
            }
        )
    return rows


def run_curnow(sequences: list[str], labels: np.ndarray, embeds=None, mode: str = "compose") -> pd.DataFrame:
    X = _stack(sequences, embeds, mode)
    y = labels.astype(int)
    rows = []
    tr, te = random_split(len(y), y)
    if split_ok(y, tr, te, min_n=20):
        rows.extend(evaluate_split(X, y, tr, te, "random"))
    clusters = hamming_clusters(sequences)
    tr, te = cluster_split(clusters)
    if split_ok(y, tr, te, min_n=20):
        rows.extend(evaluate_split(X, y, tr, te, "cluster"))
    df = pd.DataFrame(rows)
    if not df.empty:
        df["features"] = mode
    return df


def subsample_idx(y: np.ndarray, max_n: int = 4000, seed: int = 0) -> np.ndarray:
    if len(y) <= max_n:
        return np.arange(len(y))
    rng = np.random.default_rng(seed)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    n_pos = min(len(pos), max_n // 2)
    n_neg = min(len(neg), max_n - n_pos)
    take = np.concatenate(
        [rng.choice(pos, n_pos, replace=False), rng.choice(neg, n_neg, replace=False)]
    )
    rng.shuffle(take)
    return take


def run_binary(
    sequences: list[str],
    labels: np.ndarray,
    groups: list[str] | None,
    split_name: str,
    holdout: list[str] | None = None,
    min_n: int = 30,
    embeds=None,
    mode: str = "compose",
) -> pd.DataFrame:
    X = _stack(sequences, embeds, mode)
    y = labels.astype(int)
    rows = []
    tr, te = random_split(len(y), y)
    if split_ok(y, tr, te, min_n=min_n):
        rows.extend(evaluate_split(X, y, tr, te, "random"))
    if groups is not None:
        te_g = holdout
        if te_g is None:
            uniq = sorted(set(groups))
            if len(uniq) >= 2:
                te_g = [uniq[-1]]
        if te_g:
            tr, te = group_split(groups, te_g)
            if split_ok(y, tr, te, min_n=min_n):
                rows.extend(evaluate_split(X, y, tr, te, split_name))
    df = pd.DataFrame(rows)
    if not df.empty:
        df["features"] = mode
    return df


def write_results(df: pd.DataFrame, name: str) -> None:
    ensure_dirs()
    for folder in (PROCESSED, REPORTS):
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / f"{name}.csv"
        df.to_csv(path, index=False)
    (PROCESSED / f"{name}.json").write_text(df.to_json(orient="records", indent=2))
