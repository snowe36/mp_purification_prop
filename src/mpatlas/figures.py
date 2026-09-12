from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mpatlas.palette import FACE, GRID, MUTED, NEG, POS, SERIES, TEXT
from mpatlas.paths import FIGURES, ensure_dirs
from mpatlas.wrap import wrap_score

plt.rcParams.update(
    {
        "font.size": 11,
        "axes.facecolor": FACE,
        "figure.facecolor": FACE,
        "axes.edgecolor": TEXT,
        "axes.labelcolor": TEXT,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "text.color": TEXT,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _save(fig: plt.Figure, name: str, dest: Path | None = None) -> Path:
    ensure_dirs()
    path = (dest or FIGURES) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def leakage_auc(df: pd.DataFrame, dest: Path | None = None) -> Path:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    splits = list(dict.fromkeys(df["split"]))
    models = list(dict.fromkeys(df["model"]))
    x = np.arange(len(splits))
    width = 0.35 if len(models) > 1 else 0.5
    for i, model in enumerate(models):
        sub = df[df["model"] == model]
        vals = [float(sub.loc[sub["split"] == s, "roc_auc"].mean()) for s in splits]
        ax.bar(
            x + (i - 0.5) * width,
            vals,
            width=width,
            label=model,
            color=SERIES[i % len(SERIES)],
            edgecolor=TEXT,
            linewidth=0.4,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=20 if len(splits) > 4 else 0, ha="right" if len(splits) > 4 else "center")
    ax.set_ylabel("ROC-AUC")
    ax.set_ylim(0.45, 1.02)
    ax.axhline(0.5, color=GRID, lw=1)
    ax.legend(frameon=False)
    ax.set_title("Random-split skill is not the same as a center or organism holdout")
    fig.set_size_inches(8.0 if len(splits) > 4 else 6.4, 4.2)
    return _save(fig, "fig_leakage_auc.png", dest)


def funnel(counts: dict[str, int], dest: Path | None = None) -> Path:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    order = [k for k in ("selected", "cloned", "expressed", "soluble", "purified", "in pdb") if k in counts]
    if not order:
        order = list(counts)
    vals = [counts[k] for k in order]
    colors = [POS if i == 0 else SERIES[min(i, len(SERIES) - 1)] for i in range(len(order))]
    ax.barh(order[::-1], vals[::-1], color=colors[::-1], edgecolor=TEXT, linewidth=0.4)
    ax.set_xlabel("targets")
    ax.set_title("Fluorescence is early in the funnel")
    for i, v in enumerate(vals[::-1]):
        ax.text(v, i, f" {v:,}", va="center", color=MUTED, fontsize=9)
    return _save(fig, "fig_funnel.png", dest)


def wrap_panel(sequences: list[str], names: list[str] | None = None, dest: Path | None = None) -> Path:
    names = names or [f"p{i}" for i in range(len(sequences))]
    pairs = list(reversed(list(zip([wrap_score(s) for s in sequences], names, strict=True))))
    scores = [p[0] for p in pairs]
    names = [p[1] for p in pairs]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    y = np.arange(len(scores))
    colors = [POS if s.amenable else NEG for s in scores]
    ax.barh(y, [s.score for s in scores], color=colors, edgecolor=TEXT, linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("WRAP amenability (heuristic)")
    ax.set_xlim(0, 1)
    ax.set_title("B-failures still need a stabilizer path")
    return _save(fig, "fig_wrap.png", dest)


def architecture_map(df: pd.DataFrame, dest: Path | None = None) -> Path:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    if "architecture" not in df.columns:
        return _save(fig, "fig_architecture.png", dest)
    counts = df["architecture"].fillna("unknown").value_counts()
    ax.bar(counts.index.astype(str), counts.values, color=SERIES[0], edgecolor=TEXT, linewidth=0.4)
    ax.set_ylabel("proteins")
    ax.tick_params(axis="x", rotation=20)
    ax.set_title("Solved membrane proteins are not a random TM draw")
    return _save(fig, "fig_architecture.png", dest)
