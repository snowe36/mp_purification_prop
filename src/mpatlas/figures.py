from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mpatlas.palette import FACE, GRID, MUTED, NEG, PALETTE, POS, SERIES, TEXT
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


def embed_auc(df: pd.DataFrame, dest: Path | None = None) -> Path:
    sub = df[df["model"].eq("rf")].copy() if "model" in df.columns else df
    splits = list(dict.fromkeys(sub["split"]))
    feats = list(dict.fromkeys(sub["features"]))
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    x = np.arange(len(splits))
    width = 0.8 / max(len(feats), 1)
    for i, feat in enumerate(feats):
        vals = []
        for s in splits:
            hit = sub[(sub["split"] == s) & (sub["features"] == feat)]
            vals.append(float(hit["roc_auc"].mean()) if len(hit) else float("nan"))
        ax.bar(
            x + (i - (len(feats) - 1) / 2) * width,
            vals,
            width=width,
            label=feat,
            color=SERIES[i % len(SERIES)],
            edgecolor=TEXT,
            linewidth=0.4,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=20, ha="right")
    ax.set_ylabel("ROC-AUC")
    ax.set_ylim(0.45, 1.02)
    ax.axhline(0.5, color=GRID, lw=1)
    ax.legend(frameon=False, title="features")
    ax.set_title("ESM vs composition on the same leakage splits")
    return _save(fig, "fig_embed_auc.png", dest)


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


ACTION_COLORS = {
    "standard": PALETTE["sage"],
    "wrap_rescue": PALETTE["lavender"],
    "redesign": PALETTE["mustard"],
    "deprioritize": PALETTE["coral"],
}


def playbook_coefs(df: pd.DataFrame, dest: Path | None = None, n: int = 12) -> Path:
    sub = df.head(n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    colors = [POS if v else NEG for v in sub["helps_purify"]]
    ax.barh(sub["plain"], sub["coef"], color=colors, edgecolor=TEXT, linewidth=0.4)
    ax.axvline(0, color=GRID, lw=1)
    ax.set_xlabel("logreg weight (standardized)")
    ax.set_title("Sage tracks purified; coral tracks expressed-but-not-purified")
    return _save(fig, "fig_playbook_coefs.png", dest)


def playbook_actions(df: pd.DataFrame, dest: Path | None = None) -> Path:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    order = [a for a in ("standard", "wrap_rescue", "redesign", "deprioritize") if a in set(df["action"])]
    if "source" in df.columns and df["source"].nunique() > 1:
        counts = df.groupby(["source", "action"]).size().unstack(fill_value=0)
        counts = counts.reindex(columns=order, fill_value=0)
        bottom = np.zeros(len(counts))
        x = np.arange(len(counts))
        for action in order:
            vals = counts[action].to_numpy()
            ax.bar(x, vals, bottom=bottom, color=ACTION_COLORS[action], edgecolor=TEXT, linewidth=0.4, label=action)
            bottom = bottom + vals
        ax.set_xticks(x)
        ax.set_xticklabels(counts.index.astype(str))
        ax.legend(frameon=False, fontsize=8)
    else:
        counts = df["action"].value_counts().reindex(order)
        ax.bar(counts.index.astype(str), counts.values, color=[ACTION_COLORS[a] for a in counts.index], edgecolor=TEXT, linewidth=0.4)
        ax.tick_params(axis="x", rotation=15)
    ax.set_ylabel("proteins")
    ax.set_title("Recommended next step — not a purification guarantee")
    return _save(fig, "fig_playbook_actions.png", dest)


def playbook_scatter(df: pd.DataFrame, dest: Path | None = None) -> Path:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for action, sub in df.groupby("action"):
        ax.scatter(
            sub["p_purify"],
            sub["wrap_score"],
            s=18,
            alpha=0.7,
            color=ACTION_COLORS.get(str(action), SERIES[0]),
            edgecolors=TEXT,
            linewidths=0.2,
            label=str(action),
        )
    ax.axvline(0.55, color=GRID, lw=1)
    ax.set_xlabel("P(purified | expressed)")
    ax.set_ylabel("WRAP amenability")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("High purify score → standard path; else WRAP or redesign")
    return _save(fig, "fig_playbook_scatter.png", dest)


def detergents_bar(counts: dict[str, int], dest: Path | None = None, title: str | None = None) -> Path:
    items = [(k, v) for k, v in counts.items() if k and k != "None"][:12]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    if not items:
        ax.set_title(title or "GPCRdb solubilization detergents")
        return _save(fig, "fig_detergents.png", dest)
    names, vals = zip(*items, strict=True)
    ax.barh(list(names)[::-1], list(vals)[::-1], color=SERIES[0], edgecolor=TEXT, linewidth=0.4)
    ax.set_xlabel("GPCRdb solubilization rows")
    ax.set_title(title or "Starting detergent is DDM, not a purify-vs-fail score")
    return _save(fig, "fig_detergents.png", dest)
