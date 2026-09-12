"""What to try next so a membrane protein might purify.

Not buffer recipes. Not a claim that the action produces protein.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from mpatlas.detergents import prior_for
from mpatlas.features import FEATURE_NAMES, FeatureRow, featurize_seq, matrix
from mpatlas.wrap import WrapScore, wrap_score

ACTIONS = ("standard", "wrap_rescue", "redesign", "deprioritize")

PLAIN = {
    "length": "length",
    "kd_mean": "mean hydrophobicity",
    "net_charge": "net charge",
    "n_tm": "TM helix count",
    "tm_frac": "TM residue fraction",
    "longest_tm": "longest TM span",
    "longest_loop": "longest loop",
    "extra_frac": "extramembrane fraction",
    "n_sequons": "glycan sequons",
    "n_cys": "cysteine count",
    "c_term_extra": "C-terminus outside",
    "n_term_extra": "N-terminus outside",
}

PURIFY_CUT = 0.55


@dataclass(frozen=True)
class Advice:
    action: str
    p_purify: float
    wrap: float
    wrap_amenable: bool
    levers: tuple[str, ...]
    why: str


def levers_for(feat: FeatureRow) -> tuple[str, ...]:
    out: list[str] = []
    if feat.length > 900:
        out.append("shorten or split the construct")
    if feat.n_sequons >= 4:
        out.append("drop sequons or use a prokaryotic host")
    if feat.n_cys >= 8:
        out.append("reduce disulfides / skip eukaryotic folding")
    if feat.n_tm >= 1 and feat.n_term_extra < 0.5 and feat.c_term_extra < 0.5:
        out.append("add a fusion-able N- or C-terminus")
    if feat.longest_loop >= 250:
        out.append("trim a long extramembrane loop")
    return tuple(out)


def decide(p_purify: float, wrap: WrapScore, feat: FeatureRow) -> Advice:
    levers = levers_for(feat)
    if p_purify >= PURIFY_CUT:
        action = "standard"
        why = "sequence resembles TargetTrack purified (given expression)"
    elif wrap.amenable:
        action = "wrap_rescue"
        why = "standard purify looks unlikely; termini/TM belt look WRAP-amenable"
    elif levers:
        action = "redesign"
        why = "; ".join(levers)
    else:
        action = "deprioritize"
        why = "low purify score, not WRAP-amenable, no obvious construct lever"
    return Advice(action, float(p_purify), wrap.score, wrap.amenable, levers, why)


def fit_purify_logreg(sequences: list[str], y: np.ndarray) -> tuple[Pipeline, pd.DataFrame]:
    X = matrix(sequences)
    pipe = Pipeline(
        [
            ("sc", StandardScaler()),
            ("clf", LogisticRegression(max_iter=400, class_weight="balanced", random_state=0)),
        ]
    )
    pipe.fit(X, y.astype(int))
    coef = pipe.named_steps["clf"].coef_[0]
    table = pd.DataFrame(
        {
            "feature": list(FEATURE_NAMES),
            "coef": coef,
            "plain": [PLAIN.get(n, n.replace("aa_", "aa ")) for n in FEATURE_NAMES],
        }
    )
    table["helps_purify"] = table["coef"] > 0
    table["abs"] = table["coef"].abs()
    return pipe, table.sort_values("abs", ascending=False).reset_index(drop=True)


def score_sequences(pipe: Pipeline, sequences: list[str]) -> np.ndarray:
    return pipe.predict_proba(matrix(sequences))[:, 1]


def advise_one(pipe: Pipeline, seq: str) -> Advice:
    feat = featurize_seq(seq)
    p = float(pipe.predict_proba(matrix([seq]))[0, 1])
    return decide(p, wrap_score(seq, feat=feat), feat)


def advise_table(
    pipe: Pipeline,
    sequences: list[str],
    names: list[str],
    source: str = "",
    organisms: list[str | None] | None = None,
) -> pd.DataFrame:
    rows = []
    orgs = organisms or [None] * len(sequences)
    for seq, name, org in zip(sequences, names, orgs, strict=True):
        feat = featurize_seq(seq)
        p = float(pipe.predict_proba(matrix([seq]))[0, 1])
        adv = decide(p, wrap_score(seq, feat=feat), feat)
        recipe = prior_for(feat, organism=org)
        rows.append(
            {
                "accession": name,
                "source": source,
                "action": adv.action,
                "p_purify": round(adv.p_purify, 4),
                "wrap_score": round(adv.wrap, 3),
                "wrap_amenable": int(adv.wrap_amenable),
                "levers": "; ".join(adv.levers),
                "why": adv.why,
                "detergent_extract": recipe.extract,
                "detergent_polish": recipe.polish,
                "detergent_rescue": recipe.rescue,
                "detergent_avoid": recipe.avoid,
                "detergent_why": recipe.why,
            }
        )
    return pd.DataFrame(rows)
