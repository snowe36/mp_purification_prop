"""WRAP amenability heuristics. Not a design engine."""

from __future__ import annotations

from dataclasses import dataclass

from mpatlas.features import FeatureRow, featurize_seq
from mpatlas.topology import Topology, predict_topology


@dataclass(frozen=True)
class WrapScore:
    score: float
    reasons: tuple[str, ...]
    amenable: bool


def wrap_score(seq: str, topo: Topology | None = None, feat: FeatureRow | None = None) -> WrapScore:
    topo = topo or predict_topology(seq)
    feat = feat or featurize_seq(seq, topo)
    reasons: list[str] = []
    score = 0.0

    if feat.n_tm >= 1:
        score += 0.25
        reasons.append(f"n_tm={feat.n_tm}")
    else:
        reasons.append("no predicted TM belt")

    if feat.kd_mean > 0.3:
        score += 0.2
        reasons.append("hydrophobic overall")
    if feat.tm_frac >= 0.15:
        score += 0.2
        reasons.append(f"tm_frac={feat.tm_frac:.2f}")

    fusible = feat.n_term_extra > 0.5 or feat.c_term_extra > 0.5 or not feat.n_tm
    # cytoplasmic termini are also fusion-able for WRAP (cytoplasmic expression)
    n_term_kind = topo.segments[0].kind if topo.segments else "unknown"
    c_term_kind = topo.segments[-1].kind if topo.segments else "unknown"
    if n_term_kind in {"cyto", "extra"} or c_term_kind in {"cyto", "extra"}:
        fusible = True
        score += 0.15
        reasons.append(f"fusion terminus {n_term_kind}/{c_term_kind}")

    if 80 <= feat.length <= 900:
        score += 0.15
        reasons.append(f"length={feat.length}")
    elif feat.length > 900:
        reasons.append("long for single-particle WRAP-EM")
        score -= 0.05

    if feat.n_sequons >= 4:
        score -= 0.1
        reasons.append(f"glycosylation sequons={feat.n_sequons}")

    amenable = score >= 0.5 and feat.n_tm >= 1 and fusible
    return WrapScore(float(max(0.0, min(1.0, score))), tuple(reasons), amenable)
