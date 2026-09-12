from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from mpatlas.topology import KD, Topology, predict_topology

AA = tuple("ACDEFGHIKLMNPQRSTVWY")
SEQUON = re.compile(r"N[^P][ST]")


@dataclass(frozen=True)
class FeatureRow:
    length: int
    kd_mean: float
    net_charge: float
    n_tm: int
    tm_frac: float
    longest_tm: int
    longest_loop: int
    extra_frac: float
    n_sequons: int
    n_cys: int
    c_term_extra: float
    n_term_extra: float
    composition: tuple[float, ...]


FEATURE_NAMES = (
    ["length", "kd_mean", "net_charge", "n_tm", "tm_frac", "longest_tm", "longest_loop"]
    + ["extra_frac", "n_sequons", "n_cys", "c_term_extra", "n_term_extra"]
    + [f"aa_{a}" for a in AA]
)


def _charge(seq: str) -> float:
    return sum(1 if a in "KR" else -1 if a in "DE" else 0 for a in seq)


def featurize_seq(seq: str, topo: Topology | None = None) -> FeatureRow:
    seq = "".join(aa for aa in seq.upper() if aa in KD or aa in AA)
    n = max(len(seq), 1)
    topo = topo or predict_topology(seq)
    tm_len = sum(s.end - s.start + 1 for s in topo.segments if s.kind == "tm")
    extra_len = sum(s.end - s.start + 1 for s in topo.segments if s.kind == "extra")
    loops = [s.end - s.start + 1 for s in topo.segments if s.kind != "tm"]
    tms = [s.end - s.start + 1 for s in topo.segments if s.kind == "tm"]
    n_term = topo.segments[0].kind if topo.segments else "unknown"
    c_term = topo.segments[-1].kind if topo.segments else "unknown"
    sequons = 0
    for m in SEQUON.finditer(seq):
        pos = m.start() + 1
        in_tm = any(s.start <= pos <= s.end and s.kind == "tm" for s in topo.segments)
        if not in_tm:
            sequons += 1
    comp = tuple(seq.count(a) / n for a in AA)
    return FeatureRow(
        length=len(seq),
        kd_mean=sum(KD.get(a, 0.0) for a in seq) / n,
        net_charge=_charge(seq),
        n_tm=topo.n_tm,
        tm_frac=tm_len / n,
        longest_tm=max(tms) if tms else 0,
        longest_loop=max(loops) if loops else len(seq),
        extra_frac=extra_len / n,
        n_sequons=sequons,
        n_cys=seq.count("C"),
        c_term_extra=1.0 if c_term == "extra" else 0.0,
        n_term_extra=1.0 if n_term == "extra" else 0.0,
        composition=comp,
    )


def to_vector(row: FeatureRow) -> np.ndarray:
    core = [
        row.length,
        row.kd_mean,
        row.net_charge,
        row.n_tm,
        row.tm_frac,
        row.longest_tm,
        row.longest_loop,
        row.extra_frac,
        row.n_sequons,
        row.n_cys,
        row.c_term_extra,
        row.n_term_extra,
        *row.composition,
    ]
    return np.asarray(core, dtype=np.float64)


def matrix(sequences: list[str]) -> np.ndarray:
    return np.vstack([to_vector(featurize_seq(s)) for s in sequences])
