"""Kyte–Doolittle TM segments and positive-inside orientation."""

from __future__ import annotations

from dataclasses import dataclass

KD = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}

WINDOW = 19
KD_CUTOFF = 1.6
MIN_TM = 16
POSITIVE = set("KR")


@dataclass(frozen=True)
class Segment:
    start: int  # 1-indexed inclusive
    end: int
    kind: str  # tm | extra | cyto | unknown


@dataclass(frozen=True)
class Topology:
    n_tm: int
    segments: tuple[Segment, ...]
    n_in: bool
    architecture: str


def _kd_profile(seq: str) -> list[float]:
    seq = seq.upper()
    n = len(seq)
    half = WINDOW // 2
    scores = [0.0] * n
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        vals = [KD.get(seq[j], 0.0) for j in range(a, b)]
        scores[i] = sum(vals) / len(vals)
    return scores


def _tm_spans(seq: str) -> list[tuple[int, int]]:
    scores = _kd_profile(seq)
    n = len(seq)
    spans: list[tuple[int, int]] = []
    i = 0
    while i < n:
        if scores[i] < KD_CUTOFF:
            i += 1
            continue
        j = i
        while j < n and scores[j] >= KD_CUTOFF:
            j += 1
        if j - i >= MIN_TM:
            spans.append((i + 1, j))
        i = j
    return spans


def _loop_kr(seq: str, start: int, end: int) -> int:
    return sum(1 for aa in seq[start - 1 : end] if aa in POSITIVE)


def predict_topology(seq: str) -> Topology:
    """Windowed KD helices; orientation by positive-inside on odd vs even loops."""
    seq = "".join(aa for aa in seq.upper() if aa.isalpha())
    tms = _tm_spans(seq)
    n = len(seq)
    if not tms:
        architecture = "unknown"
        segs = (Segment(1, n, "unknown"),) if n else ()
        return Topology(0, segs, True, architecture)

    loops: list[tuple[int, int]] = []
    prev = 1
    for a, b in tms:
        if a > prev:
            loops.append((prev, a - 1))
        prev = b + 1
    if prev <= n:
        loops.append((prev, n))
    else:
        loops.append((n, n))

    odd_kr = 0
    even_kr = 0
    for i, (a, b) in enumerate(loops):
        kr = _loop_kr(seq, a, b)
        if i % 2 == 0:
            odd_kr += kr
        else:
            even_kr += kr
    n_in = odd_kr >= even_kr

    segs: list[Segment] = []
    prev = 1
    for k, (a, b) in enumerate(tms):
        if a > prev:
            extra_face = (k % 2 == 0 and not n_in) or (k % 2 == 1 and n_in)
            segs.append(Segment(prev, a - 1, "extra" if extra_face else "cyto"))
        segs.append(Segment(a, b, "tm"))
        prev = b + 1
    if prev <= n:
        extra_face = (len(tms) % 2 == 0 and not n_in) or (len(tms) % 2 == 1 and n_in)
        segs.append(Segment(prev, n, "extra" if extra_face else "cyto"))

    n_tm = len(tms)
    if n_tm == 1:
        architecture = "bitopic_helical"
    else:
        architecture = "polytopic_helical"
    return Topology(n_tm, tuple(segs), n_in, architecture)


def architecture_from_n_tm(n_tm: int, beta: bool = False) -> str:
    if beta:
        return "beta_barrel"
    if n_tm <= 0:
        return "unknown"
    if n_tm == 1:
        return "bitopic_helical"
    return "polytopic_helical"
