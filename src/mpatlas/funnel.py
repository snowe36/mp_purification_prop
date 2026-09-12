"""TargetTrack status attrition. Not FSEC."""

from __future__ import annotations

from mpatlas.catalog import STATUS_RANK, Record

STEPS = ("selected", "cloned", "expressed", "soluble", "purified", "in pdb")


def attrition(rows: list[Record]) -> dict[str, int]:
    ranks = [STATUS_RANK.get(r.status or "", int((r.extra or {}).get("status_rank", 0))) for r in rows]
    thresholds = {"selected": 1, "cloned": 2, "expressed": 3, "soluble": 4, "purified": 5, "in pdb": 7}
    return {step: sum(1 for x in ranks if x >= thresholds[step]) for step in STEPS}
