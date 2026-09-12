"""Published detergent screens (Lantez, Lin, Högbom, Kotov).

Committed tables only — no invented per-protein FSEC scores. Högbom Table 2
is a figure; we keep the 60-protein list, 16-detergent list, and family-level
claims from the text.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mpatlas.catalog import Record
from mpatlas.detergents import canonicalize, family_of
from mpatlas.paths import FIXTURES

HOGBOM_TARGETS = FIXTURES / "hogbom2017_targets.csv"
HOGBOM_DETS = FIXTURES / "hogbom2017_detergents.csv"
KOTOV_TARGETS = FIXTURES / "kotov2019_targets.csv"
KOTOV_DETS = FIXTURES / "kotov2019_detergents.csv"
FINDINGS = FIXTURES / "literature_findings.csv"


def _targets(path: Path, source: str, organism: str | None = None) -> list[Record]:
    if not path.exists():
        return []
    df = pd.read_csv(path)
    out: list[Record] = []
    for rec in df.to_dict("records"):
        name = str(rec.get("protein") or rec.get("accession") or "")
        org = rec.get("organism", organism)
        n_tm = rec.get("n_tm")
        extra = {k: rec[k] for k in rec if k not in {"protein", "organism"} and pd.notna(rec[k])}
        extra["kind"] = "target"
        pdb = str(rec.get("pdb_ids") or "").split(";")[0].strip()
        out.append(
            Record(
                source=source,
                question="B-conditions",
                sequence="",
                accession=name or None,
                organism=None if pd.isna(org) else str(org),
                pdb_id=pdb[:4].upper() if pdb and len(pdb) >= 4 else None,
                n_tm=None if pd.isna(n_tm) else int(n_tm),
                extra=extra,
            )
        )
    return out


def _det_rows(path: Path, source: str) -> list[Record]:
    if not path.exists():
        return []
    df = pd.read_csv(path)
    out: list[Record] = []
    for rec in df.to_dict("records"):
        abbrev = canonicalize(str(rec.get("abbrev") or rec.get("abbreviation") or ""))
        fam = rec.get("family") or family_of(abbrev)
        extra = {
            "kind": "detergent",
            "abbrev": abbrev,
            "family": None if pd.isna(fam) else str(fam),
            "chemical": rec.get("chemical", rec.get("name")),
        }
        out.append(
            Record(
                source=source,
                question="B-conditions",
                sequence="",
                accession=abbrev,
                extra=extra,
            )
        )
    return out


def _findings() -> list[Record]:
    if not FINDINGS.exists():
        return []
    df = pd.read_csv(FINDINGS)
    out: list[Record] = []
    for rec in df.to_dict("records"):
        det = canonicalize(rec.get("detergent")) if pd.notna(rec.get("detergent")) else None
        extra = {k: rec[k] for k in rec if pd.notna(rec[k])}
        extra["kind"] = "finding"
        extra["detergent"] = det
        extra["family"] = rec.get("family") or family_of(det)
        out.append(
            Record(
                source=str(rec.get("source") or "literature"),
                question="B-conditions",
                sequence="",
                accession=None if pd.isna(rec.get("protein")) else str(rec.get("protein")),
                extra=extra,
            )
        )
    return out


def load() -> list[Record]:
    return (
        _targets(HOGBOM_TARGETS, "hogbom2017", organism="Escherichia coli")
        + _det_rows(HOGBOM_DETS, "hogbom2017")
        + _targets(KOTOV_TARGETS, "kotov2019")
        + _det_rows(KOTOV_DETS, "kotov2019")
        + _findings()
    )


def screen_counts(rows: list[Record] | None = None) -> dict:
    rows = rows if rows is not None else load()
    by_src: dict[str, dict[str, int]] = {}
    for r in rows:
        d = by_src.setdefault(r.source, {"targets": 0, "detergents": 0, "findings": 0})
        kind = (r.extra or {}).get("kind")
        if kind == "target":
            d["targets"] += 1
        elif kind == "detergent":
            d["detergents"] += 1
        elif kind == "finding":
            d["findings"] += 1
    return by_src
