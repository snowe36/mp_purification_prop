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

# Kotov 2019 Table 3 (abbreviation + chemistry). Water/blank wells omitted.
KOTOV_TABLE3 = """abbrev,family,name
Z3-10,propanesulfonate,Anzergent 3-10
Z3-12,propanesulfonate,Anzergent 3-12
Z3-14,propanesulfonate,Anzergent 3-14
DMG,dimethylglycine,n-Decyl-N,N-Dimethylglycine
DOMG,dimethylglycine,n-Dodecyl-N,N-Dimethylglycine
DDAO,amine-oxide,n-Decyl-N,N-Dimethylamine-N-Oxide
UDAO,amine-oxide,n-Undecyl-N,N-Dimethylamine-Oxide
LDAO,amine-oxide,n-Dodecyl-N,N-Dimethylamine-N-Oxide
CF-3,fos-choline,Cyclofos-3
CF-4,fos-choline,Cyclofos-4
CF-5,fos-choline,Cyclofos-5
CF-6,fos-choline,Cyclofos-6
CF-7,fos-choline,Cyclofos-7
FC-12,fos-choline,Fos-Choline-12
FC-13,fos-choline,Fos-Choline-13
FC-14,fos-choline,Fos-Choline-14
FC-15,fos-choline,Fos-Choline-15
FC-16,fos-choline,Fos-Choline-16
FC-I9,fos-choline,Fos-Choline-ISO-9
FC-I11,fos-choline,Fos-Choline-ISO-11
FC-U10-11,fos-choline,Fos-Choline-UNSAT-11-10
DHPC,fos-choline,DHPC
LPC-12,lyso-PC,LysoPC-12
LPC-14,lyso-PC,LysoPC-14
CHAPS,zwitterionic,CHAPS
CHAPSO,zwitterionic,CHAPSO
LAPAO,amine-oxide,LAPAO
TRIPAO,amine-oxide,TRIPAO
T-20,peg,Tween 20
Brij-35,peg,Brij-35
TX-100,peg,Triton X-100
TX-114,peg,Triton X-114
TX-305,peg,Triton X-305
TX-405,peg,Triton X-405
NID-P40,peg,NP-40
LMNG,maltose-NG,LMNG
OGNG,glucose-NG,OGNG
DMNG,maltose-NG,DMNG
CYMAL-5-NG,maltose-NG,CYMAL-5-NG
CYMAL-6-NG,maltose-NG,CYMAL-6-NG
GDN,maltose-NG,GDN
C6E3,peg,C6E3
C6E4,peg,C6E4
C6E5,peg,C6E5
C7E5,peg,C7E5
C8E4,peg,C8E4
C8E5,peg,C8E5
C8E6,peg,C8E6
C10E5,peg,C10E5
C10E6,peg,C10E6
C10E9,peg,C10E9
C12E7,peg,C12E7
C12E8,peg,C12E8
C12E9,peg,C12E9
C12E10,peg,C12E10
C13E8,peg,C13E8
CHAP,zwitterionic,Big CHAP
CHAP-D,zwitterionic,Big CHAP Deoxy
HTG,glucoside,n-Heptyl-β-D-Thioglucopyranoside
OTG,glucoside,n-Octyl-β-D-Thioglucopyranoside
OG,glucoside,n-Octyl-β-D-Glucopyranoside
NG,glucoside,n-Nonyl-β-D-Glucopyranoside
HEGA-9,glucamide,Hega-9
HEGA-10,glucamide,Hega-10
M9,glucamide,Mega-9
M10,glucamide,Mega-10
CYMAL-3,cymal,CYMAL-3
CYMAL-4,cymal,CYMAL-4
CYMAL-5,cymal,CYMAL-5
CYMAL-6,cymal,CYMAL-6
CYMAL-7,cymal,CYMAL-7
OM,maltoside,n-Octyl-β-D-Maltopyranoside
NM,maltoside,n-Nonyl-β-D-Maltopyranoside
DM,maltoside,n-Decyl-β-D-Maltopyranoside
UDM,maltoside,n-Undecyl-β-D-Maltopyranoside
DDM,maltoside,n-Dodecyl-β-D-Maltopyranoside
TDM,maltoside,n-Tridecyl-β-D-Maltopyranoside
S-12,sucrose,Sucrose-12
"""


def _ensure_kotov_dets() -> Path:
    if not KOTOV_DETS.exists():
        KOTOV_DETS.parent.mkdir(parents=True, exist_ok=True)
        KOTOV_DETS.write_text(KOTOV_TABLE3)
    return KOTOV_DETS


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
    _ensure_kotov_dets()
    return (
        _targets(HOGBOM_TARGETS, "hogbom2017", organism="Escherichia coli")
        + _det_rows(HOGBOM_DETS, "hogbom2017")
        + _targets(KOTOV_TARGETS, "kotov2019")
        + _det_rows(_ensure_kotov_dets(), "kotov2019")
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
