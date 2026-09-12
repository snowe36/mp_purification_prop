"""Canonical detergent names and a first-line extract prior.

Priors come from GPCRdb solubilization counts plus published screens
(Lantez 2015, Lin 2016, Högbom 2017, Kotov 2019). Success-biased.
Not a purify-vs-fail model.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from mpatlas.features import FeatureRow

_STRIP = re.compile(r"[^a-z0-9]+")

# Maps a stripped token to a short name used in catalogs and figures.
ALIASES = {
    "ddm": "DDM",
    "ndodecylbdmaltoside": "DDM",
    "ndodecylbdmaltopyranoside": "DDM",
    "ndodecylalphadmaltopyranoside": "DDM",
    "dodecylmaltoside": "DDM",
    "laurylmaltoside": "DDM",
    "lmt": "DDM",
    "dm": "DM",
    "ndecylbdmaltoside": "DM",
    "ndecylbdmaltopyranoside": "DM",
    "udm": "UDM",
    "nundecylbdmaltoside": "UDM",
    "nundecylbdmaltopyranoside": "UDM",
    "nm": "NM",
    "nnonylbdmaltoside": "NM",
    "nnonylbdmaltopyranoside": "NM",
    "om": "OM",
    "lmng": "LMNG",
    "laurylmaltoseneopentylglycol": "LMNG",
    "ng310": "LMNG",
    "dmng": "DMNG",
    "decylmaltoseneopentylglycol": "DMNG",
    "gdn": "GDN",
    "gdn101": "GDN",
    "og": "OG",
    "noctylbdglucoside": "OG",
    "noctylbdglucopyranoside": "OG",
    "bog": "OG",
    "ng": "NG",
    "nnonylbdglucoside": "NG",
    "nnonylbdglucopyranoside": "NG",
    "ldao": "LDAO",
    "lda": "LDAO",
    "ndodecylnndimethylaminenoxide": "LDAO",
    "ddao": "DDAO",
    "udao": "UDAO",
    "chs": "CHS",
    "cholesterylhemisuccinate": "CHS",
    "fc12": "FC12",
    "fc120": "FC12",
    "fos12": "FC12",
    "foscholine12": "FC12",
    "tx100": "TX100",
    "tritonx100": "TX100",
    "tritonsx100": "TX100",
    "anapoex100": "TX100",
    "lapao": "LAPAO",
    "chaps": "CHAPS",
    "cymal4": "CYMAL-4",
    "cymal5": "CYMAL-5",
    "cymal6": "CYMAL-6",
    "cymal7": "CYMAL-7",
    "cm4": "CYMAL-4",
    "cm5": "CYMAL-5",
    "cm6": "CYMAL-6",
    "cm7": "CYMAL-7",
    "c8e4": "C8E4",
    "c10e5": "C10E5",
    "c12e8": "C12E8",
    "htg": "HTG",
    "digitonin": "digitonin",
    "ogng": "OGNG",
}

FAMILY = {
    "DDM": "maltoside",
    "DM": "maltoside",
    "UDM": "maltoside",
    "NM": "maltoside",
    "OM": "maltoside",
    "LMNG": "maltose-NG",
    "DMNG": "maltose-NG",
    "GDN": "maltose-NG",
    "CYMAL-4": "cymal",
    "CYMAL-5": "cymal",
    "CYMAL-6": "cymal",
    "CYMAL-7": "cymal",
    "OG": "glucoside",
    "NG": "glucoside",
    "HTG": "glucoside",
    "OGNG": "glucose-NG",
    "LDAO": "amine-oxide",
    "DDAO": "amine-oxide",
    "UDAO": "amine-oxide",
    "LAPAO": "amine-oxide",
    "FC12": "fos-choline",
    "TX100": "peg",
    "C8E4": "peg",
    "C10E5": "peg",
    "C12E8": "peg",
    "CHS": "sterol",
    "CHAPS": "zwitterionic",
    "digitonin": "sterol-glycoside",
}


def _token(name: str) -> str:
    return _STRIP.sub("", name.lower().replace("β", "b").replace("α", "a"))


def canonicalize(name: str | None) -> str | None:
    if name is None:
        return None
    raw = str(name).strip()
    if not raw or raw.lower() in {"nan", "none", "-"}:
        return None
    key = _token(raw)
    if key in ALIASES:
        return ALIASES[key]
    # GPCRdb already stores short names.
    if raw.upper() in FAMILY or raw in FAMILY:
        return raw.upper() if raw.upper() in FAMILY else raw
    if raw.upper() in ALIASES.values():
        return raw.upper()
    return raw


def family_of(name: str | None) -> str | None:
    canon = canonicalize(name)
    if not canon:
        return None
    return FAMILY.get(canon)


@dataclass(frozen=True)
class DetergentPrior:
    extract: str
    polish: str
    rescue: str
    avoid: str
    why: str


def prior_for(feat: FeatureRow, *, organism: str | None = None) -> DetergentPrior:
    """Starting recipe if the playbook action is to extract at all.

    Helical IMPs: 1% DDM extract, 0.03% DDM SEC (Lin / Kotov / Högbom FSEC buffer).
    7-TM or eukaryotic: add 0.2% CHS (GPCRdb solubilization mode).
    Rescue class is maltose-NG (Kotov transporters, GPCRdb xtal LMNG).
    Fos-choline / PEG are Kotov unfolders — Lantez still listed Fos-12 as a
    solubilization winner, which is not the same as a folded-stability hit.
    """
    n_tm = int(feat.n_tm)
    org = (organism or "").lower()
    euk = any(
        h in org
        for h in ("eukaryot", "homo sapiens", "mus musculus", "human", "insect", "sf9")
    )
    seven = 6 <= n_tm <= 8
    if seven or euk:
        extract = "1% w/v DDM + 0.2% w/v CHS"
        why = "GPCRdb GPCR mode (DDM+CHS); Kotov/Lin still polish in 0.03% DDM"
    else:
        extract = "1% w/v DDM"
        why = "Lin/Kotov/Högbom first-line helical IMP extract (DDM)"
    return DetergentPrior(
        extract=extract,
        polish="0.03% w/v DDM (IMAC/SEC)",
        rescue="LMNG / DMNG / GDN",
        avoid="fos-choline and PEG as 'stability' detergents (Kotov); glucosides weaker on IMAC (Högbom)",
        why=why,
    )
