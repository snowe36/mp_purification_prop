"""mpstruc XML: curated unique-protein view of solved membrane proteins."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from mpatlas.catalog import Record
from mpatlas.ingest.http import get
from mpatlas.paths import RAW

XML_URLS = [
    "https://blanco.biomol.uci.edu/mpstruc/listAll/mpstrucTblXml",
    "https://www.blanco.biomol.uci.edu/mpstruc/listAll/mpstrucTblXml",
]
OUT = RAW / "mpstruc.xml"


def download() -> Path:
    last = None
    for url in XML_URLS:
        try:
            get(url, OUT)
            return OUT
        except Exception as exc:  # noqa: BLE001
            last = exc
    raise RuntimeError(f"mpstruc download failed: {last}") from last


def _text(el: ET.Element | None) -> str:
    if el is None or el.text is None:
        return ""
    return el.text.strip()


def _local(tag: str) -> str:
    return tag.split("}")[-1].lower()


def _child(el: ET.Element, name: str) -> ET.Element | None:
    want = name.lower()
    for c in el:
        if _local(c.tag) == want:
            return c
    return None


def _arch(group: str) -> str:
    n = group.lower()
    if "monotopic" in n:
        return "monotopic"
    if "beta" in n:
        return "beta_barrel"
    if "bitopic" in n or "single" in n:
        return "bitopic_helical"
    if "polytopic" in n or "helical" in n or "transmembrane" in n or "alpha" in n:
        return "polytopic_helical"
    return "unknown"


def _append(out: list[Record], el: ET.Element, group: str, subgroup: str) -> None:
    pdb = _text(_child(el, "pdbCode")) or _text(_child(el, "pdb"))
    name = _text(_child(el, "name")) or pdb
    if not pdb and not name:
        return
    species = _text(_child(el, "species")) or None
    domain = _text(_child(el, "taxonomicDomain")) or None
    year_el = _child(el, "bibliography")
    year = _text(_child(year_el, "year")) if year_el is not None else ""
    out.append(
        Record(
            source="mpstruc",
            question="C",
            sequence="",
            label=1.0,
            organism=species,
            accession=name,
            pdb_id=pdb.upper()[:4] if pdb else None,
            architecture=_arch(group),
            extra={"group": group, "subgroup": subgroup, "domain": domain, "year": year},
        )
    )


def _walk(el: ET.Element, out: list[Record], group: str = "", subgroup: str = "") -> None:
    tag = _local(el.tag)
    if tag == "group":
        group = _text(_child(el, "name")) or group
        for c in el:
            _walk(c, out, group, subgroup)
        return
    if tag == "subgroup":
        subgroup = _text(_child(el, "name")) or subgroup
        for c in el:
            _walk(c, out, group, subgroup)
        return
    if tag in {"protein", "memberprotein"}:
        _append(out, el, group, subgroup)
        members = _child(el, "memberProteins")
        if members is not None:
            for c in members:
                if _local(c.tag) == "memberprotein":
                    _append(out, c, group, subgroup)
        return
    for c in el:
        _walk(c, out, group, subgroup)


def load(path: Path | None = None) -> list[Record]:
    path = path or OUT
    if not path.exists():
        return []
    tree = ET.parse(path)
    out: list[Record] = []
    _walk(tree.getroot(), out)
    seen: set[str] = set()
    uniq: list[Record] = []
    for rec in out:
        key = rec.pdb_id or rec.accession or ""
        if not key or key in seen:
            continue
        seen.add(key)
        uniq.append(rec)
    return uniq
