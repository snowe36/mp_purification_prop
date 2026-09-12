from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from mpatlas.catalog import Record
from mpatlas.ingest.http import get
from mpatlas.paths import RAW
from mpatlas.topology import architecture_from_n_tm

PDBTM_URLS = [
    "https://pdbtm.unitmp.org/data/pdbtmall.xml",
    "https://pdbtm.unitmp.org/data/pdbtmall",
]
TOPDB_URLS = [
    "https://topdb.unitmp.org/download.php?type=xml",
    "https://topdb.unitmp.org/data/topdb.xml",
]

TOPDB_OUT = RAW / "topdb.xml"
PDBTM_OUT = RAW / "pdbtm.xml"


def _try_urls(urls: list[str], dest: Path) -> Path | None:
    for url in urls:
        try:
            get(url, dest)
            if dest.exists() and dest.stat().st_size > 1000:
                return dest
        except Exception:  # noqa: BLE001
            continue
    return None


def download() -> dict[str, Path | None]:
    return {
        "topdb": _try_urls(TOPDB_URLS, TOPDB_OUT),
        "pdbtm": _try_urls(PDBTM_URLS, PDBTM_OUT),
    }


def _local(tag: str) -> str:
    return tag.split("}")[-1].lower()


def _seq_from(el: ET.Element) -> str:
    parts = []
    for child in el.iter():
        if _local(child.tag) in {"sequence", "seq"} and child.text:
            parts.append("".join(aa for aa in child.text if aa.isalpha()))
    return max(parts, key=len) if parts else ""


def load_pdbtm(path: Path | None = None) -> list[Record]:
    path = path or PDBTM_OUT
    if not path or not path.exists():
        return []
    out: list[Record] = []
    for _, el in ET.iterparse(path, events=("end",)):
        if _local(el.tag) not in {"pdbtm", "protein", "entry"}:
            continue
        pdb = el.attrib.get("ID") or el.attrib.get("id") or el.attrib.get("pdbid") or ""
        seq = _seq_from(el)
        n_tm = 0
        for child in el.iter():
            if _local(child.tag) in {"region", "membrane", "tmhelix"}:
                typ = (child.attrib.get("type") or child.attrib.get("Type") or "").lower()
                if "m" == typ or "tm" in typ or "helix" in typ:
                    n_tm += 1
        if not pdb and not seq:
            el.clear()
            continue
        out.append(
            Record(
                source="pdbtm",
                question="C",
                sequence=seq,
                label=1.0,
                pdb_id=pdb.upper()[:4] if pdb else None,
                n_tm=n_tm or None,
                architecture=architecture_from_n_tm(n_tm) if n_tm else None,
            )
        )
        el.clear()
    return out


def load_topdb(path: Path | None = None) -> list[Record]:
    path = path or TOPDB_OUT
    if not path or not path.exists():
        return []
    out: list[Record] = []
    for _, el in ET.iterparse(path, events=("end",)):
        if _local(el.tag) not in {"entry", "protein", "topdb"}:
            continue
        acc = el.attrib.get("id") or el.attrib.get("uniprot") or ""
        seq = _seq_from(el)
        if not acc and not seq:
            el.clear()
            continue
        n_tm = sum(1 for c in el.iter() if _local(c.tag) in {"tm", "membrane", "region"} and "tm" in (
            (c.attrib.get("loc") or c.attrib.get("type") or "").lower()
        ))
        out.append(
            Record(
                source="topdb",
                question="C",
                sequence=seq,
                label=None,
                accession=acc,
                n_tm=n_tm or None,
            )
        )
        el.clear()
    return out
