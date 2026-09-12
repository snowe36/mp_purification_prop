"""Load TargetTrack membrane-center trials from FASTA (Zenodo) or XML fixtures."""

from __future__ import annotations

import gzip
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path

from mpatlas.catalog import STATUS_RANK, Record, normalize_status
from mpatlas.ingest.http import get
from mpatlas.paths import RAW

ZENODO = "https://zenodo.org/records/821654/files/TargetTrack-1Jul2017.tar.gz?download=1"
TARBALL = RAW / "TargetTrack-1Jul2017.tar.gz"
EXTRACT = RAW / "targettrack"
FASTA_REL = "TargetTrack-1Jul2017/TargetTrack FASTA files/proteinTrialSeqs.fasta.gz"

MEMBRANE_CENTERS = {
    "NYCOMPS",
    "GPCR",
    "ISFI",
    "MPSBYNMR",
    "MPP",
    "CSMP",
    "MPID",
    "MPSBC",
    "ATCG3D",
    "TMPC",
    "TRANSPORTPDB",
}


def download() -> Path:
    if TARBALL.exists() and TARBALL.stat().st_size > 1_000_000:
        return TARBALL
    get(ZENODO, TARBALL)
    return TARBALL


def _fasta_path() -> Path | None:
    dest = EXTRACT / FASTA_REL
    if dest.exists():
        return dest
    if not TARBALL.exists():
        return None
    EXTRACT.mkdir(parents=True, exist_ok=True)
    with tarfile.open(TARBALL, "r:gz") as tar:
        tar.extract(tar.getmember(FASTA_REL), path=EXTRACT)
    return dest


def _local(tag: str) -> str:
    return tag.split("}")[-1].lower()


def _text(el: ET.Element | None) -> str:
    if el is None or el.text is None:
        return ""
    return el.text.strip()


def _protein_seq(raw: str) -> str:
    seq = "".join(aa for aa in raw.upper() if aa.isalpha())
    if not seq or set(seq) <= set("ATGCN"):
        return ""
    return seq


def _record(target_id: str, seq: str, status_raw: str, center: str, organism: str | None) -> Record | None:
    seq = _protein_seq(seq)
    if not seq:
        return None
    mem = center in MEMBRANE_CENTERS or "membrane" in (status_raw or "").lower()
    status = normalize_status(status_raw)
    rank = STATUS_RANK.get(status or "", 0)
    return Record(
        source="targettrack",
        question="B",
        sequence=seq,
        label=1.0 if rank >= 5 else 0.0,
        status=status,
        organism=organism,
        center=center,
        accession=target_id,
        extra={"status_rank": rank, "status_raw": status_raw, "membrane_center": mem},
    )


def _keep(rec: Record, tm_only: bool) -> bool:
    if not tm_only:
        return True
    return bool((rec.extra or {}).get("membrane_center"))


def _load_xml(path: Path, tm_only: bool) -> list[Record]:
    best: dict[tuple[str, str], Record] = {}
    for _, el in ET.iterparse(path, events=("end",)):
        if _local(el.tag) != "target":
            continue
        kids = {_local(c.tag): c for c in list(el)}
        tid = (
            _text(kids.get("targetid"))
            or _text(kids.get("target_id"))
            or el.attrib.get("id", "")
        )
        seq = _text(kids.get("sequence")) or _text(kids.get("onelettercode"))
        status_raw = _text(kids.get("status"))
        center = _text(kids.get("lab")) or _text(kids.get("center")) or "unknown"
        organism = _text(kids.get("organism")) or None
        rec = _record(tid, seq, status_raw, center, organism)
        el.clear()
        if rec is None or not _keep(rec, tm_only):
            continue
        key = (rec.accession or "", rec.sequence)
        prev = best.get(key)
        if prev is None or int((prev.extra or {}).get("status_rank", 0)) < int(
            (rec.extra or {}).get("status_rank", 0)
        ):
            best[key] = rec
    return list(best.values())


def _load_fasta(path: Path, tm_only: bool) -> list[Record]:
    open_fn = gzip.open if path.suffix == ".gz" or path.name.endswith(".gz") else open
    best: dict[tuple[str, str], Record] = {}
    header = None
    chunks: list[str] = []

    def flush() -> None:
        if not header:
            return
        parts = [x.strip() for x in header[1:].split("$")]
        if len(parts) < 4:
            return
        rec = _record(parts[0], "".join(chunks), parts[3], parts[2], None)
        if rec is None or not _keep(rec, tm_only):
            return
        key = (rec.accession or "", rec.sequence)
        prev = best.get(key)
        if prev is None or int((prev.extra or {}).get("status_rank", 0)) < int(
            (rec.extra or {}).get("status_rank", 0)
        ):
            best[key] = rec

    with open_fn(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                header = line.strip()
                chunks = []
            else:
                chunks.append(line.strip())
        flush()
    return list(best.values())


def load(path: Path | None = None, tm_only: bool = True) -> list[Record]:
    if path is not None:
        path = Path(path)
        if path.suffix.lower() == ".xml":
            return _load_xml(path, tm_only)
        return _load_fasta(path, tm_only)
    fasta = _fasta_path()
    if fasta is None or not Path(fasta).exists():
        return []
    return _load_fasta(Path(fasta), tm_only)
