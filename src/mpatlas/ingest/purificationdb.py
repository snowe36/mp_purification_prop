from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from mpatlas.catalog import Record
from mpatlas.ingest.http import get
from mpatlas.paths import RAW

OUT = RAW / "purificationdb.csv"
JSON_OUT = RAW / "purificationdb.json"

MP_DETERGENTS = (
    "ddm",
    "dm",
    "og",
    "ldao",
    "lmng",
    "gdn",
    "fc12",
    "fc-12",
    "chaps",
    "triton",
    "cholate",
    "digitonin",
    "ogng",
    "cymal",
    "lmg",
    "undecyl",
)

URLS = [
    "https://raw.githubusercontent.com/biobricks-ai/PurificationDB/main/brick/purificationdb.parquet",
    "https://raw.githubusercontent.com/biobricks-ai/PurificationDB/main/data/purificationdb.csv",
    "https://github.com/biobricks-ai/PurificationDB/raw/main/brick/purificationdb.parquet",
    "https://purificationdatabase.herokuapp.com/api/entries",
    "https://purificationdatabase.herokuapp.com/download/csv",
]


def download() -> Path | None:
    for url in URLS:
        try:
            dest = OUT if url.endswith(".csv") or "csv" in url else JSON_OUT
            if url.endswith(".parquet"):
                dest = RAW / "purificationdb.parquet"
            get(url, dest)
            if dest.exists() and dest.stat().st_size > 100:
                return dest
        except Exception:  # noqa: BLE001
            continue
    return None


def _detergent_hit(text: str) -> str | None:
    t = text.lower()
    for d in MP_DETERGENTS:
        if re.search(rf"\b{re.escape(d)}\b", t):
            return d
    return None


def load(path: Path | None = None) -> list[Record]:
    path = path or (OUT if OUT.exists() else JSON_OUT)
    parquet = RAW / "purificationdb.parquet"
    if path is None or not Path(path).exists():
        if parquet.exists():
            path = parquet
        else:
            return []
    path = Path(path)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    elif path.suffix == ".json":
        payload = json.loads(path.read_text())
        df = pd.DataFrame(payload if isinstance(payload, list) else payload.get("entries", []))
    else:
        df = pd.read_csv(path)
    if df.empty:
        return []
    cols = {c.lower(): c for c in df.columns}
    out: list[Record] = []
    for _, row in df.iterrows():
        seq = str(row.get(cols.get("sequence", ""), "") or "")
        acc = str(row.get(cols.get("uniprot", cols.get("uniprot_id", cols.get("uniprotid", ""))), "") or "")
        pdb = str(row.get(cols.get("pdb", cols.get("pdb_id", cols.get("pdbid", ""))), "") or "")
        chrom = str(row.get(cols.get("chromatography", cols.get("type", "")), "") or "SEC")
        blob = " ".join(str(v) for v in row.values if pd.notna(v))
        det = _detergent_hit(blob)
        extra = {c: row[c] for c in df.columns if c.lower() in {
            "ph", "salt", "buffer", "detergent", "additive", "nacl", "chromatography"
        }}
        extra["mp_detergent"] = det
        extra["chromatography"] = chrom
        out.append(
            Record(
                source="purificationdb",
                question="B-conditions",
                sequence="".join(aa for aa in seq if aa.isalpha()) if seq and seq != "nan" else "",
                label=None,
                accession=acc if acc not in {"", "nan"} else None,
                pdb_id=pdb.upper()[:4] if pdb not in {"", "nan"} else None,
                extra=extra,
            )
        )
    return out
