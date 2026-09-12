from __future__ import annotations

from pathlib import Path

import pandas as pd
import requests

from mpatlas.catalog import Record
from mpatlas.ingest import http as http_mod
from mpatlas.paths import RAW

OUT = RAW / "uniprot_swiss_tm.tsv"

STREAM = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?compressed=true"
    "&format=tsv"
    "&query=(reviewed:true)+AND+(ft_transmem:*)"
    "&fields=accession,id,gene_primary,organism_name,length,sequence,ft_transmem,xref_pdb"
)


def download(max_pages: int = 200) -> Path:
    """Download reviewed UniProt transmembrane entries as one compressed TSV stream."""
    del max_pages
    OUT.parent.mkdir(parents=True, exist_ok=True)
    headers = {**http_mod.HEADERS, "Accept-Encoding": "gzip"}
    with requests.get(STREAM, headers=headers, timeout=http_mod.TIMEOUT, stream=True) as r:
        r.raise_for_status()
        data = r.content
    if data[:2] == b"\x1f\x8b":
        import gzip

        text = gzip.decompress(data).decode("utf-8")
    else:
        text = data.decode("utf-8")
    OUT.write_text(text)
    return OUT


def load(path: Path | None = None) -> list[Record]:
    path = path or OUT
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t")
    seq_col = "Sequence" if "Sequence" in df.columns else None
    if seq_col is None:
        return []
    out = []
    for rec in df.to_dict("records"):
        seq = str(rec.get(seq_col, "")).strip()
        if not seq or seq == "nan":
            continue
        acc = str(rec.get("Entry", rec.get("accession", "")))
        org = str(rec.get("Organism", rec.get("organism_name", "")))
        pdb = str(rec.get("Cross-reference (PDB)", rec.get("xref_pdb", rec.get("PDB", ""))) or "")
        lineage = str(rec.get("Taxonomic lineage", rec.get("lineage", "")) or "")
        out.append(
            Record(
                source="uniprot_swiss_tm",
                question="C",
                sequence=seq,
                label=0.0,
                organism=org if org != "nan" else None,
                accession=acc,
                pdb_id=pdb.split(";")[0][:4].upper() if pdb not in {"", "nan"} else None,
                extra={"pdb_xrefs": pdb, "lineage": lineage},
            )
        )
    return out
