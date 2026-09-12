from __future__ import annotations

import pandas as pd

from mpatlas.catalog import Record
from mpatlas.paths import RAW

LABELLED = RAW / "curnow_labelled.csv"
UNLABELLED = RAW / "curnow_unlabelled.csv"

LABELLED_URL = (
    "https://raw.githubusercontent.com/Curnow-Lab-University-of-Bristol/"
    "Membrane-protein-library-ML/main/labelled_fasta.csv"
)
UNLABELLED_URL = (
    "https://raw.githubusercontent.com/Curnow-Lab-University-of-Bristol/"
    "Membrane-protein-library-ML/main/no_selected_seqs_protein.csv"
)


def download() -> None:
    from mpatlas.ingest.http import get

    if not LABELLED.exists():
        get(LABELLED_URL, LABELLED)
    if not UNLABELLED.exists():
        get(UNLABELLED_URL, UNLABELLED)


def _seq_lab_id(df: pd.DataFrame) -> tuple[str, str | None, str]:
    seq_col = "Sequence" if "Sequence" in df.columns else df.columns[min(1, len(df.columns) - 1)]
    lab_col = "Label" if "Label" in df.columns else None
    id_col = "ID" if "ID" in df.columns else df.columns[0]
    return seq_col, lab_col, id_col


def load_labelled() -> list[Record]:
    if not LABELLED.exists():
        return []
    df = pd.read_csv(LABELLED)
    seq_col, lab_col, id_col = _seq_lab_id(df)
    out = []
    for rec in df.to_dict("records"):
        seq = str(rec[seq_col]).strip()
        if not seq or seq == "nan":
            continue
        out.append(
            Record(
                source="curnow",
                question="A",
                sequence=seq,
                label=float(rec[lab_col]) if lab_col else None,
                organism="Escherichia coli",
                center="curnow_facs",
                accession=str(rec[id_col]),
                architecture="polytopic_helical",
            )
        )
    return out


def load_unlabelled() -> list[Record]:
    if not UNLABELLED.exists():
        return []
    df = pd.read_csv(UNLABELLED)
    seq_col, _, id_col = _seq_lab_id(df)
    out = []
    for rec in df.to_dict("records"):
        seq = str(rec[seq_col]).strip()
        if not seq or seq == "nan":
            continue
        out.append(
            Record(
                source="curnow_unlabelled",
                question="A",
                sequence=seq,
                label=None,
                organism="Escherichia coli",
                center="curnow_facs",
                accession=str(rec[id_col]),
                architecture="polytopic_helical",
            )
        )
    return out
