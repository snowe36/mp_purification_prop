"""Daley 2005 and Hammon 2009 GFP overexpression (question A transfer)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from mpatlas.catalog import Record
from mpatlas.paths import FIXTURES, RAW

DALEY = FIXTURES / "daley2005.csv"
HAMMON = FIXTURES / "hammon2009.csv"
UNIPROT = RAW / "uniprot_swiss_tm.tsv"


def _uniprot_maps() -> tuple[dict[str, tuple[str, str]], dict[str, tuple[str, str]]]:
    """accession -> (seq, org) and lowercased gene|ecoli -> (seq, org)."""
    by_acc: dict[str, tuple[str, str]] = {}
    by_gene: dict[str, tuple[str, str]] = {}
    if not UNIPROT.exists():
        return by_acc, by_gene
    df = pd.read_csv(UNIPROT, sep="\t")
    seq_col = "Sequence"
    acc_col = "Entry" if "Entry" in df.columns else "accession"
    gene_col = "Gene Names (primary)" if "Gene Names (primary)" in df.columns else None
    org_col = "Organism" if "Organism" in df.columns else None
    for rec in df.to_dict("records"):
        seq = str(rec.get(seq_col, "") or "")
        seq = "".join(aa for aa in seq.upper() if aa.isalpha())
        if not seq:
            continue
        acc = str(rec.get(acc_col, "") or "")
        org = str(rec.get(org_col, "") or "") if org_col else ""
        if acc and acc != "nan":
            by_acc[acc] = (seq, org)
        gene = str(rec.get(gene_col, "") or "") if gene_col else ""
        if gene and gene != "nan" and "escherichia coli" in org.lower():
            by_gene.setdefault(gene.lower(), (seq, org))
    return by_acc, by_gene


def load() -> list[Record]:
    by_acc, by_gene = _uniprot_maps()
    out: list[Record] = []
    if DALEY.exists():
        df = pd.read_csv(DALEY)
        for rec in df.to_dict("records"):
            gene = str(rec.get("gene") or "")
            hit = by_gene.get(gene.lower())
            seq, org = hit if hit else ("", "Escherichia coli")
            out.append(
                Record(
                    source="daley",
                    question="A",
                    sequence=seq,
                    label=float(rec.get("label")) if pd.notna(rec.get("label")) else None,
                    organism=org or "Escherichia coli",
                    accession=gene or None,
                    extra={
                        "gfp_ml": rec.get("gfp_ml"),
                        "cterm": rec.get("cterm"),
                        "joined": int(bool(seq)),
                    },
                )
            )
    if HAMMON.exists():
        df = pd.read_csv(HAMMON)
        for rec in df.to_dict("records"):
            acc = str(rec.get("accession") or "")
            raw = str(rec.get("sequence") or "")
            raw = "".join(aa for aa in raw.upper() if aa.isalpha())
            hit = by_acc.get(acc)
            seq, org = (raw, str(rec.get("organism") or "")) if raw else (hit if hit else ("", str(rec.get("organism") or "")))
            out.append(
                Record(
                    source="hammon",
                    question="A",
                    sequence=seq,
                    label=float(rec.get("label")) if pd.notna(rec.get("label")) else None,
                    organism=org or str(rec.get("organism") or "") or None,
                    accession=acc or str(rec.get("candidate") or "") or None,
                    extra={
                        "candidate": rec.get("candidate"),
                        "fsu": rec.get("fsu"),
                        "joined": int(bool(seq)),
                    },
                )
            )
    return out
