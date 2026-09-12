#!/usr/bin/env python
"""Fetch Hammon 2009 protein sequences from UniProt / NCBI into the fixture CSV."""

from __future__ import annotations

import re
import time
from pathlib import Path

import pandas as pd
import requests

from mpatlas.paths import FIXTURES

HAMMON = FIXTURES / "hammon2009.csv"
UA = {"User-Agent": "mp-atlas/0.1 (research)"}
UNIPROT_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$")
NCBI_RE = re.compile(r"^(NP|YP|XP|WP|AP)_")


def _parse_fasta(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    key, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if key is not None:
                out[key] = "".join(chunks)
            header = line[1:]
            key = header.split()[0]
            chunks = []
        else:
            chunks.append("".join(aa for aa in line.strip().upper() if aa.isalpha()))
    if key is not None:
        out[key] = "".join(chunks)
    return out


def _uniprot_ids(accs: list[str]) -> dict[str, str]:
    want = [a.split(".")[0] for a in accs if UNIPROT_RE.match(a.split(".")[0])]
    seqs: dict[str, str] = {}
    for i in range(0, len(want), 50):
        chunk = want[i : i + 50]
        url = "https://rest.uniprot.org/uniprotkb/accessions"
        r = requests.get(
            url,
            params={"accessions": ",".join(chunk), "format": "fasta"},
            headers=UA,
            timeout=90,
        )
        r.raise_for_status()
        parsed = _parse_fasta(r.text)
        for k, seq in parsed.items():
            acc = k.split("|")[1] if "|" in k else k.split()[0]
            seqs[acc] = seq
        time.sleep(0.3)
    return seqs


def _ncbi_ids(accs: list[str]) -> dict[str, str]:
    want = [a for a in accs if NCBI_RE.match(a)]
    seqs: dict[str, str] = {}
    for i in range(0, len(want), 40):
        chunk = want[i : i + 40]
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        r = requests.get(
            url,
            params={"db": "protein", "id": ",".join(chunk), "rettype": "fasta", "retmode": "text"},
            headers=UA,
            timeout=90,
        )
        r.raise_for_status()
        parsed = _parse_fasta(r.text)
        for k, seq in parsed.items():
            for acc in chunk:
                if acc in k or acc.split(".")[0] in k:
                    seqs[acc] = seq
                    break
            else:
                seqs[k] = seq
        time.sleep(0.4)
    return seqs


def main() -> None:
    df = pd.read_csv(HAMMON)
    accs = df["accession"].astype(str).tolist()
    uni = _uniprot_ids(accs)
    ncbi = _ncbi_ids(accs)
    seqs = []
    for acc in accs:
        seq = uni.get(acc) or uni.get(acc.split(".")[0]) or ncbi.get(acc) or ""
        seqs.append(seq)
    df["sequence"] = seqs
    df.to_csv(HAMMON, index=False)
    n = int((df["sequence"].str.len() > 10).sum())
    print(f"hammon sequences {n}/{len(df)} uniprot={len(uni)} ncbi={len(ncbi)}")


if __name__ == "__main__":
    main()
