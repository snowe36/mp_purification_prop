#!/usr/bin/env python
"""Build unique protein FASTA for ESM embedding (Swiss-Prot TM + Curnow + TOPDB + TargetTrack)."""

from __future__ import annotations

import gzip
from pathlib import Path

from mpatlas.embed import seq_key
from mpatlas.ingest import curnow, gfp, targettrack, unitmp, uniprot
from mpatlas.paths import PROCESSED, ensure_dirs

MIN_LEN = 20
MAX_LEN = 1022


def _add(store: dict[str, str], seq: str) -> None:
    seq = "".join(aa for aa in seq.upper() if aa.isalpha())
    if len(seq) < MIN_LEN:
        return
    seq = seq[:MAX_LEN]
    store[seq_key(seq)] = seq


def main() -> None:
    ensure_dirs()
    store: dict[str, str] = {}
    print("uniprot ...", flush=True)
    for r in uniprot.load():
        _add(store, r.sequence)
    print("n", len(store), flush=True)
    print("curnow ...", flush=True)
    for r in curnow.load_labelled() + curnow.load_unlabelled():
        _add(store, r.sequence)
    print("n", len(store), flush=True)
    print("topdb ...", flush=True)
    for r in unitmp.load_topdb():
        _add(store, r.sequence)
    print("n", len(store), flush=True)
    print("targettrack ...", flush=True)
    for r in targettrack.load():
        _add(store, r.sequence)
    print("n", len(store), flush=True)
    print("gfp ...", flush=True)
    for r in gfp.load():
        _add(store, r.sequence)
    print("n", len(store), flush=True)
    dest = PROCESSED / "embed_universe.fasta.gz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest, "wt") as fh:
        for k, seq in store.items():
            fh.write(f">{k}\n{seq}\n")
    meta = PROCESSED / "embed_universe.txt"
    meta.write_text(f"n={len(store)}\nmin_len={MIN_LEN}\nmax_len={MAX_LEN}\n")
    print(f"wrote {dest} n={len(store)}", flush=True)


if __name__ == "__main__":
    main()
