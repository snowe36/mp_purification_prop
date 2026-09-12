"""Natural-protein GFP overexpression cohorts (question A transfer). Optional."""

from __future__ import annotations

import pandas as pd

from mpatlas.catalog import Record
from mpatlas.paths import RAW

FILES = {
    "daley": RAW / "gfp_daley.csv",
    "hammon": RAW / "gfp_hammon.csv",
    "drew": RAW / "gfp_drew.csv",
}


def load() -> list[Record]:
    out: list[Record] = []
    for source, path in FILES.items():
        if not path.exists():
            continue
        df = pd.read_csv(path)
        cols = {c.lower(): c for c in df.columns}
        seq_col = cols.get("sequence")
        if seq_col is None:
            continue
        lab_col = cols.get("label", cols.get("expressed"))
        org_col = cols.get("organism")
        for rec in df.to_dict("records"):
            seq = "".join(aa for aa in str(rec.get(seq_col, "")).upper() if aa.isalpha())
            if not seq:
                continue
            lab = rec.get(lab_col) if lab_col else None
            org = rec.get(org_col) if org_col else None
            out.append(
                Record(
                    source=source,
                    question="A",
                    sequence=seq,
                    label=float(lab) if lab is not None and str(lab) not in {"", "nan"} else None,
                    organism=str(org) if org and str(org) != "nan" else None,
                )
            )
    return out
