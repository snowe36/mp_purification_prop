"""GPCRdb construct annotations — solubilization detergents and chromatography.

Excel from protwis/gpcrdb_data. Success-only GPCR structures. Column 'solvent type'
is the chemical (DDM, CHS); 'Solvent name' is detergent vs additive.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd

from mpatlas.catalog import Record
from mpatlas.detergents import canonicalize, family_of
from mpatlas.ingest.http import get
from mpatlas.paths import RAW

URL = (
    "https://raw.githubusercontent.com/protwis/gpcrdb_data/master/"
    "structure_data/construct_data/construct_annotations.xlsx"
)
OUT = RAW / "gpcrdb_construct_annotations.xlsx"

SOLUB = "Solubzn-Detergent"
PURIFY = "Purification-Treatment"
XTAL = "Xtal-Chemicals"


def download(dest: Path | None = None) -> Path:
    dest = dest or OUT
    get(URL, dest)
    return dest


def _drop_unnamed(df: pd.DataFrame) -> pd.DataFrame:
    keep = [c for c in df.columns if not str(c).startswith("Unnamed")]
    return df.loc[:, keep]


def _pdb(val) -> str | None:
    s = str(val or "").strip().upper()
    if len(s) >= 4 and s[:4].isalnum() and s not in {"PDB", "NAN", "NONE"}:
        return s[:4]
    return None


def _sheet(path: Path, name: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=name)
    return _drop_unnamed(df)


def _solub_rows(df: pd.DataFrame) -> list[Record]:
    out: list[Record] = []
    for rec in df.to_dict("records"):
        pdb = _pdb(rec.get("PDB"))
        chem = rec.get("solvent type (detergent/additive)")
        role = str(rec.get("Solvent name") or "").strip().lower() or None
        if not pdb or pd.isna(chem):
            continue
        name = canonicalize(str(chem))
        conc = rec.get("concentration")
        unit = rec.get("conc unit")
        acc = rec.get("Receptor entry_name (uniprot)")
        out.append(
            Record(
                source="gpcrdb",
                question="B-conditions",
                sequence="",
                accession=None if pd.isna(acc) else str(acc),
                pdb_id=pdb,
                architecture="gpcr",
                extra={
                    "sheet": "solubilization",
                    "detergent": name,
                    "family": family_of(name),
                    "role": role,
                    "conc": None if pd.isna(conc) else str(conc),
                    "conc_unit": None if pd.isna(unit) else str(unit),
                    "gpcr_entry": None if pd.isna(acc) else str(acc),
                },
            )
        )
    return out


def _purify_rows(df: pd.DataFrame) -> list[Record]:
    out: list[Record] = []
    for rec in df.to_dict("records"):
        pdb = _pdb(rec.get("PDB"))
        chrom = rec.get("Type of chromatography")
        treat = rec.get("Chemical/Enzymatic treatment")
        if not pdb:
            continue
        if (pd.isna(chrom) or not str(chrom).strip()) and (pd.isna(treat) or not str(treat).strip()):
            continue
        acc = rec.get("Receptor entry_name (uniprot)")
        out.append(
            Record(
                source="gpcrdb",
                question="B-conditions",
                sequence="",
                accession=None if pd.isna(acc) else str(acc),
                pdb_id=pdb,
                architecture="gpcr",
                extra={
                    "sheet": "purification",
                    "chromatography": None if pd.isna(chrom) else str(chrom).strip(),
                    "treatment": None if pd.isna(treat) else str(treat).strip(),
                    "remarks": None
                    if pd.isna(rec.get("REMARKS (can include detergent exchanges)"))
                    else str(rec.get("REMARKS (can include detergent exchanges)")),
                    "gpcr_entry": None if pd.isna(acc) else str(acc),
                },
            )
        )
    return out


def _xtal_rows(df: pd.DataFrame) -> list[Record]:
    out: list[Record] = []
    type_col = [c for c in df.columns if str(c).lower().startswith("chemical type")]
    type_col = type_col[0] if type_col else "Chemical Type"
    for rec in df.to_dict("records"):
        pdb = _pdb(rec.get("PDB"))
        chem = rec.get("Chemical")
        kind = rec.get(type_col) if type_col in rec else rec.get("Type")
        if not pdb or pd.isna(chem):
            continue
        name = canonicalize(str(chem))
        blob = f"{chem} {kind}".lower()
        if family_of(name) is None and "detergent" not in blob and "chs" not in blob:
            continue
        acc = rec.get("Receptor entry_name (uniprot)")
        conc = rec.get("Conc")
        unit = rec.get("Conc unit")
        out.append(
            Record(
                source="gpcrdb",
                question="B-conditions",
                sequence="",
                accession=None if pd.isna(acc) else str(acc),
                pdb_id=pdb,
                architecture="gpcr",
                extra={
                    "sheet": "xtal",
                    "detergent": name,
                    "family": family_of(name),
                    "role": None if pd.isna(kind) else str(kind),
                    "conc": None if pd.isna(conc) else str(conc),
                    "conc_unit": None if pd.isna(unit) else str(unit),
                    "gpcr_entry": None if pd.isna(acc) else str(acc),
                },
            )
        )
    return out


def load(path: Path | None = None) -> list[Record]:
    path = Path(path) if path else OUT
    if not path.exists():
        return []
    if path.suffix.lower() in {".xlsx", ".xls"}:
        solub = _sheet(path, SOLUB)
        purify = _sheet(path, PURIFY)
        xtal = _sheet(path, XTAL)
        return _solub_rows(solub) + _purify_rows(purify) + _xtal_rows(xtal)
    # tests may pass a CSV of solubilization rows
    df = pd.read_csv(path)
    return _solub_rows(df)


def summarize(rows: list[Record] | None = None) -> dict:
    rows = rows if rows is not None else load()
    solub = [r for r in rows if (r.extra or {}).get("sheet") == "solubilization"]
    dets = [r for r in solub if (r.extra or {}).get("role") == "detergent"]
    adds = [r for r in solub if (r.extra or {}).get("role") == "additive"]
    det_counts = Counter((r.extra or {}).get("detergent") for r in dets)
    add_counts = Counter((r.extra or {}).get("detergent") for r in adds)
    pdbs = {r.pdb_id for r in solub if r.pdb_id}
    combos = Counter()
    by_pdb: dict[str, list[str]] = {}
    for r in dets:
        if r.pdb_id and (r.extra or {}).get("detergent"):
            by_pdb.setdefault(r.pdb_id, []).append(str((r.extra or {})["detergent"]))
    for pdb, names in by_pdb.items():
        chs = any((r.extra or {}).get("detergent") == "CHS" for r in solub if r.pdb_id == pdb)
        key = "+".join(sorted(set(names)))
        if chs and "CHS" not in key:
            key = f"{key}+CHS"
        combos[key] += 1
    chrom = Counter(
        (r.extra or {}).get("chromatography")
        for r in rows
        if (r.extra or {}).get("sheet") == "purification" and (r.extra or {}).get("chromatography")
    )
    return {
        "n_rows": len(rows),
        "n_solubilization": len(solub),
        "n_pdb": len(pdbs),
        "detergents": dict(det_counts.most_common()),
        "additives": dict(add_counts.most_common()),
        "extract_combos": dict(combos.most_common(12)),
        "chromatography": dict(chrom.most_common(8)),
    }
