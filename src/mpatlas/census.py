from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from mpatlas.catalog import STATUS_RANK, Record
from mpatlas.funnel import attrition
from mpatlas.ingest import curnow, gfp, mpstruc, purificationdb, targettrack, uniprot, unitmp
from mpatlas.ingest.purificationdb import _detergent_hit
from mpatlas.paths import PROCESSED, REPORTS, ensure_dirs
from mpatlas.topology import predict_topology

STOP = {
    "B_min_purified": 100,
    "Bcond_min_tm": 50,
    "C_min_pos": 500,
    "C_min_bg": 2000,
}


def records_to_df(rows: list[Record]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame([r.to_dict() for r in rows])


def _domain(organism: str | None) -> str:
    if not organism:
        return "unknown"
    o = organism.lower()
    euk_hints = (
        "eukaryota",
        "homo sapiens",
        "mus musculus",
        "saccharomyces",
        "arabidopsis",
        "drosophila",
        "danio",
        "gallus",
        "rattus",
        "bos taurus",
        "caenorhabditis",
    )
    if any(h in o for h in euk_hints):
        return "eukaryote"
    if "virus" in o:
        return "virus"
    return "prokaryote"


def gather() -> dict[str, list[Record]]:
    print("loading curnow ...", flush=True)
    labelled = curnow.load_labelled()
    unlabelled = curnow.load_unlabelled()
    print("loading uniprot ...", flush=True)
    swiss = uniprot.load()
    print("loading mpstruc ...", flush=True)
    mps = mpstruc.load()
    print("loading unitmp ...", flush=True)
    pdbtm = unitmp.load_pdbtm()
    topdb = unitmp.load_topdb()
    print("loading targettrack ...", flush=True)
    tt = targettrack.load()
    print("loading purificationdb ...", flush=True)
    pdb = purificationdb.load()
    print("loading gfp ...", flush=True)
    gfp_rows = gfp.load()
    return {
        "curnow": labelled,
        "curnow_unlabelled": unlabelled,
        "gfp": gfp_rows,
        "uniprot": swiss,
        "mpstruc": mps,
        "pdbtm": pdbtm,
        "topdb": topdb,
        "targettrack": tt,
        "purificationdb": pdb,
    }


def run_census(
    bags: dict[str, list[Record]] | None = None,
    *,
    reports_dir: Path | None = None,
    processed_dir: Path | None = None,
    write_processed: bool = True,
) -> dict:
    ensure_dirs()
    reports_dir = reports_dir or REPORTS
    processed_dir = processed_dir or PROCESSED
    bags = bags or gather()
    layers = []
    for name, rows in bags.items():
        q = Counter(r.question for r in rows)
        layers.append(
            {
                "source": name,
                "n": len(rows),
                "n_with_sequence": sum(1 for r in rows if r.sequence),
                "questions": dict(q),
                "n_label_pos": sum(1 for r in rows if r.label == 1),
                "n_label_neg": sum(1 for r in rows if r.label == 0),
            }
        )

    tt = bags.get("targettrack") or []
    cloned = [r for r in tt if STATUS_RANK.get(r.status or "", 0) >= 2]
    expressed = [r for r in tt if STATUS_RANK.get(r.status or "", 0) >= 3]
    purified = [r for r in tt if STATUS_RANK.get(r.status or "", 0) >= 5]
    centers = Counter((r.center or "unknown") for r in purified)
    funnel_counts = attrition(tt)

    pdb_ids = set()
    for r in (bags.get("mpstruc") or []) + (bags.get("pdbtm") or []):
        if r.pdb_id:
            pdb_ids.add(r.pdb_id.upper()[:4])

    swiss = bags.get("uniprot") or []
    swiss_pos = 0
    for r in swiss:
        xrefs = str((r.extra or {}).get("pdb_xrefs") or r.pdb_id or "")
        hits = {p[:4].upper() for p in xrefs.replace(";", " ").split() if len(p) >= 4}
        if hits & pdb_ids:
            swiss_pos += 1
            r.label = 1.0
        r.extra = {**(r.extra or {}), "domain": _domain(r.organism)}

    pdb_rows = bags.get("purificationdb") or []
    tm_acc = {r.accession for r in swiss if r.accession}
    pdbtm_acc = {r.accession for r in (bags.get("topdb") or []) if r.accession}
    pdb_tm = []
    det_n = 0
    for r in pdb_rows:
        blob = " ".join(str(v) for v in (r.extra or {}).values())
        if r.sequence and len(pdb_rows) <= 500:
            n_tm = predict_topology(r.sequence).n_tm
        else:
            n_tm = 1 if r.accession in tm_acc or r.accession in pdbtm_acc else 0
        is_tm = n_tm >= 1 or (r.pdb_id or "").upper()[:4] in pdb_ids
        if is_tm:
            pdb_tm.append(r)
        if _detergent_hit(blob):
            det_n += 1

    labelled = bags.get("curnow") or []
    unlabelled = bags.get("curnow_unlabelled") or []
    b_ok = len(purified) >= STOP["B_min_purified"] and len(centers) >= 2
    bcond_ok = len(pdb_tm) >= STOP["Bcond_min_tm"]
    c_ok = (len(pdb_ids) >= STOP["C_min_pos"] or swiss_pos >= STOP["C_min_pos"]) and len(
        swiss
    ) >= STOP["C_min_bg"]
    c_model = swiss_pos >= 200 and (len(swiss) - swiss_pos) >= 200

    result = {
        "layers": layers,
        "targettrack": {
            "n": len(tt),
            "cloned": len(cloned),
            "expressed": len(expressed),
            "purified": len(purified),
            "centers_with_purified": dict(centers.most_common(20)),
            "funnel": funnel_counts,
            "claim_B": b_ok,
        },
        "purificationdb": {
            "n": len(pdb_rows),
            "tm_slice": len(pdb_tm),
            "detergent_rows": det_n,
            "claim_B_conditions_model": bcond_ok,
        },
        "structure": {
            "mpstruc_or_pdbtm_pdb_ids": len(pdb_ids),
            "swiss_tm": len(swiss),
            "swiss_with_structure_join": swiss_pos,
            "claim_C": c_ok,
            "claim_C_model": c_model,
        },
        "curnow": {
            "n": len(labelled),
            "n_unlabelled": len(unlabelled),
            "claim_A_local": len(labelled) >= 200,
            "claim_A_transfer": False,
            "gfp_n": len(bags.get("gfp") or []),
        },
        "stop": {
            "B": "fit" if b_ok else "census_only",
            "B-conditions": "fit" if bcond_ok else "lookup_only",
            "C": "fit" if c_ok else "census_only",
            "A": "fit_local" if labelled else "missing",
        },
    }

    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / "gate0.json").write_text(json.dumps(result, indent=2, default=str))
    (reports_dir / "gate0.md").write_text(_markdown(result))
    if write_processed:
        processed_dir.mkdir(parents=True, exist_ok=True)
        frames = {k: records_to_df(v) for k, v in bags.items() if v}
        for name, df in frames.items():
            df.to_parquet(processed_dir / f"{name}.parquet", index=False)
        pd.DataFrame(layers).to_csv(processed_dir / "census_layers.csv", index=False)
    return result


def _markdown(result: dict) -> str:
    lines = ["# Gate 0 census", ""]
    lines.append("| Source | n | with sequence |")
    lines.append("|---|---:|---:|")
    for layer in result["layers"]:
        lines.append(f"| {layer['source']} | {layer['n']} | {layer['n_with_sequence']} |")
    tt = result["targettrack"]
    lines += [
        "",
        "## B (TargetTrack)",
        f"- cloned {tt['cloned']}, expressed {tt['expressed']}, purified {tt['purified']}",
        f"- funnel: {tt.get('funnel', {})}",
        f"- claim B: **{tt['claim_B']}**",
        "",
        "## B-conditions (PurificationDB)",
        f"- n {result['purificationdb']['n']}, TM slice {result['purificationdb']['tm_slice']}, detergent rows {result['purificationdb']['detergent_rows']}",
        "",
        "## C",
        f"- PDB ids {result['structure']['mpstruc_or_pdbtm_pdb_ids']}, Swiss-Prot TM {result['structure']['swiss_tm']}, join positives {result['structure']['swiss_with_structure_join']}",
        "",
        "## Stop",
        *(f"- {k}: {v}" for k, v in result["stop"].items()),
        "",
    ]
    return "\n".join(lines)
