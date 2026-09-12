from __future__ import annotations

import json

import numpy as np
import pandas as pd

from mpatlas import figures as figmod
from mpatlas.census import run_census
from mpatlas.ingest import curnow, download_all
from mpatlas.ingest import uniprot as uniprot_mod
from mpatlas.models import run_binary, run_curnow, subsample_idx, write_results
from mpatlas.paths import FIGURES, PROCESSED, REPORTS, ensure_dirs
from mpatlas.wrap import wrap_score


def download_main() -> None:
    import sys

    skip_tt = "--skip-targettrack" in sys.argv
    status = download_all(skip_targettrack=skip_tt)
    for k, v in status.items():
        print(f"{k}: {v}", flush=True)


def census_main() -> None:
    result = run_census()
    print(result["stop"])
    print("wrote reports/gate0.md")


def _sub_frame(df: pd.DataFrame, y: np.ndarray, max_n: int = 4000) -> pd.DataFrame:
    idx = subsample_idx(y, max_n=max_n)
    return df.iloc[idx].reset_index(drop=True)


def express_main() -> None:
    ensure_dirs()
    rows = curnow.load_labelled()
    if not rows:
        raise SystemExit("no Curnow labels; run mpx-download")
    seqs = [r.sequence for r in rows]
    y = np.array([int(r.label or 0) for r in rows])
    df = run_curnow(seqs, y)
    write_results(df, "curnow_leakage")
    if not df.empty:
        figmod.leakage_auc(df, FIGURES)
        (FIGURES / "fig_curnow_auc.png").write_bytes((FIGURES / "fig_leakage_auc.png").read_bytes())
    print("A", df.to_string(index=False) if not df.empty else "no splits passed the sample floor")

    gate_path = REPORTS / "gate0.json"
    gate = json.loads(gate_path.read_text()) if gate_path.exists() else {}
    stop = gate.get("stop", {})

    tt_path = PROCESSED / "targettrack.parquet"
    if tt_path.exists() and stop.get("B") == "fit":
        tt = pd.read_parquet(tt_path)
        rank = tt["status_rank"] if "status_rank" in tt.columns else pd.Series([0] * len(tt))
        cloned = tt[rank >= 2].copy()
        if "sequence" in cloned.columns and len(cloned) > 80:
            cloned = cloned[cloned["sequence"].astype(str).str.len() > 20]
            yb = (cloned["status_rank"] >= 5).astype(int).to_numpy()
            if yb.sum() >= 30 and (yb == 0).sum() >= 30:
                cloned = _sub_frame(cloned, yb, 3000)
                yb = (cloned["status_rank"] >= 5).astype(int).to_numpy()
                centers = cloned.get("center", pd.Series(["unk"] * len(cloned))).fillna("unk").astype(str)
                holdout = sorted({c for c in centers if c != "NYCOMPS"})
                bdf = run_binary(
                    cloned["sequence"].tolist(),
                    yb,
                    centers.tolist(),
                    "center",
                    holdout=holdout or None,
                )
                write_results(bdf, "purified_leakage")
                print("B", bdf.to_string(index=False) if not bdf.empty else "skipped")
                funnel_counts = {
                    k: int(v)
                    for k, v in (gate.get("targettrack", {}).get("funnel") or {}).items()
                }
                if funnel_counts:
                    figmod.funnel(funnel_counts, FIGURES)
                fail = cloned[cloned["status_rank"] < 5].copy()
                scored = sorted(
                    (
                        (wrap_score(seq), acc, seq)
                        for seq, acc in zip(
                            fail["sequence"],
                            fail.get("accession", pd.Series([""] * len(fail))),
                            strict=False,
                        )
                    ),
                    key=lambda t: -t[0].score,
                )
                top = scored[:12]
                if top:
                    figmod.wrap_panel(
                        [seq for _, _, seq in top],
                        [str(acc)[:16] for _, acc, seq in top],
                        FIGURES,
                    )
                lines = ["score,amenable,accession,reasons"]
                for sc, acc, _seq in scored[:40]:
                    lines.append(f"{sc.score:.3f},{int(sc.amenable)},{acc},{'|'.join(sc.reasons)}")
                (PROCESSED / "wrap_b_failures.csv").write_text("\n".join(lines) + "\n")
                (REPORTS / "wrap_b_failures.csv").write_text("\n".join(lines) + "\n")

    swiss_path = PROCESSED / "uniprot.parquet"
    if swiss_path.exists() and gate.get("structure", {}).get("claim_C_model"):
        swiss = pd.read_parquet(swiss_path)
        if "sequence" in swiss.columns and "label" in swiss.columns:
            sub = swiss.dropna(subset=["sequence"]).copy()
            if sub["label"].nunique() > 1 and len(sub) > 80:
                yc = sub["label"].astype(int).to_numpy()
                sub = _sub_frame(sub, yc, 4000)
                groups = sub.get("domain", sub.get("organism", pd.Series(["unk"] * len(sub))))
                groups = groups.fillna("unk").astype(str).tolist()
                holdout = ["eukaryote"] if "eukaryote" in groups else None
                cdf = run_binary(sub["sequence"].tolist(), sub["label"].to_numpy(), groups, "organism", holdout=holdout)
                write_results(cdf, "structure_leakage")
                print("C", cdf.to_string(index=False) if not cdf.empty else "skipped")

    mp_path = PROCESSED / "mpstruc.parquet"
    if mp_path.exists():
        mp = pd.read_parquet(mp_path)
        if "architecture" in mp.columns:
            figmod.architecture_map(mp, FIGURES)

    frames = []
    for name, q in (("curnow_leakage", "A"), ("purified_leakage", "B"), ("structure_leakage", "C")):
        p = PROCESSED / f"{name}.csv"
        if p.exists():
            d = pd.read_csv(p)
            d["question"] = q
            frames.append(d)
    if frames:
        all_df = pd.concat(frames, ignore_index=True)
        all_df["split"] = all_df["question"] + " / " + all_df["split"]
        figmod.leakage_auc(all_df, FIGURES)


def panel_main() -> None:
    rows = curnow.load_labelled()[:40]
    if not rows:
        swiss = uniprot_mod.load()[:40]
        rows = swiss
    if not rows:
        raise SystemExit("no sequences")
    names = [(r.accession or f"seq{i}")[:18] for i, r in enumerate(rows)]
    figmod.wrap_panel([r.sequence for r in rows], names, FIGURES)
    scored = sorted(((wrap_score(r.sequence), r) for r in rows), key=lambda t: -t[0].score)
    print("top WRAP-amenable (heuristic)")
    for sc, r in scored[:10]:
        print(f"  {sc.score:.2f}  {r.accession}  {sc.amenable}  {', '.join(sc.reasons[:3])}")
