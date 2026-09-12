from __future__ import annotations

import json

import numpy as np
import pandas as pd

from mpatlas import figures as figmod
from mpatlas.census import run_census
from mpatlas.ingest import curnow, download_all
from mpatlas.ingest import uniprot as uniprot_mod
from mpatlas.embed import load as load_embeds
from mpatlas.ingest import gfp
from mpatlas.models import run_binary, run_curnow, run_transfer, subsample_idx, write_results
from mpatlas.paths import FIGURES, PROCESSED, REPORTS, ensure_dirs
from mpatlas.playbook import advise_table, fit_purify_logreg
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
    gpcr = result.get("gpcrdb") or {}
    dets = gpcr.get("detergents") or {}
    if dets:
        figmod.detergents_bar(dets, FIGURES)
        print("GPCRdb PDBs", gpcr.get("n_pdb"), "top detergents", list(dets.items())[:6])
    print("wrote reports/gate0.md")


def _sub_frame(df: pd.DataFrame, y: np.ndarray, max_n: int = 4000) -> pd.DataFrame:
    idx = subsample_idx(y, max_n=max_n)
    return df.iloc[idx].reset_index(drop=True)


def _modes(seqs: list[str]):
    embeds = load_embeds()
    modes = ["compose"]
    if embeds is not None and embeds.coverage(seqs) >= 0.99:
        modes = ["compose", "esm", "both"]
    elif embeds is not None:
        print(f"embeddings coverage {embeds.coverage(seqs):.3f}; compose only", flush=True)
    return embeds, modes


def _compose(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "features" not in df.columns:
        return df
    return df[df["features"].eq("compose")].copy()


def express_main() -> None:
    ensure_dirs()
    rows = curnow.load_labelled()
    if not rows:
        raise SystemExit("no Curnow labels; run mpx-download")
    seqs = [r.sequence for r in rows]
    y = np.array([int(r.label or 0) for r in rows])
    embeds, modes = _modes(seqs)
    frames = [run_curnow(seqs, y, embeds=embeds, mode=m) for m in modes]
    df = pd.concat([f for f in frames if not f.empty], ignore_index=True) if frames else pd.DataFrame()
    write_results(df, "curnow_leakage")
    plot = _compose(df)
    if not plot.empty:
        figmod.leakage_auc(plot, FIGURES)
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
                b_seqs = cloned["sequence"].tolist()
                b_embeds, b_modes = _modes(b_seqs)
                bdf = pd.concat(
                    [
                        run_binary(
                            b_seqs,
                            yb,
                            centers.tolist(),
                            "center",
                            holdout=holdout or None,
                            embeds=b_embeds,
                            mode=m,
                        )
                        for m in b_modes
                    ],
                    ignore_index=True,
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
                c_seqs = sub["sequence"].tolist()
                c_embeds, c_modes = _modes(c_seqs)
                cdf = pd.concat(
                    [
                        run_binary(
                            c_seqs,
                            sub["label"].to_numpy(),
                            groups,
                            "organism",
                            holdout=holdout,
                            embeds=c_embeds,
                            mode=m,
                        )
                        for m in c_modes
                    ],
                    ignore_index=True,
                )
                write_results(cdf, "structure_leakage")
                print("C", cdf.to_string(index=False) if not cdf.empty else "skipped")

    mp_path = PROCESSED / "mpstruc.parquet"
    if mp_path.exists():
        mp = pd.read_parquet(mp_path)
        if "architecture" in mp.columns:
            figmod.architecture_map(mp, FIGURES)

    gfp_rows = gfp.load()
    xfer_frames = []
    for src in ("daley", "hammon"):
        hit = [r for r in gfp_rows if r.source == src and r.sequence and r.label is not None]
        if len(hit) < 80:
            continue
        te_seq = [r.sequence for r in hit]
        te_y = np.array([int(r.label) for r in hit])
        x_embeds, x_modes = _modes(te_seq)
        for m in x_modes:
            try:
                xdf = run_transfer(seqs, y, te_seq, te_y, f"curnow_to_{src}", embeds=x_embeds, mode=m)
            except ValueError:
                continue
            if not xdf.empty:
                xfer_frames.append(xdf)
    if xfer_frames:
        xfer = pd.concat(xfer_frames, ignore_index=True)
        write_results(xfer, "gfp_transfer")
        print("A-transfer", xfer.to_string(index=False))

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
        figmod.leakage_auc(_compose(all_df), FIGURES)
        if "features" in all_df.columns and all_df["features"].nunique() > 1:
            figmod.embed_auc(all_df, FIGURES)


def playbook_main() -> None:
    ensure_dirs()
    tt_path = PROCESSED / "targettrack.parquet"
    if not tt_path.exists():
        raise SystemExit("no TargetTrack catalog; run mpx-census")
    tt = pd.read_parquet(tt_path)
    rank = tt["status_rank"] if "status_rank" in tt.columns else pd.Series([0] * len(tt))
    expressed = tt[(rank >= 3) & (tt["sequence"].astype(str).str.len() > 20)].copy()
    if len(expressed) < 80:
        raise SystemExit("not enough expressed TargetTrack rows")
    y = (expressed["status_rank"] >= 5).astype(int).to_numpy()
    train = _sub_frame(expressed, y, 3000)
    ytr = (train["status_rank"] >= 5).astype(int).to_numpy()
    centers = train.get("center", pd.Series(["unk"] * len(train))).fillna("unk").astype(str)
    holdout = sorted({c for c in centers if c != "NYCOMPS"})
    leak = run_binary(train["sequence"].tolist(), ytr, centers.tolist(), "center", holdout=holdout or None)
    write_results(leak, "playbook_leakage")
    print("P(purified|expressed)", leak.to_string(index=False) if not leak.empty else "skipped")

    pipe, coefs = fit_purify_logreg(train["sequence"].tolist(), ytr)
    coefs.to_csv(REPORTS / "playbook_coefs.csv", index=False)
    coefs.to_csv(PROCESSED / "playbook_coefs.csv", index=False)
    figmod.playbook_coefs(coefs, FIGURES)

    fail = expressed[expressed["status_rank"] < 5]
    if len(fail) > 800:
        fail = fail.sample(n=800, random_state=0)
    swiss_path = PROCESSED / "uniprot.parquet"
    swiss_n = 0
    frames = [
        advise_table(
            pipe,
            fail["sequence"].tolist(),
            fail.get("accession", pd.Series([f"t{i}" for i in range(len(fail))])).astype(str).tolist(),
            source="tt_expressed_not_purified",
        )
    ]
    if swiss_path.exists():
        swiss = pd.read_parquet(swiss_path).dropna(subset=["sequence"])
        swiss = swiss[swiss["sequence"].astype(str).str.len().between(60, 900)]
        if len(swiss) > 600:
            swiss = swiss.sample(n=600, random_state=0)
        swiss_n = len(swiss)
        frames.append(
            advise_table(
                pipe,
                swiss["sequence"].tolist(),
                swiss.get("accession", pd.Series([f"u{i}" for i in range(len(swiss))])).astype(str).tolist(),
                source="swiss_tm_sample",
            )
        )
    panel = pd.concat(frames, ignore_index=True)
    panel.to_csv(REPORTS / "playbook_panel.csv", index=False)
    panel.to_csv(PROCESSED / "playbook_panel.csv", index=False)
    figmod.playbook_actions(panel, FIGURES)
    figmod.playbook_scatter(panel[panel["source"] == "tt_expressed_not_purified"], FIGURES)
    gate_path = REPORTS / "gate0.json"
    if gate_path.exists():
        gpcr = json.loads(gate_path.read_text()).get("gpcrdb") or {}
        dets = gpcr.get("detergents") or {}
        if dets:
            figmod.detergents_bar(dets, FIGURES)

    print(panel.groupby(["source", "action"]).size().to_string())
    print("top wrap_rescue (B-failures)")
    wrap = panel[(panel["source"] == "tt_expressed_not_purified") & (panel["action"] == "wrap_rescue")]
    wrap = wrap.sort_values("wrap_score", ascending=False).head(8)
    for _, row in wrap.iterrows():
        print(f"  {row.wrap_score:.2f}  p={row.p_purify:.2f}  {row.accession}  {row.why[:80]}")
    print(f"scored {len(fail)} B-failures" + (f" + {swiss_n} Swiss-Prot TM" if swiss_n else ""))


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
