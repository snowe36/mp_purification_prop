"""Offline demo: Curnow fixture/labels → leakage + WRAP figures."""

from __future__ import annotations

import json

import numpy as np

from mpatlas.catalog import Record
from mpatlas.census import run_census
from mpatlas.figures import (
    funnel,
    leakage_auc,
    playbook_actions,
    playbook_coefs,
    playbook_scatter,
    wrap_panel,
)
from mpatlas.ingest.curnow import LABELLED, load_labelled
from mpatlas.models import run_curnow
from mpatlas.paths import DEMO_OUT, FIXTURES, ensure_dirs
from mpatlas.playbook import advise_table, fit_purify_logreg


def _synthetic_library(n: int = 240, seed: int = 0) -> list[Record]:
    rng = np.random.default_rng(seed)
    # Three 2-TM scaffolds. Expression tracks family, not just composition.
    cores = [
        "MIFLFVLLAIVVLITGLLMLNTSNSPYIRLIAFILLVITGLIMLK",
        "MGSPWLRLLHLILAILVLITGLIMLLNTSNSPYLRLIHFLLAILV",
        "MAKKLILIVGALLLIVGALLKNDTSPYAIAIAIAIAIVVLITGLL",
    ]
    rows = []
    for i in range(n):
        fam = i % 3
        seq = list(cores[fam] + cores[fam][::-1][:60])
        for _ in range(6):
            j = int(rng.integers(5, len(seq) - 5))
            seq[j] = rng.choice(list("AILVF"))
        label = 1.0 if fam != 2 else 0.0
        if fam == 0 and rng.random() < 0.15:
            label = 0.0
        if fam == 2 and rng.random() < 0.15:
            label = 1.0
        rows.append(
            Record(
                source="fixture",
                question="A",
                sequence="".join(seq),
                label=label,
                accession=f"fix{i}",
                organism="synthetic",
                center=f"family{fam}",
                architecture="polytopic_helical",
            )
        )
    return rows


def main() -> None:
    ensure_dirs()
    DEMO_OUT.mkdir(parents=True, exist_ok=True)
    FIXTURES.mkdir(parents=True, exist_ok=True)

    real = load_labelled() if LABELLED.exists() else []
    rows = real if len(real) >= 200 else _synthetic_library()
    source = "curnow" if rows and rows[0].source == "curnow" else "fixture"

    seqs = [r.sequence for r in rows]
    y = np.array([int(r.label or 0) for r in rows])
    df = run_curnow(seqs, y)
    df.to_csv(DEMO_OUT / "leakage.csv", index=False)
    if not df.empty:
        leakage_auc(df, DEMO_OUT)

    funnel_counts = {
        "cloned": len(rows),
        "expressed": int(y.sum()),
        "purified": max(1, int(y.sum() * 0.25)),
        "in pdb": max(1, int(y.sum() * 0.04)),
    }
    if source == "fixture":
        # Synthetic attrition so the figure exists; labeled as fixture in meta.
        funnel(funnel_counts, DEMO_OUT)
    else:
        funnel(
            {
                "labelled FACS": len(rows),
                "bright": int((y == 1).sum()),
                "dim": int((y == 0).sum()),
            },
            DEMO_OUT,
        )

    panel = rows[:: max(1, len(rows) // 12)][:12]
    wrap_panel([r.sequence for r in panel], [(r.accession or "")[:16] for r in panel], DEMO_OUT)

    bags = {
        "curnow": rows,
        "curnow_unlabelled": [],
        "gfp": [],
        "uniprot": [],
        "mpstruc": [],
        "pdbtm": [],
        "topdb": [],
        "targettrack": [],
        "purificationdb": [],
    }
    census = run_census(bags, reports_dir=DEMO_OUT, processed_dir=DEMO_OUT, write_processed=False)
    (DEMO_OUT / "gate0.json").write_text(json.dumps(census, indent=2, default=str))

    synth = _synthetic_library()
    sy = np.array([int(r.label or 0) for r in synth])
    pipe, coefs = fit_purify_logreg([r.sequence for r in synth], sy)
    panel = advise_table(
        pipe,
        [r.sequence for r in synth],
        [r.accession or f"s{i}" for i, r in enumerate(synth)],
        source="fixture",
    )
    coefs.to_csv(DEMO_OUT / "playbook_coefs.csv", index=False)
    panel.to_csv(DEMO_OUT / "playbook_panel.csv", index=False)
    playbook_coefs(coefs, DEMO_OUT)
    playbook_actions(panel, DEMO_OUT)
    playbook_scatter(panel, DEMO_OUT)

    meta = {
        "source": source,
        "n": len(rows),
        "pos": int(y.sum()),
        "leakage": df.to_dict(orient="records"),
        "note": "Fixture funnel counts are schematic unless source is curnow. Playbook fixture uses synthetic family labels, not TargetTrack.",
    }
    (DEMO_OUT / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"demo source={source} n={len(rows)} wrote {DEMO_OUT}")
    if not df.empty:
        print(df.to_string(index=False))
