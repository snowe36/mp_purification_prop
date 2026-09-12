from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
FIXTURES = DATA / "fixtures"
REPORTS = ROOT / "reports"
FIGURES = REPORTS / "figures"
DEMO = ROOT / "demo"
DEMO_OUT = DEMO / "outputs"


def ensure_dirs() -> None:
    for path in (RAW, PROCESSED, FIGURES, DEMO_OUT, FIXTURES):
        path.mkdir(parents=True, exist_ok=True)
