import numpy as np

from mpatlas.models import run_curnow, run_transfer


def _toy(n: int = 200, L: int = 40, seed: int = 0):
    rng = np.random.default_rng(seed)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    seqs = ["".join(rng.choice(aa, size=L)) for _ in range(n)]
    y = np.zeros(n, dtype=int)
    y[::2] = 1
    return seqs, y


def test_curnow_compose_returns_features():
    seqs, y = _toy()
    df = run_curnow(seqs, y, mode="compose")
    assert not df.empty
    assert set(df["features"]) == {"compose"}
    assert df["roc_auc"].between(0, 1).all()


def test_transfer_sample_floor():
    seqs, y = _toy(n=200)
    empty = run_transfer(seqs, y, seqs[:10], y[:10], "toy")
    assert empty.empty
    hit = run_transfer(seqs[:100], y[:100], seqs[100:], y[100:], "toy")
    assert not hit.empty
    assert (hit["split"] == "toy").all()
