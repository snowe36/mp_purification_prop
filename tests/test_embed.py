import hashlib
from pathlib import Path

import numpy as np

from mpatlas.embed import EmbedIndex, MAX_LEN, seq_key


def test_seq_key_truncates():
    long = "M" + "A" * (MAX_LEN + 40)
    assert seq_key(long) == seq_key(long[:MAX_LEN])
    assert seq_key("acdef") == hashlib.sha1(b"ACDEF").hexdigest()


def test_index_roundtrip(tmp_path: Path):
    seqs = ["MVLSPADKTNVKAAW", "ACDEFGHIKLMNPQRSTVWY"]
    keys = [seq_key(s) for s in seqs]
    X = np.arange(8, dtype=np.float16).reshape(2, 4)
    path = tmp_path / "e.npz"
    np.savez(path, keys=np.array(keys), X=X, model=np.array(["test"]))
    ix = EmbedIndex(path)
    assert ix.has(seqs[0])
    assert ix.has(seqs[0] + "AAAA") is False
    assert ix.coverage(seqs) == 1.0
    mat = ix.matrix(seqs)
    assert mat is not None
    assert mat.shape == (2, 4)
