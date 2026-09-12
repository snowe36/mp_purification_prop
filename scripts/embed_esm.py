#!/usr/bin/env python
"""Mean-pool ESM-2 650M over a FASTA. GPU if available. Writes npz shards then a final npz."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer

MODEL_ID = "facebook/esm2_t33_650M_UR50D"
MAX_LEN = 1022


def read_fasta(path: Path) -> list[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".gz") else open
    rows: list[tuple[str, str]] = []
    key, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if key is not None:
                    rows.append((key, "".join(chunks)))
                key = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if key is not None:
            rows.append((key, "".join(chunks)))
    return rows


def mean_pool(last: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    mask = mask.unsqueeze(-1).to(last.dtype)
    return (last * mask).sum(1) / mask.sum(1).clamp(min=1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--start", type=int, default=0)
    args = ap.parse_args()
    rows = read_fasta(Path(args.fasta))
    rows = rows[args.start :]
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"n={len(rows)} device={device} model={MODEL_ID}", flush=True)
    tok = AutoTokenizer.from_pretrained(MODEL_ID)
    model = AutoModel.from_pretrained(MODEL_ID).to(device)
    model.eval()
    keys, vecs = [], []
    with torch.no_grad():
        for i in range(0, len(rows), args.batch):
            chunk = rows[i : i + args.batch]
            seqs = [s[:MAX_LEN] for _, s in chunk]
            batch = tok(
                seqs,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=MAX_LEN + 2,
            )
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model(**batch)
            pooled = mean_pool(out.last_hidden_state, batch["attention_mask"])
            keys.extend(k for k, _ in chunk)
            vecs.append(pooled.float().cpu().numpy())
            if (i // args.batch) % 20 == 0:
                print(f"done {i + len(chunk)}/{len(rows)}", flush=True)
    X = np.concatenate(vecs, axis=0).astype(np.float16)
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest, keys=np.array(keys), X=X, model=np.array([MODEL_ID]))
    print(f"wrote {dest} shape={X.shape}", flush=True)


if __name__ == "__main__":
    main()
