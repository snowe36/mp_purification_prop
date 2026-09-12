from __future__ import annotations

import time
from pathlib import Path

import requests

TIMEOUT = 60
HEADERS = {"User-Agent": "mp-atlas/0.1 (research; +https://github.com/snowe36)"}


def get(url: str, dest: Path | None = None, retries: int = 2) -> bytes:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=HEADERS, timeout=(10, TIMEOUT), stream=True)
            r.raise_for_status()
            chunks = []
            for block in r.iter_content(1 << 20):
                chunks.append(block)
            data = b"".join(chunks)
            if dest is not None:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1 * (attempt + 1))
    raise RuntimeError(f"GET failed {url}: {last}") from last
