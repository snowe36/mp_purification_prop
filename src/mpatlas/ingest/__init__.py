from __future__ import annotations

from mpatlas.ingest import curnow, gpcrdb, mpstruc, purificationdb, screens, targettrack, uniprot, unitmp
from mpatlas.paths import ensure_dirs


def download_all(*, skip_targettrack: bool = False) -> dict[str, str]:
    ensure_dirs()
    status: dict[str, str] = {}

    try:
        print("curnow ...", flush=True)
        curnow.download()
        status["curnow"] = "ok"
    except Exception as exc:  # noqa: BLE001
        status["curnow"] = f"fail: {exc}"

    try:
        print("uniprot ...", flush=True)
        uniprot.download()
        status["uniprot"] = "ok"
    except Exception as exc:  # noqa: BLE001
        status["uniprot"] = f"fail: {exc}"

    try:
        print("mpstruc ...", flush=True)
        mpstruc.download()
        status["mpstruc"] = "ok"
    except Exception as exc:  # noqa: BLE001
        status["mpstruc"] = f"fail: {exc}"

    try:
        print("unitmp ...", flush=True)
        got = unitmp.download()
        status["unitmp"] = "ok " + ", ".join(f"{k}={v is not None}" for k, v in got.items())
    except Exception as exc:  # noqa: BLE001
        status["unitmp"] = f"fail: {exc}"

    try:
        print("purificationdb ...", flush=True)
        dest = purificationdb.download()
        status["purificationdb"] = "ok" if dest else "missing"
    except Exception as exc:  # noqa: BLE001
        status["purificationdb"] = f"fail: {exc}"

    try:
        print("gpcrdb ...", flush=True)
        gpcrdb.download()
        status["gpcrdb"] = "ok"
    except Exception as exc:  # noqa: BLE001
        status["gpcrdb"] = f"fail: {exc}"

    if skip_targettrack:
        status["targettrack"] = "skipped"
    else:
        try:
            print("targettrack ...", flush=True)
            targettrack.download()
            status["targettrack"] = "ok"
        except Exception as exc:  # noqa: BLE001
            status["targettrack"] = f"fail: {exc}"

    return status
