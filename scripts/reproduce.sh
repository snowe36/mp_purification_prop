#!/usr/bin/env bash
set -euo pipefail
export MPLCONFIGDIR="${MPLCONFIGDIR:-.mplconfig}"
mpx-demo
pytest -q
