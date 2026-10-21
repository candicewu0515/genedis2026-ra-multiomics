#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
"${PY:-python}" "$ROOT/scripts/build_submission.py"
