#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

cleanup() {
  if [[ "${KEEP_BUILD:-0}" != "1" ]]; then
    make clean
  fi
}

trap cleanup EXIT

echo "[1/2] strict tests"
make strict

echo "[2/2] ASan/UBSan tests"
make sanitize
