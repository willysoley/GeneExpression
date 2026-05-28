#!/usr/bin/env bash
set -euo pipefail

TOY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "${TOY_DIR}/.." && pwd)"

module load nextflow >/dev/null 2>&1 || true

nextflow run "${REPO_ROOT}/nf/main.nf" \
  -c "${TOY_DIR}/toy_smoke_test.config" \
  -resume
