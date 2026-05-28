#!/usr/bin/env bash
set -euo pipefail

TOY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${TOY_DIR}/.." && pwd)"
CONFIG_PATH="${GREML_TOY_CONFIG:-${TOY_DIR}/toy_smoke_test.config}"
SBATCH_CPUS="${GREML_TOY_SBATCH_CPUS:-4}"
SBATCH_MEM="${GREML_TOY_SBATCH_MEM:-64G}"

if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "ERROR: Missing toy config: ${CONFIG_PATH}" >&2
  exit 2
fi

cd "${REPO_ROOT}"

echo "Submitting toy GEUVADIS GREML run"
echo "  config: ${CONFIG_PATH}"
echo "  sbatch : -c ${SBATCH_CPUS} --mem=${SBATCH_MEM}"

GREML_FANOUT=0 GREML_ISOLATE_BY_COMBO=0 GREML_ALLOW_SHARED_RESUME=1 GREML_NEXTFLOW_CONFIG="${CONFIG_PATH}" sbatch -c "${SBATCH_CPUS}" --mem="${SBATCH_MEM}" run_greml.sh
