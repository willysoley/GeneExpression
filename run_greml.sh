#!/usr/bin/env bash
#SBATCH -J NF_GREML
#SBATCH -p mostafavilab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o nextflow_logs/nextflow_driver.out
#SBATCH -e nextflow_logs/nextflow_driver.err

set -euo pipefail

# Never derive project root from ${BASH_SOURCE[0]} in sbatch jobs:
# that path may be Slurm spool (/cm/local/apps/slurm/...).
# Defaults:
#   project root = submit dir
#   run dir      = submit dir
# Overrides:
#   GREML_PROJECT_ROOT=/abs/path/to/GeneExpression
#   GREML_RUN_DIR=/abs/path/for/results_and-logs
resolve_dir() {
  local dir="$1"
  if [[ ! -d "${dir}" ]]; then
    echo "ERROR: directory does not exist: ${dir}" >&2
    exit 2
  fi
  (cd "${dir}" && pwd -P)
}

RUN_DIR="$(resolve_dir "${GREML_RUN_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}")"
PROJECT_DIR="$(resolve_dir "${GREML_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}")"
NF_DIR="${PROJECT_DIR}/nf"

if [[ ! -f "${NF_DIR}/main.nf" || ! -f "${NF_DIR}/nextflow.config" ]]; then
  cat >&2 <<EOF
ERROR: Could not find pipeline files under PROJECT_DIR=${PROJECT_DIR}
Expected:
  ${NF_DIR}/main.nf
  ${NF_DIR}/nextflow.config

Submit from the repository root, or set:
  GREML_PROJECT_ROOT=/absolute/path/to/GeneExpression
EOF
  exit 2
fi

# Keep this stable across reruns to preserve Nextflow resume/cache behavior.
export NXF_WORK="${NXF_WORK:-/gpfs/scratch/sl8085/nextflow_work/greml_production}"
mkdir -p "${RUN_DIR}/nextflow_logs" "${RUN_DIR}/results" "${NXF_WORK}"

module load nextflow

export NXF_LOG_FILE="${RUN_DIR}/nextflow_logs/nextflow.log"

cd "${RUN_DIR}"

echo "Launching GREML workflow"
echo "Submit dir  : ${SLURM_SUBMIT_DIR:-<unset>}"
echo "Project dir : ${PROJECT_DIR}"
echo "Run dir     : ${RUN_DIR}"
echo "Work dir    : ${NXF_WORK}"
echo "Config file : ${NF_DIR}/nextflow.config"

nextflow -log "${NXF_LOG_FILE}" run "${NF_DIR}/main.nf" \
  -c "${NF_DIR}/nextflow.config" \
  -profile slurm \
  -w "${NXF_WORK}" \
  -resume \
  "$@"
