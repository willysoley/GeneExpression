#!/usr/bin/env bash
#SBATCH -J NF_GREML
#SBATCH -p mostafavilab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o nextflow_logs/nextflow_driver-%j.out
#SBATCH -e nextflow_logs/nextflow_driver-%j.err

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

sha256_file() {
  local target="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${target}" | awk '{print $1}'
    return
  fi
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${target}" | awk '{print $1}'
    return
  fi
  if command -v openssl >/dev/null 2>&1; then
    openssl dgst -sha256 "${target}" | awk '{print $NF}'
    return
  fi
  echo "ERROR: No SHA-256 tool found (need sha256sum, shasum, or openssl)." >&2
  exit 2
}

mtime_utc() {
  local target="$1"
  if stat --version >/dev/null 2>&1; then
    date -u -d "@$(stat -c '%Y' "${target}")" '+%Y-%m-%dT%H:%M:%SZ'
    return
  fi
  date -u -r "$(stat -f '%m' "${target}")" '+%Y-%m-%dT%H:%M:%SZ'
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

PREPARE_PHENO_SCRIPT="${NF_DIR}/bin/prepare_phenotypes.R"
if [[ ! -f "${PREPARE_PHENO_SCRIPT}" ]]; then
  echo "ERROR: Missing required script: ${PREPARE_PHENO_SCRIPT}" >&2
  exit 2
fi

PREPARE_SHA256="$(sha256_file "${PREPARE_PHENO_SCRIPT}")"
PREPARE_MTIME_UTC="$(mtime_utc "${PREPARE_PHENO_SCRIPT}")"
echo "prepare_phenotypes.R: ${PREPARE_PHENO_SCRIPT}"
echo "prepare_phenotypes.R mtime_utc: ${PREPARE_MTIME_UTC}"
echo "prepare_phenotypes.R sha256: ${PREPARE_SHA256}"

REQUIRED_PREP_DIAGNOSTICS=(
  "Applying mandatory GTEx-style filter"
  "Running TMM normalization on raw counts"
  "Phenotype prep completed successfully"
)
missing_diag=()
for marker in "${REQUIRED_PREP_DIAGNOSTICS[@]}"; do
  if ! grep -Fq "${marker}" "${PREPARE_PHENO_SCRIPT}"; then
    missing_diag+=("${marker}")
  fi
done
if [[ ${#missing_diag[@]} -gt 0 ]]; then
  {
    echo "ERROR: ${PREPARE_PHENO_SCRIPT} is missing required diagnostic marker(s)."
    printf '  - %s\n' "${missing_diag[@]}"
  } >&2
  exit 2
fi
echo "prepare_phenotypes.R diagnostics preflight: OK (${#REQUIRED_PREP_DIAGNOSTICS[@]} markers)"

# Keep this stable across reruns to preserve Nextflow resume/cache behavior.
export NXF_WORK="${NXF_WORK:-/gpfs/scratch/sl8085/nextflow_work/greml_production}"
mkdir -p "${RUN_DIR}/nextflow_logs" "${RUN_DIR}/results" "${NXF_WORK}"

module load nextflow

JOB_TAG="${SLURM_JOB_ID:-manual}"
export NXF_LOG_FILE="${RUN_DIR}/nextflow_logs/nextflow_${JOB_TAG}.log"

cd "${RUN_DIR}"

echo "Launching GREML workflow"
echo "Submit dir  : ${SLURM_SUBMIT_DIR:-<unset>}"
echo "Project dir : ${PROJECT_DIR}"
echo "Run dir     : ${RUN_DIR}"
echo "Work dir    : ${NXF_WORK}"
echo "Config file : ${NF_DIR}/nextflow.config"

RESUME_FLAG="${GREML_RESUME:-1}"
if [[ "${RESUME_FLAG}" == "0" ]]; then
  echo "Resume mode : disabled (GREML_RESUME=0)"
  nextflow -log "${NXF_LOG_FILE}" run "${NF_DIR}/main.nf" \
    -c "${NF_DIR}/nextflow.config" \
    -profile slurm \
    -w "${NXF_WORK}" \
    "$@"
else
  echo "Resume mode : enabled (-resume)"
  nextflow -log "${NXF_LOG_FILE}" run "${NF_DIR}/main.nf" \
    -c "${NF_DIR}/nextflow.config" \
    -profile slurm \
    -w "${NXF_WORK}" \
    -resume \
    "$@"
fi
