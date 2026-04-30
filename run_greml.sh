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
#   run root     = submit dir
# Overrides:
#   GREML_PROJECT_ROOT=/abs/path/to/GeneExpression
#   GREML_RUN_DIR=/abs/path/for/results_and-logs
# Behavior:
#   by default each parameter combination gets its own launch dir under:
#     ${GREML_RUN_DIR}/runs/<combo-label>
#   set GREML_ISOLATE_BY_COMBO=0 to keep legacy shared launch dir behavior
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

sha256_text() {
  local text="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    printf '%s' "${text}" | sha256sum | awk '{print $1}'
    return
  fi
  if command -v shasum >/dev/null 2>&1; then
    printf '%s' "${text}" | shasum -a 256 | awk '{print $1}'
    return
  fi
  if command -v openssl >/dev/null 2>&1; then
    printf '%s' "${text}" | openssl dgst -sha256 | awk '{print $NF}'
    return
  fi
  echo "ERROR: No SHA-256 tool found (need sha256sum, shasum, or openssl)." >&2
  exit 2
}

to_lower() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]'
}

sanitize_slug() {
  local raw="$1"
  raw="$(to_lower "${raw}")"
  raw="$(printf '%s' "${raw}" | sed -E 's/[^a-z0-9._-]+/_/g; s/_+/_/g; s/^_+//; s/_+$//')"
  if [[ -z "${raw}" ]]; then
    raw="na"
  fi
  printf '%s' "${raw}"
}

normalize_bool() {
  local raw
  raw="$(to_lower "$1")"
  case "${raw}" in
    1|true|t|yes|y|on) printf 'true' ;;
    0|false|f|no|n|off) printf 'false' ;;
    *)
      echo "ERROR: Invalid boolean value '${1}'." >&2
      exit 2
      ;;
  esac
}

# Parse pipeline parameter values from CLI args.
# Supports --key=value and --key value, and treats --flag as true.
extract_param_value() {
  local default="$1"
  shift
  local value="${default}"
  local i arg key next

  for ((i=0; i<${#CLI_ARGS[@]}; i++)); do
    arg="${CLI_ARGS[$i]}"
    for key in "$@"; do
      if [[ "${arg}" == "--${key}" ]]; then
        if (( i + 1 < ${#CLI_ARGS[@]} )); then
          next="${CLI_ARGS[$((i+1))]}"
          if [[ "${next}" == --* ]]; then
            value="true"
          else
            value="${next}"
            i=$((i + 1))
          fi
        else
          value="true"
        fi
        break
      fi
      if [[ "${arg}" == "--${key}="* ]]; then
        value="${arg#*=}"
        break
      fi
    done
  done

  printf '%s' "${value}"
}

canonicalize_cli_args() {
  local -a normalized=()
  local -a sorted=()
  local i arg next line

  for ((i=0; i<${#CLI_ARGS[@]}; i++)); do
    arg="${CLI_ARGS[$i]}"
    if [[ "${arg}" == --*=* ]]; then
      normalized+=("${arg}")
      continue
    fi
    if [[ "${arg}" == --* ]]; then
      if (( i + 1 < ${#CLI_ARGS[@]} )) && [[ "${CLI_ARGS[$((i+1))]}" != --* ]]; then
        normalized+=("${arg}=${CLI_ARGS[$((i+1))]}")
        i=$((i + 1))
      else
        normalized+=("${arg}=true")
      fi
      continue
    fi
    normalized+=("ARG:${arg}")
  done

  if [[ ${#normalized[@]} -eq 0 ]]; then
    return 0
  fi

  while IFS= read -r line; do
    sorted+=("${line}")
  done < <(printf '%s\n' "${normalized[@]}" | LC_ALL=C sort)

  printf '%s\n' "${sorted[@]}"
}

split_csv_to_array() {
  local csv="$1"
  local -n out_arr="$2"
  local item
  local -a cleaned=()
  out_arr=()
  IFS=',' read -r -a out_arr <<< "${csv}"
  for item in "${out_arr[@]}"; do
    item="$(printf '%s' "${item}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
    if [[ -n "${item}" ]]; then
      cleaned+=("${item}")
    fi
  done
  out_arr=("${cleaned[@]}")
}

mtime_utc() {
  local target="$1"
  if stat --version >/dev/null 2>&1; then
    date -u -d "@$(stat -c '%Y' "${target}")" '+%Y-%m-%dT%H:%M:%SZ'
    return
  fi
  date -u -r "$(stat -f '%m' "${target}")" '+%Y-%m-%dT%H:%M:%SZ'
}

RUN_ROOT="$(resolve_dir "${GREML_RUN_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}")"
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

CLI_ARGS=("$@")

FANOUT_CHILD="${GREML_FANOUT_CHILD:-0}"
if [[ "${FANOUT_CHILD}" == "1" ]]; then
  fanout_enabled="false"
else
  fanout_mode="$(to_lower "${GREML_FANOUT:-auto}")"
  case "${fanout_mode}" in
    auto)
      if [[ ${#CLI_ARGS[@]} -eq 0 ]]; then
        fanout_enabled="true"
      else
        fanout_enabled="false"
      fi
      ;;
    1|true|t|yes|y|on)
      fanout_enabled="true"
      ;;
    0|false|f|no|n|off)
      fanout_enabled="false"
      ;;
    *)
      echo "ERROR: GREML_FANOUT must be one of: auto,true,false,1,0 (got '${GREML_FANOUT}')." >&2
      exit 2
      ;;
  esac
fi

if [[ "${fanout_enabled}" == "true" ]]; then
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: GREML_FANOUT is enabled but sbatch is not available in PATH." >&2
    exit 2
  fi

  fanout_sleep_sec="${GREML_FANOUT_INTERVAL_SEC:-1}"
  if ! [[ "${fanout_sleep_sec}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: GREML_FANOUT_INTERVAL_SEC must be a non-negative integer (got '${fanout_sleep_sec}')." >&2
    exit 2
  fi

  fanout_snp_sets_raw="${GREML_FANOUT_SNP_SETS:-hm3_no_mhc,all_snps}"
  fanout_expr_raw="${GREML_FANOUT_EXPRESSION_SOURCES:-tpm,tmm}"
  fanout_norm_raw="${GREML_FANOUT_NORMALIZATION_TYPES:-irnt,inverse_normal,raw}"

  split_csv_to_array "${fanout_snp_sets_raw}" fanout_snp_sets
  split_csv_to_array "${fanout_expr_raw}" fanout_expr_sources
  split_csv_to_array "${fanout_norm_raw}" fanout_norm_types

  if [[ ${#fanout_snp_sets[@]} -eq 0 || ${#fanout_expr_sources[@]} -eq 0 || ${#fanout_norm_types[@]} -eq 0 ]]; then
    echo "ERROR: One or more GREML_FANOUT_* lists were empty." >&2
    exit 2
  fi

  mkdir -p "${RUN_ROOT}/nextflow_logs"
  echo "Fanout mode: submitting matrix runs from one driver job"
  echo "  SNP sets            : ${fanout_snp_sets[*]}"
  echo "  Expression sources  : ${fanout_expr_sources[*]}"
  echo "  Normalization types : ${fanout_norm_types[*]}"
  echo "  Submit spacing      : ${fanout_sleep_sec}s"

  submitted_n=0
  for snp_item in "${fanout_snp_sets[@]}"; do
    snp_item_norm="$(sanitize_slug "${snp_item}")"
    case "${snp_item_norm}" in
      hm3|hm3_no_mhc|hm3_no_hla|true|t|1|yes|y|on)
        hm3_bool="true"
        snp_label="hm3_no_mhc"
        ;;
      all|all_snps|false|f|0|no|n|off)
        hm3_bool="false"
        snp_label="all_snps"
        ;;
      *)
        echo "ERROR: Unknown SNP set token in GREML_FANOUT_SNP_SETS: '${snp_item}'" >&2
        echo "Allowed tokens include: hm3_no_mhc, all_snps" >&2
        exit 2
        ;;
    esac

    for expr_item in "${fanout_expr_sources[@]}"; do
      expr_norm="$(to_lower "${expr_item}")"
      if [[ "${expr_norm}" != "tpm" && "${expr_norm}" != "tmm" ]]; then
        echo "ERROR: Invalid expression source in GREML_FANOUT_EXPRESSION_SOURCES: '${expr_item}'" >&2
        exit 2
      fi

      for norm_item in "${fanout_norm_types[@]}"; do
        norm_norm="$(to_lower "${norm_item}")"
        case "${norm_norm}" in
          irnt|inverse_normal|raw) ;;
          ukb_irnt) norm_norm="irnt" ;;
          tmm_only) norm_norm="raw" ;;
          *)
            echo "ERROR: Invalid normalization type in GREML_FANOUT_NORMALIZATION_TYPES: '${norm_item}'" >&2
            exit 2
            ;;
        esac

        combo_label="$(sanitize_slug "${snp_label}_${expr_norm}_${norm_norm}")"
        child_cmd=(
          sbatch
          --job-name "GREML_${combo_label}"
          --output "${RUN_ROOT}/nextflow_logs/nextflow_driver-${combo_label}-%j.out"
          --error "${RUN_ROOT}/nextflow_logs/nextflow_driver-${combo_label}-%j.err"
          --export "ALL,GREML_FANOUT=0,GREML_FANOUT_CHILD=1"
          "${PROJECT_DIR}/run_greml.sh"
          "${CLI_ARGS[@]}"
          --use_hm3_no_hla "${hm3_bool}"
          --expression_source "${expr_norm}"
          --normalization_type "${norm_norm}"
        )
        submit_out="$("${child_cmd[@]}")"
        submitted_n=$((submitted_n + 1))
        echo "[fanout:${submitted_n}] ${combo_label} -> ${submit_out}"

        if (( fanout_sleep_sec > 0 )); then
          sleep "${fanout_sleep_sec}"
        fi
      done
    done
  done

  echo "Fanout submission complete: ${submitted_n} jobs."
  exit 0
fi

expression_source="$(to_lower "$(extract_param_value "tmm" expression_source expression-source expressionSource)")"
normalization_type="$(to_lower "$(extract_param_value "irnt" normalization_type normalization-type normalizationType)")"
use_hm3_no_hla="$(normalize_bool "$(extract_param_value "true" use_hm3_no_hla use-hm3-no-hla useHm3NoHla)")"
peer_nk="$(extract_param_value "auto" peer_nk peer-nk peerNk)"
peer_max_genes="$(extract_param_value "0" peer_max_genes peer-max-genes peerMaxGenes)"
n_pcs="$(extract_param_value "5" n_pcs n-pcs nPcs)"
if [[ "${normalization_type}" == "ukb_irnt" ]]; then
  normalization_type="irnt"
fi
if [[ "${normalization_type}" == "tmm_only" ]]; then
  normalization_type="raw"
fi
if [[ "${use_hm3_no_hla}" == "true" ]]; then
  snp_mode="hm3_no_mhc"
else
  snp_mode="all_snps"
fi
combo_label_default="$(sanitize_slug "${snp_mode}_${expression_source}_${normalization_type}_peer${peer_nk}_pmg${peer_max_genes}_npc${n_pcs}")"
combo_hash_material="$(canonicalize_cli_args)"
if [[ -z "${combo_hash_material}" ]]; then
  combo_hash_material="<defaults>"
fi
combo_hash="$(sha256_text "${combo_hash_material}")"
combo_hash_short="${combo_hash:0:12}"

ISOLATE_BY_COMBO="${GREML_ISOLATE_BY_COMBO:-1}"
if [[ "${ISOLATE_BY_COMBO}" == "0" ]]; then
  RUN_DIR="${RUN_ROOT}"
  RUN_KEY="shared_launchdir"
else
  if [[ "${ISOLATE_BY_COMBO}" != "1" ]]; then
    echo "ERROR: GREML_ISOLATE_BY_COMBO must be 0 or 1 (got '${ISOLATE_BY_COMBO}')." >&2
    exit 2
  fi
  combo_label="$(sanitize_slug "${GREML_COMBO_LABEL:-${combo_label_default}}")"
  RUN_KEY="${combo_label}"
  RUN_DIR="${RUN_ROOT}/runs/${RUN_KEY}"
fi

# Keep work/cache stable per RUN_KEY to preserve resume semantics while avoiding
# cross-combo contention in a shared scratch root.
NXF_WORK_ROOT="${NXF_WORK_ROOT:-/gpfs/scratch/sl8085/nextflow_work/greml_production}"
export NXF_WORK="${NXF_WORK:-${NXF_WORK_ROOT}/${RUN_KEY}}"
mkdir -p "${RUN_DIR}/.nextflow" "${RUN_DIR}/nextflow_logs" "${RUN_DIR}/results" "${NXF_WORK}"

module load nextflow

JOB_TAG="${SLURM_JOB_ID:-manual}"
export NXF_LOG_FILE="${RUN_DIR}/nextflow_logs/nextflow_${JOB_TAG}.log"

RESUME_FLAG="${GREML_RESUME:-1}"
if [[ "${ISOLATE_BY_COMBO}" == "0" && "${RESUME_FLAG}" != "0" && "${GREML_ALLOW_SHARED_RESUME:-0}" != "1" ]]; then
  echo "ERROR: GREML_ISOLATE_BY_COMBO=0 with resume enabled is unsafe and can trigger Nextflow LOCK collisions." >&2
  echo "Set GREML_ISOLATE_BY_COMBO=1 (recommended), or GREML_RESUME=0 for a fresh non-resume run." >&2
  exit 2
fi

# Serialize launches per RUN_DIR to prevent concurrent -resume collisions.
LOCK_WAIT_SEC="${GREML_RUN_LOCK_WAIT_SEC:-0}"
if ! [[ "${LOCK_WAIT_SEC}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: GREML_RUN_LOCK_WAIT_SEC must be a non-negative integer (got '${LOCK_WAIT_SEC}')." >&2
  exit 2
fi

LAUNCH_LOCK_FILE="${RUN_DIR}/.nextflow/.launch.lock"
LOCK_ACQUIRED="0"
LOCK_IMPL=""
LOCK_DIR=""

cleanup_launch_lock() {
  if [[ "${LOCK_IMPL}" == "mkdir" && "${LOCK_ACQUIRED}" == "1" && -n "${LOCK_DIR}" ]]; then
    rm -f "${LOCK_DIR}/owner.txt" 2>/dev/null || true
    rmdir "${LOCK_DIR}" 2>/dev/null || true
  fi
}
trap cleanup_launch_lock EXIT

if command -v flock >/dev/null 2>&1; then
  LOCK_IMPL="flock"
  exec 9>"${LAUNCH_LOCK_FILE}"
  if (( LOCK_WAIT_SEC > 0 )); then
    flock -w "${LOCK_WAIT_SEC}" 9 || {
      echo "ERROR: timed out waiting ${LOCK_WAIT_SEC}s for launch lock ${LAUNCH_LOCK_FILE}." >&2
      exit 75
    }
  else
    flock -n 9 || {
      echo "ERROR: another launch is active for RUN_KEY=${RUN_KEY}; refusing concurrent start in ${RUN_DIR}." >&2
      echo "Set GREML_RUN_LOCK_WAIT_SEC to wait instead of failing fast." >&2
      exit 75
    }
  fi
  LOCK_ACQUIRED="1"
else
  LOCK_IMPL="mkdir"
  LOCK_DIR="${RUN_DIR}/.nextflow/.launch.lock.d"
  START_EPOCH="$(date +%s)"
  while true; do
    if mkdir "${LOCK_DIR}" 2>/dev/null; then
      LOCK_ACQUIRED="1"
      {
        echo "host=$(hostname)"
        echo "pid=$$"
        echo "run_key=${RUN_KEY}"
        echo "started_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
      } > "${LOCK_DIR}/owner.txt"
      break
    fi
    if (( LOCK_WAIT_SEC == 0 )); then
      echo "ERROR: another launch is active for RUN_KEY=${RUN_KEY}; refusing concurrent start in ${RUN_DIR}." >&2
      echo "Set GREML_RUN_LOCK_WAIT_SEC to wait instead of failing fast." >&2
      exit 75
    fi
    NOW_EPOCH="$(date +%s)"
    if (( NOW_EPOCH - START_EPOCH >= LOCK_WAIT_SEC )); then
      echo "ERROR: timed out waiting ${LOCK_WAIT_SEC}s for launch lock ${LOCK_DIR}." >&2
      exit 75
    fi
    sleep 1
  done
fi

cd "${RUN_DIR}"

echo "Launching GREML workflow"
echo "Submit dir  : ${SLURM_SUBMIT_DIR:-<unset>}"
echo "Project dir : ${PROJECT_DIR}"
echo "Run root    : ${RUN_ROOT}"
echo "Run key     : ${RUN_KEY}"
echo "Combo label : ${combo_label_default}"
echo "Combo hash  : ${combo_hash_short}"
echo "Isolation   : $([[ "${ISOLATE_BY_COMBO}" == "1" ]] && echo "enabled" || echo "disabled")"
echo "Run dir     : ${RUN_DIR}"
echo "Work dir    : ${NXF_WORK}"
echo "Config file : ${NF_DIR}/nextflow.config"
echo "Launch lock : ${LOCK_IMPL} (wait=${LOCK_WAIT_SEC}s)"

NEXTFLOW_CMD=(
  nextflow -log "${NXF_LOG_FILE}" run "${NF_DIR}/main.nf"
  -c "${NF_DIR}/nextflow.config"
  -profile slurm
  -w "${NXF_WORK}"
)

if [[ "${RESUME_FLAG}" == "0" ]]; then
  echo "Resume mode : disabled (GREML_RESUME=0)"
else
  RESUME_SESSION="${GREML_RESUME_SESSION:-}"
  if [[ -n "${RESUME_SESSION}" ]]; then
    echo "Resume mode : enabled (-resume ${RESUME_SESSION})"
    NEXTFLOW_CMD+=(-resume "${RESUME_SESSION}")
  else
    echo "Resume mode : enabled (-resume)"
    NEXTFLOW_CMD+=(-resume)
  fi
fi

"${NEXTFLOW_CMD[@]}" "$@"
