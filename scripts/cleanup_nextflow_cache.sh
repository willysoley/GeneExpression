#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  cleanup_nextflow_cache.sh [--run-dir DIR] [--nxf-work DIR] [--delete-results] [--execute]

Description:
  Cleans Nextflow cache/log/work artifacts for this GREML workflow.
  Default mode is DRY RUN (prints what would be removed).

Options:
  --run-dir DIR       Run directory containing .nextflow/, nextflow_logs/, results/ (default: $PWD)
  --nxf-work DIR      Nextflow work directory to remove (default: $NXF_WORK or config default)
  --delete-results    Also remove GREML result subdirs: results/data, results/pca, results/summary
  --execute           Actually delete files/directories
  -h, --help          Show this help
EOF
}

RUN_DIR="$(pwd -P)"
NXF_WORK_DIR="${NXF_WORK:-/gpfs/scratch/sl8085/nextflow_work/greml_production}"
DELETE_RESULTS=0
EXECUTE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir)
      RUN_DIR="$2"
      shift 2
      ;;
    --nxf-work)
      NXF_WORK_DIR="$2"
      shift 2
      ;;
    --delete-results)
      DELETE_RESULTS=1
      shift
      ;;
    --execute)
      EXECUTE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

RUN_DIR="$(cd "$RUN_DIR" && pwd -P)"

targets=()
for p in \
  "${RUN_DIR}/.nextflow" \
  "${RUN_DIR}/nextflow_logs" \
  "${RUN_DIR}/work" \
  "${NXF_WORK_DIR}"; do
  if [[ -e "$p" ]]; then
    targets+=("$p")
  fi
done

if [[ $DELETE_RESULTS -eq 1 ]]; then
  for p in \
    "${RUN_DIR}/results/data" \
    "${RUN_DIR}/results/pca" \
    "${RUN_DIR}/results/summary"; do
    if [[ -e "$p" ]]; then
      targets+=("$p")
    fi
  done
fi

mapfile -t file_targets < <(find "${RUN_DIR}" -maxdepth 1 -type f \
  \( -name '.nextflow.log*' -o -name '.nextflow.history*' -o -name '.nxf*' -o -name 'trace*.txt' -o -name 'timeline*.html' -o -name 'report*.html' -o -name 'dag*.svg' \) \
  2>/dev/null | sort)

if [[ ${#file_targets[@]} -gt 0 ]]; then
  targets+=("${file_targets[@]}")
fi

if [[ ${#targets[@]} -eq 0 ]]; then
  echo "No cache/log/work targets found."
  exit 0
fi

echo "Run dir     : ${RUN_DIR}"
echo "NXF work dir: ${NXF_WORK_DIR}"
echo "Mode        : $([[ $EXECUTE -eq 1 ]] && echo 'EXECUTE' || echo 'DRY RUN')"
echo ""
echo "Targets:"
printf '  %s\n' "${targets[@]}"

if [[ $EXECUTE -eq 0 ]]; then
  echo ""
  echo "Dry run only. Re-run with --execute to delete."
  exit 0
fi

echo ""
echo "Deleting targets..."
for t in "${targets[@]}"; do
  rm -rf "$t"
done
echo "Cleanup complete."
