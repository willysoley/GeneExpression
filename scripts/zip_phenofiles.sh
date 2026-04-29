#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  zip_phenofiles.sh [--data-dir DIR] [--out FILE]

Description:
  Archives phenotype artifacts into a single tar.gz bundle.
  Included patterns:
    *.phenotypes.tsv
    *.qcovar
    *.gene_index_map.txt
    *.filtered_gene_ids.txt
EOF
}

DATA_DIR="results/data"
OUT_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir)
      DATA_DIR="$2"
      shift 2
      ;;
    --out)
      OUT_FILE="$2"
      shift 2
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

DATA_DIR="$(cd "$DATA_DIR" && pwd -P)"
if [[ ! -d "$DATA_DIR" ]]; then
  echo "Data directory not found: $DATA_DIR" >&2
  exit 1
fi

if [[ -z "$OUT_FILE" ]]; then
  OUT_FILE="${DATA_DIR}/phenofiles_$(date +%Y%m%d_%H%M%S).tar.gz"
fi

tmp_list="$(mktemp)"
trap 'rm -f "$tmp_list"' EXIT

find "$DATA_DIR" -maxdepth 1 -type f \
  \( -name '*.phenotypes.tsv' -o -name '*.qcovar' -o -name '*.gene_index_map.txt' -o -name '*.filtered_gene_ids.txt' \) \
  -print | sort > "$tmp_list"

if [[ ! -s "$tmp_list" ]]; then
  echo "No phenotype files found in: $DATA_DIR" >&2
  exit 1
fi

echo "Creating archive: $OUT_FILE"
tar -czf "$OUT_FILE" -T "$tmp_list"
echo "Done."
