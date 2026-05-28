#!/usr/bin/env bash
set -euo pipefail

TOY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${TOY_DIR}/results/data"
SUMMARY_FILE="${TOY_DIR}/results/summary/final_heritability_summary.tsv"

PHENO_FILE=$(find "${DATA_DIR}" -maxdepth 1 -name '*.phenotypes.tsv' | head -1)
QCOVAR_FILE=$(find "${DATA_DIR}" -maxdepth 1 -name '*.qcovar' | head -1)
MAP_FILE=$(find "${DATA_DIR}" -maxdepth 1 -name '*.gene_index_map.txt' | head -1)

for file in "${PHENO_FILE}" "${QCOVAR_FILE}" "${MAP_FILE}" "${SUMMARY_FILE}"; do
  test -n "${file}"
  test -s "${file}"
done

echo "summary"
wc -l "${SUMMARY_FILE}"
head -10 "${SUMMARY_FILE}"

echo
echo "phenotype"
head -5 "${PHENO_FILE}"

echo
echo "qcovar"
head -5 "${QCOVAR_FILE}"
