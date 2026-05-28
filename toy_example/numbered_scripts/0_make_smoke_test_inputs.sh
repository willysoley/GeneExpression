#!/usr/bin/env bash
set -euo pipefail

TOY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${TOY_DIR}/smoke_test_data"
SOURCE_SDRF="/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260106_0_test_GREML/E-GEUV-1.sdrf.txt"

mkdir -p "${DATA_DIR}"
cp "${SOURCE_SDRF}" "${DATA_DIR}/E-GEUV-1.sdrf.txt"

TOY_DIR="${TOY_DIR}" python3 - <<'PY'
from pathlib import Path
from collections import OrderedDict
import csv
import os

toy_dir = Path(os.environ["TOY_DIR"])
data_dir = toy_dir / "smoke_test_data"
sdrf_file = data_dir / "E-GEUV-1.sdrf.txt"
fam_file = Path("/gpfs/data/mostafavilab/sool/analysis/GeneExpression/20260106_0_test_GREML/geno_chr22_keep.fam")
out_tpm = data_dir / "toy_tpm.tsv"
out_counts = data_dir / "toy_counts.tsv"

with fam_file.open() as handle:
    fam_ids = [line.split()[0] for line in handle if line.strip()]

with sdrf_file.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    mapping = OrderedDict()
    for row in reader:
        ancestry = row.get("Characteristics[ancestry category]", "")
        if ancestry not in {"British", "Finnish", "Tuscan", "Utah"}:
            continue
        source_name = row.get("Source Name", "").strip()
        ena_run = row.get("Comment[ENA_RUN]", "").strip()
        if source_name and ena_run and source_name not in mapping:
            mapping[source_name] = ena_run

samples = []
missing = []
for fid in fam_ids:
    if fid in mapping:
        samples.append((fid, mapping[fid]))
    else:
        missing.append(fid)

if missing:
    raise SystemExit("Missing SDRF mapping for FAM IDs: " + ", ".join(missing[:10]))

sample_cols = [f"{ena_run}_1" for _, ena_run in samples]

genes = ["TOY_GENE_A", "TOY_GENE_B", "TOY_GENE_C", "TOY_GENE_D"]
with out_tpm.open("w") as tpm, out_counts.open("w") as counts:
    header = "gene_id\t" + "\t".join(sample_cols) + "\n"
    tpm.write(header)
    counts.write(header)

    for gene_index, gene in enumerate(genes):
        tpm_values = []
        count_values = []
        for sample_index in range(len(sample_cols)):
            if gene_index == 0:
                base = sample_index + 1
            elif gene_index == 1:
                base = len(sample_cols) - sample_index
            elif gene_index == 2:
                base = 3 + (sample_index % 4)
            else:
                base = 2 + ((sample_index + gene_index) % 3)
            tpm_values.append(f"{base:.2f}")
            count_values.append(str(base * 10 + gene_index))
        tpm.write(gene + "\t" + "\t".join(tpm_values) + "\n")
        counts.write(gene + "\t" + "\t".join(count_values) + "\n")

print(f"Wrote {out_tpm} and {out_counts} with {len(genes)} genes and {len(sample_cols)} samples")
PY
