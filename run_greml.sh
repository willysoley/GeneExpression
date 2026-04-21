#!/usr/bin/env bash
#SBATCH -J NF_GREML
#SBATCH -p mostafavilab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o nextflow_logs/nextflow_driver.out
#SBATCH -e nextflow_logs/nextflow_driver.err

set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NF_DIR="${BASE}/nf"

# Keep this stable across reruns to preserve Nextflow resume/cache behavior.
export NXF_WORK="${NXF_WORK:-/gpfs/scratch/sl8085/nextflow_work/greml_production}"
mkdir -p "${BASE}/nextflow_logs" "${BASE}/results" "${NXF_WORK}"

module load nextflow

export NXF_LOG_FILE="${BASE}/nextflow_logs/nextflow.log"

cd "${BASE}"

echo "Launching GREML workflow"
echo "Project dir : ${BASE}"
echo "Work dir    : ${NXF_WORK}"
echo "Config file : ${NF_DIR}/nextflow.config"

nextflow run "${NF_DIR}/main.nf" \
  -c "${NF_DIR}/nextflow.config" \
  -profile slurm \
  -resume \
  "$@"
