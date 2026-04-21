#!/usr/bin/env bash
#SBATCH -J NF_GREML
#SBATCH -p mostafavilab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o nextflow_logs/nextflow_driver.out
#SBATCH -e nextflow_logs/nextflow_driver.err

set -euo pipefail

# Resolve project directory from script location, but write logs/results to the
# directory where sbatch was submitted from.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"
NF_DIR="${PROJECT_DIR}/nf"

# Keep this stable across reruns to preserve Nextflow resume/cache behavior.
export NXF_WORK="${NXF_WORK:-/gpfs/scratch/sl8085/nextflow_work/greml_production}"
mkdir -p "${RUN_DIR}/nextflow_logs" "${RUN_DIR}/results" "${NXF_WORK}"

module load nextflow

export NXF_LOG_FILE="${RUN_DIR}/nextflow_logs/nextflow.log"

cd "${RUN_DIR}"

echo "Launching GREML workflow"
echo "Project dir : ${PROJECT_DIR}"
echo "Run dir     : ${RUN_DIR}"
echo "Work dir    : ${NXF_WORK}"
echo "Config file : ${NF_DIR}/nextflow.config"

nextflow run "${NF_DIR}/main.nf" \
  -c "${NF_DIR}/nextflow.config" \
  -profile slurm \
  -resume \
  "$@"
