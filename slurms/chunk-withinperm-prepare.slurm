#!/bin/bash
#SBATCH --job-name=prepchunk50k
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --account=uoa04097
#SBATCH --output=logs/prepchunk50k_%j.out
#SBATCH --error=logs/prepchunk50k_%j.err

set -euo pipefail

BUNDLE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p "${BUNDLE_DIR}/logs"
cd "${BUNDLE_DIR}"

module purge || true
module load R || true

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

Rscript figure-scripts/prepare-chunk-withinperm-base-50k.R

echo "Done base prep in ${BUNDLE_DIR}"
