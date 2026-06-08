#!/bin/bash
#SBATCH --job-name=chunkwperm50k
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --account=uoa04097
#SBATCH --output=logs/chunkwperm50k_%A_%a.out
#SBATCH --error=logs/chunkwperm50k_%A_%a.err

set -euo pipefail

# NOTE: this script calls figure-scripts/chunk-withinperm-worker-50k.R which is
# not included in this repository. It is retained here for reference only.
# The local alternative is to run figure-scripts/analyse-chunk-withinperm-nulls-50k.R
# directly after extracting the Zenodo hpc-nulls.tar.gz, which re-runs the
# permutations locally from data/chunk_withinperm_base/chunk_withinperm_base_50k.rds.
#
# Usage (if worker script is available on HPC):
#   sbatch --array=1-200 --export=ALL,N_PERM=100 slurms/chunk-withinperm-run.slurm
# Mapping:
#   1..N_PERM          -> chunk250kb_withinperm iter 1..N_PERM
#   N_PERM+1..2*N_PERM -> chunk500kb_withinperm iter 1..N_PERM

N_PERM="${N_PERM:-100}"
BUNDLE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
BASE_FILE="${BUNDLE_DIR}/results/chunk-withinperm-nulls-50k/hpc_base/chunk_withinperm_base_50k.rds"

mkdir -p "${BUNDLE_DIR}/logs"
cd "${BUNDLE_DIR}"

if [[ ! -f "${BASE_FILE}" ]]; then
  echo "Missing base file: ${BASE_FILE}"
  echo "Run prepare job first:"
  echo "  sbatch slurms/chunk-withinperm-prepare.slurm"
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID}"
MAX_TASK=$((2 * N_PERM))
if (( TASK_ID < 1 || TASK_ID > MAX_TASK )); then
  echo "Task ${TASK_ID} outside expected range 1..${MAX_TASK}; exiting."
  exit 0
fi

if (( TASK_ID <= N_PERM )); then
  SCHEME="chunk250kb_withinperm"
  ITER="${TASK_ID}"
else
  SCHEME="chunk500kb_withinperm"
  ITER=$((TASK_ID - N_PERM))
fi

module purge || true
module load R || true

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "Task ${TASK_ID}/${MAX_TASK}: ${SCHEME} iter ${ITER}"
Rscript figure-scripts/chunk-withinperm-worker-50k.R "${SCHEME}" "${ITER}"
echo "Done ${SCHEME} iter ${ITER}"
