#!/bin/bash --login
#SBATCH --job-name=stable_k
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=run/stable_k_%j.out
#SBATCH --account=bioinformaticscore

# Clustering stability sweep over k, via src/run_stable_k.R.
#
# Needs results/${NO_DBL_DATE}_no_doublets/merged.sce.Rds, written by
# 01_all_cells_clustering_annotation_no_doublets.Rmd. Nothing else.
#
# Cost is 5 values of k times N_ITERATIONS graph builds, each on ~80% of about
# 61,000 cells, plus the pairwise ARI comparisons. Memory is modest because the
# script drops fastMNN's reconstructed assay and keeps only the corrected
# coordinates, so 64G is generous. Time is the uncertain one: the 8 hour
# request is an estimate, not a measurement. Read the seff output at the end
# and trim both on a rerun.
#
# To time it before committing, run once with N_ITERATIONS=3 and K_RANGE=30.
#
# run_stable_k.R installs anything missing into R_LIBS_PROJECT. Compute nodes
# frequently have no outbound network, so do the install once from a login node
# before the first submission:
#
#   cd "${PROJECT_DIR}"
#   module load R/4.4.1-gfbf-2023b
#   INSTALL_ONLY=1 Rscript src/run_stable_k.R
#
# It is a no-op once everything is present, so it is safe to repeat.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"
export PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"

# Where anything the module does not provide gets installed. Left unset, the
# script falls back to R_LIBS_USER and then to a versioned directory under your
# home, taking the first location it can actually write to and reporting which.
# Set this only to override that.
# export R_LIBS_PROJECT="/path/you/can/write/to"
export NO_DBL_DATE="20260813"
export K_RANGE="10,20,30,40,50"
export N_ITERATIONS="10"
export PROP="0.8"
export SEED="42"

mkdir -p "${PROJECT_DIR}/run"

# fail loudly rather than leaving a half-written results directory
set -euo pipefail

MERGED="${PROJECT_DIR}/results/${NO_DBL_DATE}_no_doublets/merged.sce.Rds"
if [[ ! -f "${MERGED}" ]]; then
  echo "missing ${MERGED}" >&2
  echo "Run knit_01_no_doublets.sh first." >&2
  exit 1
fi

echo "host:    $(hostname)"
echo "started: $(date)"
echo "object:  ${MERGED}"
echo "k:       ${K_RANGE}  iterations: ${N_ITERATIONS}  prop: ${PROP}"

# BLAS held to one thread so it does not oversubscribe the allocation
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript "${PROJECT_DIR}/src/run_stable_k.R"

echo "finished: $(date)"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
