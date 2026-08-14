#!/bin/bash --login
#SBATCH --job-name=knit_12
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=run/knit_12_%j.out
#SBATCH --account=bioinformaticscore

# Knits 12_doublets_versus_intermediates.Rmd.
#
# Memory is the binding constraint rather than CPU. all.sce, merged.sce and
# combined.sce are loaded together, and the four scDblFinder runs each build a
# kNN graph over real plus artificial cells. 128G is generous; check the
# seff output below and trim on a rerun.
#
# Time is dominated by the scDblFinder runs on the first knit. Afterwards the
# calls are cached in the results Robjects directory and a reknit is minutes,
# so the 8 hour request only applies until that cache exists.

module purge
module load R-bundle-Bioconductor        # adjust to the exact module on HPCC

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RMD="${PROJECT_DIR}/src/12_doublets_versus_intermediates.Rmd"
OUT_DIR="${PROJECT_DIR}/html"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

# fail loudly rather than producing a partial html
set -euo pipefail

echo "host:    $(hostname)"
echo "started: $(date)"
echo "Rmd:     ${RMD}"

# one worker per allocated core, and BLAS held to one thread so the two do not
# oversubscribe each other
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript -e "
  rmarkdown::render(
    input       = '${RMD}',
    output_dir  = '${OUT_DIR}',
    knit_root_dir = '${PROJECT_DIR}',
    params      = list(project_dir = '${PROJECT_DIR}'),
    envir       = new.env()
  )
"

echo "finished: $(date)"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
