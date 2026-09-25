#!/bin/bash --login
#SBATCH --job-name=atlas_notech
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=run/atlas_notech_%j.out
#SBATCH --account=bioinformaticscore

# Knits 13b_atlas_no_technical_genes.Rmd, which is 13 with the mitochondrial,
# rRNA and ribosomal protein genes removed from the objects before anything
# runs. It is a comparison against 13, not a replacement for it.
#
# Nothing is shared but the inputs, so the two runs cannot overwrite each other:
#
#   results/<date>_atlas                       13
#   results/<date>_atlas_no_technical_genes    13b, including its own Robjects
#
# 13 keeps its cached mnn.sce, merged.sce and combined.sce in <date>_no_doublets.
# 13b writes its own into its results directory, because they describe a
# different gene set. all.sce and the doublet statistics are shared inputs and
# are read, never written.
#
# To run 17 on this atlas instead of the default one:
#
#   ATLAS_SUFFIX=_no_technical_genes sbatch src/knit_17_tscan.sh
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="${RESULTS_DATE:-20260813}"
K="${K:-30}"
OUT_DIR="${PROJECT_DIR}/html"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

ALL_SCE="${PROJECT_DIR}/results/${RESULTS_DATE}_no_doublets/all.sce.Rds"
[[ -f "${ALL_SCE}" ]] || {
  echo "missing ${ALL_SCE}" >&2
  echo "Knit 13_atlas_and_stage_correspondence.Rmd first; 13b reuses its objects." >&2
  exit 1; }

GENE_DESC="${PROJECT_DIR}/results/${RESULTS_DATE}_no_doublets/gene_descritption.csv"
[[ -f "${GENE_DESC}" ]] || {
  echo "missing ${GENE_DESC}" >&2
  echo "13b reads it to decide which genes are technical." >&2
  exit 1; }

Rscript -e "
  need = c('rmarkdown','tidyverse','SingleCellExperiment','scran','scater',
           'scuttle','batchelor','BiocSingular','scDblFinder','DT','patchwork')
  miss = Filter(function(p) !requireNamespace(p, quietly = TRUE), need)
  if (length(miss) > 0) {
    message('missing R packages: ', paste(miss, collapse = ', ')); quit(status = 1)
  }
"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "host:    $(hostname)"
echo "started: $(date)"
echo "k:       ${K}"

Rscript -e "
  rmarkdown::render(
    input         = '${PROJECT_DIR}/src/13b_atlas_no_technical_genes.Rmd',
    output_dir    = '${OUT_DIR}',
    output_file   = '13b_atlas_no_technical_genes.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir  = '${PROJECT_DIR}',
                         results_date = '${RESULTS_DATE}',
                         k            = ${K}),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/13b_atlas_no_technical_genes.html"
echo "results:  ${PROJECT_DIR}/results/${RESULTS_DATE}_atlas_no_technical_genes"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
