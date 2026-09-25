#!/bin/bash --login
#SBATCH --job-name=tscan_all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=run/tscan_%j.out
#SBATCH --account=bioinformaticscore

# Knits 17_germ_pseudotime_tscan.Rmd for one lineage, by default "all". The html
# goes to html/ and the results to results/<date>_pseudotime_tscan_<lineage>.
#
# The all lineage is the analysis being reported: clusters 5 and 15 to 18, with
# 12_transformer excluded, the subclusters annotated from a previous run, the
# cells annotated as somatic dropped, and the ordering rooted at the primordial
# germ cell subcluster.
#
# Inputs:
#   - results/<date>_atlas/cell_cluster_labels.csv, from 13
#   - results/<date>_no_doublets/all.sce.Rds and gene_descritption.csv
#   - PREV_NODES, the per-cell node assignment of the annotated run
#   - data/<marker_file> and data/<extra_marker_file>, for the identity scoring
#
# PREV_NODES carries the annotation and the somatic-cell drops. Re-running this
# document overwrites pseudotime_per_cell.csv, so it must point at a copy taken
# from the annotated run:
#
#   cp results/<date>_pseudotime_tscan_all/pseudotime_per_cell.csv \
#      results/<date>_pseudotime_tscan_all/pseudotime_per_cell_annotated.csv
#
# To run with the somatic cells kept, for the with-and-without comparison the
# methods need, clear it:
#
#   PREV_NODES= sbatch src/knit_17_tscan.sh
#
# Other lineages:
#
#   LINEAGE=female sbatch src/knit_17_tscan.sh
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="${RESULTS_DATE:-20260813}"
LINEAGE="${LINEAGE:-all}"
ATLAS_SUFFIX="${ATLAS_SUFFIX:-}"
OUT_DIR="${PROJECT_DIR}/html"

# unset rather than absent means "run without the annotation and the drops"
if [[ -z "${PREV_NODES+x}" ]]; then
  PREV_NODES="${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_tscan_all/pseudotime_per_cell_annotated.csv"
fi

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

LABELS="${PROJECT_DIR}/results/${RESULTS_DATE}_atlas${ATLAS_SUFFIX}/cell_cluster_labels.csv"
[[ -f "${LABELS}" ]] || {
  echo "missing ${LABELS}" >&2
  echo "Knit 13_atlas_and_stage_correspondence.Rmd first." >&2
  exit 1; }

SCE="${PROJECT_DIR}/results/${RESULTS_DATE}_no_doublets/all.sce.Rds"
[[ -f "${SCE}" ]] || {
  echo "missing ${SCE}" >&2
  echo "Check RESULTS_DATE." >&2
  exit 1; }

if [[ -z "${PREV_NODES}" ]]; then
  echo "note: PREV_NODES is empty."
  echo "      Subclusters will not be annotated and no cells will be dropped."
  echo "      This is the run to use for the with-and-without comparison."
elif [[ ! -f "${PREV_NODES}" ]]; then
  echo "note: ${PREV_NODES} not found."
  echo "      17 will report this and run without the annotation or the drops."
  echo "      Copy the annotated run's pseudotime_per_cell.csv there first."
fi

# fail before the knit rather than partway through
Rscript -e "
  need = c('rmarkdown','tidyverse','SingleCellExperiment','scran','scater',
           'batchelor','BiocSingular','TSCAN','igraph','DT','patchwork','purrr')
  miss = Filter(function(p) !requireNamespace(p, quietly = TRUE), need)
  if (length(miss) > 0) {
    message('missing R packages: ', paste(miss, collapse = ', ')); quit(status = 1)
  }
"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "host:       $(hostname)"
echo "started:    $(date)"
echo "lineage:    ${LINEAGE}"
echo "atlas:      ${RESULTS_DATE}_atlas${ATLAS_SUFFIX}"
echo "prev nodes: ${PREV_NODES:-<none>}"

Rscript -e "
  rmarkdown::render(
    input         = '${PROJECT_DIR}/src/17_germ_pseudotime_tscan.Rmd',
    output_dir    = '${OUT_DIR}',
    output_file   = '17_germ_pseudotime_tscan_${LINEAGE}.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir     = '${PROJECT_DIR}',
                         results_date    = '${RESULTS_DATE}',
                         lineage         = '${LINEAGE}',
                         atlas_suffix    = '${ATLAS_SUFFIX}',
                         prev_nodes_file = '${PREV_NODES}'),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/17_germ_pseudotime_tscan_${LINEAGE}.html"
echo "results:  ${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_tscan_${LINEAGE}${ATLAS_SUFFIX}"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
