#!/bin/bash --login
#SBATCH --job-name=monocle_all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=run/monocle_%j.out
#SBATCH --account=bioinformaticscore

# Knits 18_germ_pseudotime_monocle3.Rmd for one lineage, by default "all":
# clusters 5 and 15 to 18, with 12_transformer excluded.
#
# By default no cell is dropped, so 5.1 is kept. monocle3 learns a principal
# graph over the cells rather than a tree over cluster centroids, so it does not
# need the subclusters and is a fair place to look at 5.1 before deciding.
#
# To drop the cells annotated as somatic, matching what 17 currently does:
#
#   DROP="Migrating neural crest cells" sbatch src/knit_18_monocle.sh
#   DROP="Migrating neural crest cells;Granulosa" sbatch src/knit_18_monocle.sh
#
# DROP needs PREV_NODES, which names the run whose subclusters were annotated.
# 17 and 18 only order the same cells when their DROP settings agree, so check
# drop_annotations in 17 before comparing the two.
#
# If 17 has already run for this lineage, 18 compares its ordering with the
# TSCAN one.
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="${RESULTS_DATE:-20260813}"
LINEAGE="${LINEAGE:-all}"
ATLAS_SUFFIX="${ATLAS_SUFFIX:-}"
DROP="${DROP:-}"
CORES="${SLURM_CPUS_PER_TASK:-8}"
OUT_DIR="${PROJECT_DIR}/html"

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
[[ -f "${SCE}" ]] || { echo "missing ${SCE}; check RESULTS_DATE." >&2; exit 1; }

if [[ -n "${DROP}" ]]; then
  if [[ ! -f "${PREV_NODES}" ]]; then
    echo "DROP is set but ${PREV_NODES} was not found." >&2
    echo "The cells to drop are identified by barcode from that file." >&2
    exit 1
  fi
  echo "dropping: ${DROP}"
else
  echo "note: no cell is dropped, so 5.1 is kept."
fi

Rscript -e "
  need = c('rmarkdown','tidyverse','SingleCellExperiment','scran','scater',
           'batchelor','BiocSingular','monocle3','DT','patchwork','ggraph','tidygraph')
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
    input         = '${PROJECT_DIR}/src/18_germ_pseudotime_monocle3.Rmd',
    output_dir    = '${OUT_DIR}',
    output_file   = '18_germ_pseudotime_monocle3_${LINEAGE}.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir      = '${PROJECT_DIR}',
                         results_date     = '${RESULTS_DATE}',
                         lineage          = '${LINEAGE}',
                         atlas_suffix     = '${ATLAS_SUFFIX}',
                         prev_nodes_file  = '${PREV_NODES}',
                         drop_annotations = '${DROP}',
                         cores            = ${CORES}),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/18_germ_pseudotime_monocle3_${LINEAGE}.html"
echo "results:  ${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_monocle3_${LINEAGE}${ATLAS_SUFFIX}"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
