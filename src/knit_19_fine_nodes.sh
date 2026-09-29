#!/bin/bash --login
#SBATCH --job-name=fine_nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=run/fine_nodes_%j.out
#SBATCH --account=bioinformaticscore

# Knits 19_germ_pseudotime_fine_nodes.Rmd, by default for the "all" lineage:
# clusters 5 and 15 to 18, 12_transformer excluded, nothing removed.
#
# Every cluster is subclustered at sub_k = 10 with a floor of 30 cells, so the
# MST is fitted to many small waypoints rather than a dozen large ones. Node
# distances use dist.method = "mnn", which measures how close two nodes are at
# their boundaries rather than at their centers.
#
# The collaborators' annotation is carried per cell, by barcode, from
# PREV_NODES. It is never attached to a node, so nothing here depends on the
# subclustering being stable or on any node carrying a label.
#
# Knobs worth turning:
#
#   SUB_K=20 sbatch src/knit_19_fine_nodes.sh          fewer, larger nodes
#   MIN_NODE=50 sbatch src/knit_19_fine_nodes.sh       raise the size floor
#   DIST=simple sbatch src/knit_19_fine_nodes.sh       Euclidean centroids, much faster
#   OUTGROUP=TRUE sbatch src/knit_19_fine_nodes.sh     let unrelated nodes detach
#
# OUTGROUP=TRUE is the one to try for the 5.1 question. It reroutes any edge
# longer than 1.5x the median through an artificial outgroup, splitting the tree
# into components. A node that belongs nowhere then detaches itself instead of
# being removed by hand.
#
# MNN runs a cell-based neighbor search for every pair of nodes, so its cost
# grows with the node count. If the job runs long or runs out of memory, raise
# SUB_K to cut the number of nodes, or fall back to DIST=simple.
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="${RESULTS_DATE:-20260813}"
LINEAGE="${LINEAGE:-all}"
ATLAS_SUFFIX="${ATLAS_SUFFIX:-}"
SUB_K="${SUB_K:-10}"
MIN_NODE="${MIN_NODE:-30}"
DIST="${DIST:-mnn}"
MNN_K="${MNN_K:-50}"
OUTGROUP="${OUTGROUP:-FALSE}"
OUT_DIR="${PROJECT_DIR}/html"

if [[ -z "${PREV_NODES+x}" ]]; then
  PREV_NODES="${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_tscan_all/pseudotime_per_cell_annotated.csv"
fi

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

case "${DIST}" in
  simple|scaled.full|scaled.diag|slingshot|mnn) ;;
  *) echo "DIST must be one of simple, scaled.full, scaled.diag, slingshot, mnn" >&2; exit 1 ;;
esac
case "${OUTGROUP}" in
  TRUE|FALSE) ;;
  *) echo "OUTGROUP must be TRUE or FALSE" >&2; exit 1 ;;
esac

LABELS="${PROJECT_DIR}/results/${RESULTS_DATE}_atlas${ATLAS_SUFFIX}/cell_cluster_labels.csv"
[[ -f "${LABELS}" ]] || {
  echo "missing ${LABELS}" >&2
  echo "Knit 13_atlas_and_stage_correspondence.Rmd first." >&2
  exit 1; }

SCE="${PROJECT_DIR}/results/${RESULTS_DATE}_no_doublets/all.sce.Rds"
[[ -f "${SCE}" ]] || { echo "missing ${SCE}; check RESULTS_DATE." >&2; exit 1; }

if [[ ! -f "${PREV_NODES}" ]]; then
  echo "note: ${PREV_NODES} not found."
  echo "      The trajectory will be built, but with no annotation to read it against."
fi

Rscript -e "
  need = c('rmarkdown','tidyverse','SingleCellExperiment','scran','scater',
           'batchelor','BiocSingular','BiocNeighbors','TSCAN','igraph','DT','patchwork',
           'ggraph','tidygraph','ggrepel')
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
echo "sub_k:      ${SUB_K}   min node: ${MIN_NODE}"
echo "distance:   ${DIST}    mnn.k: ${MNN_K}    outgroup: ${OUTGROUP}"
echo "prev nodes: ${PREV_NODES:-<none>}"

Rscript -e "
  rmarkdown::render(
    input         = '${PROJECT_DIR}/src/19_germ_pseudotime_fine_nodes.Rmd',
    output_dir    = '${OUT_DIR}',
    output_file   = '19_germ_pseudotime_fine_nodes_${LINEAGE}.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir     = '${PROJECT_DIR}',
                         results_date    = '${RESULTS_DATE}',
                         atlas_suffix    = '${ATLAS_SUFFIX}',
                         lineage         = '${LINEAGE}',
                         sub_k           = ${SUB_K},
                         min_node_cells  = ${MIN_NODE},
                         dist_method     = '${DIST}',
                         mnn_k           = ${MNN_K},
                         outgroup        = ${OUTGROUP},
                         prev_nodes_file = '${PREV_NODES}'),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/19_germ_pseudotime_fine_nodes_${LINEAGE}.html"
echo "results:  ${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_fine_${LINEAGE}${ATLAS_SUFFIX}"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
