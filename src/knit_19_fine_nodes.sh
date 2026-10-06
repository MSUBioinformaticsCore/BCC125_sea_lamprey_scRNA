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
DIST="${DIST:-simple}"
MNN_K="${MNN_K:-50}"
OUTGROUP="${OUTGROUP:-FALSE}"
EDGE_BOOT="${EDGE_BOOT:-100}"
EDGE_BOOT_CELLS="${EDGE_BOOT_CELLS:-100}"
GRAPH_SETS="${GRAPH_SETS:-Primordial germ cells,Migrating germ cells,Spermatocytes,Oocyte}"
# nodes to remove, by node id, after the evidence section reports what each
# node expresses. Empty removes nothing.
DROP="${DROP:-}"
EXTRA="${EXTRA:-neural_crest_markers_2026-09-25.csv,immediate_early_markers_2026-10-06.csv}"
EVIDENCE="${EVIDENCE:-Dissociation response}"
# A run that changes the tree keeps its own html and results directory, so two
# settings can sit side by side, so an earlier run is never overwritten.
# Built from DROP and OUTGROUP unless RUN_SUFFIX is set directly.
if [[ -z "${RUN_SUFFIX+x}" ]]; then
  RUN_SUFFIX=""
  [[ -n "${DROP}" ]]              && RUN_SUFFIX="${RUN_SUFFIX}_dropped"
  [[ "${OUTGROUP}" == "TRUE" ]]   && RUN_SUFFIX="${RUN_SUFFIX}_outgroup"
fi

OUT_DIR="${PROJECT_DIR}/html"

if [[ -z "${PREV_NODES+x}" ]]; then
  PREV_NODES="${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_tscan_all/pseudotime_per_cell_annotated.csv"
fi

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

# Nothing already on disk is overwritten without being asked for. Set
# RUN_SUFFIX to keep both runs, or FORCE=1 to replace the earlier one.
HTML_OUT="${OUT_DIR}/19_germ_pseudotime_fine_nodes_${LINEAGE}${RUN_SUFFIX}.html"
RES_OUT="${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_fine_${LINEAGE}${ATLAS_SUFFIX}${RUN_SUFFIX}"
if [[ "${FORCE:-0}" != "1" ]]; then
  for p in "${HTML_OUT}" "${RES_OUT}"; do
    if [[ -e "${p}" ]]; then
      echo "refusing to overwrite ${p}" >&2
      echo "RUN_SUFFIX=_something to keep both, or FORCE=1 to replace it." >&2
      exit 1
    fi
  done
fi

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
echo "dropping:   ${DROP:-<nothing; read the Removal section, then set DROP>}"

Rscript -e "
  rmarkdown::render(
    input         = '${PROJECT_DIR}/src/19_germ_pseudotime_fine_nodes.Rmd',
    output_dir    = '${OUT_DIR}',
    output_file   = '19_germ_pseudotime_fine_nodes_${LINEAGE}${RUN_SUFFIX}.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir     = '${PROJECT_DIR}',
                         results_date    = '${RESULTS_DATE}',
                         atlas_suffix    = '${ATLAS_SUFFIX}',
                         run_suffix      = '${RUN_SUFFIX}',
                         lineage         = '${LINEAGE}',
                         sub_k           = ${SUB_K},
                         min_node_cells  = ${MIN_NODE},
                         dist_method     = '${DIST}',
                         mnn_k           = ${MNN_K},
                         outgroup        = ${OUTGROUP},
                         edge_boot       = ${EDGE_BOOT},
                         edge_boot_cells = ${EDGE_BOOT_CELLS},
                         graph_marker_sets = '${GRAPH_SETS}',
                         drop_nodes        = '${DROP}',
                         extra_marker_file = '${EXTRA}',
                         removal_evidence_sets = '${EVIDENCE}',
                         prev_nodes_file = '${PREV_NODES}'),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/19_germ_pseudotime_fine_nodes_${LINEAGE}${RUN_SUFFIX}.html"
echo "results:  ${PROJECT_DIR}/results/${RESULTS_DATE}_pseudotime_fine_${LINEAGE}${ATLAS_SUFFIX}${RUN_SUFFIX}"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
