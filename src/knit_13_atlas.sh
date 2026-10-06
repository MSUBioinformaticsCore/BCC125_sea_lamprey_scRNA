#!/bin/bash --login
#SBATCH --job-name=atlas
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=run/atlas_%j.out
#SBATCH --account=bioinformaticscore

# Knits 13_atlas_and_stage_correspondence.Rmd: the all-cells clustering, the
# stage correspondence, the split labelling and the marker work that every
# downstream document reads.
#
# It writes cell_cluster_labels.csv, which 17 and 19 take their clusters from,
# so a rerun that changes the labelling changes what those documents see.
#
# Objects are cached in <date>_no_doublets and reused across runs. Raising K
# changes the clustering and the cache signature, so that run recomputes.
#
# Knobs:
#
#   K=20 sbatch src/knit_13_atlas.sh              a different clustering k
#   IEG= sbatch src/knit_13_atlas.sh              skip the dissociation check
#   FORCE=1 sbatch src/knit_13_atlas.sh           overwrite an existing run
#
# 13b is the same document with the mitochondrial, rRNA and ribosomal protein
# genes stripped, as a comparison:
#
#   sbatch src/knit_13b_no_technical_genes.sh
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="${RESULTS_DATE:-20260813}"
K="${K:-30}"
MARKERS="${MARKERS:-canonical_markers_S1_only_2026-08-14.csv}"
IEG="${IEG-immediate_early_markers_2026-10-06.csv}"
OUT_DIR="${PROJECT_DIR}/html"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

# Nothing already on disk is overwritten without being asked for. The results
# directory is the one every downstream document reads, so losing it silently
# would be worse than a failed submission.
HTML_OUT="${OUT_DIR}/13_atlas_and_stage_correspondence.html"
RES_OUT="${PROJECT_DIR}/results/${RESULTS_DATE}_atlas"
if [[ "${FORCE:-0}" != "1" ]]; then
  for p in "${HTML_OUT}" "${RES_OUT}"; do
    if [[ -e "${p}" ]]; then
      echo "refusing to overwrite ${p}" >&2
      echo "FORCE=1 to replace it. Everything downstream reads this directory." >&2
      exit 1
    fi
  done
fi

ALL_SCE="${PROJECT_DIR}/results/${RESULTS_DATE}_no_doublets/all.sce.Rds"
[[ -f "${ALL_SCE}" ]] || {
  echo "missing ${ALL_SCE}" >&2
  echo "Run src/knit_01_no_doublets.sh first; 13 reuses its objects." >&2
  exit 1; }

GENE_DESC="${PROJECT_DIR}/results/${RESULTS_DATE}_no_doublets/gene_descritption.csv"
[[ -f "${GENE_DESC}" ]] || { echo "missing ${GENE_DESC}" >&2; exit 1; }

for f in "${MARKERS}" ${IEG:+"${IEG}"}; do
  [[ -f "${PROJECT_DIR}/data/${f}" ]] || {
    echo "missing ${PROJECT_DIR}/data/${f}" >&2; exit 1; }
done

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

echo "host:       $(hostname)"
echo "started:    $(date)"
echo "k:          ${K}"
echo "markers:    ${MARKERS}"
echo "dissoc set: ${IEG:-<skipped>}"

Rscript -e "
  rmarkdown::render(
    input         = '${PROJECT_DIR}/src/13_atlas_and_stage_correspondence.Rmd',
    output_dir    = '${OUT_DIR}',
    output_file   = '13_atlas_and_stage_correspondence.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir     = '${PROJECT_DIR}',
                         results_date    = '${RESULTS_DATE}',
                         k               = ${K},
                         marker_file     = '${MARKERS}',
                         ieg_marker_file = '${IEG}'),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${HTML_OUT}"
echo "results:  ${RES_OUT}"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
