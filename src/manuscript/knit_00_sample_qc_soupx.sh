#!/bin/bash --login
#SBATCH --job-name=00_soupx_qc
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=run/00_soupx_qc_%A_%a.out
#SBATCH --account=bioinformaticscore

# Knits src/manuscript/00_sample_qc_soupx.Rmd once per library, as a SLURM
# array: task 0 is UND1, 1 is UND2, 2 is TF and 3 is TM. Submit from the
# project directory:
#
#   sbatch src/manuscript/knit_00_sample_qc_soupx.sh
#
# One library only, for example TM after choosing its estimate:
#
#   RHO_SOURCE=spermatocyte sbatch --array=3 src/manuscript/knit_00_sample_qc_soupx.sh
#
# Each task writes results/soupx/<sample>/, including <sample>_soupx.Rds, which
# 01_atlas_soupx.Rmd reads, and html/00_sample_qc_soupx_<sample>.html.
#
# RHO_SOURCE picks which contamination estimate is applied: auto (the default,
# SoupX's autoEstCont) or spermatocyte. Both are computed and reported either
# way. If the chosen one cannot be made for a library, that task stops; resubmit
# that task with the other RHO_SOURCE.
#
# Inputs, per library:
#   - data/<sample>/raw_matrix/, the unfiltered matrix with the empty droplets
#   - data/<sample>/<sample>_cds_mito_filter.Rds, the cells the per-sample QC kept
#   - data/gene_descritption.csv and data/canonical_markers_S1_only_2026-08-14.csv
#
# Time and memory are estimates. The raw matrix holds every barcode, so memory
# is the one to watch. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RMD="${PROJECT_DIR}/src/manuscript/00_sample_qc_soupx.Rmd"
OUT_DIR="${PROJECT_DIR}/html"
RHO_SOURCE="${RHO_SOURCE:-auto}"

SAMPLES=(UND1_C00004336 UND2_C00004362 TF_C00004337 TM_C00004374)
SAMPLE="${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

case "${RHO_SOURCE}" in
  auto|spermatocyte) ;;
  *) echo "RHO_SOURCE must be auto or spermatocyte, not ${RHO_SOURCE}" >&2
     exit 1 ;;
esac

# Preflight, so a missing input fails in a second rather than after the raw
# matrix has been read.
need=(
  "${RMD}"
  "${PROJECT_DIR}/data/${SAMPLE}/raw_matrix/matrix.mtx.gz"
  "${PROJECT_DIR}/data/${SAMPLE}/raw_matrix/barcodes.tsv.gz"
  "${PROJECT_DIR}/data/${SAMPLE}/raw_matrix/features.tsv.gz"
  "${PROJECT_DIR}/data/${SAMPLE}/${SAMPLE}_cds_mito_filter.Rds"
  "${PROJECT_DIR}/data/gene_descritption.csv"
  "${PROJECT_DIR}/data/canonical_markers_S1_only_2026-08-14.csv"
)
missing=0
for f in "${need[@]}"; do
  [[ -f "${f}" ]] || { echo "missing ${f}" >&2; missing=1; }
done
[[ "${missing}" -eq 0 ]] || exit 1

Rscript -e "
  need = c('rmarkdown','tidyverse','Matrix','SingleCellExperiment','monocle3',
           'scater','scran','SoupX','DT','igraph')
  miss = Filter(function(p) !requireNamespace(p, quietly = TRUE), need)
  if (length(miss) > 0) {
    message('missing R packages: ', paste(miss, collapse = ', '))
    if ('SoupX' %in% miss)
      message('SoupX is on CRAN: install.packages(\"SoupX\") into your R library.')
    quit(status = 1)
  }
"

echo "host:       $(hostname)"
echo "started:    $(date)"
echo "sample:     ${SAMPLE}"
echo "rho source: ${RHO_SOURCE}"

# the intermediates go to a per-task directory, so the four tasks knitting the
# same Rmd at once do not overwrite each other's files
WORK_DIR="${PROJECT_DIR}/run/knit_00_soupx_${SAMPLE}"
mkdir -p "${WORK_DIR}"

Rscript -e "
  rmarkdown::render(
    input             = '${RMD}',
    output_dir        = '${OUT_DIR}',
    output_file       = '00_sample_qc_soupx_${SAMPLE}.html',
    intermediates_dir = '${WORK_DIR}',
    knit_root_dir     = '${PROJECT_DIR}',
    params            = list(project_dir = '${PROJECT_DIR}',
                             sampleID    = '${SAMPLE}',
                             rho_source  = '${RHO_SOURCE}'),
    envir             = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/00_sample_qc_soupx_${SAMPLE}.html"
echo "results:  ${PROJECT_DIR}/results/soupx/${SAMPLE}"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
