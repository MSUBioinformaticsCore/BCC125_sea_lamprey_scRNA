#!/bin/bash --login
#SBATCH --job-name=pseudotime
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=run/pseudotime_%j.out
#SBATCH --account=bioinformaticscore

# Knits 17_germ_pseudotime_tscan.Rmd and 18_germ_pseudotime_monocle3.Rmd per
# lineage. Each html goes to html/ and each result to its own directory:
#
#   results/<date>_pseudotime_tscan_female     results/<date>_pseudotime_monocle3_female
#   results/<date>_pseudotime_tscan_male       results/<date>_pseudotime_monocle3_male
#   results/<date>_pseudotime_tscan_all        results/<date>_pseudotime_monocle3_all
#
# The "all" lineage orders every germ cluster in one tree or graph.
#
# Needs results/<date>_atlas/cell_cluster_labels.csv from 13, and all.sce.Rds and
# gene_descritption.csv in results/<date>_no_doublets.
#
# TSCAN runs before monocle3 for each lineage, because the monocle3 document
# compares its ordering with the TSCAN one when that file exists.
#
# Run a subset by setting LINEAGES or METHODS:
#
#   LINEAGES=all METHODS=monocle3 sbatch src/knit_17_18_pseudotime.sh
#   LINEAGES=male sbatch src/knit_17_18_pseudotime.sh
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="${RESULTS_DATE:-20260813}"
LINEAGES="${LINEAGES:-female male all}"
METHODS="${METHODS:-tscan monocle3}"
OUT_DIR="${PROJECT_DIR}/html"
CORES="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

LABELS="${PROJECT_DIR}/results/${RESULTS_DATE}_atlas/cell_cluster_labels.csv"
[[ -f "${LABELS}" ]] || {
  echo "missing ${LABELS}" >&2
  echo "Knit 13_atlas_and_stage_correspondence.Rmd first." >&2
  exit 1; }

# fail before the first knit rather than partway through
NEEDED="c('rmarkdown','TSCAN','monocle3','scatterpie','ggraph','tidygraph','batchelor','BiocSingular')"
Rscript -e "
  miss = Filter(function(p) !requireNamespace(p, quietly = TRUE), ${NEEDED})
  if (length(miss) > 0) { message('missing R packages: ', paste(miss, collapse = ', ')); quit(status = 1) }
"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "host:     $(hostname)"
echo "started:  $(date)"
echo "lineages: ${LINEAGES}"
echo "methods:  ${METHODS}"

for lineage in ${LINEAGES}; do
  for method in ${METHODS}; do
    case "${method}" in
      tscan)
        RMD="${PROJECT_DIR}/src/17_germ_pseudotime_tscan.Rmd"
        EXTRA=""
        ;;
      monocle3)
        RMD="${PROJECT_DIR}/src/18_germ_pseudotime_monocle3.Rmd"
        EXTRA=", cores = ${CORES}"
        ;;
      *)
        echo "unknown method '${method}'; use tscan or monocle3" >&2
        exit 1
        ;;
    esac

    OUT_FILE="$(basename "${RMD}" .Rmd)_${lineage}.html"
    echo
    echo "knitting ${OUT_FILE}: $(date)"

    Rscript -e "
      rmarkdown::render(
        input         = '${RMD}',
        output_dir    = '${OUT_DIR}',
        output_file   = '${OUT_FILE}',
        knit_root_dir = '${PROJECT_DIR}',
        params        = list(project_dir  = '${PROJECT_DIR}',
                             results_date = '${RESULTS_DATE}',
                             lineage      = '${lineage}'${EXTRA}),
        envir         = new.env()
      )
    "
  done
done

echo
echo "finished: $(date)"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
