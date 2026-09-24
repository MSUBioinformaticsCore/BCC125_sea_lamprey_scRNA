#!/bin/bash --login
#SBATCH --job-name=bulk_validation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=run/bulk_validation_%j.out
#SBATCH --account=bioinformaticscore

# Knits 16_bulk_validation.Rmd, the comparison against the Yasmin et al. 2022
# bulk series (PRJNA749754). The html goes to html/ and the results to
# results/<date>_bulk_validation.
#
# Inputs:
#   - the nf-core/rnaseq gene matrix for the Yasmin reads, under NFCORE_DIR
#   - results/<date>_no_doublets/all.sce.Rds, for the pseudobulk
#   - data/yasmin_PRJNA749754_sample_groups.csv, the per-animal records
#   - data/<marker_file>, optional; the composition section skips without it
#
# The slow step is the label permutation: six contrasts times N_PERM edgeR
# refits, about 1,200 fits at the default 200. Those fits run in parallel across
# the allocated cores, so the wall time scales with cpus-per-task. The draws are
# generated in one stream before the fits fan out, so the result is identical
# whatever the core count. Everything else in the document is serial. It is
# cached in
# results/<date>_bulk_validation/Robjects/label_permutation_null.Rds and reused
# on a rerun unless the gene set, the group sizes, the trend scores, the count
# filter or N_PERM change.
#
# For a quick trial, cut the permutations and the bootstrap:
#
#   N_PERM=20 N_BOOT=200 sbatch src/knit_16_bulk_validation.sh
#
# To force the permutation to recompute, delete that Robjects file first.
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="${RESULTS_DATE:-20260813}"
NFCORE_DIR="${NFCORE_DIR:-/mnt/research/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA/results/yasmin_nfcore}"
MARKER_FILE="${MARKER_FILE:-canonical_markers_S1_only_2026-08-14.csv}"
# Which group assignment to use. group_suppT1labels follows Supplementary
# Table 1 and the main text's description of mid males; group_textcounts
# reproduces the group sizes printed in the text. They differ for two males.
GROUP_COL="${GROUP_COL:-group_suppT1labels}"
N_PERM="${N_PERM:-200}"
N_BOOT="${N_BOOT:-2000}"
N_RAND="${N_RAND:-2000}"
MIN_PB="${MIN_PB:-10}"
CORES="${SLURM_CPUS_PER_TASK:-1}"
OUT_DIR="${PROJECT_DIR}/html"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

# Preflight. Each of these would otherwise fail partway through the knit, after
# the pseudobulk has already been built.

SCE="${PROJECT_DIR}/results/${RESULTS_DATE}_no_doublets/all.sce.Rds"
[[ -f "${SCE}" ]] || {
  echo "missing ${SCE}" >&2
  echo "The pseudobulk is summed from this object. Check RESULTS_DATE." >&2
  exit 1; }

SHEET="${PROJECT_DIR}/data/yasmin_PRJNA749754_sample_groups.csv"
[[ -f "${SHEET}" ]] || {
  echo "missing ${SHEET}" >&2
  echo "This holds the per-animal lengths, life stages and group labels." >&2
  exit 1; }

# The document accepts the run's top level or its star_salmon subdirectory, and
# prefers the length-scaled matrix. Check the same candidates here so a wrong
# NFCORE_DIR fails in a second rather than after the knit starts.
FOUND=""
for rel in star_salmon/salmon.merged.gene_counts_length_scaled.tsv \
           salmon/salmon.merged.gene_counts_length_scaled.tsv \
           star_salmon/salmon.merged.gene_counts.tsv \
           salmon/salmon.merged.gene_counts.tsv \
           salmon.merged.gene_counts_length_scaled.tsv \
           salmon.merged.gene_counts.tsv; do
  if [[ -f "${NFCORE_DIR}/${rel}" ]]; then FOUND="${rel}"; break; fi
done

[[ -n "${FOUND}" ]] || {
  echo "no salmon gene matrix under ${NFCORE_DIR}" >&2
  echo "Point NFCORE_DIR at the nf-core/rnaseq output for the Yasmin data," >&2
  echo "either the run's top level or its star_salmon subdirectory." >&2
  echo "Found there instead:" >&2
  ls -1 "${NFCORE_DIR}" 2>/dev/null | head -20 >&2 || true
  exit 1; }

case "${FOUND}" in
  *length_scaled*) ;;
  *) echo "note: only the unscaled matrix is present (${FOUND})."
     echo "      The transcript-length correction will be absent." ;;
esac

MARKERS="${PROJECT_DIR}/data/${MARKER_FILE}"
[[ -f "${MARKERS}" ]] || \
  echo "note: ${MARKERS} not found; the composition section will be skipped."

# fail before the knit rather than partway through
Rscript -e "
  need = c('rmarkdown','tidyverse','SingleCellExperiment','edgeR','patchwork',
           'ggrepel','DT','Matrix')
  miss = Filter(function(p) !requireNamespace(p, quietly = TRUE), need)
  if (length(miss) > 0) {
    message('missing R packages: ', paste(miss, collapse = ', ')); quit(status = 1)
  }
"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "host:        $(hostname)"
echo "started:     $(date)"
echo "counts:      ${NFCORE_DIR}/${FOUND}"
echo "permutations: ${N_PERM} on ${CORES} core(s)"

Rscript -e "
  rmarkdown::render(
    input         = '${PROJECT_DIR}/src/16_bulk_validation.Rmd',
    output_dir    = '${OUT_DIR}',
    output_file   = '16_bulk_validation.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir          = '${PROJECT_DIR}',
                         results_date         = '${RESULTS_DATE}',
                         nfcore_dir           = '${NFCORE_DIR}',
                         marker_file          = '${MARKER_FILE}',
                         group_column         = '${GROUP_COL}',
                         n_perm_labels        = ${N_PERM},
                         n_boot               = ${N_BOOT},
                         n_rand_sets          = ${N_RAND},
                         min_pseudobulk_count = ${MIN_PB},
                         cores                = ${CORES}),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/16_bulk_validation.html"
echo "results:  ${PROJECT_DIR}/results/${RESULTS_DATE}_bulk_validation"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
