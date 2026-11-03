#!/bin/bash --login
#SBATCH --job-name=01_atlas_soupx
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=run/01_atlas_soupx_%j.out
#SBATCH --account=bioinformaticscore

# Knits src/manuscript/01_atlas_soupx.Rmd, the atlas rebuilt on SoupX-corrected
# counts. Submit from the project directory once all four
# knit_00_sample_qc_soupx.sh tasks have finished:
#
#   sbatch src/manuscript/knit_01_atlas_soupx.sh
#
# Or chain it on the array so it starts only if every task succeeded:
#
#   jid=$(sbatch --parsable src/manuscript/knit_00_sample_qc_soupx.sh)
#   sbatch --dependency=afterok:${jid} src/manuscript/knit_01_atlas_soupx.sh
#
# The html goes to html/01_atlas_soupx.html and the results to
# results/20260813_atlas_soupx. Its cached objects go to
# results/20260813_no_doublets_soupx, apart from the uncorrected ones, so the
# uncorrected atlas and everything that reads it are untouched.
#
# Resources follow the uncorrected atlas's knit script. Read the seff output and
# trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="20260813"
GENOME_DIR="/mnt/research/bioinformaticsCore/shared/Genomes/Petromyzon_marinus"
RMD="${PROJECT_DIR}/src/manuscript/01_atlas_soupx.Rmd"
OUT_DIR="${PROJECT_DIR}/html"
RESULTS="${PROJECT_DIR}/results"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

# Preflight. Each of these is a hard stop inside the document.
need=(
  "${RMD}|the document itself"
  "${PROJECT_DIR}/src/manuscript/scRNA_seq_functions_manuscript.R|the shared functions"
  "${RESULTS}/soupx/UND1_C00004336/UND1_C00004336_soupx.Rds|corrected counts; knit_00_sample_qc_soupx.sh task 0"
  "${RESULTS}/soupx/UND2_C00004362/UND2_C00004362_soupx.Rds|corrected counts; knit_00_sample_qc_soupx.sh task 1"
  "${RESULTS}/soupx/TF_C00004337/TF_C00004337_soupx.Rds|corrected counts; knit_00_sample_qc_soupx.sh task 2"
  "${RESULTS}/soupx/TM_C00004374/TM_C00004374_soupx.Rds|corrected counts; knit_00_sample_qc_soupx.sh task 3"
  "${RESULTS}/${RESULTS_DATE}_atlas/cell_cluster_labels.csv|the uncorrected labels the cell type names are carried from; written by 01_atlas"
  "${RESULTS}/${RESULTS_DATE}_no_doublets/gene_descritption.csv|gene names for the tables and labels"
  "${RESULTS}/${RESULTS_DATE}_stable_k/stability_by_k.csv|the k stability sweep, from src/slurm_stable_k.sh"
  "${PROJECT_DIR}/data/canonical_markers_S1_only_2026-08-14.csv|the curated markers"
  "${PROJECT_DIR}/data/immediate_early_markers_2026-10-06.csv|the immediate-early genes"
  "${GENOME_DIR}/GCF_010993605.1_kPetMar1.pri_gene_ontology.gaf|GO annotation"
  "${GENOME_DIR}/GOs_names.txt|GO term names"
)
missing=0
for entry in "${need[@]}"; do
  f="${entry%%|*}"
  if [[ ! -f "${f}" ]]; then
    echo "missing ${f}" >&2
    echo "  (${entry#*|})" >&2
    missing=1
  fi
done
[[ "${missing}" -eq 0 ]] || exit 1

Rscript -e "
  need = c('rmarkdown','tidyverse','SingleCellExperiment','scater','scran',
           'patchwork','RColorBrewer','DT','batchelor','BiocSingular',
           'scDblFinder','igraph')
  miss = Filter(function(p) !requireNamespace(p, quietly = TRUE), need)
  if (length(miss) > 0) {
    message('missing R packages: ', paste(miss, collapse = ', ')); quit(status = 1)
  }
"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "host:     $(hostname)"
echo "started:  $(date)"

Rscript -e "
  rmarkdown::render(
    input         = '${RMD}',
    output_dir    = '${OUT_DIR}',
    output_file   = '01_atlas_soupx.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir = '${PROJECT_DIR}'),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/01_atlas_soupx.html"
echo "results:  ${RESULTS}/${RESULTS_DATE}_atlas_soupx"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
