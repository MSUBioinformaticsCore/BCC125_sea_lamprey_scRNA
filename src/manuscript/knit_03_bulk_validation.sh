#!/bin/bash --login
#SBATCH --job-name=03_bulk_validation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=run/03_bulk_validation_%j.out
#SBATCH --account=bioinformaticscore

# Knits src/manuscript/03_bulk_validation.Rmd, the comparison against the
# Yasmin et al. 2022 bulk series (PRJNA749754). Submit from the project
# directory, after 01_atlas and 02_germ_pseudotime have run:
#
#   sbatch src/manuscript/knit_03_bulk_validation.sh
#
# The html goes to html/ and the results to
# results/20260813_bulk_validation_manuscript.
#
# The slow step is the label permutation: four contrasts times 200 edgeR refits.
# The document runs them on 8 cores (perm_cores in its settings chunk), which is
# why cpus-per-task is 8. Fewer cores than that oversubscribes the node rather
# than failing. The result is cached in
# results/20260813_bulk_validation_manuscript/Robjects/label_permutation_null.Rds
# and reused on a rerun when the inputs it was built on are unchanged. When a
# cached file is passed over, the log names the inputs that differ.
#
# To force the permutation to recompute, delete that Robjects file first.
#
# Time and memory are estimates. Read the seff output and trim them on a rerun.

module purge
module load R/4.4.1-gfbf-2023b

export R_LIBS_SITE="/opt/software-current/2023.06/x86_64/generic/software/R-bundle-CRAN/2024.06-foss-2023b"

PROJECT_DIR="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
RESULTS_DATE="20260813"
NFCORE_DIR="/mnt/research/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA/results/yasmin_nfcore"
GENOME_DIR="/mnt/research/bioinformaticsCore/shared/Genomes/Petromyzon_marinus"
RMD="${PROJECT_DIR}/src/manuscript/03_bulk_validation.Rmd"
OUT_DIR="${PROJECT_DIR}/html"
RESULTS="${PROJECT_DIR}/results"

mkdir -p "${PROJECT_DIR}/run" "${OUT_DIR}"

set -euo pipefail

if [[ "${SLURM_CPUS_PER_TASK:-1}" -lt 8 ]]; then
  echo "note: ${SLURM_CPUS_PER_TASK:-1} core(s) allocated, but the permutation"
  echo "      runs on 8. It will be slower than the time requested assumes."
fi

# Preflight. Every one of these is a hard stop inside the document, so a missing
# file fails here in a second rather than after the pseudobulk and the
# permutation have run.
need=(
  "${RMD}|the document itself"
  "${RESULTS}/${RESULTS_DATE}_no_doublets/all.sce.Rds|the pseudobulk is summed from this; written by 01_atlas"
  "${RESULTS}/${RESULTS_DATE}_no_doublets/gene_descritption.csv|gene names for the tables and labels"
  "${PROJECT_DIR}/data/yasmin_PRJNA749754_sample_groups.csv|the per-animal lengths, life stages and group labels"
  "${RESULTS}/${RESULTS_DATE}_atlas/cluster_marker_genes.csv|atlas markers; written by 01_atlas"
  "${RESULTS}/${RESULTS_DATE}_germ_pseudotime_manuscript/node_marker_genes.csv|fine node markers; written by 02_germ_pseudotime"
  "${RESULTS}/${RESULTS_DATE}_germ_pseudotime_manuscript/pseudotime_genes.csv|path genes; written by 02_germ_pseudotime"
  "${RESULTS}/${RESULTS_DATE}_germ_pseudotime_manuscript/pseudotime_per_cell.csv|per-cell pseudotime; written by 02_germ_pseudotime"
  "${RESULTS}/${RESULTS_DATE}_germ_pseudotime_manuscript/node_graph_layout.csv|tree layout for Figure 7; written by 02_germ_pseudotime"
  "${RESULTS}/${RESULTS_DATE}_germ_pseudotime_manuscript/mst_edges.csv|tree edges for Figure 7; written by 02_germ_pseudotime"
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

# The document accepts the run's top level or its star_salmon subdirectory and
# prefers the length-scaled matrix. Check the same candidates here.
FOUND=""
for rel in star_salmon/salmon.merged.gene_counts_length_scaled.tsv \
           salmon/salmon.merged.gene_counts_length_scaled.tsv \
           star_salmon/salmon.merged.gene_counts.tsv \
           salmon/salmon.merged.gene_counts.tsv \
           salmon.merged.gene_counts_length_scaled.tsv \
           salmon.merged.gene_counts.tsv; do
  if [[ -f "${NFCORE_DIR}/${rel}" ]]; then FOUND="${rel}"; break; fi
done
if [[ -z "${FOUND}" ]]; then
  echo "no salmon gene matrix under ${NFCORE_DIR}" >&2
  missing=1
fi

[[ "${missing}" -eq 0 ]] || exit 1

case "${FOUND}" in
  *length_scaled*) ;;
  *) echo "note: only the unscaled matrix is present (${FOUND})."
     echo "      The transcript-length correction will be absent." ;;
esac

Rscript -e "
  need = c('rmarkdown','tidyverse','SingleCellExperiment','edgeR','patchwork',
           'ggrepel','DT','Matrix','fgsea','BiocParallel')
  miss = Filter(function(p) !requireNamespace(p, quietly = TRUE), need)
  if (length(miss) > 0) {
    message('missing R packages: ', paste(miss, collapse = ', ')); quit(status = 1)
  }
"

# one thread per process, since the parallelism is across the permutations
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PERM_CACHE="${RESULTS}/${RESULTS_DATE}_bulk_validation_manuscript/Robjects/label_permutation_null.Rds"

echo "host:        $(hostname)"
echo "started:     $(date)"
echo "counts:      ${NFCORE_DIR}/${FOUND}"
echo "cores:       ${SLURM_CPUS_PER_TASK:-1}"
if [[ -f "${PERM_CACHE}" ]]; then
  echo "permutation: cached file present, reused if its inputs match"
else
  echo "permutation: no cached file, computing"
fi

Rscript -e "
  rmarkdown::render(
    input         = '${RMD}',
    output_dir    = '${OUT_DIR}',
    output_file   = '03_bulk_validation.html',
    knit_root_dir = '${PROJECT_DIR}',
    params        = list(project_dir = '${PROJECT_DIR}'),
    envir         = new.env()
  )
"

echo
echo "finished: $(date)"
echo "html:     ${OUT_DIR}/03_bulk_validation.html"
echo "results:  ${RESULTS}/${RESULTS_DATE}_bulk_validation_manuscript"
echo
echo "resource use for this job:"
seff "${SLURM_JOB_ID}" || true
