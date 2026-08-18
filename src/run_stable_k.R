#!/usr/bin/env Rscript

# Clustering stability sweep over k for the doublet-filtered atlas.
#
# Answers whether the k = 30 used in
# 01_all_cells_clustering_annotation_no_doublets.Rmd sits on a plateau or on a
# slope. find_stable_k clusters an independent subsample of cells at each
# iteration and compares pairs of iterations on the cells they share, so the
# ARI it reports is stability against perturbation of the data rather than
# repetition of a deterministic call.
#
# Configured by environment variables, set in slurm_stable_k.sh.
#
# Missing packages are installed into a project library. Compute nodes often
# have no outbound network, in which case the install fails and the job stops
# with instructions. To do the install once, from a login node:
#
#   INSTALL_ONLY=1 Rscript src/run_stable_k.R
#
# after which the batch job finds everything already present and installs
# nothing.

env <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

project_dir  <- env("PROJECT_DIR",
                    "/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA")
no_dbl_date  <- env("NO_DBL_DATE", "20260813")
k_range      <- as.integer(strsplit(env("K_RANGE", "10,20,30,40,50"), ",")[[1]])
n_iterations <- as.integer(env("N_ITERATIONS", "10"))
prop         <- as.numeric(env("PROP", "0.8"))
seed         <- as.integer(env("SEED", "42"))
install_only <- nzchar(Sys.getenv("INSTALL_ONLY", unset = ""))

# ---- packages ------------------------------------------------------------

# A writable library, since the module's is read only. Candidates are tried in
# order and the first one that can actually be written to wins, rather than
# assuming any particular location is writable: research space is often mounted
# read only or quota limited, and a home directory is not.
#
# Set R_LIBS_PROJECT to force a specific path.
lib_candidates <- c(
  Sys.getenv("R_LIBS_PROJECT", unset = ""),
  Sys.getenv("R_LIBS_USER",    unset = ""),
  file.path(path.expand("~"), "R", paste0(R.version$platform, "-library"),
            paste(R.version$major, sub("[.].*", "", R.version$minor), sep = ".")))
lib_candidates <- unique(lib_candidates[nzchar(lib_candidates)])

is_writable <- function(d) {
  if (!dir.exists(d))
    suppressWarnings(dir.create(d, recursive = TRUE, showWarnings = FALSE))
  dir.exists(d) && file.access(d, mode = 2) == 0
}

lib_path <- NULL
for (cand in lib_candidates) {
  if (is_writable(cand)) { lib_path <- cand; break }
  message("not writable, skipping: ", cand)
}

if (is.null(lib_path))
  stop("no writable library path. Tried:\n  ",
       paste(lib_candidates, collapse = "\n  "),
       "\nSet R_LIBS_PROJECT to somewhere you can write and resubmit.")

message("library path: ", lib_path)
.libPaths(c(lib_path, .libPaths()))

cran_repo <- env("CRAN_REPO", "https://cloud.r-project.org")
n_cpus    <- max(1L, as.integer(env("SLURM_CPUS_PER_TASK", "1")))

cran_pkgs <- c("mclust", "igraph", "ggplot2")
bioc_pkgs <- c("SingleCellExperiment", "scran")

have    <- function(p) requireNamespace(p, quietly = TRUE)
missing <- function(pkgs) pkgs[!vapply(pkgs, have, logical(1))]

install_cran <- function(pkgs) {
  if (!length(pkgs)) return(invisible(NULL))
  message("installing from CRAN into ", lib_path, ": ", paste(pkgs, collapse = ", "))
  try(install.packages(pkgs, lib = lib_path, repos = cran_repo, Ncpus = n_cpus),
      silent = FALSE)
}

install_bioc <- function(pkgs) {
  if (!length(pkgs)) return(invisible(NULL))
  install_cran(missing("BiocManager"))
  if (!have("BiocManager")) return(invisible(NULL))
  message("installing from Bioconductor into ", lib_path, ": ",
          paste(pkgs, collapse = ", "))
  try(BiocManager::install(pkgs, lib = lib_path, ask = FALSE, update = FALSE,
                           Ncpus = n_cpus), silent = FALSE)
}

install_cran(missing(cran_pkgs))
install_bioc(missing(bioc_pkgs))

still <- missing(c(cran_pkgs, bioc_pkgs))
if (length(still)) {
  stop("could not load or install: ", paste(still, collapse = ", "), "\n",
       "If this ran on a compute node, it most likely has no outbound network. ",
       "Install once from a login node with:\n",
       "  INSTALL_ONLY=1 R_LIBS_PROJECT=", lib_path,
       " Rscript ", file.path(project_dir, "src", "run_stable_k.R"), "\n",
       "then resubmit the batch job.")
}

if (install_only) {
  message("all packages present in ", lib_path, "; INSTALL_ONLY set, stopping here.")
  print(sessionInfo())
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
  library(igraph)
  library(ggplot2)
})

source(file.path(project_dir, "src", "scRNA_seq_functions.R"))

# ---- fail early, before spending the allocation ---------------------------

no_dbl_dir  <- file.path(project_dir, "results", paste0(no_dbl_date, "_no_doublets"))
merged_file <- file.path(no_dbl_dir, "merged.sce.Rds")
if (!file.exists(merged_file))
  stop("missing ", merged_file, ". Run 01_all_cells_clustering_annotation_no_doublets.Rmd first.")

results_dir <- file.path(project_dir, "results", paste0(no_dbl_date, "_stable_k"))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

message("reading ", merged_file)
merged.sce <- readRDS(merged_file)

# The clustering being interrogated was built on the fastMNN corrected space,
# so the sweep must use that space. find_stable_k defaults to "PCA", which this
# object does not carry at all.
use_dimred <- "corrected"
if (!use_dimred %in% reducedDimNames(merged.sce))
  stop("reducedDim '", use_dimred, "' not found. Present: ",
       paste(reducedDimNames(merged.sce), collapse = ", "))

# ---- strip to the coordinates -------------------------------------------

# fastMNN output carries a reconstructed assay. Nothing here reads it, and
# find_stable_k subsets the object once per iteration, so it is dropped rather
# than carried through every subset.
lean <- SingleCellExperiment(assays = list(), colData = colData(merged.sce))
reducedDim(lean, use_dimred) <- reducedDim(merged.sce, use_dimred)
colnames(lean) <- colnames(merged.sce)
rm(merged.sce); invisible(gc())

message(sprintf("%d cells, %d dimensions in '%s'",
                ncol(lean), ncol(reducedDim(lean, use_dimred)), use_dimred))
message(sprintf("k: %s | iterations: %d | prop: %.2f",
                paste(k_range, collapse = ", "), n_iterations, prop))

# ---- run -----------------------------------------------------------------

t0 <- Sys.time()
stab <- find_stable_k(lean,
                      k_range      = k_range,
                      n_iterations = n_iterations,
                      use.dimred   = use_dimred,
                      prop         = prop,
                      seed         = seed)
message("elapsed: ", format(round(difftime(Sys.time(), t0), 1)))

saveRDS(stab, file = file.path(results_dir, "stable_k.Rds"))

# ---- summarise -----------------------------------------------------------

viz <- plot_stability_metrics(stab)
df  <- viz$data
df$n_iterations <- n_iterations
df$prop         <- prop
df$mean_shared  <- vapply(stab, function(x) x$mean_shared, numeric(1))

write.csv(df, file = file.path(results_dir, "stability_by_k.csv"), row.names = FALSE)

ggsave(file.path(results_dir, "stability_n_clusters.png"), viz$p1, width = 6, height = 4)
ggsave(file.path(results_dir, "stability_ari.png"),        viz$p2, width = 6, height = 4)
ggsave(file.path(results_dir, "stability_cv.png"),         viz$p3, width = 6, height = 4)

cat("\n=== stability by k ===\n")
print(df[order(df$k), c("k", "mean_n_clusters", "sd_n_clusters",
                        "cv_n_clusters", "mean_ari", "min_ari", "mean_shared")],
      row.names = FALSE)

best <- df$k[which.max(df$mean_ari)]
cat("\nhighest mean ARI at k =", best, "\n")
if (30 %in% df$k)
  cat("k = 30 was used for the published clustering; its mean ARI is ",
      signif(df$mean_ari[df$k == 30], 3), "\n", sep = "")

cat("\nRead the ARI curve for a plateau rather than a peak. A k sitting on a\n",
    "flat region means the choice is not load bearing. A sharp drop either\n",
    "side means the cluster count is a knife edge and should be reported.\n", sep = "")

cat("\nwritten to ", results_dir, "\n", sep = "")
sessionInfo()
