# Claude Code prompt: standalone doublet QC document

Supersedes `prompt_doublet_exploration.md` for this workflow. That file assumed modifying the
per-sample QC Rmds; this one adds a separate document instead and leaves existing code untouched.

Paste everything below the line.

---

Write a new R Markdown document, `src/04_doublet_qc.Rmd`, that adds doublet scores to the
existing merged objects for this project. Do not modify any existing script or Rmd, and do not
rerun clustering, integration or trajectory analysis.

Read `src/01_all_cells_clustering_annotation.Rmd` first to match its conventions: YAML params
block, `results_dir` and `Robject_dir` construction, chunk options, `code_folding: hide`,
DT::datatable for tables, and the `eval=FALSE` pattern used for expensive chunks with
pre-computed objects loaded for rendering.

## Existing objects

Loaded in `01_all_cells_clustering_annotation.Rmd` from
`params$project_dir/results/20260312/Robjects`:

- `all.sce` is a **named list** of per-sample SingleCellExperiments with raw counts,
  log-normalized values and PCA. Each carries a `sample` column in colData. This is the correct
  input for doublet scoring: unintegrated, per library.
- `merged.sce` is the `fastMNN` output. Carries `colLabels` with 19 Louvain clusters, a `batch`
  column identifying the sample, and UMAP and TSNE reduced dimensions.
- `combined.sce` is the uncorrected counts, built by `do.call(cbind, all.sce2)`.
- `merged.cds` is a monocle3 cell_data_set built from `counts(combined.sce)` with
  `cell_metadata = merged.sce@colData`.

## What the document should do

### 1. Score doublets per sample

Run `scDblFinder::scDblFinder` on each element of `all.sce` separately. Do not pass `samples =`,
since each list element is already one library. Set a seed. Put this in an `eval=FALSE` chunk and
save the result, following the existing pattern.

Collect the results into a single data frame with columns for sample, cell barcode,
`scDblFinder.score` and `scDblFinder.class`.

### 2. Join to merged.sce, carefully

**Verify the join before trusting it.** Cell barcodes may be duplicated across samples, so match
on the sample and barcode pair rather than barcode alone. Check whether `colnames(merged.sce)`
are unique and whether they correspond to `colnames` of the concatenated `all.sce` elements in
the same order.

Report explicitly:

- whether `colnames(merged.sce)` are unique
- the number of cells matched, and the number unmatched in either direction
- whether the row order after joining is identical to the column order of `merged.sce`

If any cell fails to match, stop and report rather than proceeding with a partial or misaligned
join. A silently misaligned doublet score is worse than no doublet score.

Then add `scDblFinder.score` and `scDblFinder.class` to `colData(merged.sce)`, and to
`colData(merged.cds)` so monocle3 plotting can use them.

### 3. Save to a new results directory

Save updated `merged.sce` and `merged.cds` to a **new dated results directory**, following the
existing convention. Do not overwrite the objects in `results/20260312/Robjects`.

### 4. Diagnostics, no filtering

**Do not remove any cells. Do not subset any object.** The purpose is to characterize and report.
Two features of this dataset resemble doublets but are real biology:

1. The germ cell clusters form a differentiation trajectory, so transcriptional intermediates
   genuinely resemble a mixture of the states either side of them.
2. Oocytes have high RNA content, which inflates library size, the other feature doublet callers
   key on. Elevated library size in an oocyte cluster is expected.

Do not use cluster names as evidence about physical cell state. Labels here derive from marker
expression and GO enrichment, so a name referring to cytokinesis or cysts describes a
transcriptional program, not an observed morphology.

Produce:

- Overall percent of cells called doublet, per sample and across the dataset.
- A per-cluster table: cluster, n cells, median doublet score, percent called doublet, median
  UMI count (`sum`), median genes detected (`detected`). Sorted by percent called doublet.
  Render with `datatable` and write to CSV in `results_dir`.
- The discriminating test: identify any cluster combining a high doublet rate **with** elevated
  library size, which indicates true multiplets, versus clusters with elevated scores but
  ordinary UMI counts, which is more consistent with transitional states. State which pattern
  each flagged cluster shows.
- For flagged clusters, the fraction of cells co-expressing germ cell markers (DDX4, PIWIL1,
  SYCP3) and somatic markers (FOXL2, S100A1, SOX18, RBPMS) simultaneously. Genuine doublets
  co-express incompatible lineages; genuine intermediates shift within one. Check that these
  symbols exist in `rownames` before using them, and report any that do not.
- Figures: doublet score on the existing TSNE, and a violin of score per cluster.

### 5. Text for the manuscript

End the document with a short section containing two or three sentences suitable for the Methods
and two or three for the reviewer response, reporting what was done, the numbers found, and that
no cells were removed, with the reasoning.

## Style

- Match the conventions in `01_all_cells_clustering_annotation.Rmd`.
- US English spelling. Never use em dashes; use commas, semicolons, colons or separate sentences.
- Be concise and direct. Flag anything you cannot verify from the objects rather than assuming it.
- If the results contradict the expectation above, say so plainly. A cluster that genuinely looks
  like multiplets is a real finding and I need to know.

## Follow-on, not part of this document

The germ cell subset used in `02_sex_cell_subanalysis.Rmd` and `03_trajectory_analysis.Rmd` will
need the same scores. Note at the end of the document how they should be propagated, but do not
modify those files.
