# Claude Code prompt: test the female germ cells for continuous structure

Context: after doublet removal, the female germ cells do not resolve into discrete subclusters at
`k = 8`, so the TSCAN trajectory used previously cannot be reconstructed. The question is whether
they nonetheless form a continuum, since Louvain clustering finds communities and a smooth
gradient has none. Failing to cluster is weak evidence against a trajectory.

Paste everything below the line.

---

Write a new R Markdown document, `src/06_female_germ_cell_continuum.Rmd`, testing whether the
doublet-filtered female germ cells form a continuous transcriptional gradient. Do not modify any
existing script or Rmd.

Read `src/02_sex_cell_subanalysis.Rmd` and `src/03_trajectory_analysis.Rmd` first, to find object
names and to match conventions: the YAML params block, `results_dir` and `Robject_dir`
construction, chunk options, `code_folding: hide`, `DT::datatable` for tables, and the pattern of
putting expensive steps in `eval=FALSE` chunks with pre-computed objects loaded for rendering.

## Background

Roughly 1,700 female germ cells survive doublet filtering. Subclustering with
`buildSNNGraph(k = 8)` followed by Louvain does not split them. This is not a cell number problem
and not an over-smoothing problem; `k = 8` is already a small neighborhood.

**Do not attempt to force discrete clusters by lowering `k` or raising resolution.** If a
continuum is present, clustering is the wrong tool, and tuning parameters until clusters appear
would be selecting a parameter to obtain a desired result.

## Analyses

### 1. Confirm clustering genuinely fails

Run `buildSNNGraph` plus Louvain across a range of `k`, for example 5, 8, 10, 15, 20, 30, and
report the number of clusters at each. Include a cluster stability assessment, for example
bootstrap resampling with the adjusted Rand index between runs, or `bluster::clusterSweep` if
convenient. The purpose is to document that the absence of substructure is robust rather than a
single-parameter outcome.

### 2. Test for a continuous gradient

Compute PCA and a diffusion map on the female germ cells alone. Report the variance explained by
the leading components.

Then establish what the leading component actually tracks:

- Correlate DC1 and PC1 with total counts (`sum`), genes detected (`detected`), and
  mitochondrial or ribosomal fraction if available. **A first component that tracks detection rate
  is a technical artifact, not a trajectory.** Report these correlations before anything else.
- Correlate DC1 and PC1 with library of origin, to check whether the gradient is a sample effect.
- Plot expression of germ cell and oocyte markers against DC1, using **uncorrected** normalized
  log-expression: DDX4, PIWIL1, SYCP3, ALPL, PRDM14, and the female-associated genes used in this
  project such as CUZD1, S100A1, ZP genes and AEP1. Confirm each symbol exists in `rownames` and
  report any that do not.
- For each marker, fit a monotonic trend against DC1 and report direction and effect size.

### 3. Decide, and state the decision

State plainly which of the following the data support:

- **A continuum exists.** DC1 tracks maturation markers and not detection rate. In this case
  compute diffusion pseudotime or fit a principal curve, order the cells, and repeat the marker
  trend analysis against that ordering. Report the ordering with cells colored by library.
- **No structure.** The leading components track technical covariates, or markers are flat. In
  this case say so and stop. Do not fit a trajectory to noise.

A pseudotime can always be computed from a blob of cells, so the marker and technical-covariate
checks decide this, not the existence of an ordering.

### 4. Compare to the pre-filtering result

Report how many female germ cells were present before doublet filtering and how many survive, and
whether the pre-filtering cells did resolve into subclusters at the same `k`. If they did and the
filtered cells do not, state whether that difference is attributable to cell number or to the
removal of specific populations.

## Outputs

- `SuppTable_female_continuum_diagnostics.csv` with, per leading component, the variance explained
  and the correlations with technical covariates and library.
- `SuppTable_female_marker_trends.csv` with the marker trend direction, effect size and
  significance against the ordering.
- Figures: cluster count against `k`; DC1 against total counts and against genes detected; marker
  expression against DC1 colored by library.
- A closing section reporting the decision from section 3 in prose, plus two or three sentences
  suitable for the Methods describing what was done and what was found.

## Style

- Match the conventions in the existing Rmds.
- US English spelling. Never use em dashes; use commas, semicolons, colons or separate sentences.
- Be concise and direct. Flag anything you cannot verify from the objects rather than assuming it.
- If the leading component correlates with detection rate more strongly than with any marker, say
  so plainly and conclude that no trajectory is supported. That is a real result and I need to
  know it rather than have it worked around.
