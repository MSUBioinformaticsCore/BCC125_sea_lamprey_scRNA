# Claude Code prompt: male germ cells as a positive control for the continuum method

Companion to `prompt_female_germ_cell_continuum.md`. The female germ cells do not subcluster after
doublet filtering, so a continuum method is being considered for them. The male germ cells do
subcluster and have an existing TSCAN trajectory, which makes them a positive control: if the
continuum method reproduces the TSCAN ordering in males, applying it to females is justified.

Paste everything below the line.

---

Write a new R Markdown document, `src/07_male_germ_cell_continuum.Rmd`, applying the same analysis
used for the female germ cells in `src/06_female_germ_cell_continuum.Rmd` to the doublet-filtered
male germ cells, and using the male cells to validate the method.

Read `src/06_female_germ_cell_continuum.Rmd` first and **use identical methods, parameters and
functions throughout**. The two documents must be directly comparable. Where a choice was made for
the female analysis, make the same choice here rather than optimizing separately. Also read
`src/02_sex_cell_subanalysis.Rmd` and `src/03_trajectory_analysis.Rmd` for object names and
conventions.

## Why this document exists

The male germ cells resolve into discrete subclusters and have an existing TSCAN pseudotime. The
female germ cells do not resolve. If a continuum method applied to the male cells recovers
substantially the same ordering as TSCAN, the method is validated on data where the answer is
known, and its application to the female cells is defensible. If it does not, the female result
cannot be trusted either, and that needs to be known before anything is written.

## Analyses

### 1. Clustering, for comparability

Run `buildSNNGraph` plus Louvain across the same range of `k` used in the female document, and
report the number of clusters at each, with the same stability assessment. This documents that
clustering succeeds here and fails there under identical settings, which is the contrast the
Methods will rest on.

### 2. Continuum diagnostics, identical to the female analysis

Compute PCA and a diffusion map on the male germ cells alone, using the same parameters.

Report before anything else, exactly as in the female document:

- correlation of DC1 and PC1 with total counts (`sum`), genes detected (`detected`), and
  mitochondrial or ribosomal fraction if available
- correlation of DC1 and PC1 with library of origin

If DC1 tracks detection rate in the male cells too, then the female result is uninterpretable
regardless of what it shows, and this needs stating prominently.

### 3. Marker trends

Plot expression against DC1 on **uncorrected** normalized log-expression, using male germ cell and
spermatocyte markers: BOLL, PIWIL1, SYCP3, DDX4, and any others used to annotate the male
subclusters in this project. Confirm each symbol exists in `rownames` and report any that do not.
Fit monotonic trends and report direction and effect size, using the same procedure as the female
document.

### 4. The validation: continuum ordering against TSCAN

This is the point of the document.

Compute diffusion pseudotime or fit a principal curve to the male germ cells, using the same
method chosen for the female analysis. Then compare that ordering to the existing TSCAN pseudotime
for the male path:

- Spearman correlation between the two orderings, overall and within each male subcluster.
- A scatter plot of the two orderings, with cells colored by subcluster and by library.
- Whether the subclusters appear in the same order along both.

State plainly whether the continuum method reproduces the TSCAN ordering. A strong correlation
means the method works on this dataset. A weak one means it does not, and the female result must
be treated as uninterpretable rather than as evidence of no structure.

### 5. Direct comparison to the female result

A short section tabulating the two analyses side by side: number of cells, whether clustering
resolves subpopulations, variance explained by leading components, correlations with technical
covariates, and marker trend results. This is what the Methods and the reviewer response will draw
on.

## Outputs

- `SuppTable_male_continuum_diagnostics.csv` and `SuppTable_male_marker_trends.csv`, matching the
  female tables in structure.
- `SuppTable_male_female_comparison.csv` for section 5.
- Figures: cluster count against `k`; DC1 against total counts and genes detected; marker
  expression against DC1 colored by library; continuum ordering against TSCAN pseudotime.
- A closing section reporting whether the continuum method is validated in the male cells, and
  what that implies for the female analysis, plus two or three sentences suitable for the Methods.

## Style

- Match the conventions in the existing Rmds.
- US English spelling. Never use em dashes; use commas, semicolons, colons or separate sentences.
- Be concise and direct. Flag anything you cannot verify from the objects rather than assuming it.
- If the continuum method does not reproduce the TSCAN ordering in the male cells, say so plainly.
  That would invalidate the approach for the female cells and I need to know it rather than have
  it worked around.
