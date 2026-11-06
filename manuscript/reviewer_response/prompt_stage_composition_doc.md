# Claude Code prompt: shared populations and trajectory support

Answers reviewer 2 major comment 2, that developmental stage is confounded with collection time
and technical batch, and that fastMNN cannot by itself distinguish a genuine stage effect from a
fully confounded batch effect.

Paste everything below the line.

---

Write a new R Markdown document, `src/05_stage_composition.Rmd`, assessing what evidence exists
that the germ cell trajectory reflects development rather than the collection batch. Do not modify
any existing script or Rmd, and do not rerun clustering, integration or trajectory inference. Load
the existing objects.

Read `src/01_all_cells_clustering_annotation.Rmd` and `src/03_trajectory_analysis.Rmd` first, to
find the object names and to match conventions: the YAML params block, `results_dir` and
`Robject_dir` construction, chunk options, `code_folding: hide`, `DT::datatable` for tables, and
the pattern of putting expensive steps in `eval=FALSE` chunks with pre-computed objects loaded for
rendering.

## Why this is needed, and what does not count as evidence

The two larval libraries were collected in July 2023 and the two transformer libraries in
November 2023, so developmental stage is confounded with collection date, season and technical
batch.

**An important constraint on the interpretation.** Some of this confound is inherent to the
biology rather than to the design. If primordial germ cells are present only in larvae and
differentiated germ cells only in transformers, then the early positions in the trajectory will
be larval and the late positions transformer no matter how the experiment is run. That
observation is therefore uninformative: a genuine developmental trajectory and a corrected batch
effect both predict it. **Do not treat composition at the ends of the trajectory as evidence
either way.**

What does discriminate is whether the integration is aligning populations that genuinely exist in
both stages, and whether the ordering is supported by evidence that does not depend on library of
origin. Those are the two things this document should establish.

## Sample and stage mapping

Four libraries. Confirm the exact identifiers in the objects rather than assuming:

- Undifferentiated.1 and Undifferentiated.2 are larvae, collected July 2023
- Female and Male are transformers, collected November 2023

Derive a `stage` variable with values "larva" and "transformer", and keep the library identifier
as a separate variable, since the male and female branches are single libraries and that matters
for interpretation.

## Analyses

### 1. Which populations are shared between stages

The primary question. For every cluster in the full dataset and every germ cell subcluster,
report the number and percentage of cells contributed by each library and by each stage. Then
classify each population:

- **shared**, at least 10% of cells from each stage
- **larva-dominated**, over 90% larval
- **transformer-dominated**, over 90% transformer

Report how many populations fall into each class, and list the shared ones by name. Distinguish
somatic from germ cell populations in this summary, because a germ cell state present in both
stages is much stronger evidence than a shared somatic population.

For each shared population, also report whether its marker genes are consistent across stages:
compute markers separately within larval cells and within transformer cells, and report the
overlap in the top 50. A shared cluster whose markers agree across stages is being aligned on
genuine common biology rather than being merged by the correction.

### 2. Marker behavior along pseudotime

Evidence for the ordering that does not depend on library of origin. Using **uncorrected**
normalized log-expression, not the fastMNN output, plot expression of established germ cell
markers against pseudotime for each branch: DDX4, PIWIL1, SYCP3, BOLL, and any others used for
annotation in this project. Confirm each symbol exists in `rownames` and report any that do not.

For each marker and branch, fit a simple monotonic trend and report whether expression increases,
decreases or is flat, with the effect size. Markers behaving as expected across a developmental
progression support the ordering; markers that are flat or move against expectation should be
reported too.

Colour points by library, so it is visible whether any apparent trend is carried entirely by one
library.

### 3. Composition along pseudotime, reported without interpretation

For each branch, bin pseudotime into deciles and report the proportion of cells from each stage
and each library per decile. Report this as a table and a stacked bar figure.

Present this descriptively. Do not use it to argue for or against the trajectory, for the reason
given above. Its purpose is transparency, since the reviewer asked how the libraries are
distributed and the answer should be available.

### 4. Node and subcluster composition

Number and percentage of cells from each library at each node of the minimum spanning tree, and
in each germ cell subcluster. Report every node and subcluster, including any that look anomalous,
and note any dominated by a single library. Do not omit or aggregate these away.

### 5. Figures

- Shared, larva-dominated and transformer-dominated populations, as a labelled bar chart.
- Marker expression against pseudotime, faceted by marker and branch, coloured by library.
- Stage and library proportion by pseudotime decile, faceted by branch.
- Node composition, as stacked bars ordered by pseudotime.

## Outputs

- `SuppTable_population_sharing.csv`, `SuppTable_marker_trends.csv`,
  `SuppTable_pseudotime_composition.csv`, `SuppTable_node_composition.csv`, in `results_dir`.
- The figures above, saved to `results_dir`.
- A short section at the end reporting the key numbers in prose, plus two or three sentences
  suitable for a reviewer response, stating how many populations are shared between stages,
  whether any germ cell populations are among them, and whether marker behavior along pseudotime
  supports the inferred ordering.

## Style

- Match the conventions in the existing Rmds.
- US English spelling. Never use em dashes; use commas, semicolons, colons or separate sentences.
- Be concise and direct. Flag anything you cannot verify from the objects rather than assuming it.
- If few or no germ cell populations are shared between stages, or if markers do not behave
  monotonically along pseudotime, say so plainly. Those are real results and I need to know,
  because they determine whether the trajectory claim can stand.
