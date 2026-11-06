# Claude Code prompt: doublet diagnostics

Two prompts. Run the first against the per-sample QC Rmds, the second once the scores have
propagated to the integrated object.

The work splits because `scDblFinder` should run per sample on unintegrated data, where the
neighborhood structure it simulates from is intact, but the interpretation depends on clusters,
which only exist after integration.

---

## Prompt 1: modify the per-sample QC Rmds

Paste everything below the line.

---

I have per-sample R Markdown QC documents for a single-cell RNA-seq experiment. Read them first
to learn the existing structure, object names and conventions, then modify each to add doublet
scoring. Keep the changes minimal and consistent with the existing style.

### Context

Sea lamprey gonad, PIP-seq (Fluent BioSciences PIPseq T20 3' v4.0, loaded at 4e6 cells/ml).
Four libraries: two larval replicates, one female transformer, one male transformer. Current QC
uses `emptyDrops` from DropletUtils. Downstream steps use scran, scater and batchelor, ending in
fastMNN integration and Louvain clustering.

A reviewer asked how doublets were identified or removed. No doublet detection was performed.
I am adding it as a diagnostic for the revision.

### The critical constraint

**Compute doublet scores. Do not remove any cells. Do not add any filtering step.**

Simulated-doublet methods assume a doublet resembles a mixture of two distinct clusters. Three
features of this dataset share that signature but are real biology:

1. The germ cell data are an explicit differentiation trajectory, so transcriptional
   intermediates are the object of study and genuinely resemble a mixture of the states on
   either side of them.
2. Oocytes have high RNA content, which inflates library size, the other feature doublet callers
   use. Elevated library size in an oocyte cluster is expected and is not evidence of multiplets.

A third possibility, worth testing rather than assuming: germ cells develop in cysts, and if
cysts were incompletely dissociated they would appear as multiplets. This is a hypothesis, not
an established feature of this dataset. Dissociation used papain with repeated trituration and
two passes through a 40 micron strainer, which would break up or exclude most cysts, so do not
assert it. Report it only if the data support it.

The goal is to characterize and report, not to clean.

### What to add to each Rmd

1. A new section after the existing QC filtering, running `scDblFinder` on that sample's
   filtered SCE. Run it per sample, without a `samples =` argument, since each document handles
   one library.

2. Store `scDblFinder.score` and `scDblFinder.class` in `colData`. **Verify they persist into
   whatever object the document saves for downstream use.** This matters: `fastMNN` returns a new
   SingleCellExperiment carrying only the corrected values, so colData will not survive
   integration automatically. If the saved object is later rebuilt, note in the document what
   needs to be carried forward and by what key, usually the cell barcode in `colnames`.

3. A short summary in each document: number and percent of cells called doublet, median and
   interquartile range of the score, and the score distribution as a histogram or violin.

4. A table of median UMI count and median genes detected, split by `scDblFinder.class`. If cells
   called doublet do not have elevated library size, that is informative and worth stating.

5. A sentence in the document stating that no cells were removed on the basis of these scores,
   with the reasoning above in brief.

### Style

- Match the existing Rmd conventions: chunk naming, caching, figure sizing, output format.
- US English spelling. Never use em dashes; use commas, semicolons, colons or separate sentences.
- Do not restructure or reorganize the existing documents beyond adding this section.
- If the documents differ from one another in structure, adapt rather than forcing uniformity,
  and tell me where they diverge.

---

## Prompt 2: cluster-level interpretation, after integration

Run this once the scores are present in the integrated object.

Paste everything below the line.

---

Doublet scores from `scDblFinder` are now in the integrated SingleCellExperiment for a lamprey
gonad single-cell dataset. Write a new analysis script interpreting them at cluster level. Do not
remove any cells.

There are 19 clusters in the full dataset and a germ cell subclustering of 15 subclusters used
for a TSCAN pseudotime trajectory.

### Analysis

1. Confirm the scores actually propagated through `fastMNN` and were not silently dropped or
   misaligned. Check that `colnames` match between the pre-integration and integrated objects
   before trusting any of what follows. Report if they do not.

2. Per-cluster diagnostic table: cluster, n cells, median doublet score, percent called doublet,
   median UMI count, median genes detected. Sorted by percent called doublet. Produce for the
   full dataset and for the germ cell subset.

3. The discriminating test: does any cluster combine a high doublet rate **with** elevated
   library size? That combination indicates true multiplets. A high score with ordinary UMI
   counts is more consistent with a transitional state the simulator cannot distinguish from a
   mixture. State which pattern each flagged cluster shows.

   Treat the oocyte cluster separately, since high RNA content there is expected biology rather
   than evidence of multiplets. Do not use cluster names as evidence of physical cell state:
   cluster labels here derive from marker expression and GO enrichment, so a label referring to
   cytokinesis or to cysts describes a transcriptional program, not an observed morphology.

4. Within flagged clusters, report the fraction of cells co-expressing germ cell markers
   (DDX4, PIWIL1, SYCP3) and somatic markers (FOXL2, S100A1, SOX18, RBPMS) simultaneously.
   Genuine doublets co-express incompatible lineages; genuine intermediates shift within one.

5. For the germ cell subset, plot doublet score against pseudotime. If elevated scores
   concentrate at intermediate pseudotime, that supports a transitional interpretation and is
   worth a supplementary panel.

### Outputs

- `SuppTable_doublet_diagnostics.csv` for both the full dataset and the germ cell subset.
- Figures: doublet score per cluster (violin), and doublet score against pseudotime.
- A short markdown summary with the numbers, naming which clusters were flagged and whether each
  fits the multiplet or the intermediate pattern.
- Two or three sentences for the Methods, and two or three for the reviewer response, stating
  what was done, what was found, and that cells were retained.

### Style

- R, matching the conventions in this repo.
- US English spelling. Never use em dashes.
- Be concise and direct. Flag anything you cannot verify from the data rather than assuming it.
- If the numbers contradict the expectation above, say so plainly. If a cluster does look like
  genuine multiplets, that is a real finding and I need to know.
