# Claude Code update: benchmark and stabilize the marker overlap comparison

Paste everything below the line as a follow-up.

---

The marker overlap result needs a benchmark and a more stable metric before it can be
interpreted. Currently we have: across 6 shared germ cell populations, top 50 markers computed
within larval cells and within transformer cells overlap by 2% to 42%, median 38%.

That number is far above chance, since two random top 50 lists drawn from roughly 22,000 genes
would overlap by about 0.1 genes. But it cannot be interpreted without knowing how much overlap
is achievable at all, because top 50 marker lists are unstable even between random halves of
identical data. Please add the following.

## 1. Split-half control, to establish the ceiling

Within a single stage, split the cells of the shared clusters at random into two halves, compute
markers on each half, and calculate top 50 overlap exactly as done for the cross-stage
comparison. Use the same cluster set, the same `findMarkers` arguments, and the same overlap
function.

Requirements:

- Run within larval cells and within transformer cells separately.
- Repeat over at least 10 random splits and report the median and range, since a single split is
  noisy.
- Report per cluster as well as overall.

This gives the ceiling. If the within-stage split-half median is close to 38%, the cross-stage
agreement is essentially as high as the data allow, and the clusters should be described as
consistent between stages. If the split-half median is much higher, say 70% or more, then the
cross-stage figure indicates real divergence. Report the comparison explicitly rather than
leaving the reader to judge 38% in isolation.

## 2. Rank correlation across all genes

Top 50 overlap discards most of the information and depends on an arbitrary cutoff. Add, for
each shared cluster, the Spearman correlation between the larval and transformer test statistics
across all tested genes, using whichever statistic `findMarkers` returns for the relevant
comparison, for example the summary log fold change.

Report this alongside the overlap percentages. It is far more stable to differences in cell
number, and it may give a cleaner result. Compute the same quantity for the split-half control so
the two are comparable.

## 3. Identify and explain the 2% case

One of the six shared populations shows essentially no marker agreement between stages. Do not
leave it inside a range.

- Name the population.
- Report its cell count in each stage.
- Report its split-half overlap within each stage.

If it has very few cells in one stage, low overlap is most likely a power artifact rather than a
biological difference, and should be described that way. If it has ample cells in both stages and
still shows no agreement, that is a substantive finding and the cluster's identity is questionable.
State which of these applies.

## 4. Cell counts throughout

Report cells per stage for every shared population wherever overlap is reported. Overlap
percentages are not interpretable without them, and any population where one stage contributes
only a small number of cells should be flagged.

## Outputs

- Update `SuppTable_population_sharing.csv` to include, per cluster: cells per stage, cross-stage
  top 50 overlap, within-stage split-half overlap for each stage, and cross-stage rank
  correlation.
- A figure comparing cross-stage overlap against the split-half ceiling, per cluster.
- Revise the summary text to state the cross-stage figure relative to the ceiling, not on its own.

## Style

US English spelling. Never use em dashes. Be concise and direct. If the split-half control shows
that the ceiling is low, say so plainly: that would mean marker list overlap is not a usable
measure here and we should rely on the rank correlation instead.
