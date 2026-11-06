# Claude Code update: fix the comparison background in stage_split_markers

Paste everything below the line as a follow-up.

---

There is a methodological problem with `stage_split_markers` that biases the very comparison it
is meant to support. Please fix it.

## The problem

The function selects clusters with at least `min_cells` **within each stage separately**. So in
larvae, a cluster is tested against whichever clusters passed the threshold in larvae, and in
transformers it is tested against a possibly different set. The comparison backgrounds are not
the same.

That matters because the purpose of the function is to compare marker overlap between stages. A
gene can appear as a marker for a cluster in one stage simply because a competing cluster was
absent from the background there, not because the cluster is genuinely defined differently. With
the current code, low marker overlap between stages is partly an artifact of unequal backgrounds.

The asymmetry is expected to be substantial here: the four meiotic male clusters and the
differentiated oocyte cluster should be sparse or absent in larvae, so they will drop out of the
larval background while remaining in the transformer background.

## The fix

Compute the eligible cluster set as the **intersection** across stages, so that both `findMarkers`
calls use an identical set of clusters, then subset each stage to that common set before testing.

Report explicitly:

- the clusters eligible in each stage individually
- the intersection actually used
- for every cluster excluded, which stage it failed in and its cell count in each stage

The excluded clusters are informative in their own right, since a cluster present in one stage and
absent in the other is a stage-restricted population. Present that as a table rather than letting
it disappear silently.

## Also check

Confirm the `findMarkers` arguments match those used in
`src/01_all_cells_clustering_annotation.Rmd`, in particular `test.type`, `direction` and
`pval.type`. The stage-split comparison is internally consistent either way, but if the settings
differ from the original analysis, any statement relating these markers to the ones reported in
the manuscript is not valid. If they differ, either align them or state the difference in the
document.

## Reporting

Update the marker overlap summary so it is clear that overlap was computed on a common cluster
background. If the overlap changes materially after the fix, report both the before and after
values, since that difference is itself worth knowing.

## Style

US English spelling. Never use em dashes. Be concise and direct. If the intersection leaves too
few clusters for the comparison to be meaningful, say so plainly rather than proceeding.
