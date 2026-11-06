# Analyses required by the reviewer comments

Scope: only comments requiring computation. Text edits, antibody validation and figure
work are excluded. Ordered by value, not by reviewer numbering.

Decision taken: no new sequencing. Analyses therefore serve to characterize what the existing
four libraries can and cannot support, rather than to add power.

---

## 1. Cell counts per cluster per sample

**Answers:** R2 minor 5 directly. Also supplies evidence for R2 major 1 (sex confound), R2 major
2 (stage confound) and R1's replication concern. Highest value per unit effort in this list.

R2 asks specifically whether clusters 16 to 19 are dominated by the single male library. Answer
that with the table rather than in prose.

```r
library(SingleCellExperiment); library(dplyr); library(tidyr)

tab <- as.data.frame(table(cluster = colLabels(sce), sample = sce$Sample)) %>%
  pivot_wider(names_from = sample, values_from = Freq)

tab <- tab %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  rowwise() %>%
  mutate(across(where(is.numeric) & !c(total),
                ~ .x, .names = "n_{.col}")) %>%
  ungroup()

# proportions alongside counts
prop <- as.data.frame(prop.table(table(colLabels(sce), sce$Sample), margin = 1)) %>%
  setNames(c("cluster", "sample", "prop")) %>%
  pivot_wider(names_from = sample, values_from = prop)

write.csv(tab,  "SuppTable_cluster_by_sample_counts.csv", row.names = FALSE)
write.csv(prop, "SuppTable_cluster_by_sample_props.csv", row.names = FALSE)

# the specific claim R2 asks about
prop %>% filter(cluster %in% c("16","17","18","19"))
```

Repeat for the 15 germ cell subclusters, since the trajectory claims rest on those.

**Expect this to confirm the reviewers' suspicion** for clusters 16 to 19. Report it plainly;
the table is the concession, and arguing around it will not work.

---

## 2. Stage composition along the trajectory

**Answers:** R2 major 2.

Use `pseudotime_stage_composition.R` in this folder. The key result is already visible in
Figure 7A: the PGC root and the intermediate nodes contain cells from all four libraries, while
terminal nodes are library-specific. Extract the per-node percentages and quote them.

Also report, for each MST node, the count and proportion from each library, as a supplementary
table. This is the same table as task 1 restricted to germ cell subclusters, so produce once and
cite twice.

---

## 3. Verify that "transcript isoforms" are distinct loci

**Answers:** R2 major 4. Reviewer explicitly asks for verification "using genomic coordinates and
gene models".

Gene series to check: `wnt5*`, `hykk*`, `spopl*`, `cuzd1*`, `zp*`, `s100a1*`, `aep1*`, `cdh4*`,
plus anything else described in the text as an isoform, transcript or variant.

```bash
# GFF for the assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/993/605/GCF_010993605.1_kPetMar1.pri/GCF_010993605.1_kPetMar1.pri_genomic.gff.gz

# coordinates for a list of LOC IDs
zcat GCF_010993605.1_kPetMar1.pri_genomic.gff.gz | awk -F'\t' '$3=="gene"' \
  | grep -F -f loc_ids.txt \
  | sed -E 's/.*ID=gene-([^;]+).*/\1\t&/' \
  | cut -f1,2,5,6,8 > gene_coords.tsv
```

```r
# tidier version
library(rtracklayer)
gff <- import("GCF_010993605.1_kPetMar1.pri_genomic.gff.gz")
genes <- gff[gff$type == "gene"]

ids <- c("LOC116937338","LOC116949348","LOC116937873","LOC116938414",   # wnt5
         "LOC116957691","LOC116957686","LOC116957689","LOC116957684",   # hykk
         "LOC116957694","LOC116957693","LOC116957687","LOC116956112")   # spopl, cdh4, cuzd1

res <- as.data.frame(genes[genes$gene %in% ids])[, c("gene","seqnames","start","end","strand")]
res
```

**What to report:** seqname, start, end and strand for each. Distinct non-overlapping intervals,
especially on different scaffolds, settle it. Include as a supplementary table, since the
reviewer asked for verification rather than assertion.

**Then revise the terminology everywhere**, including the abstract. Suggested wording is in
`manuscript_review_notes.md` section 3.

---

## 4. GO thresholds: adjusted versus raw p

**Answers:** R2 minor 6.

Methods state Benjamini-Hochberg correction; Results filter at p < 0.01. Reviewer wants to know
which. As written, FDR is computed and then discarded.

```r
# re-tabulate existing GO results at an explicit FDR threshold
go_results %>%
  mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)

# how many terms survive under each rule
c(raw_p01 = sum(go_results$pvalue < 0.01),
  fdr05   = sum(p.adjust(go_results$pvalue, "BH") < 0.05),
  fdr01   = sum(p.adjust(go_results$pvalue, "BH") < 0.01))
```

**Decide and be consistent.** If the headline terms survive FDR correction, switch to FDR and say
so. If they do not, that is important to know before responding, and the honest move is to report
raw p explicitly and describe the analysis as exploratory. Do not leave the Methods and Results
describing different procedures.

Note also that the GO universe was all assayed genes (21,967) while only 14,859 of 17,580
protein-coding genes carry any annotation. Consider reporting results with the annotated set as
universe as a sensitivity check.

---

## 5. Direct search for DMRT1 in the assembly

**Answers:** R2 minor 7. Reviewer correctly notes that absence from the annotation does not
demonstrate absence from the genome.

Three approaches, cheapest first.

**a. Search your existing annotation for the DM domain.** If any lamprey gene carries a DM
DNA-binding domain, it is a *dmrt*-family member regardless of its symbol. This uses data you
already have.

```bash
# find DM-domain InterPro accessions among your annotations
# VERIFY the accession at interpro.org before relying on it; IPR001275 is the DM DNA-binding domain
cat GCF_010993605.1_kPetMar1.pri_gene_ontology.gaf | grep -v '^[!#]' \
  | grep 'IPR001275' | cut -f2,3 | sort -u
```

**b. tblastn a DMRT1 protein against the genome.**

```bash
makeblastdb -in GCF_010993605.1_kPetMar1.pri_genomic.fna -dbtype nucl -out kpetmar1
tblastn -query dmrt1_proteins.faa -db kpetmar1 -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length evalue bitscore sstart send' \
  -out dmrt1_hits.tsv
```

Use several query proteins: human DMRT1, mouse Dmrt1, zebrafish dmrt1, and a chicken or
amphibian sequence. Report best hits with coordinates, identity and E-value.

**c. Synteny**, if hits are found: check whether flanking genes correspond to the DMRT1
neighborhood in gnathostomes. Only worth doing if (b) gives candidates.

**Report the outcome either way.** A negative result stated as "no DM-domain-containing sequence
was recovered by tblastn at E < 1e-5 using queries from four species" is a real finding and much
stronger than "not annotated".

---

## 6. Genes driving the neurogenesis enrichment

**Answers:** R1's fourth major comment.

```r
# genes behind GO:0050767 regulation of neurogenesis in the PGC vs cyst comparison
neuro <- gaf %>% filter(GOID == "GO:0050767") %>% pull(Gene)
hits  <- intersect(neuro, de_genes_pgc_vs_cyst$gene)

de_genes_pgc_vs_cyst %>%
  filter(gene %in% hits) %>%
  arrange(FDR) %>%
  select(gene, logFC, FDR)

# what else are these genes annotated with? tests the pleiotropy explanation R1 proposes
gaf %>% filter(Gene %in% hits) %>% count(GOID, sort = TRUE)
```

R1 suggests the enrichment may reflect pleiotropic functions in migration, polarity,
cytoskeletal organization or signaling rather than neurogenesis. The second query tests that
directly: if the same genes carry cytoskeletal or migration terms, say so. This also connects to
the family-level annotation argument already made for reviewer 2 major 3, since a GO term name
reflects the system in which the function was characterized.

Do the same for GO:0003341 cilium movement and GO:0060271 cilium assembly, which are more
plausible in germ cells and worth discussing rather than leaving unexplained.

---

## 7. HCR probe specificity against wnt5 paralogs

**Answers:** R2 major 5, the computational half. The antibody half is not yours.

If probe sequences are available from the vendor, align them against the four wnt5-family
transcripts to check for cross-hybridization.

```bash
# transcript sequences for the four loci, then check probe matches
blastn -task blastn-short -query hcr_probes.fa -subject wnt5_family_transcripts.fa \
  -outfmt '6 qseqid sseqid pident length mismatch evalue' -evalue 1000
```

Report the number of mismatches between each probe and each non-target paralog. If probes cannot
distinguish paralogs, say so and qualify the corresponding claims.

---

## 8. QC reporting, and the doublet gap

**Answers:** R2 minor 1.

The reviewer asks for per-sample viability, QC thresholds, and how doublets were identified or
removed. The Methods describe `emptyDrops` but **no doublet detection is mentioned anywhere**.

Assemble a per-sample QC table: cells before and after filtering, median UMIs, median genes
detected, mitochondrial percentage if computed, `emptyDrops` FDR threshold, and any other cutoffs.

```r
qc <- data.frame(
  sample        = names(sce_list),
  droplets      = sapply(sce_list, ncol),
  cells_kept    = sapply(sce_list_filt, ncol),
  median_umi    = sapply(sce_list_filt, function(x) median(colSums(counts(x)))),
  median_genes  = sapply(sce_list_filt, function(x) median(colSums(counts(x) > 0)))
)
write.csv(qc, "SuppTable_QC_per_sample.csv", row.names = FALSE)
```

### Doublets: diagnose, do not filter

Simulated-doublet methods assume a doublet resembles a mixture of two distinct clusters. In a
differentiation continuum, genuine intermediates have the same signature, so the method's
positive signal and the biology are confounded. This dataset is a particularly poor case for
naive filtering:

- The germ cell data are explicitly a trajectory, so intermediates are the object of study.
- The manuscript describes germ cell cysts and a cytokinetic spermatocyte cluster. Germ cells in
  cysts are connected by intercellular bridges, and cells in cytokinesis are physically two
  nuclei in shared cytoplasm. These are not artifacts.
- Oocytes carry high RNA content, inflating library size, which is the other feature doublet
  callers use.

Filtering could remove the trajectory intermediates and the cytokinetic cluster, which is to say
the result. **Run detection as a diagnostic and report it; do not subset.**

```r
library(scDblFinder)
sce <- scDblFinder(sce, samples = sce$Sample)   # scores only, no filtering

as.data.frame(colData(sce)) %>%
  group_by(label) %>%
  summarise(median_score = median(scDblFinder.score),
            pct_called   = mean(scDblFinder.class == "doublet") * 100,
            n = n()) %>%
  arrange(desc(median_score))
```

**The more informative check, and the one to report.** A genuine doublet cluster expresses
markers of two incompatible lineages simultaneously and has elevated library size. A genuine
intermediate shows a graded shift within one lineage.

```r
as.data.frame(colData(sce)) %>%
  group_by(label) %>%
  summarise(median_umi = median(sum), median_genes = median(detected), n = n())

# germ versus somatic markers in the same cells
plotExpression(sce, features = c("DDX4","PIWIL1","SYCP3","FOXL2","S100A1","SOX18"), x = "label")
```

If no cluster combines elevated library size with simultaneous expression of incompatible
lineage markers, report that with numbers. It answers the reviewer's question more convincingly
than a doublet score, and it justifies retaining the cells.

Also report the platform's expected multiplet rate at the loading concentration used
(4 x 10^6 cells/ml, PIPseq T20 v4.0) rather than estimating it, and state explicitly that
cells were not removed on the basis of doublet scores, with the reasoning above.

---

## 9. Confirm the two-sample t-test is two-sided

**Answers:** R1 minor, page 9 line 36.

Check the call in the analysis code. `findMarkers` in scran uses two-sided tests by default
unless `direction=` is set. State the direction explicitly in the Methods.

---

## Not analysis, but blocked on you

- **Sample ID discrepancy** (R2 minor 3): Methods say C00003262, Results say C00004362. Resolve
  against GEO GSE334273.
- **Larval replicate definition** (R2 minor 2): whether each larval library is one animal or a
  pool, and whether dissociation and library prep were independent. You know the pipeline; the
  wet-lab detail may need Chung-Davidson.
- **Batch structure** (R2 major 2): already drafted in `reviewer_response_CONSOLIDATED.md`.
