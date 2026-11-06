# Resolving the DMRT1 question

## RESULT (checked 2026-08-03)

**Assembly coordinates match.** kPetMar1.pri chr24 is NC_046092.1, length 14,075,897 bp. The
reported DMRT1 interval, 14,041,939 to 14,064,356, falls inside it, ending about 11.5 kb from the
chromosome end.

**Release 100 annotates no gene over that interval.** The nearest annotated genes are:

| Gene | GeneID | Coordinates | Strand |
|---|---|---|---|
| TRAF6 | 116945680 | 13,975,597 to 13,987,716 | minus |
| TTC17 | 116945679 | 13,987,888 to 14,040,537 | plus |

TTC17 ends 1,402 bp before the DMRT1 locus begins. The final ~35 kb of chr24, from 14,040,537 to
the chromosome end, carries no gene model in Release 100.

**Consequence.** Reads from the DMRT1 locus were not assigned to any feature during
quantification, so the existing count matrices contain no information about it. Its expression
cannot be reported without re-quantification. Proceed to step 4 below.

**Likely explanation, usable in the manuscript.** The locus is subtelomeric. Transcript evidence
is typically sparser and assembly quality lower at chromosome ends, which is a plausible reason
Gnomon did not build a model there while a PASA-based re-annotation with targeted evidence did.

**Synteny check now available, and explicitly requested by reviewer 2.** The flanking genes are
TRAF6 and TTC17. Determine whether either is a DMRT1 neighbor in gnathostomes. Agreement
corroborates the assignment; disagreement is worth reporting too. This answers the "where
feasible, synteny" part of the comment without any new sequencing.

---

Reviewer 2 minor comment 7 asks whether a *dmrt1* homolog was searched for directly in the
assembly. Separately, Jeramiah Smith's lab emailed to say they annotated one.

- Paper: Timoshevskiy et al., PNAS, PMID 40504158, doi 10.1073/pnas.2421883122
- Annotation: `GSE280458_PM.gene_annot_2023.exons_renum.alias.gff3.gz` at GEO GSE280458
- Reported DMRT1 location: **chr24:14,041,939-14,064,356**

**Do not re-quantify before running the check below.** The likely outcome is that Release 100
already has a gene model over this interval, under a LOC identifier, in which case your existing
count matrices already contain these reads and nothing needs rerunning.

---

## Step 1: is the annotation on the same assembly?

Everything downstream depends on this. Their "chr24" must correspond to kPetMar1.pri chr24.
The Smith lab has published assemblies including germline-specific chromosomes, which
kPetMar1.pri does not fully contain, so this is not safe to assume.

```bash
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE280458&format=file&file=GSE280458%5FPM%2Egene%5Fannot%5F2023%2Eexons%5Frenum%2Ealias%2Egff3%2Egz" \
  -O PM.gene_annot_2023.gff3.gz

# header, and the seqid naming convention
zcat PM.gene_annot_2023.gff3.gz | head -30
zcat PM.gene_annot_2023.gff3.gz | awk -F'\t' '!/^#/ {print $1}' | sort -u | head -30

# sequence lengths, if declared, are the decisive check against kPetMar1
zcat PM.gene_annot_2023.gff3.gz | grep '^##sequence-region' | head -30
```

Compare chromosome lengths against the kPetMar1.pri assembly report:

```bash
wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/993/605/GCF_010993605.1_kPetMar1.pri/GCF_010993605.1_kPetMar1.pri_assembly_report.txt
grep -v '^#' GCF_010993605.1_kPetMar1.pri_assembly_report.txt | awk -F'\t' '{print $1, $3, $7, $9}' | head -30
```

If lengths match, the coordinate systems are the same and you can proceed directly. If they do
not, skip to step 4.

## Step 2: find their DMRT1 entry

```bash
zcat PM.gene_annot_2023.gff3.gz \
  | awk -F'\t' '$1=="chr24" && $4 < 14064356 && $5 > 14041939' | head -20

# the alias field is why this file is useful; see what they call it
zcat PM.gene_annot_2023.gff3.gz | grep -i 'dmrt' | head
```

Record the gene ID, exact coordinates and strand.

## Step 3: does Release 100 have a gene there?

This is the decisive question for whether re-analysis is needed.

```bash
# RefSeq accession for chr24 in kPetMar1.pri
CHR24=$(grep -P '^24\t' GCF_010993605.1_kPetMar1.pri_assembly_report.txt | cut -f7)
echo "chr24 = $CHR24"

zcat GCF_010993605.1_kPetMar1.pri_genomic.gff.gz \
  | awk -F'\t' -v c="$CHR24" '$1==c && $3=="gene" && $4 < 14064356 && $5 > 14041939'
```

**If a gene model is returned:**

Extract its LOC identifier. That gene is in your count matrix already. Check it directly:

```r
loc <- "LOC1169XXXXX"     # from the awk output above
loc %in% rownames(merged.sce)

plotExpression(merged.sce, features = loc, x = "label")
plotTSNE(merged.sce, colour_by = loc)

# is it a cluster marker anywhere?
lapply(scatter_markers, function(x) x[loc, c("summary.logFC","FDR")])
```

Then simply report it. No re-quantification, no re-clustering, no figure regeneration beyond one
added panel. Revise the manuscript text from "DMRT1 is not annotated" to a statement naming the
LOC gene, citing Timoshevskiy et al. for the homology assignment, and describing its expression.

**If no gene model is returned:**

Reads from that interval were never assigned to any feature, so your matrices contain no
information about it and its expression cannot be reported without re-quantification. Go to
step 4.

## Step 4: if the locus is unannotated in Release 100

Options, cheapest first.

**a. Correct the statement and stop.** Honest, adequate for the reviewer, costs nothing. Draft:

> "A *dmrt1* homolog is not present in NCBI *Petromyzon marinus* Annotation Release 100, the
> annotation against which our data were quantified. Following submission, a *dmrt1* homolog was
> reported in an independent annotation of this assembly at chr24 (NC_046092.1):14,041,939 to
> 14,064,356 (Timoshevskiy et al. 2025). We confirmed that this interval carries no gene model in
> Release 100: the nearest annotated genes are *TTC17*, ending at 14,040,537, and *TRAF6*, and the
> terminal ~35 kb of chromosome 24 is unannotated in that release. Reads originating from this
> locus were therefore not assigned to any feature during quantification, and we are unable to
> report its expression in the present dataset. We note this as a limitation, and that the
> subtelomeric position of the locus is a plausible reason for its absence from the automated
> annotation. [Add synteny result here.]"

This is stronger than the current text because it reports a search rather than an absence, gives
coordinates, and explains the discrepancy.

**b. Targeted re-quantification of one locus.** Only if you have access to alignments or raw
reads per cell barcode. Count reads overlapping the interval per barcode and add a single row to
the matrix. Feasible but depends on what the vendor returned; if you only have gene-level
matrices, this is not possible without realignment.

**c. Full re-quantification against the 2023 annotation.** Not warranted for one gene. It would
invalidate every existing cluster, marker list, GO result and figure, and the annotation change
would itself need justifying to reviewers mid-revision.

Recommendation: (a), unless (b) turns out to be cheap.

---

## Worth doing regardless of the outcome

**Reply to Smith.** They volunteered useful information about your preprint. A citation costs
nothing and they may have more relevant annotation. Their group produced the germline chromosome
work you already cite (Timoshevskaya et al. 2023, Timoshevskiy et al. 2025), so this is a natural
connection.

**Consider their alias file as corroboration for gene naming.** Reviewers questioned the Wang et
al. symbol assignments, which are unidirectional top-BLAST-hit calls. If the 2023 annotation
carries independent gene names, agreement between the two for your key marker genes is orthogonal
support, obtainable without changing any analysis:

```bash
# their names for loci overlapping your marker genes
zcat PM.gene_annot_2023.gff3.gz | awk -F'\t' '$3=="gene"' | grep -iE 'ddx4|piwil1|sycp3|foxl2|sox18|rbpms|s100a1|boll|cxcr4'
```

Do not switch annotations. Use it only as a cross-check, and report agreement or disagreement
per gene.

**This strengthens the paper independently of the reviewers.** DMRT1 is the canonical vertebrate
sex-determination gene, and the manuscript argues lamprey occupies a transitional state between
environmental and genetic sex determination. Being able to describe *dmrt1* expression across
germ and somatic clusters, even if the pattern is unremarkable, is a substantive addition to that
argument. The current text has to say the gene could not be assessed.
