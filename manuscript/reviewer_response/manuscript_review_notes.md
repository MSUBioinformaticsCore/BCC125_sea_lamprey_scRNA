# Manuscript read-through: issues found

My read of `PIP-seq_SL_Sex_v6_text_SLH_genes_italic.docx`. Verify each before acting; I am
working from the text alone, without the figures or supplementary tables.

## Tier 1: fix before anyone else reads it

### 1. Words are missing from at least five sentences

Almost certainly lost when gene names were italicized. Every instance I found:

- **Abstract:** "female germ cells differentially expressed various transcript isoforms of **,**
  CUZD1, S100A1, and ZP genes" — a gene name is missing after "of".
- **Abstract:** "genes with changing **[ ]** along the female differentiation pathway" and the
  same construction again for the male pathway. Probably "expression".
- **Results para 48:** "The top 50 genes exhibiting dynamic **[ ]** along the male pathway".
- **Results para 62:** "Notably, **[ ]** was expressed in 10% of fibroblasts originating from the
  male sample" — the subject is missing. Context suggests SOX9.
- **Results para 62:** "DMRT1 is not annotated in the sea lamprey **[ ]**" — probably "genome".
- **Results para 52:** "SPOPL.2" appears without its LOC number while every sibling gene has one.

Diff v6 against v5 to catch any I missed. This is the single most damaging category because it
is visible on a first read and suggests the manuscript was not proofed.

### 2. NR5A1 is called the androgen receptor, and the paper contradicts itself

Results para 60: "pre-granulosa and granulosa cells differentially expressed CYP17A1 and
**androgen receptor NR5A1**".

NR5A1 encodes steroidogenic factor 1 (SF-1), a nuclear receptor. The androgen receptor is AR
(NR3C4). Para 62 then correctly says "SF-1, the protein encoded by NR5A1 gene". Fix para 60.

Also check **CYP11b2.1**: 11β-hydroxylase is CYP11B1; CYP11B2 is aldosterone synthase. Confirm
which you mean.

### 3. The "transcript isoform" error is pervasive and changes the biology

Distinct LOC numbers are distinct loci, not transcripts of one gene. Instances found:

| Location | Text | Reality |
|---|---|---|
| Abstract | "various transcript isoforms of CUZD1, S100A1, and ZP genes" | distinct genes |
| Para 44 | "multiple CUZD1 transcript isoforms... CUZD1.15: LOC116956112" | distinct genes |
| Para 49 | "WNT5 transcript variants" (LOC116937338, LOC116949348, LOC116937873, LOC116938414) | four loci |
| Para 52 | HYKK, HYKK.5, HYKK.6, HYKK.4 (four LOC numbers) | four loci |
| Para 52 | SPOPL, SPOPL.1, SPOPL.2, CDH4.3 | distinct loci |
| Para 53 | "some transcripts of these genes were enriched in oogonia... other transcripts... in somatic cells and male germ cells" | distinct paralogs |

**Para 53 is the one that matters most.** "Different transcripts of the same gene show different
cell-type specificity" is a claim about isoform-level regulation and is genuinely interesting.
"Different paralogs show different cell-type specificity" is a much more ordinary observation.
As written you are claiming the former and your data show the latter. The same applies to the
WNT5 and CUZD1 results, which are presented as isoform switching.

Suggested replacements: "WNT5-family genes", "several CUZD1-related genes", "four HYKK-related
loci". Where you want to retain emphasis, "paralogs" is defensible if you can support duplication;
"related genes" requires nothing.

## Tier 2: preempt, because reviewers will raise these

### 4. One library per sex

Female (C00004337) and male (C00004374) transformers are single libraries, each pooled from
three animals. Every sex-specific claim, including the divergent male and female trajectories
that the title and abstract are built on, rests on n=1 versus n=1. Pooling reduces individual
variation but provides no replication, so library-level effects cannot be separated from sex.

Reviewer 3 approached this obliquely via batch. Someone will ask directly. State it in the
Limitations yourself, alongside the point that both libraries were processed together so the
comparison is at least within-batch.

### 5. The 40 µm filtration almost certainly biases against oocytes

The cell suspension passes a 40 µm strainer twice, once after dissociation and once before
PIP-seq loading. Any cell larger than 40 µm is excluded. Lamprey oocytes in newly transformed
juveniles are considerably larger than this, so the "oocyte" cluster can only contain the
smallest early oocytes, or fragments.

This matters for interpretation: the female differentiation trajectory terminates in whatever
oocyte-like cells passed a 40 µm filter, not in oocytes generally. Report the size cutoff as a
limitation and describe the terminal state accordingly. Papain-based dissociation with a neuron
isolation kit will also draw a question about its suitability for gonad; a sentence justifying
the choice would help.

### 6. Two homology claims in para 63 are not sound

- "SRCAP transcript (LOC116950250), with **67% identity with PIE-1 mRNA**". SRCAP and PIE-1 are
  not homologs. mRNA-level percent identity to a *C. elegans* gene is not evidence of functional
  equivalence, and at 67% over an unstated length it is not evidence of homology at all.
- "**CDK9 transcript, with 70% identity with p-TEFb**". P-TEFb is a complex; CDK9 is its
  catalytic subunit. A transcript cannot have percent identity to a protein complex.

Either drop these comparisons or replace them with proper reciprocal-BLAST or phylogenetic
support and state alignment length, coverage and E-value. As written they are the weakest
sentences in the paper and they sit in a paragraph making an interesting evolutionary argument.

### 7. Percent-of-cells-expressing comparisons

Para 62 compares "10% of fibroblasts" against "4-5% of oogonia/oocytes". Detection rates in
single-cell data depend on sequencing depth and cell size, both of which differ systematically
across these populations. Either normalize for detection, or present these as observations
rather than as contrasts.

## Tier 3: framing

### 8. Causal language exceeds the evidence

The title says the work "reveals genes involved in germ cell development and sexual
differentiation". The data show differential expression. This is reviewer 2's point generalized
beyond the GO terms, and it recurs: "may be important for", "play important roles during germ
cell development", "suggesting potential interactions". Reviewer 2 has already signalled they
are reading for this. A pass converting "involved in" to "associated with" and "reveals" to
"identifies" would cost nothing and remove the target.

### 9. The abstract undersells the validation

Nine markers validated by immunohistochemistry, immunofluorescence and HCR RNA-FISH in intact
tissue is the strongest thing in the paper and the answer to two of the three reviewer comments.
The abstract does not mention it. Given reviewers form their view from the abstract, one
sentence there would change how the rest is read.

### 10. Cluster number is unjustified

Nineteen clusters from Louvain at k = 30, with no statement of how resolution was chosen or
whether cluster stability was assessed. Given that several biological claims rest on specific
subcluster identities, a sentence on resolution choice, or a stability check, would close an
easy line of attack.

## What is strong, and worth protecting

- The orthogonal validation, as above. Unusual and convincing.
- Honest reporting of negatives: SOX9 barely detected, DMRT1 not annotated, steroid manipulation
  historically unsuccessful. This reads as careful rather than weak, and it supports the
  transitional GSD/ESD argument better than a positive result would.
- The chromosome 81 germline-specific analysis, which connects the atlas to the programmed DNA
  loss literature and is a genuinely novel use of the data.
- The resource value. A first gonadal single-cell atlas for a cyclostome is worth publishing
  independently of the mechanistic claims, and the framing could lean on that more confidently.
