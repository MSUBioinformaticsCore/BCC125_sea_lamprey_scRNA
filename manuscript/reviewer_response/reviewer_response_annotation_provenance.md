# Response to reviewer: annotation database, orthology mapping, and GO evidence

> ## ⚠ BLOCKED: cell-type half is unverified
>
> **The GO half of this document is fully verified against the GAF and ready to use.**
>
> **The cell-type half is not.** Every statement in this document about marker-gene provenance
> (that markers were "established directly in lamprey", that no non-lamprey marker set was
> applied, that no atlas label transfer occurred) is currently **unconfirmed**. It originates
> from an early working assumption, not from the analysis record.
>
> Do not send any of that text until the marker provenance is confirmed. If markers in fact
> came from zebrafish/mouse/human sets, the reviewer's comment is correct as written and the
> response needs to concede considerably more than it currently does.
>
> **First step is cheap:** the collaborator cited the marker source in the supplemental. Read
> that citation before escalating; it will likely say whether the markers are lamprey-derived
> or borrowed, which is the only fact currently blocking this section. Only if it is ambiguous
> does this need to go back to her.
>
> Every claim below marked **[UNVERIFIED]** depends on that answer.

## The one sentence this reduces to

Everything in the GO half of this document is an elaboration of a single claim. If the response
says nothing else, it should say this:

> An annotation on a lamprey protein is not a function copied from a counterpart gene in any
> particular species; it is the assertion that a function is conserved across a protein family,
> together with the evidence that the lamprey protein belongs to that family. Family membership
> is the whole of the evidence that the function applies in lamprey.

Note "any particular species", not "human". The experimentally characterized members behind an
ancestral node annotation may be in mouse, *Drosophila*, yeast, zebrafish or any other organism
in the family; naming human specifically is both inaccurate and needlessly anthropocentric.

Everything else here supports that sentence: the versions establish *which* family definitions
and *which* curation state, the `with/from` field establishes *which* family each gene was
assigned to, and the IEA status establishes that no lamprey member of any family has been
characterized directly.

## The request

> Cell-type labels and GO interpretations rely substantially on cross-species annotation.
> Please report the annotation database and version, orthology-mapping procedure and, where
> available, the evidence supporting the transferred GO annotations.

**This is a documentation request, not a challenge to the analysis.** The opening statement is
correct: the work does rely substantially on cross-species annotation. There is nothing to
rebut. What is being asked for is three specific things, and the response should simply supply
them.

| Asked for | Status |
|---|---|
| **Annotation database and version** | Fully answerable. Two provenance events: structural annotation and GO annotation, with different dates and versions. |
| **Orthology-mapping procedure** | Answerable for gene symbols, where a genuine orthology procedure exists (NCBI EGAP, and Wang et al. for LOC genes). For GO there is no orthology step; describe the family-assignment and tree-placement procedure instead, and say plainly that orthology is not used. |
| **Evidence supporting the annotations** | Fully answerable, and the strongest part: every annotation carries a method reference and a populated `with/from` field. |

Two notes on tone. First, do not spend the response explaining that the reviewer used
"transferred" imprecisely. Terms do move from characterized proteins in other species to
uncharacterized ones; "transferred" is reasonable shorthand, and correcting it reads as
evasion. Describe the mechanism accurately and let that stand as the clarification.

Second, the only substantive gap is the cell-type half, where the marker provenance is
genuinely unknown (see the banner above). That is the part requiring work, not argument.

### Gene symbols: two sources, both cross-species, both documentable

This is the part of the reviewer's question you *can* answer precisely, and it is where the
phrase "orthology-mapping procedure" actually applies. Symbols in your data have two origins.

**1. NCBI symbols (used where they exist).** These are not arbitrary. EGAP assigns them by an
explicit orthology procedure:

- Predicted genes are evaluated for orthology against a reference species using **pairwise
  protein alignment combined with local synteny**. Orthology is assigned **only where there is
  a clear 1:1 relationship** between loci; microsynteny is used specifically to resolve
  closely related paralogs. Where an ortholog is identified, the symbol and gene name are
  transferred from the reference species.
- Where no ortholog can be determined, the gene is named after the most similar SwissProt
  protein with a **`-like` suffix**, explicitly marking the assignment as putative.
- Where neither is possible, the gene receives `LOC` + GeneID.

Two things follow that are worth stating in the response. First, this is a synteny-aware,
1:1-constrained procedure, substantially more conservative than best-BLAST-hit naming, and
designed to avoid exactly the paralog-mislabeling failure the reviewer is likely worried
about. Second, **the `-like` suffix is a built-in confidence flag**: a symbol ending in
`-like` denotes similarity only, not orthology. Worth reporting how many of your marker genes
carry `-like` symbols; it is a concrete, quantitative answer to "what evidence supports
this?".

```bash
# how many genes in your marker/DE lists carry putative "-like" names vs. asserted orthologs
grep -c -- '-like' gene_descriptions.csv
```

**2. Wang et al. 2024 symbols (used for LOC genes otherwise).** Their procedure, from the
published Methods:

- LOC gene IDs → RefSeq protein IDs (bioDBnet); sequences retrieved with `efetch`.
- `blastp` (BLAST+) against `nr`, restricted with `-entrez_query "Homo sapiens [organism]"`.
- **The single top-aligned human RefSeq protein** was taken and converted to its human gene
  symbol. 8,984 LOC genes annotated this way, written in UPPERCASE.
- 274 genes with no significant human match were BLASTed against other species and annotated
  in lowercase.
- Genes sharing the same top protein were disambiguated by appending a suffix.
- Net effect: symbol coverage rose from 35.86% to 68.92%.

**Assess this honestly: it is unidirectional top-BLAST-hit naming.** State the limitations
yourself rather than letting the reviewer find them:

- **No reciprocal best-hit check.** A top hit in one direction does not establish orthology.
- **No synteny and no phylogenetic analysis.** This is the material contrast with the NCBI
  symbols, which require 1:1 correspondence supported by local synteny.
- **No reported E-value, identity or coverage threshold.** "Significant match" is undefined.
- **Lamprey is close to the worst case for this method.** The high GC content of lamprey
  coding sequence makes lamprey paralogs more similar to one another than to their gnathostome
  orthologs, so a top human hit is systematically vulnerable to picking the wrong family
  member. For some families (e.g. KCNA) no clear lamprey orthology can be established at all.
- **The suffix convention marks the weakest cases.** Multiple lamprey genes mapping to one
  human protein is unresolved paralogy by definition; any suffixed symbol should be treated as
  low-confidence.

None of this makes the annotation unusable; Wang et al. used it for exactly this purpose and
were transparent about it. But it is **homology, not orthology**, and it cannot bear
identity claims on its own.

### Two confidence tiers: separate them

Your symbols come from two procedures of very different strength, and you can tell them apart:

| Tier | Source | Procedure | Strength |
|---|---|---|---|
| **A** | NCBI (~35.9% of genes) | Protein alignment + local synteny, 1:1 only; `-like` suffix flags similarity-only | Reasonably strong |
| **B** | Wang et al. (fills to ~68.9%) | Unidirectional top blastp hit vs. human; no reciprocal, no synteny | Weak, suggestive only |

Reporting this split is a strong move: it shows you understand the annotation you used, and it
answers "orthology-mapping procedure" precisely for both layers. Case suffices as the flag:
Wang UPPERCASE = human-derived, lowercase = other species, suffixed = unresolved paralogy.

### The substantive fix: validate the marker genes specifically

The marker genes driving cell-type calls are probably a manageable set (tens, not thousands).
Rather than defending Tier B annotation in general, re-examine those specific genes:

1. For each marker gene underlying a cluster label, record which tier its symbol came from.
2. For Tier B markers, do a reciprocal check: blastp the lamprey protein against human, then
   the top human protein back against lamprey, and confirm it returns the original gene.
   Alternatively check whether NCBI assigns a 1:1 ortholog for that gene.
3. Report the outcome, and relabel any cluster whose identity rests on a marker that fails.

This converts the weakest part of your response into its strongest: *"we re-examined the N
marker genes underlying our cell-type assignments; X were confirmed by reciprocal best hit or
NCBI 1:1 orthology, Y are supported by sequence similarity alone, and we have qualified those
labels accordingly."* That is a concrete, verifiable answer, and it is a few hours of work.

### Where cross-species information does and does not enter

Four places it could enter, and what is true of each here:

1. **Gene identity / quantification**: kPetMar1 gene IDs throughout. No cross-species input.
2. **GO enrichment**: no cross-species input. The mapping is joined on the GAF's gene symbol
   field (`V3`), which carries **NCBI** symbols, against the NCBI annotation used for
   quantification. The join key is simply an identifier for the same gene, and the GO terms
   themselves were determined from protein sequence, so the choice of key does not affect which
   terms a gene receives. Crucially, the Wang et al. homology-derived symbols appear nowhere in
   this mapping and cannot have influenced enrichment results. Word the claim in terms of *which
   symbol vocabulary* was used rather than "identifiers not names", which is imprecise and
   invites a pointless correction.

   The one mechanical check worth running: confirm every detected gene joins, i.e. that no Wang
   symbols were substituted into the expression matrix rownames at any point.
   `sum(!expr_genes %in% gaf$Gene)` should be 0 or near it.
3. **Cluster identification / cell-type labeling**: **cross-species information did
   contribute here.** Homology-derived symbols for LOC genes (Wang et al. 2024) assisted in
   recognizing and interpreting marker genes when assigning cluster identities.
4. **Gene names in text and figures**: descriptive, from the same homology mapping.

**Do not write "no analysis step used gene symbols."** It is false for cluster identification,
where homology-derived symbols did contribute, and it is the wrong axis anyway. The
distinction that matters is not identifiers versus symbols but **which symbol vocabulary**:

> GO annotations were matched to genes using the NCBI annotation's own gene symbols, as
> distributed in the GAF. The homology-derived symbols of Wang et al. were not used in this
> mapping and could not have influenced enrichment results. Those symbols did contribute to
> cluster identification, as described below.

That is accurate, defensible against anyone reading the code, and makes the point cleanly. It
also pre-empts the objection that the mapping "used gene names", which is technically true but
irrelevant, since the terms were determined from protein sequence and the symbol serves only
as a key.

**The distinction that matters for the cell-type half:** whether homology-derived symbols were
used to *generate candidate identities that were then confirmed against lamprey-specific
marker evidence*, or whether any cluster's identity rests on the symbol alone. The first is
routine and easily defended; the second is what the reviewer is worried about. Establish which
is true for every cluster before answering, and if any cluster's label rests substantially on
homology-derived names without independent lamprey support, say so and soften that label.

Residual honesty point: lamprey orthology is genuinely hard: high GC content makes lamprey
paralogs more similar to each other than to their gnathostome orthologs, and for some families
no clear orthology can be established at all. Describe the symbols as annotations of homology,
not as orthology assertions. Do not claim they are verified orthologs.

### The GO half: do not overclaim

It is tempting to say "no cross-species inference occurred." **That claim is false and a
knowledgeable reviewer will catch it.** The GO terms in the GAF are assigned by matching
lamprey proteins to protein signatures, and then mapping those signatures to GO terms, but
the *signature→GO mappings themselves are built from experimental evidence in other species*.
The file cites exactly two methods (confirmed from column 6, see Results below):

- **InterPro2GO**: PMID:22301074 (Burge et al. 2012), corresponding to GO_REF:0000002.
  An InterPro curator assigns GO terms to each InterPro entry based on the shared
  characteristics of the matched protein set, explicitly evaluating sequence alignments and
  **the experimental literature** for those proteins, overwhelmingly from model organisms.
  Any protein matching the signature inherits the terms. 20,570 annotations (30.2%).
- **TreeGrafter**: PMID:30032202 (Tang et al. 2019), corresponding to GO_REF:0000118.
  The query protein is assigned to its best-matching PANTHER family, **grafted onto that
  family's phylogenetic tree**, and then inherits GO annotations propagated down from the
  appropriate ancestral nodes. Those ancestral annotations were themselves established from
  experimental evidence in extant members of the family, i.e. in other species.
  47,575 annotations (69.8%). Runs as part of InterProScan.

#### There is no cross-species-free step in GO annotation

Worth stating plainly, because it settles how much of the reviewer's premise to accept. Family
assignment itself is cross-species: a PANTHER family HMM is a statistical model built from an
alignment of that family's members across reference genomes, and the pre-computed MSA and tree
topology are likewise built from other species' sequences. Assigning a lamprey protein to a
family asserts that its sequence resembles an aggregate model of other species' proteins.

| Step | Source of information |
|---|---|
| Family model the protein is scored against | Sequences from other species |
| Tree it is grafted onto | Sequences from other species |
| Ancestral node annotations inherited | Experiments in other species |
| InterPro2GO entry-to-term mappings | Literature on other species |

The only lamprey-specific inputs are the predicted protein sequence and the expression data.
Note that the family assignment is not a third lamprey-specific input: it is the *output* of
comparing the lamprey sequence against cross-species models, so describing it as lamprey-derived
smuggles the cross-species step back in as though it were native. Even the protein sequence is
not wholly independent, since Gnomon uses protein alignments from other species as evidence
when building gene models, though lamprey genomic sequence and RNA-seq evidence dominate.

Do not look for a step to characterize as species-independent; there isn't one. What remains
precisely true, and worth stating, is that the *unit* of inference is a protein family or a
phylogenetic node, never a designated 1:1 ortholog, and that no annotation was copied from a
named gene in another species.

Two lamprey-specific consequences: a divergent protein may fail to match any family model,
contributing to the 2,721 protein-coding genes with no terms; and it may match the wrong
family, which the last-common-ancestor fallback mitigates but does not eliminate.

**Worth checking:** whether *Petromyzon marinus* is among PANTHER 17.0's reference genomes. If
it is not, the models your proteins are scored against contain no lamprey sequence whatsoever,
which is a sharper and more honest version of this point.

#### PANTHER appears in both routes. They are not the same thing.

This is the easiest part of the provenance to garble, so state it carefully.

| | InterPro2GO (20,570 annots) | TreeGrafter (47,575 annots) |
|---|---|---|
| PANTHER's role | One of 11 signature databases feeding InterPro entries | Source of the phylogenetic trees and PAINT node annotations |
| Mechanism | Flat lookup in the `interpro2go` file | Placement of the query into a family tree |
| Resolution | Every protein matching an entry gets identical terms | Terms depend on where in the tree the protein falls |
| Can it separate paralogs? | No | Yes, in principle |

The same PANTHER family definitions underlie both. What differs is how GO terms are extracted:
a static entry-to-term table in one case, phylogenetic position in the other.

**What `interpro2go` actually is.** A plain text file, one line per entry-to-term pair, e.g.:

```
InterPro:IPR000023 Phosphofructokinase domain > GO:6-phosphofructokinase activity ; GO:0003872
InterPro:IPR000023 Phosphofructokinase domain > GO:glycolytic process ; GO:0006096
InterPro:IPR000025 Melatonin receptor family > GO:melatonin receptor activity ; GO:0008502
```

InterProScan reports which entries a protein matched; the pipeline looks those accessions up
in this file and emits the corresponding terms. The curation is in building the table, not in
applying it.

**The binding constraint: a term must be true of *every* protein the entry matches.** Because
the mapping is applied automatically to the entire matched set, often thousands of proteins
across all sequenced organisms, a curator can only assign terms they will assert for all of
them, including members nobody has studied. It is the intersection of what is true across the
family, not the union. Three consequences:

- **Specificity is capped by entry breadth.** `IPR000023 Phosphofructokinase domain` can carry
  `6-phosphofructokinase activity` because every member catalyses that reaction, but could not
  carry a tissue- or context-specific process term. A narrow entry such as `IPR000025 Melatonin
  receptor family` can carry `melatonin receptor activity`; a broad class-A GPCR entry could
  not.
- **Heterogeneous entries get nothing.** Where the matched proteins are too functionally
  diverse for any term to hold across them, the entry has no GO mapping. This is common for
  homologous superfamily entries and contributes to the 2,721 protein-coding genes in your
  annotation with no terms.
- **The judgment comes from characterized members**, overwhelmingly human, mouse, *Drosophila*
  and yeast. Extending it to a lamprey protein is an assumption of functional conservation
  across the family. This is the accurate form of the reviewer's cross-species concern: not
  that terms were copied from a human ortholog, but that each term asserts family-wide
  conservation, with family membership as the sole evidence that it holds in lamprey. This is
  the accurate description of the evidence, and it is what the reviewer asked for.

**Signature scope matters and is recorded.** InterPro types each entry as Family, Domain,
Homologous Superfamily, Repeat, or Site. PANTHER signatures are family-level, matching whole
proteins rather than domains. Since PANTHER supplies 48,431 of your `with/from` entries against
InterPro's 31,670, the majority of your annotations rest on whole-protein family membership,
not on domain composition. If asked whether these are "just domain predictions", the answer is
no.

#### How TreeGrafter actually places a sequence (from Tang et al. 2019, Methods 2.2)

1. Query scored against the PANTHER **family** HMM library with `hmmscan` (HMMER3) to select
   the best-matching tree and produce a query-to-HMM alignment.
2. That alignment is used to insert the query into the family's pre-computed multiple sequence
   alignment (chosen over re-running MAFFT, which was the bottleneck).
3. **RAxML v8 in parsimony mode** grafts the query onto the reference tree. Where several graft
   points are equally parsimonious, the **last common ancestor** is used as consensus.
4. Annotations are traced root-to-graft-point: ancestral **gain**-of-function annotations are
   inherited, ancestral **loss**-of-function annotations block propagation.

Points to use in the response:

- **The LCA fallback is conservative.** An ambiguously-placed divergent sequence grafts deeper
  and inherits more general terms rather than a confidently wrong specific one. This is a
  direct answer to "how reliable is this for an organism as divergent as lamprey".
- **Placement uses parsimony, not maximum likelihood**, traded for speed. Divergent sequences
  are where this is weakest. Do not oversell placement accuracy.
- **Coverage is gated by PAINT curation.** At publication only ~4,650 of 14,710 PANTHER
  families had curated GO annotations. Genes in uncurated families receive nothing by this
  route, which partly explains the 2,721 protein-coding genes with no terms in your file.
- **TreeGrafter terms are more specific, not less.** The paper reports more annotations per
  protein and more specific terms than InterPro2GO where the two are comparable. Do not
  describe the 70% TreeGrafter fraction as shallow or generic; that characterization applies
  to the InterPro2GO fraction.
- **Validation taxa were horse, anole, mosquito, poplar, fungus, archaeon, bacterium.** Broad,
  but nothing near cyclostomes. If pressed on accuracy in lamprey, concede that no
  lamprey-proximate benchmark exists rather than citing their figures as if they transfer.

**Correction, verified against the GO_REF:0000118 page:** PAN-GO is not irrelevant here, it is
*upstream*. TreeGrafter propagates from ancestral nodes "which were manually annotated using
the PAN-GO method (see GO_REF:0000033)". So the ancestral annotations your lamprey genes
inherit were curator-assigned; what is automated is the placement of the lamprey protein into
the tree and the propagation step. This is better for you than the flat "all automated"
framing: the functional assignments at the nodes are manually curated, and the electronic step
is the grafting. State it that way.

Note this is TreeGrafter, **not** PAINT/PAN-GO (GO_REF:0000033) as the direct reference. The distinction matters if
challenged: PAINT annotations are curator-reviewed and carry evidence code IBA, whereas
TreeGrafter performs the tree placement and propagation automatically, hence IEA. The
underlying tree annotations derive from the same curated PANTHER resource, but no curator
reviewed the lamprey assignment itself. Do not cite GO_REF:0000033; it is not in your file.

So the honest statement is: **no ortholog is identified for any lamprey gene, but the
functional knowledge being applied is derived from experimental work in other species,
propagated at the level of protein family or phylogenetic node rather than ortholog pair.**

This is a *better* rebuttal than a flat denial. It answers the question the reviewer actually
asked ("what evidence supports these?") with a real, traceable answer: the reference field
identifies the method for every annotation and the `with/from` field is **100% populated**
with the specific PANTHER/InterPro/UniProtKB entry the term came from. Reframe from "you are
mistaken" to "here is the precise mechanism, and here is the complete evidence trail."

The one thing you genuinely owe them either way: these annotations are **entirely IEA
(68,145/68,145)**: no lamprey gene in this set has been experimentally characterized. Say
so up front rather than letting them discover it.

---

## Numbers to fill in before sending

Run these against your GAF and paste the real values into the placeholders below.

```bash
GAF=GCF_010993605.1_kPetMar1.pri_gene_ontology.gaf

# 1. Header. This file has BOTH "!" metadata lines and a "#" column-name line,
#    so filter on '^[!#]' for the data rows (done below).
grep '^!' $GAF

# 2. Evidence-code distribution (col 7)
cat $GAF | grep -v '^[!#]' | cut -f7 | sort | uniq -c | sort -rn

# 3. GO_REF / reference IDs (col 6): identifies the assignment method
cat $GAF | grep -v '^[!#]' | cut -f6 | sort | uniq -c | sort -rn

# 4. Assigned-by (col 15)
cat $GAF | grep -v '^[!#]' | cut -f15 | sort | uniq -c | sort -rn

# 4b. with/from (col 8): THE EVIDENCE TRAIL the reviewer is asking for.
#     Expect InterPro IDs (IPR...) and/or PANTHER node IDs (PTN...). Report which.
cat $GAF | grep -v '^[!#]' | cut -f8 | grep -v '^$' | tr '|' '\n' \
  | sed 's/:.*//' | sort | uniq -c | sort -rn | head
# fraction of annotations that carry any with/from support at all
cat $GAF | grep -v '^[!#]' | awk -F'\t' '{n++; if($8!="") w++} END {print w"/"n, w/n}'

# 5. Coverage: annotated genes, and split by aspect (P/F/C)
cat $GAF | grep -v '^[!#]' | cut -f3 | sort -u | wc -l
cat $GAF | grep -v '^[!#]' | cut -f9 | sort | uniq -c

# 6. How many of YOUR detected/expressed genes have >=1 GO term
#    (this is the number that defines your enrichment universe)
```

```r
# R version of #6: adjust object names.
# This file has "!" metadata lines AND a "#" column-name line. readr's `comment`
# takes a single string, so filter both out explicitly rather than relying on it.
lines <- readr::read_lines(GAF)
gaf <- readr::read_tsv(I(lines[!grepl("^[!#]", lines)]), col_names = c(
  "DB","GeneID","Symbol","Qualifier","GO_ID","Reference","Evidence_Code",
  "With_From","Aspect","Gene_Name","Gene_Synonym","Type","Taxon","Date",
  "Assigned_By","Annot_Ext","Gene_Product_Form_ID"))
annotated <- unique(gaf$GeneID)   # GeneID (col 2), not Symbol
length(intersect(expressed_genes, annotated)) / length(expressed_genes)
```

### Results from this file

Everything below is confirmed from the file itself; no placeholders remain except your own
dataset's gene counts.

**From the header:**

```
!gaf-version: 2.2
!generated-by: NCBI
!date-generated: 2023-12-11
!interProScan version: 5.64-96.0
!PANTHER version: 17.0
!GO version: 2023-07-27
```

| Field | Value |
|---|---|
| GAF generated | **2023-12-11** by NCBI (GAF format 2.2) |
| InterProScan | **5.64-96.0** |
| PANTHER | **17.0** |
| GO ontology release | **2023-07-27** |
| Total annotations | 68,145 |
| Evidence code (col 7) | **IEA, 68,145 (100%)**. No experimental, no manually curated. |
| Reference (col 6) | PMID:30032202 (TreeGrafter) 47,575 · PMID:22301074 (InterPro2GO) 20,570 |
| Assigned_By (col 15) | RefSeq (100%) |
| With/From (col 8) | PANTHER 48,431 · InterPro 31,670 · UniProtKB 1,893; **68,145/68,145 (100%) populated** |
| Distinct genes annotated | 14,859 (~85% of the 17,580 protein-coding genes in Release 100) |
| Aspect (col 9) | F 25,246 · P 24,269 · C 18,630 |

**The date gap is worth stating explicitly.** The structural annotation is Release 100
(submitted 2020-03-20, pipeline v8.3); the GO annotation was generated 2023-12-11 against GO
release 2023-07-27. These are two different provenance events on the same assembly, and
reporting a single "annotation version" would conflate them. Report both.

Two things this confirms, and one still open:

- **Confirmed:** the file cites no GO_REF IDs at all; it cites the two method papers
  directly. Cite the PMIDs. Mentioning the equivalent GO_REF (0000002 / 0000118) alongside is
  fine and helpful, but the PMIDs are what is actually in the data.
- **Confirmed:** every annotation carries a `with/from` entry. This is the strongest single
  fact in your response: the reviewer asked for "evidence supporting the transferred
  annotations" and you can point to a complete, per-annotation source record.
- **Still open:** whether the PANTHER `with/from` entries are ancestral node IDs (`PTN…`) or
  family IDs (`PTHR…`). TreeGrafter should emit node IDs, which are more specific and better
  for you. Check before writing the sentence about node-level traceability:

```bash
cat $GAF | grep -v '^[!#]' | cut -f8 | tr '|' '\n' | grep '^PANTHER' \
  | sed 's/PANTHER:\(....\).*/\1/' | sort | uniq -c
```

**Confirmed:** RefSeq *Petromyzon marinus* Annotation Release 100 was submitted 20 Mar 2020
(pipeline software v8.3), but the GAF header shows the GO annotation was generated
2023-12-11, consistent with GO terms having been added to the Eukaryotic Genome Annotation
Pipeline output in late 2023. Use the header date and versions for the GO provenance, and the
Release 100 date for the structural annotation. Do not merge them.

---

## Draft methods text

> **Genome and gene annotation.** Reads were quantified against the *Petromyzon marinus*
> reference assembly kPetMar1.pri (GenBank/RefSeq accession GCF_010993605.1, Vertebrate
> Genomes Project, assembly date 2020-02-24), using the structural gene annotation from NCBI
> *Petromyzon marinus* Annotation Release 100 (submitted 2020-03-20; NCBI Eukaryotic Genome
> Annotation Pipeline software version 8.3), which comprises 22,167 genes and pseudogenes,
> of which 17,580 are protein-coding.
>
> **Gene Ontology annotation.** GO annotations were taken directly from the GO Annotation
> File distributed with that assembly
> (`GCF_010993605.1_kPetMar1.pri_gene_ontology.gaf`; GAF format 2.2, generated by NCBI on
> 2023-12-11 using InterProScan v5.64-96.0 and PANTHER v17.0, against Gene Ontology release
> 2023-07-27), comprising 68,145 annotations across 14,859 genes (25,246 molecular function,
> 24,269 biological process, 18,630 cellular component). We note that the GO annotation
> (2023-12-11) postdates the structural annotation of this assembly (Release 100, submitted
> 2020-03-20) and was generated in a separate pipeline pass. These annotations are generated
> by NCBI, which runs InterProScan on the annotated lamprey proteins, comparing each against a
> library of statistical models of known protein families and domains, and assigns GO terms by
> two methods: TreeGrafter, which assigns each protein to its best-matching PANTHER family,
> grafts it onto that family's phylogenetic tree and propagates GO annotations from the
> appropriate ancestral nodes (47,575 annotations, 69.8%; Tang et al. 2019, PMID:30032202);
> and InterPro2GO, in which curator-assigned mappings between InterPro entries and GO terms
> are applied to any protein matching the corresponding entry (20,570 annotations, 30.2%;
> Burge et al. 2012, PMID:22301074). Both references are those cited in the GAF itself. All
> annotations are attributed to RefSeq.
>
> **Orthology mapping.** No orthology-mapping step is involved in producing these GO
> annotations: no ortholog is identified for any lamprey gene, and no annotation is copied from
> a designated ortholog in a reference species. Cross-species information is nonetheless
> present throughout the procedure. The family models against which each lamprey protein is
> scored, and the phylogenetic trees onto which it is grafted, are constructed from the
> sequences of other species; the ancestral node annotations inherited by a grafted protein
> were established from experimental evidence in extant members of the family, predominantly in
> model organisms; and the InterPro entry-to-term mappings are assigned by curators from the
> experimental literature on the proteins matching each entry. The only lamprey-specific input
> is the predicted protein sequence itself; the assignment of that sequence to a family is the
> result of comparing it against these cross-species models. What distinguishes
> this procedure from ortholog-based transfer is therefore the unit at which functional
> knowledge is applied, namely a protein family or phylogenetic node, rather than a designated
> counterpart gene in a particular reference species.
>
> Accordingly, all 68,145 annotations carry the evidence code IEA (Inferred from Electronic
> Annotation) with respect to *P. marinus*; none is manually curated for this species. The
> supporting inference for each annotation is recorded in the `with/from` field, which is
> populated for 100% of annotations (PANTHER, 48,431; InterPro, 31,670; UniProtKB, 1,893),
> allowing every term to be traced to the specific family, node or entry from which it
> derives. Of the [N] genes detected in our dataset, [N] ([X]%) carried at least one GO term;
> enrichment tests used this annotated subset as the background universe rather than all
> detected genes.
>
> **Cell-type assignment.** Cluster identities were assigned using marker genes established
> in the lamprey literature ([CITATIONS]) and, where available, published lamprey in situ
> hybridization data ([CITATIONS]). Marker genes were matched to annotated lamprey gene models
> by [gene symbol / RefSeq accession, state which]. No marker sets from other species were
> used, and no computational label transfer from a non-lamprey reference atlas was performed.

---

## Second reviewer point: mammalian-sounding GO terms and provisional cell types

> Enrichment for terms such as "binding of sperm to zona pellucida" and "positive regulation of
> acrosome reaction" does not by itself demonstrate that the identified lamprey genes perform
> the same roles characterized in other vertebrates. Likewise, the specific pre-granulosa,
> granulosa and follicular subtype assignments should be presented as marker-supported,
> provisional transcriptomic identities unless independently validated.

### What is actually being asked

This is an **epistemic** objection, not a claim that the terms are inapplicable. The reviewer
is not saying the enrichment is spurious; they are saying enrichment alone cannot
*demonstrate* the genes perform those roles. That is correct and cheap to concede.

Note also that the second sentence hands you the exact wording they will accept:
"marker-supported, provisional transcriptomic identities". Use their phrase. When a reviewer
specifies an acceptable claim level, adopting it verbatim closes the point.

### Do not simply concede the terms are meaningless in lamprey

Both named processes have real lamprey biology behind them. Check these against the primary
literature before using them, but the picture from a literature search is:

- **ZP-domain proteins are present in the lamprey egg envelope.** Egg coats across metazoans
  (zona pellucida, vitelline envelope, chorion) are built from ZP-module proteins. Sea lamprey
  ZP proteins show homology to the ZP1/ZP2 subtypes of jawed vertebrates. So the gene family is
  genuinely present and genuinely functions in the egg coat; what is mammal-specific is the
  *name* "zona pellucida", not the underlying biology.
- **Lamprey sperm have an acrosome and undergo an acrosome reaction**, described as exocytosis
  of the acrosomal vesicle and extrusion of an acrosomal filament, with the sperm digesting its
  way through the egg envelopes. Sea lamprey egg morphology and fertilization events were
  described by Kobayashi and Yamamoto (1994). So this term is not merely a mammalian import.

This matters: the honest response is not "these terms are artefacts of mammalian
nomenclature." It is that the terms name processes with genuine lamprey counterparts, that GO
vocabulary is inherited from the model systems where the processes were characterized, and
that enrichment nonetheless does not establish that these particular lamprey genes carry out
these particular roles. That is a stronger and more accurate position than conceding the terms
are meaningless.

### Can GO enrichment support the cell-type assignments?

**Effectively no, and the reason is non-independence rather than invalidity.** The GO
annotations come from cross-species protein family assignment. The cell-type labels come from
markers recognized through cross-species gene symbols. Both trace back to the same underlying
source. Citing GO enrichment as corroboration of a cell-type call restates the original
evidence in different vocabulary rather than adding to it.

Keep the enrichment, but change what it is claimed to show:

| Do not claim | Claim instead |
|---|---|
| Enrichment confirms this cluster is granulosa | This cluster's transcriptional program is enriched for gene families associated with follicular somatic cell function in other vertebrates |
| These genes mediate sperm binding | These genes belong to families that mediate sperm binding where characterized; their function in lamprey is untested |
| GO supports the identity assignment | Identity is marker-supported and provisional; GO describes the program, not the identity |

**Independent validation would be:** in situ hybridization or immunostaining showing the marker
localized to the expected cell population in lamprey gonad; histological correlation of
clusters with morphologically identified cell types; or functional data. If you have any of
this, it belongs in the response. If not, adopt the provisional framing and say so.

### Concrete strengthening move

For each flagged term, report the genes driving it and their family assignment. Something like:
"the enrichment of *binding of sperm to zona pellucida* is driven by N genes of the ZP-module
family, which forms the egg envelope across chordates including lamprey; the GO term name
reflects the mammalian structure in which the function was characterized." That converts a
suspicious-looking term into a demonstration that you understand the annotation, and it directly
answers what the reviewer is probing.

```r
# genes behind a flagged term
gaf |> dplyr::filter(GO_ID == "GO:0007339") |> dplyr::select(GeneID, Symbol, With_From)
```

---

## Draft point-by-point response

> **Reviewer comment:** *Cell-type labels and GO interpretations rely substantially on
> cross-species annotation. Please report the annotation database and version,
> orthology-mapping procedure and, where available, the evidence supporting the transferred
> GO annotations.*
>
> We agree, and we should have documented this in the original submission. Our work does rely
> substantially on cross-species annotation, and the Methods did not say where it enters or on
> what evidence it rests. We have revised accordingly, and set out below the annotation
> sources and versions, the mapping procedures used at each step, and the evidence available
> for each class of annotation. Cross-species information enters at three points, which differ
> in kind and in the weight they bear: gene nomenclature, cell-type labeling, and GO
> annotation.
>
> **Cell-type labels.** No computational label transfer from a non-lamprey reference atlas was
> performed, and no marker set defined in zebrafish, mouse or human was applied to our data;
> cluster identities rest on markers established directly in lamprey ([CITATIONS]). The
> reviewer is nonetheless correct that cross-species information contributed, and we have
> made this explicit in the revised Methods. Because a majority of genes in the kPetMar1
> annotation carry only LOC identifiers, we used the homology-based gene symbols of Wang et
> al. 2024 to recognize candidate marker genes when interpreting clusters. Candidate
> identities arrived at in this way were then [confirmed against lamprey-specific marker
> evidence, STATE HOW]. We have added the full marker list, its supporting citations, and the
> gene ID to symbol mapping as Supplementary Tables [N].
>
> **Gene identifiers and gene names.** GO annotations were matched to genes using the gene
> symbols of the NCBI annotation itself, as distributed in the GAF, against the same NCBI
> annotation used for quantification. The homology-derived symbols of Wang et al. 2024 were not
> used in this mapping, and could not have influenced the enrichment results; the GO terms
> associated with each gene were determined from its predicted protein sequence, independently
> of any symbol assigned to it. Gene symbols did, however, contribute to cluster
> identification, and we describe their provenance and its limitations here.
>
> Symbols derive from two sources. Where the NCBI annotation assigns a symbol, we use it; these
> are assigned by the NCBI Eukaryotic Genome Annotation Pipeline on the basis of pairwise
> protein alignment together with local synteny, and only where a clear 1:1 relationship
> between loci can be established, with microsynteny used to resolve closely related paralogs.
> Genes named from SwissProt similarity alone carry a `-like` suffix denoting a putative
> assignment. Where no NCBI symbol exists, we use the annotation of Wang et al. 2024, who
> assigned symbols to 8,984 LOC-numbered genes by retrieving each predicted protein and taking
> the single top-scoring `blastp` alignment against human sequences in `nr`, together with 274
> further genes matched against other species; genes sharing a top-matching protein were
> disambiguated by suffix. This raised symbol coverage of the reference transcriptome from
> 35.86% to 68.92%.
>
> We note explicitly that the latter procedure establishes **homology, not orthology**: it is
> unidirectional, without a reciprocal best-hit requirement, and without synteny or
> phylogenetic corroboration. This limitation is of particular consequence in lamprey, where
> the compositional properties of coding sequence can make paralogs more similar to one another
> than to their gnathostome orthologs, and where orthology cannot be confidently resolved for
> some gene families. Symbols carrying a disambiguating suffix denote unresolved
> one-to-many correspondence and should be regarded as provisional. Accordingly we treat these
> symbols as descriptive rather than as assertions of gene identity, and [describe the
> reciprocal validation of marker genes and its outcome here, see checklist]. The full mapping
> from kPetMar1 gene IDs to symbols and their source is provided as Supplementary Table [N]
> (`gene_descriptions.csv`).
>
> **GO annotations.** These derive from the GO Annotation File released with RefSeq assembly
> GCF_010993605.1 (*Petromyzon marinus* Annotation Release 100), generated by NCBI on
> 2023-12-11 using InterProScan v5.64-96.0 and PANTHER v17.0 against GO release 2023-07-27,
> and comprising 68,145 annotations over 14,859 genes. NCBI generates this file by
> running InterProScan on the annotated lamprey proteins and assigning GO terms by two
> distinct routes, both of which draw on PANTHER protein family definitions but differ in how
> terms are obtained. Under the InterPro2GO route, InterProScan reports which InterPro entries
> a protein matches, and the GO terms curators have mapped to those entries are applied
> directly; every protein matching a given entry therefore receives the same terms. Under the
> TreeGrafter route, the protein is instead placed within the phylogenetic tree of its
> best-matching PANTHER family, and inherits terms from the annotated ancestral nodes above its
> position, so that terms depend on where within the family the protein falls.
>
> We agree with the reviewer that cross-species information underlies these annotations, and we
> describe below both where it enters and what form it takes. No orthology-mapping step is
> involved: no ortholog is identified for any lamprey gene, and no annotation is copied from a
> designated ortholog in a reference species.
>
> Cross-species information is nonetheless present at every stage of the procedure, and we do
> not wish to imply otherwise. The family models against which each lamprey protein is scored,
> and the phylogenetic trees onto which it is grafted, are constructed from the sequences of
> other species. The ancestral node annotations a grafted protein inherits were established
> from experimental evidence in extant members of the family, predominantly in model organisms
> (Tang et al. 2019, PMID:30032202). The InterPro entry-to-term mappings are assigned by
> curators on the basis of the experimental literature for the proteins matching each entry
> (Burge et al. 2012, PMID:22301074). The only lamprey-specific input is the predicted protein
> sequence itself; the assignment of that sequence to a family is the result of comparing it
> against these cross-species models, not an independent lamprey-derived observation.
>
> What distinguishes this from ortholog-based transfer is therefore not the absence of
> cross-species information but the unit at which it is applied: functional knowledge is
> inherited at the level of a protein family or a phylogenetic node, rather than from a
> designated counterpart gene in a particular reference species. Terms are inherited from the nearest annotated ancestral node(s) above the graft
> point rather than from the family as a whole, with curated loss-of-function events blocking
> propagation. The specificity of an annotation therefore depends on where the sequence is
> placed: where placement is ambiguous, the last common ancestor of the candidate positions is
> used, yielding more general terms rather than confidently specific ones.
>
> **Supporting evidence.** The evidence trail is recorded in the GAF itself. The reference
> field identifies the assignment method for every annotation: TreeGrafter (PMID:30032202)
> for 47,575 annotations and InterPro2GO (PMID:22301074) for 20,570. The `with/from`
> field, which is populated for 100% of annotations (68,145/68,145), records the specific
> PANTHER entry (48,431), InterPro entry (31,670) or UniProtKB accession (1,893) from which
> each term was derived. Every annotation used in this study is therefore individually
> traceable to its source. We have deposited the exact GAF used, together with the per-gene
> annotation table, as Supplementary Data [N].
>
> **What these annotations assert.** It may help to state what an IEA annotation of this kind
> claims. A GO term mapped to an InterPro entry must be valid for every protein that entry
> matches, since it is applied automatically to all of them; curators therefore assign only
> terms they are prepared to assert across the whole matched set, which is why term specificity
> is bounded by the breadth of the entry and why many entries carry no GO mapping at all. The
> equivalent constraint in the phylogenetic route is that a term is asserted at an ancestral
> node and inherited by its descendants unless a curated loss event intervenes. In both cases
> the judgment is made from whichever members of the family have been experimentally
> characterized, predominantly in model organisms such as human, mouse, *Drosophila* and yeast.
> An annotation on a lamprey protein is therefore not a function copied from a counterpart gene
> in any particular species; it is the assertion that a function is conserved across a protein
> family, together with the evidence that the lamprey protein belongs to that family. Family
> membership is the whole of the evidence that the function applies in lamprey.
>
> **On evidence quality.** We now state the principal limitation explicitly. With respect to
> *Petromyzon marinus*, all 68,145 annotations carry the evidence code IEA: no lamprey gene in
> this set has been experimentally characterized, and none has been manually curated for this
> species. We note in particular that the PANTHER-derived annotations come from TreeGrafter,
> which places sequences onto pre-annotated family trees automatically; unlike curator-reviewed
> phylogenetic annotations (evidence code IBA), no curator has assessed the placement of any
> lamprey protein. Annotation of this kind is also biased toward well-characterized protein
> families, so genes that cannot be assigned to a recognized family receive no terms and
> are absent from the annotated set; 14,859 of 17,580 protein-coding genes (84.5%) carry at
> least one term. We have added the following to the Discussion/Limitations:
>
> > "GO annotations for *P. marinus* are entirely electronic (IEA), derived from assignment of
> > proteins to known protein families rather than from experimental characterization of
> > lamprey genes. The
> > functional assignments they encode originate from experimental work in other, mostly
> > model, species and are inherited at the level of protein family or phylogenetic node. Each
> > annotation asserts that a function is conserved throughout a protein family and that the
> > lamprey protein belongs to that family; family membership is the entirety of the evidence
> > that the function applies in lamprey, and the assertion of conservation spans the
> > approximately 500 million years separating cyclostomes from the species in which these
> > functions were characterized. Enrichment results should therefore be read as hypotheses
> > grounded in protein family membership rather than as evidence of validated function in
> > lamprey, and are expected to be conservatively biased toward conserved, well-characterized
> > protein families. Functions that have diverged in the cyclostome lineage, and
> > lamprey-specific genes that cannot be assigned to a recognized family, are systematically
> > invisible to this analysis."
>
> To ensure our enrichment tests were not additionally biased by incomplete annotation, we
> restricted the background universe to detected genes carrying at least one GO term ([N] of
> [N]; [X]%) rather than to all detected genes. [If you did not do this, do it now and report
> whether the conclusions change; this is the single most likely follow-up request.]

---

## Before you send

- [ ] Remaining placeholders are only `[N]`/`[X]` for **your own dataset's** gene counts and
      `[CITATIONS]` for the lamprey marker papers. All annotation-side values are now filled
      from the file.
- [x] ~~Check col 6 before naming GO_REF IDs.~~ Done: the file cites **PMIDs, not GO_REFs**:
      PMID:30032202 (TreeGrafter) and PMID:22301074 (InterPro2GO). Cite the PMIDs. The
      equivalent GO_REFs are 0000118 and 0000002; **GO_REF:0000033 (PAINT) is not used here**
      and must not be cited.
- [ ] Do not write "no cross-species inference was performed." Cross-species experimental
      evidence underlies these mappings; only the *ortholog-pairing step* is absent. Claiming
      otherwise is falsifiable in one click and would undermine the rest of the response.
- [ ] Report the two provenance dates separately: structural annotation 2020-03-20 (Release
      100, pipeline v8.3), GO annotation 2023-12-11 (InterProScan 5.64-96.0, PANTHER 17.0, GO
      release 2023-07-27). Conflating them into one "annotation version" is exactly the kind
      of imprecision this reviewer is probing for.
- [ ] Verify the GO-enrichment claim, since the response now rests on it: confirm the join used
      GAF column 2 (GeneID), not column 3 (Symbol). This is the one clean "no cross-species
      input" statement in the response; do not make it without checking.
- [ ] **Go cluster by cluster.** For each, record whether its identity is supported by
      lamprey-specific marker evidence independent of the homology-derived symbol, or whether
      the symbol is doing the work. Any cluster in the second category should be labeled more
      tentatively in the revision. Doing this yourself is far better than having the reviewer
      ask which clusters fall into which category.
- [ ] **Tier every marker gene's symbol** (NCBI vs. Wang) and reciprocally validate the Tier B
      ones. This is the highest-value item in this document: it is bounded work on a small gene
      set, and it converts the weakest part of the response into a concrete result. Flag any
      marker with a suffixed or lowercase Wang symbol as provisional.
- [ ] Count `-like` symbols among your marker genes and report the number. It is a free,
      quantitative confidence statement using NCBI's own flag.
- [ ] Sanity-check `gene_descriptions.csv` for many-to-one symbol assignments. Even though
      nothing computes on them, duplicated symbols in figures or text will read as an error:
      `dups <- names(which(table(gd$symbol) > 1))`
- [ ] Optional but strong: note that GO release 2023-07-27 is now ~2.5 years old. If you want
      to pre-empt a follow-up, state that you used the annotation as distributed with the
      assembly for reproducibility rather than re-running InterProScan against a current GO
      release. That is a defensible choice; just make it explicitly rather than by omission.
- [ ] Confirm the GAF header actually names an InterProScan version. If it does not, soften to
      "InterProScan-based protein signature analysis as implemented in the NCBI Eukaryotic
      Genome Annotation Pipeline" and cite the NCBI Insights announcement rather than naming a
      version you cannot verify.
- [ ] Confirm your enrichment background universe. If it was all detected genes rather than
      all *annotated* detected genes, re-run before responding; do not promise a re-analysis
      you have not done.
- [ ] Assemble the lamprey marker-gene table with citations as a supplementary table. This is
      the concrete artifact that resolves the cell-type half of the comment.

## References

- NCBI *Petromyzon marinus* Annotation Release 100: https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Petromyzon_marinus/100/
- Assembly GCF_010993605.1 (kPetMar1.pri): https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_010993605.1/
- "Gene Ontology Terms for NCBI RefSeq Eukaryotic Genomes", NCBI Insights, 15 Nov 2023: https://ncbiinsights.ncbi.nlm.nih.gov/2023/11/15/gene-ontology-terms-refseq-eukaryotic-genomes/
- GO evidence code IEA: https://wiki.geneontology.org/index.php/Inferred_from_Electronic_Annotation_(IEA)
- GO reference (GO_REF) collection: https://geneontology.org/gorefs.html
**Cited in the GAF itself (cite these):**

- **PMID:30032202**, Tang, Li, Thomas & Mi (2019), "TreeGrafter: phylogenetic tree-based annotation of proteins with Gene Ontology terms and other annotations", *Bioinformatics*: https://pubmed.ncbi.nlm.nih.gov/30032202/
- **PMID:22301074**, Burge et al. (2012), "Manual GO annotation of predictive protein signatures: the InterPro approach to GO curation", *Database*: https://academic.oup.com/database/article/doi/10.1093/database/bar068/430535

**Equivalent GO references (optional supporting citations):**

- GO_REF:0000118, TreeGrafter-generated GO annotations: https://geneontology.org/GO_REF/0000118
- GO_REF:0000002, InterPro2GO: https://geneontology.org/GO_REF/0000002.html
- InterPro2GO mapping file: https://www.ebi.ac.uk/GOA/InterPro2GO

**Not applicable to this file, do not cite:**

- GO_REF:0000033 (PAINT/PAN-GO) is **not** cited in your GAF's reference column, but it is the
  upstream provenance of the PANTHER ancestral-node annotations that TreeGrafter propagates.
  Cite it as background if useful; do not present it as the reference for your annotations.

Superseded note, retained so the reasoning is auditable: an earlier draft of this document
said GO_REF:0000033 was simply inapplicable. That was too strong, see the correction above.
Superficially similar to
  TreeGrafter but curator-reviewed and coded IBA; absent from this GAF.
