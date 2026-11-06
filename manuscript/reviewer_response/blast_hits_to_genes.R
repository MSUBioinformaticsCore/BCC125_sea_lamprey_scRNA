# Map BLAST hit coordinates to annotated genes, without bedtools.
# Companion to dmrt1_search.sh step 3. Works for any tblastn output against kPetMar1.

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

## ---- inputs -----------------------------------------------------------------
hits_file <- "dmrt1_vs_genome.tsv"   # tblastn outfmt 6 with sstart send as cols 7,8
gff_file  <- "GCF_010993605.1_kPetMar1.pri_genomic.gff.gz"
flank     <- 5000                    # also report genes within this distance

## ---- read hits --------------------------------------------------------------
hits <- read.delim(hits_file, header = FALSE,
                   col.names = c("query","seqname","pident","length",
                                 "evalue","bitscore","sstart","send"))

# tblastn gives sstart > send on the minus strand
hits <- hits %>%
  mutate(start  = pmin(sstart, send),
         end    = pmax(sstart, send),
         strand = ifelse(sstart < send, "+", "-"))

hits_gr <- GRanges(hits$seqname, IRanges(hits$start, hits$end),
                   strand = hits$strand)
mcols(hits_gr) <- hits[, c("query","pident","length","evalue","bitscore")]

## ---- read gene models -------------------------------------------------------
gff   <- import(gff_file)
genes <- gff[gff$type == "gene"]

## ---- overlaps ---------------------------------------------------------------
# Ignore strand: tblastn frame does not have to match the annotated gene strand.
ov <- findOverlaps(hits_gr, genes, ignore.strand = TRUE)

overlapping <- data.frame(
  query    = hits_gr$query[queryHits(ov)],
  evalue   = hits_gr$evalue[queryHits(ov)],
  pident   = hits_gr$pident[queryHits(ov)],
  seqname  = as.character(seqnames(hits_gr))[queryHits(ov)],
  hit_start= start(hits_gr)[queryHits(ov)],
  hit_end  = end(hits_gr)[queryHits(ov)],
  gene     = genes$gene[subjectHits(ov)],
  gene_id  = genes$Dbxref[subjectHits(ov)] |> sapply(function(x) paste(x, collapse = ";")),
  gene_biotype = genes$gene_biotype[subjectHits(ov)]
) %>% arrange(evalue)

print(overlapping, row.names = FALSE)
write.csv(overlapping, "blast_hits_overlapping_genes.csv", row.names = FALSE)

## ---- hits in intergenic sequence --------------------------------------------
# The interesting case: a real hit with no gene model over it means the locus exists
# but was never annotated, which is exactly what the reviewer is asking about.
intergenic <- hits_gr[-unique(queryHits(ov))]

if (length(intergenic) > 0) {
  nearest_idx <- nearest(intergenic, genes, ignore.strand = TRUE)
  d <- distance(intergenic, genes[nearest_idx], ignore.strand = TRUE)

  unannotated <- data.frame(
    query   = intergenic$query,
    evalue  = intergenic$evalue,
    pident  = intergenic$pident,
    seqname = as.character(seqnames(intergenic)),
    start   = start(intergenic),
    end     = end(intergenic),
    nearest_gene = genes$gene[nearest_idx],
    distance_bp  = d
  ) %>% arrange(evalue)

  cat("\nHits with no overlapping gene model:\n")
  print(unannotated, row.names = FALSE)
  write.csv(unannotated, "blast_hits_unannotated.csv", row.names = FALSE)
} else {
  cat("\nAll hits fall within annotated genes.\n")
}

## ---- link to your expression data -------------------------------------------
# Any gene recovered here can be checked against the single-cell data directly.
# The GAF and the expression matrix are keyed on gene symbol, and genes.bed carries
# the same symbols, so this joins without extra mapping.
#
# hit_genes <- unique(overlapping$gene)
# intersect(hit_genes, rownames(sce))
# plotExpression(sce, features = intersect(hit_genes, rownames(sce)), x = "label")
