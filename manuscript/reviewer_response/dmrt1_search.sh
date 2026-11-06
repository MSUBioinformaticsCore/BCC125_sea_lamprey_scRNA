#!/usr/bin/env bash
# Search kPetMar1.pri for a DMRT1 homolog, and map any hits back to annotated genes.
#
# Answers reviewer 2 minor comment 7: "Absence from the current annotation does not
# necessarily demonstrate absence from the genome."
#
# Strategy, cheapest first:
#   Step 1  blastp against the annotated proteins. Fast, and hits carry gene IDs.
#   Step 2  tblastn against the genome. Catches loci that were never annotated.
#   Step 3  map genomic hit coordinates back to the GFF to see whether they fall in,
#           near, or outside any annotated gene.
#
# Requires: BLAST+ and EDirect. bedtools for step 3 (an R alternative is given).
#   conda install -c bioconda blast entrez-direct bedtools

set -euo pipefail

BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/993/605/GCF_010993605.1_kPetMar1.pri
PFX=GCF_010993605.1_kPetMar1.pri

## ---------------------------------------------------------------- reference files
# Skip any you already have.
wget -nc $BASE/${PFX}_protein.faa.gz
wget -nc $BASE/${PFX}_genomic.fna.gz
wget -nc $BASE/${PFX}_genomic.gff.gz
wget -nc $BASE/${PFX}_feature_table.txt.gz     # maps XP_ protein IDs to GeneIDs
gunzip -kf ${PFX}_protein.faa.gz ${PFX}_genomic.fna.gz

## ---------------------------------------------------------------- query proteins
# Several species, so a single divergent query does not decide the result.
esearch -db protein -query '(DMRT1[Gene Name]) AND (Homo sapiens[Organism] OR Mus musculus[Organism] OR Danio rerio[Organism] OR Gallus gallus[Organism] OR Xenopus tropicalis[Organism]) AND refseq[filter]' \
  | efetch -format fasta > dmrt1_queries.faa

# POSITIVE CONTROL. dmrt2 and dmrt3 are usually present in lamprey. If the pipeline
# recovers those but not dmrt1, the search demonstrably works and the dmrt1 result
# is meaningful. Without this control a negative result proves nothing.
esearch -db protein -query '(DMRT2[Gene Name] OR DMRT3[Gene Name]) AND (Homo sapiens[Organism] OR Danio rerio[Organism]) AND refseq[filter]' \
  | efetch -format fasta > dmrt23_queries.faa

grep -c '^>' dmrt1_queries.faa dmrt23_queries.faa

## ---------------------------------------------------------------- step 1: proteins
makeblastdb -in ${PFX}_protein.faa -dbtype prot -out kpetmar1_prot

for q in dmrt1 dmrt23; do
  blastp -query ${q}_queries.faa -db kpetmar1_prot -evalue 1e-3 -max_target_seqs 20 \
    -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
    -out ${q}_vs_proteins.tsv
  echo "== ${q}: top protein hits =="
  sort -k7,7g ${q}_vs_proteins.tsv | head -20
done

# Map the XP_ hits to GeneIDs and symbols using the feature table.
zcat ${PFX}_feature_table.txt.gz | awk -F'\t' 'NR==1 || $1=="CDS"' \
  | cut -f11,15,16,17 > prot2gene.tsv     # product_accession, symbol, GeneID, name
echo "== dmrt1 protein hits mapped to genes =="
cut -f2 dmrt1_vs_proteins.tsv | sort -u | grep -F -f - prot2gene.tsv | sort -u

## ---------------------------------------------------------------- step 2: genome
makeblastdb -in ${PFX}_genomic.fna -dbtype nucl -out kpetmar1_gen

tblastn -query dmrt1_queries.faa -db kpetmar1_gen -evalue 1e-5 -max_target_seqs 50 \
  -outfmt '6 qseqid sseqid pident length evalue bitscore sstart send' \
  -num_threads 4 -out dmrt1_vs_genome.tsv

echo "== dmrt1 genomic hits =="
sort -k5,5g dmrt1_vs_genome.tsv | head -30

## ---------------------------------------------------------------- step 3: hits to genes
# Convert hits to BED. tblastn reports sstart > send on the minus strand, so order them.
awk -F'\t' 'BEGIN{OFS="\t"}
  { s=($7<$8)?$7:$8; e=($7<$8)?$8:$7; print $2, s-1, e, $1, $6, ($7<$8)?"+":"-" }' \
  dmrt1_vs_genome.tsv | sort -k1,1 -k2,2n > dmrt1_hits.bed

zcat ${PFX}_genomic.gff.gz | awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {
    match($9, /ID=gene-[^;]+/); id=substr($9, RSTART+8, RLENGTH-8);
    print $1, $4-1, $5, id, ".", $7 }' | sort -k1,1 -k2,2n > genes.bed

# Genes overlapping a hit
bedtools intersect -a dmrt1_hits.bed -b genes.bed -wa -wb > dmrt1_hits_in_genes.tsv
echo "== hits falling inside annotated genes =="
cut -f4,10 dmrt1_hits_in_genes.tsv | sort -u

# Nearest annotated gene for hits that fall outside any gene. An unannotated locus is
# the interesting outcome: a real hit in intergenic sequence.
bedtools intersect -a dmrt1_hits.bed -b genes.bed -v > dmrt1_hits_intergenic.bed
bedtools closest -a dmrt1_hits_intergenic.bed -b genes.bed -d > dmrt1_hits_nearest_gene.tsv
echo "== hits outside annotated genes, with nearest gene and distance =="
cat dmrt1_hits_nearest_gene.tsv

echo
echo "Reporting checklist:"
echo " - state the E-value threshold and the number and source of query sequences"
echo " - report the dmrt2/dmrt3 positive control outcome"
echo " - for any hit, give scaffold, coordinates, percent identity and E-value"
echo " - say whether hits fall inside an annotated gene, and if so which"
echo " - a clean negative is a real result and stronger than 'not annotated'"
