# make sample sheet for nfcore-rna-seq
library(tidyverse)
library(readxl)

meta = read_excel("/mnt/research/bioinformaticsCore/projects/FowkesGajan/BCC276_FowkesGajan_RNAseq/data/Comparisons for RNA Seq 7-2026.xlsx",
                  range = "A1:C25")

meta$Replicate = rep(c(1:3), 8)
meta$Timepoint = gsub(" ", "", meta$Timepoint)

meta = 
  meta %>%
  mutate(SampleName = paste0("Tx.", `Tx`, "_", Timepoint, "_rep", Replicate),
         Sample = as.character(Sample),
         SampleGroup = paste0("Tx.", `Tx`, "_", Timepoint))

write.table(meta, file = "/mnt/research/bioinformaticsCore/projects/FowkesGajan/BCC276_FowkesGajan_RNAseq/data/meta.csv", 
            sep = ",", 
            quote=F,
            row.names = F)

data_dir1 = "/mnt/research/bioinformaticsCore/projects/FowkesGajan/BCC276_FowkesGajan_RNAseq/data/fastq/20260720_mRNASeq_PE150"

R1 = list.files(data_dir1, 
                pattern = "_R1.fastq.gz", 
                full.names = TRUE)

R2 = list.files(data_dir1, 
                pattern = "_R2.fastq.gz", 
                full.names = TRUE)

sample = gsub("_R1.fastq.gz", "", basename(R1))
  
samp1 = data.frame(Sample = sample,
                  fastq_1 = R1,
                  fastq_2 = R2,
                  strandedness = "auto") %>%
  left_join(meta)


samp_final = data.frame(sample = samp1$SampleName,
                        fastq_1 = samp1$fastq_1,
                        fastq_2 = samp1$fastq_2,
                        strandedness = "auto")

write.table(samp_final, file = "/mnt/research/bioinformaticsCore/projects/FowkesGajan/BCC276_FowkesGajan_RNAseq/data/sample_sheet.csv", 
            sep = ",", 
            quote=F,
            row.names = F)
