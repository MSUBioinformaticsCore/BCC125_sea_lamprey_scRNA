# make sample sheet for nfcore-rna-seq
library(tidyverse)

project_dir = "/mnt/research/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
meta = read.csv(paste0(project_dir, "/data/SraRunTable.csv"))
meta_sample = 
  meta %>% 
  select(Run, Sample.Name, sex) %>%
  mutate(Run.Sample = paste0(Run, "_", Sample.Name))

data_dir1 = "/mnt/scratch/hickeys6/yasmin_et_at_PRJNA749754/fastq"

R1 = list.files(data_dir1, 
                pattern = "_1.fastq", 
                full.names = TRUE)

R2 = list.files(data_dir1, 
                pattern = "_2.fastq", 
                full.names = TRUE)

run = gsub("_1.fastq", "", basename(R1))
  
samp1 = data.frame(Run = run,
                  fastq_1 = R1,
                  fastq_2 = R2,
                  strandedness = "auto") %>%
  left_join(meta_sample)


samp_final = data.frame(sample = samp1$Run.Sample,
                        fastq_1 = samp1$fastq_1,
                        fastq_2 = samp1$fastq_2,
                        strandedness = "auto")

write.table(samp_final, file = paste0(project_dir, "/data/sample_sheet.csv"), 
            sep = ",", 
            quote=F,
            row.names = F)
