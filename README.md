# Gonadal PIP-seq reveals genes involved in germ cell development and sexual differentiation in sea lamprey (*Petromyzon marinus*)

**Yu-Wen Chung-Davidson**<sup>1\*</sup>, **Stephanie L. Hickey**<sup>2\*</sup>, Matthew P. Bernard<sup>3</sup>, Sreeja Sarasamma<sup>1</sup>, Ciro Amato<sup>4</sup>, Martin Estermann<sup>4</sup>, Jacob Kimmel<sup>1</sup>, Lily Gorman<sup>5</sup>, Tyler Buchinger<sup>1</sup>, Braden A. Idalski<sup>6</sup>, Nicholas S. Johnson<sup>6</sup>, Humphrey Hung-Chang Yao<sup>4</sup>, Weiming Li<sup>1</sup>

<sup>\*</sup> Co-first authors

<sup>1</sup> Department of Fisheries & Wildlife, Michigan State University, East Lansing, MI 48824
<sup>2</sup> Bioinformatics Core Facility, Michigan State University, East Lansing, MI 48824
<sup>3</sup> Department of Pharmacology & Toxicology, Michigan State University, East Lansing, MI 48824
<sup>4</sup> Reproductive Developmental Biology Group, National Institute of Environmental Health Sciences, Research Triangle Park, NC 27709
<sup>5</sup> Program of Environmental Biology & Plant Biology, Michigan State University, East Lansing, MI 48824
<sup>6</sup> Hammond Bay Biological Station, Great Lakes Science Center, US Geological Survey, Millersburg, MI 49759

---

## Analysis Reports and Code
Stephanie Hickey, PhD — MSU Bioinformatics Core

Raw and processed data are deposited in GEO: [GSE334273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334273)

### PIP-seq Preprocessing & QC

- [Female transformer sample QC](https://msubioinformaticscore.github.io/BCC125_sea_lamprey_scRNA/html/TF_sample_qc.html)
- [Male transformer sample QC](https://msubioinformaticscore.github.io/BCC125_sea_lamprey_scRNA/html/TM_sample_qc.html)
- [Undetermined larvae sample 1 QC](https://msubioinformaticscore.github.io/BCC125_sea_lamprey_scRNA/html/UND1_sample_qc.html)
- [Undetermined larvae sample 2 QC](https://msubioinformaticscore.github.io/BCC125_sea_lamprey_scRNA/html/UND2_sample_qc.html)

### PIP-seq Analysis

- [1. Clustering & annotation (all cells)](https://msubioinformaticscore.github.io/BCC125_sea_lamprey_scRNA/html/01_atlas.html)
  - Figure 1 : `Fig1_tsne.png`
  - Supplementary Table S1 : `median_stats.csv`
  - Supplementary Table S2 : `cluster_marker_genes_unfiltered.csv`
  - Supplementary Table S4 : `cluster_GO_overlap_pvalues.csv`
- [2. Germ cell trajectory analysis](https://msubioinformaticscore.github.io/BCC125_sea_lamprey_scRNA/html/02_germ_pseudotime.html) 
  - Figure 6 : `figure6.png`
  - Figure 8 : `figure8.png`
  - Figure 9 : `figure9.png`
  - Supplementary Table S6 : `node_marker_genes.csv`
  - Supplementary Table S7 : `pseudotime_genes.csv`
  - Supplementary Table S8 : `pseudotime_GO_pvalues.csv`
- [3. Comparisons with Yasmin et al](https://msubioinformaticscore.github.io/BCC125_sea_lamprey_scRNA/html/03_bulk_validation.html)
  - Figure 5 : `bulk_validation_summary.png`
  - Figure 7 : `fine_node_summary.png`
  - Supplementary Table S5 : `bulk_vs_pseudobulk_per_gene.csv`
  
### Gene Expression Figures

- [Gene expression figures]()
  - Figure 10 : `Fig10_Chr81_genes.png`
  - Supplementary Figure S5 : `S5_PIWIL1.png`
  - Supplementary Figure S6 : `S6_SYCP3.png`
  - Meta data for all cells (deposited on GEO): `cell_metadata.tsv.gz`

### Yasmin et al reprocessing

The bulk RNA-seq from Yasmin et al (SRA BioProject [PRJNA749754](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA749754)) was downloaded and requantified against the same genome as the PIP-seq data, GCF_010993605.1 (kPetMar1.pri). The scripts are in `src/manuscript` and are run in this order.

1. **`download_SRA.sh`** (SLURM) downloads every run listed in `data/SRR_Acc_List.txt`. It uses SRA Toolkit 3.0.10: `prefetch` fetches the runs, and `fasterq-dump --split-3` converts each one to paired FASTQ files on scratch.
2. **`makeSampleSheet.R`** writes `data/sample_sheet.csv`, the nf-core/rnaseq input. Each row is one run, with its two FASTQ files and strandedness set to `auto`, so the pipeline infers it. Each sample is named `<Run>_<Sample.Name>` from `data/SraRunTable.csv`.
3. **`run_rnaseq.sh`** (SLURM) runs nf-core/rnaseq 3.26.0 with Singularity. It uses the kPetMar1.pri genome FASTA and GTF, and writes to `results/yasmin_nfcore`.
   - Settings that differ from the pipeline defaults:
     - Salmon is run with `--seqBias --gcBias`.
     - The featureCounts biotype QC groups on `transcript_biotype`.
     - **`slurm.config.bioinformaticscore.sh`** is a Nextflow config that submits each pipeline process as its own SLURM job under the `bioinformaticscore` account.
     - **`custom_multiqc.yml`** sets the section order, the modules and the table columns of the MultiQC report.

`03_bulk_validation.Rmd` reads `salmon.merged.gene_counts_length_scaled.tsv` from `results/yasmin_nfcore` and the sample groups from `data/yasmin_PRJNA749754_sample_groups.csv`.
