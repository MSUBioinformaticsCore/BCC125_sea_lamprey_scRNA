#!/bin/bash --login
#SBATCH --job-name=rnaseq
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --cpus-per-task=1
#SBATCH --output=run/rnaseq-%j.out

# Load Nextflow
module purge
module load Nextflow

# Set the paths to the genome files
PROJECT="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA"
WORK="/mnt/scratch/hickeys6/yasmin_et_at_PRJNA749754"
GENOME_DIR="/mnt/research/bioinformaticsCore/shared/Genomes/Petromyzon_marinus" 
FASTA="$GENOME_DIR/GCF_010993605.1_kPetMar1.pri_genomic.fna" 
GTF="$GENOME_DIR/GCF_010993605.1_kPetMar1.pri_genomic.gtf" 

# Define the samplesheet, outdir, workdir, and config
SAMPLESHEET="$PROJECT/data/sample_sheet.csv" # Example path to sample sheet
OUTDIR="$PROJECT/results/yasmin_nfcore" # Example path to results directory
WORKDIR="$WORK/work"
# Clone the pipeline onto scratch rather than using ~/.nextflow/assets:
# singularity doesn't bind $HOME here, so projectDir/bin never mounts into the
# container and bin/ scripts (e.g. deseq2_qc.r) fail with "command not found".
SOFTWARE="$WORK/nfcore_rnaseq_3.26.0"

mkdir -p $WORKDIR
mkdir -p $OUTDIR
[ -d $SOFTWARE ] || nextflow clone nf-core/rnaseq -r 3.26.0 $SOFTWARE

# Run the RNA-seq analysis
nextflow run $SOFTWARE -resume -profile singularity -work-dir $WORKDIR \
--input $SAMPLESHEET \
--outdir $OUTDIR \
--fasta $FASTA \
--gtf $GTF \
--extra_salmon_quant_args '--seqBias --gcBias' \
--skip_bigwig \
--skip_stringtie \
-c $PROJECT/src/slurm.config.bioinformaticscore.sh \
--multiqc_config $PROJECT/src/custom_multiqc.yml

