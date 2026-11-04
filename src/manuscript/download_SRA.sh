#!/bin/bash --login
#SBATCH --job-name=sra_download       
#SBATCH --nodes=1                     
#SBATCH --cpus-per-task=8             
#SBATCH --mem=32G                     
#SBATCH --time=04:00:00              
#SBATCH --output=run/sra_job_%j.out       
#SBATCH --account=bioinformaticscore

# 1. Load the SRA Toolkit module (update the name to match your HPCC exact module path)
module purge
module load SRA-Toolkit/3.0.10-gompi-2023a

# 2. Define directories/files
OUTPUT_DIR="/mnt/scratch/hickeys6/yasmin_et_at_PRJNA749754"
SRR_ACC_LIST="/mnt/ufs18/rs-013/bioinformaticsCore/projects/chong_davidson/BCC125_sea_lamprey_scRNA/data/SRR_Acc_List.txt"
mkdir -p "$OUTPUT_DIR"

echo "Starting SRA batch download..."

# 3. Prefetch the compressed SRA files using the list file
prefetch --option-file $SRR_ACC_LIST

echo "Prefetch complete. Starting conversion to FASTQ format..."

# 4. Multi-threaded conversion using local node temp space ($TMPDIR) for peak speed
cat $SRR_ACC_LIST | xargs -n 1 -P 4 fasterq-dump \
    --threads 2 \
    --split-3 \
    --temp "$TMPDIR" \
    --outdir "$OUTPUT_DIR"

echo "FASTQ extraction completed successfully. Files are located in $OUTPUT_DIR"
