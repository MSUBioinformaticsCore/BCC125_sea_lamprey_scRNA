#!/bin/bash
#SBATCH --job-name=sra_download       # Job name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --cpus-per-task=8             # Match fasterq-dump threading
#SBATCH --mem=32G                     # Memory required for caching/processing
#SBATCH --time=04:00:00               # Time limit (Hrs:Min:Sec)
#SBATCH --output=sra_job_%j.out       # Standard output log

# 1. Load the SRA Toolkit module (update the name to match your HPCC exact module path)
module purge
module load SRA-Toolkit/3.0.10-gompi-2023a

# 2. Define directories/files
OUTPUT_DIR=""
$SRR_ACC_LIST = 
mkdir -p "$OUTPUT_DIR"

echo "Starting SRA batch download..."

# 3. Prefetch the compressed SRA files using the list file
prefetch --option-file $SRR_ACC_LIST

echo "Prefetch complete. Starting conversion to FASTQ format..."

# 4. Multi-threaded conversion using local node temp space ($TMPDIR) for peak speed
cat $SRR_ACC_LIST | xargs fasterq-dump \
    --threads 8 \
    --split-files \
    --temp "$TMPDIR" \
    --outdir "$OUTPUT_DIR"

echo "FASTQ extraction completed successfully. Files are located in $OUTPUT_DIR"
