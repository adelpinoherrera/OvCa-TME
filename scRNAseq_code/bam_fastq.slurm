#!/bin/bash
#SBATCH --job-name=bam_fastq_T89         # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adelpinoherrera@ufl.edu   # Your email for notifications
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=8                   # Number of CPU cores per task
#SBATCH --mem=6gb                           # Memory request
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec
#SBATCH --output=bam_fastq_%j.log    # Standard output and error log

#load necessary modules
module load cellranger

# Define file locations
BAM_URL="https://sra-pub-src-1.s3.amazonaws.com/SRR12244041/T89_possorted_genome_bam.bam.1"
BAM_FILE="/blue/ferrallm/00_data/single-cell/OvCa/GSE154600/T89_bam/T89_possorted_genome_bam.bam.1"
FASTQ_DIR="/blue/ferrallm/00_data/single-cell/OvCa/GSE154600/T89_bam/fastq"

# Download the .bam.1 file
echo "Downloading BAM file..."
wget $BAM_URL -O $BAM_FILE

# Check if the BAM file was downloaded successfully
if [ -f "$BAM_FILE" ]; then
    echo "BAM file downloaded successfully."
else
    echo "Failed to download BAM file."
    exit 1
fi

# Convert BAM to FASTQ using bamtofastq
echo "Converting BAM to FASTQ..."
bamtofastq $BAM_FILE $FASTQ_DIR

# Completion message
echo "BAM to FASTQ conversion complete."
