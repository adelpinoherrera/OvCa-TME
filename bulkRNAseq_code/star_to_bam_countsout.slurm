#!/bin/bash
#SBATCH --job-name=star_bam_counts_6273         # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adelpinoherrera@ufl.edu   # Your email for notifications
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=8                   # Number of CPU cores per task
#SBATCH --mem=64gb                           # Memory request
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec
#SBATCH --output=star_bam_counts_%j.log    # Standard output and error log

# Load necessary modules (modify based on your cluster)
module load star/2.7.11b    # Use appropriate version

# Set paths
GENOME_DIR="/blue/ferrallm/00_data/single-cell/reference/STAR_index_2.7.11b"
FASTQ_DIR="/blue/ferrallm/00_data/RNAseq/ovarian-cancer/IRB202200692_202501/GE7899-MFerrall/12_6273/raw"
OUTPUT_DIR="/blue/ferrallm/00_data/RNAseq/ovarian-cancer/IRB202200692_202501/GE7899-MFerrall/12_6273/out_counts"
SAMPLE_NAME="6273"

# Define multiple FASTQ files for the same sample (comma-separated)
FASTQ_R1="$FASTQ_DIR/6273_S12_L007_R1_001.fastq.gz,$FASTQ_DIR/6273_S12_L009_R1_001.fastq.gz,$FASTQ_DIR/6273_S41_L008_R1_001.fastq.gz"
FASTQ_R2="$FASTQ_DIR/6273_S12_L007_R2_001.fastq.gz,$FASTQ_DIR/6273_S12_L009_R2_001.fastq.gz,$FASTQ_DIR/6273_S41_L008_R2_001.fastq.gz"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run STAR alignment
echo "Running STAR alignment for sample $SAMPLE_NAME..."
STAR --runThreadN $SLURM_CPUS_ON_NODE \
     --genomeDir $GENOME_DIR \
     --readFilesIn $FASTQ_R1 $FASTQ_R2 \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --twopassMode Basic \
     --quantMode GeneCounts \
     --outSAMattrIHstart 0 \
     --outReadsUnmapped Fastx \
     --outFileNamePrefix $OUTPUT_DIR/${SAMPLE_NAME}_

echo "STAR alignment complete!"
