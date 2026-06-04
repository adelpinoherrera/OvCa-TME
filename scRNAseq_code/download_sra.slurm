#!/bin/bash
#SBATCH --job-name=download_SRA_Cancer1           # Job name
#SBATCH --mail-type=ALL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adelpinoherrera@ufl.edu    # Your email address for notifications
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --mem=4gb                         # Memory request
#SBATCH --time=24:00:00                   # Time limit (12 hours)
#SBATCH --output=download_SRA_%j.log          # Standard output and error log

# Load necessary modules (if any)
# If you've installed SRA Toolkit in a custom path, add it to PATH
export PATH=/blue/ferrallm/03_shared_tools/SRA_Toolkit/bin:$PATH

# Define the SRA accession number and output directory
OUTPUT_DIR=/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379608_Cancer1
LOG_FILE=sra_log.txt

# Run prefetch command
prefetch SRX12379608 --max-size 0 --output-directory "$OUTPUT_DIR" > "$OUTPUT_DIR/$LOG_FILE" 2>&1
