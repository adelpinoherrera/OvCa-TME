#!/bin/bash
#SBATCH --job-name=sra_to_fastq_Cancer1              # Job name
#SBATCH --mail-type=ALL                      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adelpinoherrera@ufl.edu       # Where to send mail
#SBATCH --ntasks=1                           # Run on a single CPU
#SBATCH --cpus-per-task=4                    # Number of CPU cores to use (adjust as needed)
#SBATCH --mem=8gb                            # Job memory request
#SBATCH --time=48:00:00                      # Time limit hrs:min:sec
#SBATCH --output=sra_to_fastq_%j.log         # Standard output and error log

#load compression modules
module load gcc/5.2.0
module load pigz/2.4

# Add SRA Toolkit installed in shared tools to PATH
export PATH=/blue/ferrallm/03_shared_tools/SRA_Toolkit/bin:$PATH

# Convert sra files to fastq

for subfolder in /blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379608_Cancer1/SRR*; do
    # Get all .sra files in the subfolder
    sra_files=($(ls "$subfolder"/*.sra | sort))
    last_sra_file="${sra_files[-1]}"

    # Convert each .sra file to paired-end FASTQ
    for sra_file in "${sra_files[@]}"; do
        if [[ "$sra_file" == "$last_sra_file" ]]; then
            # If it's the last file, run the conversion and log the output
            fasterq-dump --split-files "$sra_file" --outdir "$subfolder" > "${subfolder}/fastq_log.txt" 2>&1
        else
            # For other files, just run the conversion
            fasterq-dump --split-files "$sra_file" --outdir "$subfolder" > "${subfolder}/fastq_log.txt" 2>&1
        fi
    done
done

#compress fastq files to fastq.gz

for subfolder in /blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379608_Cancer1/SRR*; do
    echo "Processing $subfolder"
    
    # Compress all FASTQ files in the current subfolder using pigz
    pigz -p 8 "$subfolder"/*.fastq
    
done

echo "Compression complete!"
