#!/bin/bash
#SBATCH --job-name=cellranger_1_PB_multi          # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adelpinoherrera@ufl.edu           # Where to send mail
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --cpus-per-task=16                       # Number of CPU cores per task (adjust as needed)
#SBATCH --mem=132gb                              # Memory request (adjust based on data size)
#SBATCH --time=200:00:00                         # Time limit hrs:min:sec (adjust as needed)
#SBATCH --output=cellranger_multi_%j.log               # Standard output and error log

pwd; hostname; date

# Load CellRanger module and version
module load cellranger/8.0.0

# Run CellRanger Count

cellranger multi \
  --id=cellranger_multi_tutorial \
  --localcores=16 \
  --csv=/blue/ferrallm/00_data/single-cell/OvCa/GE7897-MFerral/UF_1/UF_1_multi_config_8.0.csv 
