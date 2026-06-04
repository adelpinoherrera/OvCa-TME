# OvCa-TME
Systems-Level Mapping of the Tumor Microenvironment in Ovarian Cancer

This repository contains files to download, preprocess and analyze mIF, bulk RNA- and scRNA-seq datasets. This page is divided in different sections based on the modality used to collect the data 

## mIF
- mIF-SpiderPlot.R: contains the code to create a spider plot based on different cell proportions. Comments contain information on how the proportions were scaled to be plotted on the same plot

## scRNA-seq
### Scripts needed to download and align each sample
- bam_fast.slurm: bash script to download scRNAseq sequencing samples that are in a .bam format and convert them to compressed fastq files .fastq.gz needed for cellranger (should get 2 .fastq.gz files representing forward and reverse reads)
- download_sra.slurm: bash script to download scRNAseq sequencing samples that are in .sra formats
- sra_to_compressed_fastq.slurm: bash script to convert .sra files to compressed fastq files .fastq.gz needed for cellranger (should get 2 .fastq.gz files representing forward and reverse reads)
- cellranger_GEX.slurm: bash script to align .fastq.gz files to the human genome (GRCh38) using cellranger
### Scripts to process all the samples in conjunction
