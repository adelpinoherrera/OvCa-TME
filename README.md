# OvCa-TME
Systems-Level Mapping of the Tumor Microenvironment in Ovarian Cancer

This repository contains files to download, preprocess and analyze mIF, bulk RNA- and scRNA-seq datasets. This page is divided in different sections based on the modality used to collect the data 

## mIF
- mIF-SpiderPlot.R: R script containing the code to create a spider plot based on different cell proportions. Comments contain information on how the proportions were scaled to be plotted on the same plot

## bulk RNA-seq
### Scripts needed to align each sample 

## scRNA-seq
### Scripts needed to download and align each sample
- bam_fast.slurm: bash script to download scRNAseq sequencing samples that are in a .bam format and convert them to compressed fastq files .fastq.gz needed for cellranger (should get 2 .fastq.gz files representing forward and reverse reads)
- download_sra.slurm: bash script to download scRNAseq sequencing samples that are in .sra formats
- sra_to_compressed_fastq.slurm: bash script to convert .sra files to compressed fastq files .fastq.gz needed for cellranger (should get 2 .fastq.gz files representing forward and reverse reads)
- cellranger_GEX.slurm: bash script to align .fastq.gz files to the human genome (GRCh38) using cellranger
### Scripts to process all the samples in conjunction
- individual_rds-scRNAseq.R: R script to create individual Seurat files for each sample from cellranger outputs. The singleR reference was also merged into a single Seurat object using this script 
- scRNAseq.R: R script containing all the analysis done on the scRNAseq samples. Briefly:
  -  Pre-processing, quality control and normalization
  -  Formatting singleR reference and run singleR for cell type classification for the entire object
  -  Sensitivity classification: re-assignment of sensitive and resistant cells and assignment of unlabeled cancer cells. In this section, the sankey plots per stage were also created
  -  Extract data for COMET: resistant vs all
  -  Pathway analysis: resistant vs all
  -  Pseudotime analysis of clusters 2 (epithelial cluster) and 3 (fibroblast cluster) using Monocle


