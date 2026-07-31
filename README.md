# OvCa-TME
Systems-Level Mapping of the Tumor Microenvironment in Ovarian Cancer

This repository contains files to download, preprocess and analyze mIF, bulk RNA- and scRNA-seq datasets. This page is divided in different sections based on the modality used to collect the data 

## mIF_code
- mIF-SpiderPlot.R: R script containing the code to create a spider plot based on different cell proportions. Comments contain information on how the proportions were scaled to be plotted on the same plot

## bulkRNA-seq_code
### Scripts needed to align each sample 
- star_to_bam_countsout.slurm: bash script to align .fastq.gz files to the human genome (GRCh38) using STAR
### Scripts to process all the samples in conjunction
- bulkRNAseq.R: R script containing all the analysis done on the bulkRNAseq samples. Briefly:
  - Pre-processing samples
  - Normalizing samples and differential gene expression: cancerous vs non-cancerous
  - Prepare matix for CIBERSORTx: prepare query and reference to upload to CIBERSORTx
- TCGA-OV_bulkRNAseq.R: R script containing all the analysis done on the TCGA-OV bulkRNAseq samples. Briefly:
  - Download and pre-processing samples: download TCGA-OC dataset and extracting metadata of interest. The metadata was exported to excel to find treatment-free intervals accurately 
  - Normalizing samples 
  - Score samples with the COMET signature

## scRNA-seq_code
### Publicly available datasets
#### Scripts needed to download and align each sample
- bam_fast.slurm: bash script to download scRNAseq sequencing samples that are in a .bam format and convert them to compressed fastq files .fastq.gz needed for cellranger (should get 2 .fastq.gz files representing forward and reverse reads)
- download_sra.slurm: bash script to download scRNAseq sequencing samples that are in .sra formats
- sra_to_compressed_fastq.slurm: bash script to convert .sra files to compressed fastq files .fastq.gz needed for cellranger (should get 2 .fastq.gz files representing forward and reverse reads)
- cellranger_GEX.slurm: bash script to align .fastq.gz files to the human genome (GRCh38) using cellranger counts
### Samples processed in-house 
#### Scripts needed to align each sample
- cellranger_multi.slurm: bash script to align .fastq.gz files to the human genome (GRCh38) using cellranger multi which requires a .csv file as an input indicating where the gene expression, BCR and TCR fastq raw files are located 
- UF_1_multi_config_8.0.csv: .csv files indicating the location of gene expression, BCR and TCR fastq raw files needed as an input for cellranger multi
### Scripts to process all the samples in conjunction
- individual_rds-scRNAseq.R: R script to create individual Seurat files for each sample from cellranger outputs. The singleR reference was also merged into a single Seurat object using this script 
- scRNAseq.R: R script containing all the analysis done on the scRNAseq samples. Briefly:
  -  Pre-processing, quality control and normalization
  -  Formatting singleR reference and run singleR for cell type classification for the entire object
  -  Sensitivity classification: re-assignment of sensitive and resistant cells and assignment of unlabeled cancer cells. In this section, the sankey plots per stage were also created
  -  Extract data for COMET: resistant vs all
  -  Pathway analysis: resistant vs all
  -  Pseudotime analysis of clusters 2 (epithelial cluster) and 3 (fibroblast cluster) using Monocle

## pathway_analysis_code
- pathway_analysis_CvsN_bulkandscRNAseq.R: R script for pathway analysis for upregulated and downregulated Hallmark pathways in bulk RNA- and scRNA-seq datasets and shared pathway identification

## math_modeling_code 
- Mathematical-Modeling.ipynb: Julia script for parameterization of logistic growth models for sensitive and resistant cells under standard and M2 TAM conditions
