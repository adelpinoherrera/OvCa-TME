#This script shows the preprocessing, quality control, normalization and downstream analysis for the scRNAseq cohort used in the study 

#Samples were downloaded from publicly available datasets (GSE184880, GSE154600 and GSE181955) and merged with 4 samples processed in-house. FASTQ files were downloaded from all samples and ran through cellranger version 8.0.0 for alignment to the human genome

#Features, barcodes and matrix files from each sample were used to create Seurat objects (.rds objects), code shown in individual_rds-scRNAseq.R These are the first inputs for this script

#Load libraries
library('ggplot2')
library('SeuratObject')
library('Seurat')
library('patchwork')
library('forcats')
library('tidyverse') 
library('data.table')
library('future')
plan("sequential")
options(future.globals.maxSize = Inf)
library('sctransform')
library('harmony')
library('clustree')
library('plyr')
library('hdf5r')
library('devtools')
#only 1 time BiocManager::install("SingleR")
library('SingleR')
library(BiocParallel)
#only 1 time BiocManager::install("ArrayExpress")
library(ArrayExpress)
library('celldex')
library('BiocManager') 
library('monocle3')
#only 1 time BiocManager::install('glmGamPoi')
library('glmGamPoi')
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(msigdbr)
library(ggrepel)

#set up output directories
outdir <- '/results/scRNAseq/'
dir.create(outdir)
dir <- '/data/scRNAseq/'

#BiocParallel::register(BiocParallel::MulticoreParam(workers = 8)) 
#options(future.globals.maxSize = Inf) #if not SCTransform won't run, 50 * 1024^3

######################################################################
###Preprocessing, quality control and normalization 
######################################################################
outdir1 <- '/results/scRNAseq/preprocessing/'
dir.create(outdir1)
###Read seurat objects 
T59 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-T59_GSE154600-Unprocessed-RDS_2025-04-14.rds'))
T61 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-T61_GSE154600-Unprocessed-RDS_2025-04-14.rds'))
T76 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-T76_GSE154600-Unprocessed-RDS_2025-04-14.rds'))
T89 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-T89_GSE154600-Unprocessed-RDS_2025-04-14.rds'))
T90 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-T90_GSE154600-Unprocessed-RDS_2025-04-14.rds'))
Normal1 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Normal1_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Normal2 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Normal2_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Normal3 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Normal3_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Normal4 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Normal4_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Normal5 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Normal5_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Cancer1 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Cancer1_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Cancer2 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Cancer2_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Cancer3 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Cancer3_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Cancer4 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Cancer4_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Cancer5 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Cancer5_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Cancer6 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Cancer6_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
Cancer7 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-Cancer7_GSE184880-Unprocessed-RDS_2025-04-14.rds'))
N1_normal <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-N1_normal_GSE181955-Unprocessed-RDS_2025-04-14.rds'))
OMT1_CD45 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-OMT1_CD45+_GSE181955-Unprocessed-RDS_2025-04-14.rds'))
OMT3_CD45 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-OMT3_CD45+_GSE181955-Unprocessed-RDS_2025-04-14.rds'))
T1 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-T1_GSE181955-Unprocessed-RDS_2025-04-14.rds'))
T6 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-T6_GSE181955-Unprocessed-RDS_2025-04-14.rds'))
UF_1 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-UF_1_multi-Unprocessed-RDS_2025-04-14.rds'))
UF_2 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-UF_2_multi-Unprocessed-RDS_2025-04-14.rds'))
UF_3 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-UF_3_multi-Unprocessed-RDS_2025-04-14.rds'))
UF_4 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-UF_4_multi-Unprocessed-RDS_2025-04-14.rds'))

###Merge objects
ovca <- merge(x=T59, y=list(T61,T76,T89,T90,Normal1,Normal2,Normal3,Normal4,Normal5,Cancer1,Cancer2,Cancer3,
                            Cancer4,Cancer5,Cancer6,Cancer7,N1_normal,OMT1_CD45,OMT3_CD45,T1,T6,UF_1,UF_2,UF_3,UF_4), merge.data=TRUE, project="ovca_3GSEand4CORE")
#saveRDS(ovca, paste(outdir1,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-Unprocessed-RDS.rds",sep=""))

###Join layers
ovca <- JoinLayers(ovca, assay="RNA", new.layer="counts")
saveRDS(ovca, paste(outdir1,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-JoinedLayers-Unprocessed.rds",sep=""))

###Check dimensions of the dataset
dim(ovca)
#  36601 192118, almost 200,000 cells prior to filtering 

summary(ovca@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4    1220    2291    2947    4073   15510 

####FILTERING AND REMOVING MT- GENES####
## Determining how upper cutoff for nFeatures compares
nFeatUpper_ovca <- mean(ovca@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(ovca@meta.data$nFeature_RNA, na.rm=TRUE) #7588.164
nFeatLower_ovca <- 450 #can filter out less cells later if it got rid off too many immune cells

#make a violin plot to see how many cells we should keep
pdf(paste(outdir1, "OvCa_3GSEand4CORE-ViolinPlot_nFeatureRNA", ".pdf",sep=""), width = 11, height = 6)
vlnplot1 <- VlnPlot(ovca, features = "nFeature_RNA", group.by = "orig.ident")
print(vlnplot1)
dev.off()

#store mitochondrial percentage in object meta data
ovca <- PercentageFeatureSet(ovca, pattern = "^MT-", col.name = "percent.mt")
perMitoUpper_ovca <- 25 #from the plot, keeping the cells from 0 to 25 seems like a good option 

#make a violin plot to see how many cells we should keep
pdf(paste(outdir1, "OvCa_3GSEand4CORE-ViolinPlot_Percent.mito", ".pdf",sep=""), width = 11, height = 6)
vlnplot2 <- VlnPlot(ovca, features = "percent.mt", group.by = "orig.ident")
print(vlnplot2)
dev.off()

###Visualize all QC metrics
pdf(paste(outdir1,"OvCa_3GSEand4CORE-VlnPlot_nFeature+nCount+PercentMito",".pdf",sep=""), width = 11, height = 6)
v <- VlnPlot(ovca, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(v)
dev.off()

pdf(paste(outdir1, "OvCa_3GSEand4CORE-FeatureScatter_nCountvMito_nCountvnFeature", ".pdf",sep=""), width = 20, height = 6)
plot1 <- FeatureScatter(ovca, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ovca, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cPlot <- plot1 + plot2
print(cPlot)
dev.off()

## original dataset dimensions
dim(ovca)
#[1]  36601 192118

## filtered based on approach
dim(subset(ovca, subset = nFeature_RNA > nFeatLower_ovca & nFeature_RNA < nFeatUpper_ovca & percent.mt < perMitoUpper_ovca))
#36601 165296

#filter out the cells based on features and mitochondrial content
ovca <-subset(ovca, subset = nFeature_RNA > nFeatLower_ovca & nFeature_RNA < nFeatUpper_ovca & percent.mt < perMitoUpper_ovca)
saveRDS(ovca, paste(outdir1,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-JoinedLayers-PostFiltering.rds", sep=""))

####NORMALIZE DATASET USING SCT TRANSFORM####
# run sctransform to normalize the samples based on percent.mt: substitute normalizedata, scaledata, findvariablefeatures
ovca <- SCTransform(ovca, conserve.memory = TRUE, verbose = FALSE)
dim(ovca)
#33612 165296
#save normalized objects 
saveRDS(ovca, paste(outdir1,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-JoinedLayers-PostSCT.rds",sep=""))

####RUN PCA####
#when running PCA after SCT, there's no need to specify VariableFeatures but it can be specified if needed 
ovca <- RunPCA(ovca, features=VariableFeatures(object=ovca))
# PC_ 1 
# Positive:  CD74, C1QC, IL1B, C1QB, SPP1, HLA-DRA, C1QA, APOE, CCL3, HLA-DRB1 
# CCL4L2, HLA-DPA1, CCL4, CCL3L1, CXCL8, RNASE1, HLA-DQA1, HLA-DRB5, IFI30, LYZ 
# FTL, CCL20, HLA-DPB1, TYROBP, PLAUR, G0S2, SOD2, CD83, SLC16A10, NFKBIA 
# Negative:  COL1A1, COL1A2, COL3A1, SPARC, TAGLN, IGFBP7, LUM, DCN, ACTA2, MMP11 
# CALD1, C11orf96, PLA2G2A, FN1, COL6A3, COL6A2, TIMP3, KCNIP4, CTHRC1, COL6A1 
# POSTN, RARRES2, MGP, COL11A1, MEG8, SPARCL1, COL5A2, PTGDS, IGFBP2, TPM2 
# PC_ 2 
# Positive:  COL1A1, COL1A2, COL3A1, LUM, FN1, SPARC, TAGLN, DCN, IGFBP7, MMP11 
# APOE, ACTA2, C1QC, C1QB, TIMP1, PLA2G2A, CALD1, CTHRC1, COL6A3, IL1B 
# POSTN, C11orf96, SPP1, CD74, COL6A2, COL11A1, C1QA, TIMP3, SFRP2, COL6A1 
# Negative:  GNLY, CCL5, NKG7, SKAP1, GZMA, CD247, TOX, PTPRC, GZMB, CD96 
# CD69, PRKCH, THEMIS, ARHGAP15, PLCG2, CD2, FYN, LTB, PARP8, AOAH 
# IL7R, LINC01934, AL136456.1, PTPN22, PDE3B, CTSW, CD52, KLRB1, GZMH, KLRD1 
# PC_ 3 
# Positive:  COL1A1, GNLY, COL1A2, CCL5, COL3A1, CCL4, NKG7, LUM, MMP11, FN1 
# GZMB, SPARC, GZMA, CD247, PTPRC, POSTN, CTHRC1, DCN, SKAP1, TOX 
# CD96, CD69, COL11A1, SFRP2, AOAH, CCL4L2, PRKCH, THEMIS, KLRD1, CTSW 
# Negative:  RHEX, WFDC2, AC024230.1, MECOM, CLU, SLPI, CLDN4, ERBB4, THSD4, KRT8 
# MMP7, KRT18, LCN2, NPAS3, KRT17, ELF3, TACSTD2, KRT7, MUC16, CRYAB 
# CD24, NRG3, TM4SF1, KRT19, CHODL, CLDN3, ITGB8, CP, KIAA1217, IFIT2 
# PC_ 4 
# Positive:  IGHG3, IGHG1, IGKC, IGHG4, IGHGP, IGLC2, JCHAIN, MZB1, IGLC3, IGLC1 
# DERL3, IGHG2, LUM, IGHM, CD79A, FN1, MMP11, SSR4, PTGDS, SFRP2 
# IFNG-AS1, IGHA1, COL1A1, HES1, POU2AF1, XBP1, CTHRC1, IGKV4-1, COL11A1, COMP 
# Negative:  GNLY, IGFBP7, VWF, C11orf96, KCNIP4, PLA2G2A, SPARCL1, LDB2, MGP, CCL4 
# MEG8, C7, ACSM3, ANO2, CALCRL, CCL5, MEG3, STAR, COL4A1, SERPINE2 
# EMCN, FLT1, CCL4L2, NKG7, PTPRB, DLC1, RBMS3, ADGRL4, NAV3, FBXL7 
# PC_ 5 
# Positive:  VWF, IGFBP7, IGHG3, IGHG1, IGHG4, SPARCL1, LDB2, IGKC, ANO2, C11orf96 
# CALCRL, KCNIP4, MGP, PLA2G2A, IGHGP, JCHAIN, EMCN, FLT1, CCL14, C7 
# ADGRL4, ACSM3, PTPRB, MEG8, COL4A1, MZB1, IGLC2, STAR, CLDN5, INSR 
# Negative:  COL1A1, MMP11, COL1A2, LUM, WFDC2, COL3A1, RHEX, GNLY, SLPI, FN1 
# CTHRC1, POSTN, CLU, CLDN4, AC024230.1, COL11A1, PTGDS, SFRP2, MMP7, KRT17 
# KRT8, LCN2, CRYAB, IFIT2, KRT18, TACSTD2, KIF26B, KRT7, ELF3, VCAN 

pdf(paste(outdir1, "OvCa_3GSEand4CORE-JoinedLayers-PCA_DimPlot",".pdf",sep=""), width = 11, height = 6)
DimPlot(ovca, reduction = "pca", pt.size = 1)
dev.off()

# Determine dimensionality of dataset
pdf(paste(outdir1,"OvCa_3GSEand4CORE-JoinedLayers-PCA_ElbowPlot",".pdf",sep=""), width = 11, height = 6)
ElbowPlot(ovca, ndims = 50) #from this plot, we can see that 40 PCs will capture all the deviation in the dataset
dev.off()
saveRDS(ovca, paste(outdir1,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-JoinedLayers-PostPCA.rds", sep=""))

####RUN HARMONY#### 
ovcaH <- RunHarmony(object=ovca, group.by.vars="orig.ident") #can do other groupings if necessary
#initially had reduction="pca" it didn't run with that

# Transposing data matrix
# Initializing state using k-means centroids initialization
# Harmony 1/10
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         Harmony 2/10
#       0%   10   20   30   40   50   60   70   80   90   100%
#         [----|----|----|----|----|----|----|----|----|----|
#            **************************************************|
#            0%   10   20   30   40   50   60   70   80   90   100%
#            [----|----|----|----|----|----|----|----|----|----|
#               **************************************************|
#               Harmony 3/10
#             0%   10   20   30   40   50   60   70   80   90   100%
#               [----|----|----|----|----|----|----|----|----|----|
#                  **************************************************|
#                  0%   10   20   30   40   50   60   70   80   90   100%
#                  [----|----|----|----|----|----|----|----|----|----|
#                     **************************************************|
#                     Harmony 4/10
#                   0%   10   20   30   40   50   60   70   80   90   100%
#                     [----|----|----|----|----|----|----|----|----|----|
#                        **************************************************|
#                        0%   10   20   30   40   50   60   70   80   90   100%
#                        [----|----|----|----|----|----|----|----|----|----|
#                           **************************************************|
#                           Harmony 5/10
#                         0%   10   20   30   40   50   60   70   80   90   100%
#                           [----|----|----|----|----|----|----|----|----|----|
#                              **************************************************|
#                              0%   10   20   30   40   50   60   70   80   90   100%
#                              [----|----|----|----|----|----|----|----|----|----|
#                                 **************************************************|
#                                 Harmony 6/10
#                               0%   10   20   30   40   50   60   70   80   90   100%
#                                 [----|----|----|----|----|----|----|----|----|----|
#                                    **************************************************|
#                                    0%   10   20   30   40   50   60   70   80   90   100%
#                                    [----|----|----|----|----|----|----|----|----|----|
#                                       **************************************************|
#                                       Harmony 7/10
#                                     0%   10   20   30   40   50   60   70   80   90   100%
#                                       [----|----|----|----|----|----|----|----|----|----|
#                                          **************************************************|
#                                          0%   10   20   30   40   50   60   70   80   90   100%
#                                          [----|----|----|----|----|----|----|----|----|----|
#                                             **************************************************|
#                                             Harmony converged after 7 iterations
#                                           Warning message:
#                                             Quick-TRANSfer stage steps exceeded maximum (= 8264800) 


pdf(paste(outdir1, "OvCa_3GSEand4CORE-JoinedLayers-PCA_DimPlot_PostHarmony",".pdf",sep=""), width = 11, height = 6)
DimPlot(ovcaH, reduction = "harmony", pt.size = 1)
dev.off()
saveRDS(ovcaH, paste(dir,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-JoinedLayers-PostHarmony.rds", sep=""))

####FIND NEIGHBORS AND CLUSTERS####
#for the HARMONY dataset
ovcaH <- FindNeighbors(ovcaH, dims=1:50, reduction = "harmony")
ovcaH <- FindClusters(ovcaH, resolution=seq(0.025,0.3, by=0.025)) #25 clusters
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
# 
# Number of nodes: 165296
# Number of edges: 6068314
# 
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.9930
#    Number of communities: 8
#    Elapsed time: 104 seconds
#    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#    
#    Number of nodes: 165296
#    Number of edges: 6068314
#    
#    Running Louvain algorithm...
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         Maximum modularity in 10 random starts: 0.9875
#       Number of communities: 10
#       Elapsed time: 80 seconds
#       Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#       
#       Number of nodes: 165296
#       Number of edges: 6068314
#       
#       Running Louvain algorithm...
#       0%   10   20   30   40   50   60   70   80   90   100%
#         [----|----|----|----|----|----|----|----|----|----|
#            **************************************************|
#            Maximum modularity in 10 random starts: 0.9825
#          Number of communities: 13
#          Elapsed time: 94 seconds
#          Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#          
#          Number of nodes: 165296
#          Number of edges: 6068314
#          
#          Running Louvain algorithm...
#          0%   10   20   30   40   50   60   70   80   90   100%
#            [----|----|----|----|----|----|----|----|----|----|
#               **************************************************|
#               Maximum modularity in 10 random starts: 0.9778
#             Number of communities: 15
#             Elapsed time: 90 seconds
#             Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#             
#             Number of nodes: 165296
#             Number of edges: 6068314
#             
#             Running Louvain algorithm...
#             0%   10   20   30   40   50   60   70   80   90   100%
#               [----|----|----|----|----|----|----|----|----|----|
#                  **************************************************|
#                  Maximum modularity in 10 random starts: 0.9741
#                Number of communities: 17
#                Elapsed time: 90 seconds
#                Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                
#                Number of nodes: 165296
#                Number of edges: 6068314
#                
#                Running Louvain algorithm...
#                0%   10   20   30   40   50   60   70   80   90   100%
#                  [----|----|----|----|----|----|----|----|----|----|
#                     **************************************************|
#                     Maximum modularity in 10 random starts: 0.9708
#                   Number of communities: 19
#                   Elapsed time: 87 seconds
#                   Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                   
#                   Number of nodes: 165296
#                   Number of edges: 6068314
#                   
#                   Running Louvain algorithm...
#                   0%   10   20   30   40   50   60   70   80   90   100%
#                     [----|----|----|----|----|----|----|----|----|----|
#                        **************************************************|
#                        Maximum modularity in 10 random starts: 0.9677
#                      Number of communities: 20
#                      Elapsed time: 90 seconds
#                      Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                      
#                      Number of nodes: 165296
#                      Number of edges: 6068314
#                      
#                      Running Louvain algorithm...
#                      0%   10   20   30   40   50   60   70   80   90   100%
#                        [----|----|----|----|----|----|----|----|----|----|
#                           **************************************************|
#                           Maximum modularity in 10 random starts: 0.9647
#                         Number of communities: 21
#                         Elapsed time: 101 seconds
#                         Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                         
#                         Number of nodes: 165296
#                         Number of edges: 6068314
#                         
#                         Running Louvain algorithm...
#                         0%   10   20   30   40   50   60   70   80   90   100%
#                           [----|----|----|----|----|----|----|----|----|----|
#                              **************************************************|
#                              Maximum modularity in 10 random starts: 0.9627
#                            Number of communities: 22
#                            Elapsed time: 90 seconds
#                            Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                            
#                            Number of nodes: 165296
#                            Number of edges: 6068314
#                            
#                            Running Louvain algorithm...
#                            0%   10   20   30   40   50   60   70   80   90   100%
#                              [----|----|----|----|----|----|----|----|----|----|
#                                 **************************************************|
#                                 Maximum modularity in 10 random starts: 0.9602
#                               Number of communities: 24
#                               Elapsed time: 100 seconds
#                               Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                               
#                               Number of nodes: 165296
#                               Number of edges: 6068314
#                               
#                               Running Louvain algorithm...
#                               0%   10   20   30   40   50   60   70   80   90   100%
#                                 [----|----|----|----|----|----|----|----|----|----|
#                                    **************************************************|
#                                    Maximum modularity in 10 random starts: 0.9577
#                                  Number of communities: 25
#                                  Elapsed time: 108 seconds
#                                  Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                                  
#                                  Number of nodes: 165296
#                                  Number of edges: 6068314
#                                  
#                                  Running Louvain algorithm...
#                                  0%   10   20   30   40   50   60   70   80   90   100%
#                                    [----|----|----|----|----|----|----|----|----|----|
#                                       **************************************************|
#                                       Maximum modularity in 10 random starts: 0.9557
#                                     Number of communities: 26
#                                     Elapsed time: 117 seconds


#create the clustree plot to pick a resolution 
pdf(paste0(outdir1, 'OvCa_3GSEand4CORE-JoinedLayers-clustree_withHarmony.pdf'), width = 9, height = 7)
clustree(ovcaH, prefix = "SCT_snn_res.")
dev.off()

resUseH <- 0.025  #8 clusters
ovcaH$clusterResolution_0.025 <- as.factor(as.numeric(as.character(ovcaH$SCT_snn_res.0.025)))

saveRDS(ovcaH, paste(outdir1,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-JoinedLayers-PostNeighborANDClustering_withHarmony.rds", sep=""))

####RUN UMAP####
#for the HARMONY dataset 
ovcaH <- RunUMAP(ovcaH, dims = 1:50, reduction = "harmony")
# 15:41:59 UMAP embedding parameters a = 0.9922 b = 1.112
# 15:41:59 Read 192118 rows and found 50 numeric columns
# 15:41:59 Using Annoy for neighbor search, n_neighbors = 30
# 15:41:59 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      15:43:02 Writing NN index file to temp file /scratch/local/65141394/Rtmptjp0K8/file16ea2450ee8d7a
#    15:43:02 Searching Annoy index using 1 thread, search_k = 3000
#    15:45:00 Annoy recall = 100%
#    15:45:02 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
#    15:45:12 Initializing from normalized Laplacian + noise (using RSpectra)
#    15:46:22 Commencing optimization for 200 epochs, with 9175296 positive edges
#    15:46:22 Using rng type: pcg
#    Using method 'umap'
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         15:48:34 Optimization finished


#divide by clusters based on the resolution chosen 
pdf(paste0(outdir1, 'UMAP_res=', resUseH, '-OvCa_3GSEand4CORE-JoinedLayers-UMAP_clusters_withHarmony.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap", group.by = paste0("clusterResolution_", resUseH))
dev.off()
#divive by samples
pdf(paste0(outdir1, 'OvCa_3GSEand4CORE-JoinedLayers-UMAP_origIdent_withHarmony.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap",group.by = 'orig.ident')
dev.off()
#divide by diagnosis: normal or HGSOC
pdf(paste0(outdir1, 'OvCa_3GSEand4CORE-JoinedLayers-UMAPbyDiagnosis-normalGrey.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap",group.by = 'diagnosis', cols = c("normal" = "#a9a9a9", "HGSOC"= "#00CED1"))
dev.off()

#divide by tissue type:normal, omentum, primary 
pdf(paste0(outdir1, 'OvCa_3GSEand4CORE-JoinedLayers-UMAP_tissueType_withHarmony.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap",group.by = 'tissue')
dev.off()
saveRDS(ovcaH, paste(outdir1,"OvCa-scRNAseq-Pilot-ovca_3GSEand4CORE-JoinedLayers-PostUMAP_withHarmony.rds", sep=""))

#get proportion of normal vs HGSOC per cluster
df <- data.frame("cluster_0.025"=ovcaH$clusterResolution_0.025,"diagnosis"=ovcaH$diagnosis)
df$cluster_0.025 <- paste("cluster", df$cluster_0.025)
split.df <- split(df,df$cluster_0.025)

#use this 
for (i in 1:length(split.df)){
  freq_table <- as.data.frame(table(split.df[[i]]))
  this_step <- freq_table$Freq/sum(freq_table$Freq)
  this_step.df <- as.data.frame(t(as.data.frame(this_step)))
  colnames(this_step.df) <- freq_table$diagnosis
  rownames(this_step.df) <- unique(freq_table$cluster_0.025)
  if (i>1){
    all.df <- rbind.fill(all.df,this_step.df)
  }else{
    all.df <- this_step.df
  }
}

#make NA values 0 
all.df[is.na(all.df)] <-0

#rename the rows with the right names 
row.names(all.df)<- names(split.df)

#write file
write.csv(all.df, '/results/scRNAseq/preprocessing/OvCa_3GSEand4CORE-JoinedLayers-cellDiagnosis_perCluster_resolution0.025_withHarmony.csv')

##################################################################################################
###Formatting singleR reference and run singleR for cell type classification for the entire object 
##################################################################################################
#define a new output directory for the singleR results 
outdir2 <- '/results/scRNAseq/singleR_for_cellType_classification/'
dir.create(outdir2)
####REFERENCE FORMATTING####
##GSE173682
reference_GSE173682 <- readRDS(paste0(dir,'OvCa-scRNAseq-Pilot-reference-GSE173682-JoinedLayers-unprocessed_2025-04-21.rds'))
dim(reference_GSE173682)#33538 15120
reference_GSE17_cell_names <- colnames(reference_GSE173682) #to see the intersection 
 
#import the metadata for the reference too 
#import the metadata which was located in TableS2 from the paper, only the HGSOC sheet 
metadata_GSE173682 <- read.csv(paste0(dir,'TableS2_Multi-omic_sc_landscape_human_gynecologic_malignancies.csv'))
metadata_GSE17_cell_names <- metadata_GSE173682[['Barcode']]
length(metadata_GSE173682[['SingleR.ovar']]) #13646

####Subset of the reference to only include the cells that are part of the metadata
#find the intersect between the cell names from the reference and the metadata 
common_cells_GSE17 <- intersect(reference_GSE17_cell_names, metadata_GSE17_cell_names) #13646

#subset the seurat object to only contain those cells 
subset_reference_GSE17 <- subset(reference_GSE173682, cells = common_cells_GSE17) #looks good and the cells kept their order 
rownames(subset_reference_GSE17) #gene names
colnames(subset_reference_GSE17) #cell names in this format AAACCCAAGACCAAAT-1_1

#see if the metadata and the seurat object match 
matches <- metadata_GSE173682$Barcode == colnames(subset_reference_GSE17) 
sum(matches) #13646 TRUE values so the cell names match 

saveRDS(subset_reference_GSE17, paste(outdir2,"OvCa-scRNAseq-Pilot-subset-reference-GSE173682-JoinedLayers-unprocessed.rds"))

###Process the reference with SCT
subset_reference_GSE17.norm <- SCTransform(subset_reference_GSE17, verbose = FALSE)
  
saveRDS(subset_reference_GSE17.norm, paste(outdir2,"OvCa-scRNAseq-Pilot-subset-reference-GSE173682-JoinedLayers-PostSCT.rds"))

###Run PCA
subset_reference_GSE17.norm <- RunPCA(subset_reference_GSE17.norm, features = VariableFeatures(object = subset_reference_GSE17.norm))
# PC_ 1 
# Positive:  KRT18, SLPI, WFDC2, MT1G, KRT19, S100A9, KRT8, KRT7, PI3, MAL2 
# GSTP1, ACTG1, SPINT2, RPL8, KRT6A, MT-CO2, MT1H, S100A6, CD9, FOLR1 
# PNOC, TACSTD2, KRT17, KLK6, LCN2, MT-ND4, CD24, S100A13, CLDN7, RPL30 
# Negative:  COL1A1, COL1A2, COL3A1, DCN, SPARC, LUM, C11orf96, VIM, COL6A2, TNFAIP6 
# COL6A1, IGFBP7, FN1, POSTN, CTSL, CTSK, CTHRC1, COL6A3, MMP11, SFRP2 
# MMP1, PDPN, MT2A, SERPINF1, COL5A2, PMP22, TIMP1, PCOLCE, LGALS1, SPARCL1 
# PC_ 2 
# Positive:  COL1A1, COL1A2, COL3A1, DCN, SPARC, LUM, C11orf96, COL6A2, COL6A1, IGFBP7 
# POSTN, CTSK, CTHRC1, COL6A3, MMP11, RARRES2, SFRP2, FN1, SERPINF1, MMP1 
# COL5A2, AEBP1, PCOLCE, KRT18, MT2A, SPARCL1, MMP2, CALD1, S100A6, SLPI 
# Negative:  FCER1G, HLA-DRA, FTL, APOE, SRGN, LAPTM5, C1QB, CD74, TYROBP, EREG 
# THBS1, C1QA, HLA-DRB1, C15orf48, SLC16A10, SAMSN1, CXCL5, CTSL, CXCL3, HLA-DPA1 
# C1QC, MS4A7, HLA-DPB1, PLA2G7, SPP1, GPR183, IL1B, C5AR1, KYNU, HMOX1 
# PC_ 3 
# Positive:  S100A9, MT1G, MT1E, PNOC, MT1H, CXCL1, RPL8, KLK6, S100A6, MT1M 
# RPL30, COX6C, S100A4, LY6E, MT1X, RPS18, RPS24, CD9, RIPK2, MAL2 
# TFPI2, KRT23, ATP1B1, PABPC1, PTX3, VDAC2, SDF2L1, BCAT1, LAPTM4B, FAT1 
# Negative:  WFDC2, PI3, IFI27, SLPI, PAGE2, UCHL1, FOLR1, MT-CO2, KRT6A, TUBA4A 
# CTAG2, CD74, BIRC3, CLU, AREG, GIPC1, GSTP1, HMOX1, MT-ATP6, GAL 
# MCUB, EIF4EBP1, UPK1B, RHOD, SFN, PAGE5, MT-CO1, EIF3G, MT-CYB, LCN2 
# PC_ 4 
# Positive:  LUM, CTSK, MMP11, MMP1, RARRES2, SFRP2, PDPN, POSTN, DCN, CTHRC1 
# WFDC2, SLPI, IFI27, FN1, TNFAIP6, TMEM158, MMP2, AEBP1, GREM1, RSAD2 
# COL1A1, SFRP4, COL3A1, CD74, COL6A1, PI3, HLA-DRA, COL11A1, FBLN1, S100A6 
# Negative:  IGFBP7, MT2A, COL4A1, MGP, C11orf96, TIMP1, HMOX1, CXCL8, COL4A2, RGS5 
# IL1RL1, AKR1C1, SERPINE1, G0S2, GNG11, IGFBP3, BGN, EDNRB, ANGPT2, COLEC11 
# IQCG, STC1, IER3, SLC7A2, TBX3, NPTX2, CXCL2, IL33, C7, PDK4 
# PC_ 5 
# Positive:  CXCL8, CXCL3, CXCL1, CCL20, IL1B, CXCL2, CXCL5, MT2A, CCL3, IER3 
# TIMP1, G0S2, C11orf96, SERPINB2, SLPI, EREG, WFDC2, CTSL, TNFAIP6, BCL2A1 
# IL6, CCL3L1, PLAUR, CCL4, S100A9, AKR1C1, C15orf48, PTGS2, KRT18, PI3 
# Negative:  CD7, CXCR4, TXNIP, PTPRC, GNLY, ARHGDIB, MALAT1, IL7R, ANGPT2, B2M 
# CST7, CD247, TRBC2, HLA-A, RGS1, TSC22D3, A2M, ESM1, CD3D, STK4 
# SERPINE1, SMAP2, FXYD5, HLA-B, CYTIP, SRGN, LEPROTL1, CD2, TCF4, ADGRL4 

#delete the numbers before the cell type on the metadata table 
metadata_GSE173682$cell.type
# [1] "9-Fibroblast"             "3-Epithelial cell"        "4-Fibroblast"             "17-B cell"               
# [5] "3-Epithelial cell"        "15-Empty/Epithelial cell" "4-Fibroblast"             "6-Macrophage"            
# [9] "3-Epithelial cell"        "2-Epithelial cell"        "5-T cell"                 "3-Epithelial cell"       
# [13] "9-Fibroblast"             "6-Macrophage"             "4-Fibroblast"             "2-Epithelial cell" 

metadata_GSE173682$cell.type <- gsub("^[0-9]+-", "", metadata_GSE173682$cell.type)
#check on the number of cells per cell type in the metadata table
table(metadata_GSE173682$cell.type)
# B cell Empty/Epithelial cell      Endothelial cell       Epithelial cell            Fibroblast            Macrophage 
# 192                   308                   334                  5743                  3597                  1947 
# NK cell                T cell 
# 418                  1107 

#re-order the cell names to have the same order as the seurat object
metadata_GSE173682.reorder <- metadata_GSE173682 %>% arrange(match(Barcode, colnames(subset_reference_GSE17.norm)))
#include the celltype to the seurat object 
subset_reference_GSE17.norm@meta.data[["cell.type"]] <- metadata_GSE173682.reorder$cell.type
table(subset_reference_GSE17.norm$cell.type)
# B cell Empty/Epithelial cell      Endothelial cell       Epithelial cell            Fibroblast            Macrophage 
# 192                   308                   334                  5743                  3597                  1947 
# NK cell                T cell 
# 418                  1107
saveRDS(subset_reference_GSE17.norm, paste(outdir2,"OvCa-scRNAseq-Pilot-subset-reference-GSE173682-JoinedLayers-PostPCA-withCellTypes.rds"))
#downsample to have the sample number of cells for each cell type, if not singleR is going to favor the label for the most amount of cells
Idents(subset_reference_GSE17.norm) <- "cell.type"
subset_reference_GSE17.norm.balanced <- subset(subset_reference_GSE17.norm, downsample = 192)
table(subset_reference_GSE17.norm.balanced$cell.type)
# B cell Empty/Epithelial cell      Endothelial cell       Epithelial cell            Fibroblast            Macrophage 
# 192                   192                   192                   192                   192                   192 
# NK cell                T cell 
# 192                   192 


#all markers
Idents(subset_reference_GSE17.norm) <- subset_reference_GSE17.norm$cell.type #by cell name but group by cell type 
markers <- FindAllMarkers(subset_reference_GSE17.norm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste(outdir2,'AllMarkers_referenceGSE173682.csv'))

#create a heatmap 
top_markers <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

pdf(paste0(outdir2, 'Heatmap_top10genesPerCellType_referenceGSE173682.pdf'), width = 12, height = 12)
DoHeatmap(subset_reference_GSE17.norm, features = top_markers$gene)
# Warning message:
#   In DoHeatmap(subset_reference_GSE17.norm, features = top_markers$gene) :
#   The following features were omitted as they were not found in the scale.data slot for the SCT assay: CLEC1A, ROBO4, IGHV3-23, IGLV3-1, TNNI3, PRSS8, CDH1, CLDN6
dev.off()

####RUN SINGLER#### 
#both downsampled and not downsampled gave similar results so keep the full dataset to have more information
pred.sample_ALL <- SingleR(test = ovcaH@assays[["SCT"]]@data, 
                           ref = subset_reference_GSE17.norm@assays[["SCT"]]@data, 
                           labels = metadata_GSE173682[['cell.type']], 
                           assay.type.ref = "data",
                           BPPARAM = MulticoreParam(workers=4))

#save output as csv so we can re-load this cell type assignment at any time 
out.cellType.csv <- data.frame("Cell Name" = pred.sample_ALL@rownames, "Orig.ident" = ovcaH$orig.ident, "Predicted cell type" = pred.sample_ALL@listData[["labels"]],
                               "Prune cell type" = pred.sample_ALL@listData[["pruned.labels"]], "Predicted Score" = pred.sample_ALL@listData[["scores"]],
                               "Delta value" = pred.sample_ALL$delta.next)
table(out.cellType.csv$Predicted.cell.type)
# B cell Empty/Epithelial cell      Endothelial cell       Epithelial cell            Fibroblast            Macrophage 
# 2652                  2312                  4214                 16659                 26685                 24150 
# NK cell                T cell 
# 14538                 74086 

write.csv(out.cellType.csv, paste0(outdir2, "ALL_cellType_GSE173682reference_SCTnorm.csv", row.names = F))

###add cell type to seurat object 
ovcaH$predicted.cellType <- out.cellType.csv$Predicted.cell.type
ovcaH$pruned.cellType <- out.cellType.csv$Prune.cell.type
ovcaH$predicted.score.bcell <- out.cellType.csv$Predicted.Score.B.cell
ovcaH$predicted.score.emptycell <- out.cellType.csv$Predicted.Score.Empty.Epithelial.cell
ovcaH$predicted.score.endocell <- out.cellType.csv$Predicted.Score.Endothelial.cell
ovcaH$predicted.score.epicell <- out.cellType.csv$Predicted.Score.Epithelial.cell
ovcaH$predicted.score.fibroblast <- out.cellType.csv$Predicted.Score.Fibroblast
ovcaH$predicted.score.macrophage <- out.cellType.csv$Predicted.Score.Macrophage
ovcaH$predicted.score.NKcell <- out.cellType.csv$Predicted.Score.NK.cell
ovcaH$predicted.score.Tcell <- out.cellType.csv$Predicted.Score.T.cell
ovcaH$delta.value.cellType <- out.cellType.csv$Delta.value

saveRDS(ovcaH, paste(outdir2, "OvCa-scRNAseq-Pilot_ovca-JoinedLayers_ALL-PostUMAPmodel_withHarmony_withcellTypes.rds", sep = ""))

#save certain plots 
pdf(paste0(outdir2, 'OvCa-ALL_UMAPbyCellType.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap",group.by = 'predicted.cellType')
dev.off()

pdf(paste0(outdir2, 'OvCa-ALL_UMAPbycellType-GraphPadcolorsNEWImmune_noleyend.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap", group.by = "predicted.cellType", cols = c("B cell" = "#DE8BF9", "Endothelial cell" = "#A00000", 
                                                                                "Epithelial cell" = "#00CED1","Fibroblast" = "#C5944E", 
                                                                                "Macrophage" = "#410257", "NK cell" = "#79069A",
                                                                                "T cell"= "#B856D7")) + theme(legend.position = "none")
dev.off()

pdf(paste0(outdir2, 'OvCa-ALL_DotPlotByDiagnosis_CD163-MS4A1-CD8A-ACTA2-FCGR3A-WFDC2-NEW.pdf'),width = 4.5, height = 6.5)
DotPlot(ovcaH, features = c("WFDC2","ACTA2","CD163","FCGR3A","CD8A","MS4A1"), group.by = "diagnosis", dot.scale = 15, cols = c("lightgrey", "blue")) + coord_flip()
dev.off()

pdf(paste0(outdir2, 'OvCa-ALL_DotPlotbyCanonicalMarkersperCellType-reorder.pdf'),width = 10, height = 5)
ovcaH$predicted.cellType <- factor(ovcaH$predicted.cellType, levels = c("T cell","Endothelial cell","B cell","Fibroblast","Empty/Epithelial cell","Epithelial cell","Macrophage","NK cell"))
DotPlot(ovcaH, group.by = "predicted.cellType", features = c("FCGR3A","CD163","WFDC2","ACTA2","MS4A1","ESM1","CD8A"),dot.scale = 12)
dev.off()

#save tables 
counts.orig.ident <- as.data.frame(table(ovcaH$orig.ident, ovcaH$predicted.cellType))
write.csv(counts.orig.ident, paste(outdir2, 'Allcounts_perOrig.ident_perPredicted.CellType.csv'))

counts.diagnosis <- as.data.frame(table(ovcaH$diagnosis, ovcaH$predicted.cellType))
write.csv(counts.diagnosis, paste(outdir2, 'Allcounts_perDiagnosis_perPredicted.CellType.csv'))

counts.clusters <- as.data.frame(table(ovcaH$clusterResolution_0.025, ovcaH$predicted.cellType))
write.csv(counts.clusters, paste(outdir2, 'Allcounts_perCluster_perPredicted.CellType.csv'))

#################################################################################
###Sensitivity classification: re-assignment of sensitive and resistant cells and 
#assignment of unlabeled cancer cells
#################################################################################
####Sensitivity re-labeling with correct subsetting 
outdir3 <- '/results/singleR_for_sensitivity_classfication/'
dir.create(outdir3)
###rename all the refractory cells as resistant in a new response slot 
table(ovcaH$response)
# na refractory  resistant  sensitive 
# 112824      19279      23468       9725 
ovcaH$response.2 <- ovcaH$response
table(ovcaH$response.2)
# na refractory  resistant  sensitive 
# 112824      19279      23468       9725

ovcaH$response.2[ovcaH$response.2 == 'refractory'] <- 'resistant'
table(ovcaH$response.2)
# na       resistant sensitive 
# 112824     42747      9725 

###subset the object to only include the resistant and sensitive samples 
ovcaH.SandR <- subset(ovcaH, response.2 %in% c('resistant', 'sensitive'))
table(ovcaH.SandR$response.2)
# resistant sensitive 
# 42747      9725
#save object
saveRDS(ovcaH.SandR, paste(outdir3, "OvCa-scRNAseq-Pilot_ovca-JoinedLayers_ONLYsensitiveANDresistant-PostUMAPmodel_withHarmony.rds", sep = ""))

Idents(ovcaH.SandR) <- "response.2"
ovcaH.SandR.reference <- subset(ovcaH.SandR, downsample = 4862) #total sensitive 9725 and half it is 4862
table(ovcaH.SandR.reference$response.2)
#resistant sensitive 
#4862      4862 
saveRDS(ovcaH.SandR.reference, paste(outdir3, "reference-JoinedLayers_9724sensitiveANDresistant-withEmbeddings.rds", sep = ""))

#now do the query, use the entire object and eliminate the cells with the barcodes from the reference 
reference_barcodes <- colnames(ovcaH.SandR.reference)
query_barcodes <- !(colnames(ovcaH.SandR) %in% reference_barcodes)
table(query_barcodes)
# FALSE  TRUE 
# 9724 42748

ovcaH.SandR.query <- ovcaH.SandR[ ,query_barcodes] 
table(ovcaH.SandR.query$response.2)
# resistant sensitive, total 42748
# 37885      4863 

saveRDS(ovcaH.SandR.query, paste(outdir3, "query-JoinedLayers_42748sensitiveANDresistant-withEmbeddings.rds", sep = ""))

####RUN SINGLER TO RE-LABEL HALF OF THE ALREADY LABELED SENSITIVE AND RESISTANT CELLS - SET 1
metadata_reference.SandR <- cbind(colnames(ovcaH.SandR.reference), ovcaH.SandR.reference$response.2)
colnames(metadata_reference.SandR) <- c("Barcode", "response.2")
metadata_reference.SandR <- as.data.frame(metadata_reference.SandR)

pred.sample_query.SandR.joined.SCT.set1 <- SingleR(test = ovcaH.SandR.query@assays[["SCT"]]@data, 
                                                   ref = ovcaH.SandR.reference@assays[["SCT"]]@data, 
                                                   labels = metadata_reference.SandR[['response.2']], 
                                                   assay.type.ref = "data",
                                                   BPPARAM = MulticoreParam(workers=4))
out.query.SandR.set1.csv <- data.frame("Cell Name" = pred.sample_query.SandR.joined.SCT.set1@rownames, "Orig.ident" = ovcaH.SandR.query$orig.ident, "Predicted response" = pred.sample_query.SandR.joined.SCT.set1@listData[["labels"]], "Prune response" = pred.sample_query.SandR.joined.SCT.set1@listData[["pruned.labels"]], "Response.2" = ovcaH.SandR.query$response.2, "Predicted Score" = pred.sample_query.SandR.joined.SCT.set1@listData[["scores"]], "Delta value" = pred.sample_query.SandR.joined.SCT.set1$delta.next)
table(out.query.SandR.set1.csv$Orig.ident, out.query.SandR.set1.csv$Predicted.response)
#                   resistant sensitive
# T59_GSE154600     11700      3398
# T61_GSE154600     15396      1696
# T76_GSE154600      4234      1461
# T89_GSE154600       157      2503
# T90_GSE154600        89      2114

write.csv(out.query.SandR.set1.csv, paste(outdir3,'singleRresults_querySet1_labelsANDscores.csv', row.names = F))

####SWITCH THE QUERY AND REFERENCE FROM SET 1 TO RE-ASSIGNED THE OTHER ALREADY LABELED SENSITIVE AND RESISTANT CELLS - SET 2

table(ovcaH.SandR.query$response.2)
# resistant sensitive 
# 37885      4863 

table(ovcaH.SandR.reference$response.2)
# resistant sensitive 
# 4862      4862

#get all the barcodes
SandR.barcodes <- colnames(ovcaH.SandR) #52472

#get barcodes from reference and query 
reference.set1.barcodes <- colnames(ovcaH.SandR.reference) # 9724 cells

query.set1.barcodes <- colnames(ovcaH.SandR.query) #42748

#find the indices to keep for the query and the reference 
reference.set2.barcodes <- !(SandR.barcodes %in% reference.set1.barcodes) #52472-9724 = 42748
table(reference.set2.barcodes)
# FALSE  TRUE 
# 9724 42748 
ovcaH.SandR.reference.set2 <- ovcaH.SandR[ ,reference.set2.barcodes] 
table(ovcaH.SandR.reference.set2$response.2)
# resistant sensitive, total 42748
# 37885      4863 

#downsample the number of resistant cells to only have 4863
set.seed(42)
Idents(ovcaH.SandR.reference.set2) <- "response.2"
ovcaH.SandR.reference.set2.balanced <- subset(ovcaH.SandR.reference.set2, downsample= 4863)

table(ovcaH.SandR.reference.set2.balanced$response.2)
# resistant sensitive 
# 4863      4863 

#QUERY should be just the old references
ovcaH.SandR.query.set2 <- ovcaH.SandR.reference
table(ovcaH.SandR.query.set2$response.2)
# resistant sensitive 
# 4862      4862 
reference.set2.barcodes <- colnames(ovcaH.SandR.reference.set2.balanced)
query.set2.barcodes <- colnames(ovcaH.SandR.query.set2)

intersection.reference <- intersect(reference.set1.barcodes, reference.set2.barcodes) #empty
intersection.query <- intersect(query.set1.barcodes, query.set2.barcodes) #empty

intersection.reference.object <- intersect(colnames(reference.set1.barcodes), colnames(reference.set2.barcodes)) #NULL
intersection.query.object <- intersect(colnames(query.set1.barcodes), colnames(query.set2.barcodes)) #NULL
#no matches between the old dataset and the new one, different cells in each object 

intersection.reference.query <- intersect(colnames(ovcaH.SandR.reference.set2.balanced), colnames(ovcaH.SandR.query.set2)) #empty
#no matches between the reference and the query either 

###save the seurat object set2
saveRDS(ovcaH.SandR.reference.set2.balanced, paste(outdir3, "referenceSet2-JoinedLayers_9726sensitiveANDresistant-withEmbeddingsANDHarmonyReductions.rds", sep = ""))
saveRDS(ovcaH.SandR.query.set2, paste(outdir3, "querySet2-JoinedLayers_9724sensitiveANDresistant-withEmbeddingsANDHarmonyReductions.rds", sep = "")) 


####run singleR for the new set of SandR cells - new subsetting
###make a table with the metadata 
metadata_reference.SandR.set2 <- data.frame("Barcode" = colnames(ovcaH.SandR.reference.set2.balanced), "response.2" = ovcaH.SandR.reference.set2.balanced$response.2)

pred.sample_query.SandR.joined.SCT.set2 <- SingleR(test = ovcaH.SandR.query.set2@assays[["SCT"]]@data, 
                                                   ref = ovcaH.SandR.reference.set2.balanced@assays[["SCT"]]@data, 
                                                   labels = metadata_reference.SandR.set2[['response.2']], 
                                                   assay.type.ref = "data",
                                                   BPPARAM = MulticoreParam(workers=4))
out.query.SandR.set2.csv <- data.frame("Cell Name" = pred.sample_query.SandR.joined.SCT.set2@rownames, "Orig.ident" = ovcaH.SandR.query.set2$orig.ident, "Predicted response" = pred.sample_query.SandR.joined.SCT.set2@listData[["labels"]], "Prune response" = pred.sample_query.SandR.joined.SCT.set2@listData[["pruned.labels"]], "Response.2" = ovcaH.SandR.query.set2$response.2, "Predicted Score" = pred.sample_query.SandR.joined.SCT.set2@listData[["scores"]], "Delta value" = pred.sample_query.SandR.joined.SCT.set2$delta.next)
table(out.query.SandR.set2.csv$Orig.ident, out.query.SandR.set2.csv$Predicted.response)
#      SET2          resistant sensitive
# T59_GSE154600      1453       475
# T61_GSE154600      2004       183
# T76_GSE154600       518       229
# T89_GSE154600       135      2496
# T90_GSE154600       88       2143

#        SET1       resistant sensitive
# T59_GSE154600     11700      3398
# T61_GSE154600     15396      1696
# T76_GSE154600      4234      1461
# T89_GSE154600       157      2503
# T90_GSE154600        89      2114

#          OG       resistant sensitive
# T59_GSE154600     17026         0
# T61_GSE154600     19279         0
# T76_GSE154600      6442         0
# T89_GSE154600         0      5291
# T90_GSE154600         0      4434

write.csv(out.query.SandR.set2.csv, paste(outdir3,'singleRresults_querySet2_labelsANDscores.csv', row.names = F))

###add Predicted.response, Prune.response, Predicted.Score.resistant, Predicted.Score.sensitive, 
###Delta.value to the seurat object > ovcaH.SandR
intersection.query.table <- intersect(out.query.SandR.set1.csv$Cell.Name, out.query.SandR.set2.csv$Cell.Name) #empty, no cellnames shared 

##join both the set1 and set2 tables 
predicted.set1.set2.table <- rbind(out.query.SandR.set1.csv, out.query.SandR.set2.csv)
table(predicted.set1.set2.table$Orig.ident,predicted.set1.set2.table$Predicted.response)
#                   resistant sensitive
# T59_GSE154600     13153      3873
# T61_GSE154600     17400      1879
# T76_GSE154600      4752      1690
# T89_GSE154600       292      4999
# T90_GSE154600       177      4257

##re-order the Cell.Name on the table to match the colnames from the ovcaH.SandR
predicted.set1.set2.table.reorder <- predicted.set1.set2.table %>% arrange(match(Cell.Name, colnames(ovcaH.SandR)))
write.csv(predicted.set1.set2.table.reorder, paste(outdir3,'singleRresults_ALLcellsSandRreordered_labelsANDscores.csv', row.names = F))

##add meta.data variables to the seurat object 
ovcaH.SandR@meta.data[["predicted.response"]] <- predicted.set1.set2.table.reorder$Predicted.response
ovcaH.SandR@meta.data[["prune.response"]] <- predicted.set1.set2.table.reorder$Prune.response
ovcaH.SandR@meta.data[["predicted.score.resistant"]] <- predicted.set1.set2.table.reorder$Predicted.Score.resistant
ovcaH.SandR@meta.data[["predicted.score.sensitive"]] <- predicted.set1.set2.table.reorder$Predicted.Score.sensitive
ovcaH.SandR@meta.data[["delta.value"]] <- predicted.set1.set2.table.reorder$Delta.value

saveRDS(ovcaH.SandR, paste(outdir3, "ovca-JoinedLayers_ONLYsensitiveANDresistant_withHarmony_withPredictedResponse.rds", sep = ""))

####RUN SINGLER FOR THE UNLABELED CELLS IN THE OVCAH OBJECT USING THE OVCAH.SANDR AS REFERENCE
table(ovcaH$response)
# na refractory  resistant  sensitive 
# 112824      19279      23468       9725 
#subset the cells not assigned to sensitivity
ovcaH.noresponse <- subset(ovcaH, response %in% c('na'))
dim(ovcaH.noresponse) #33612 112824
saveRDS(ovcaH.noresponse, paste(outdir3, "ovca-JoinedLayers_noResponse-PostUMAPmodel_withHarmony.rds", sep = ""))

###create a balance object from the SandR reference to run singleR 
dim(ovcaH.SandR) #33612 52472
table(ovcaH.SandR$predicted.response)
# resistant sensitive 
# 35774     16698 
#downsample resistant cells to have the same number of the sensitive cells
Idents(ovcaH.SandR) <- "predicted.response"
ovcaH.SandR.balanced <- subset(ovcaH.SandR, downsample = 16698)
table(ovca.SandR.balanced$predicted.response)
# resistant sensitive 
# 16698     16698 

###start singleR 
metadata_reference.SandR.ALL <- data.frame("Barcode" = colnames(ovcaH.SandR.balanced), "predicted.response" = ovcaH.SandR.balanced$predicted.response)
#predicted.response and prune.response look the same for this dataset

pred.sample_query.noresponse<- SingleR(test = ovcaH.noresponse@assays[["SCT"]]@data, 
                                       ref = ovcaH.SandR.balanced@assays[["SCT"]]@data, 
                                       labels = metadata_reference.SandR.ALL[['predicted.response']], 
                                       assay.type.ref = "data",
                                       BPPARAM = MulticoreParam(workers=4))
out.query.noresponse.csv <- data.frame("Cell Name" = pred.sample_query.noresponse@rownames, "Orig.ident" = ovcaH.noresponse$orig.ident, "Diagnosis" = ovcaH.noresponse$diagnosis, "Predicted response" = pred.sample_query.noresponse@listData[["labels"]], "Prune response" = pred.sample_query.noresponse@listData[["pruned.labels"]], "Predicted Score" = pred.sample_query.noresponse@listData[["scores"]], "Delta value" = pred.sample_query.noresponse$delta.next)
table(out.query.noresponse.csv$Orig.ident, out.query.noresponse.csv$Predicted.response)
#                        resistant sensitive
# Cancer1_GSE184880          979      6187
# Cancer2_GSE184880         1010      1685
# Cancer3_GSE184880          255      4292
# Cancer4_GSE184880          722      1519
# Cancer5_GSE184880          440      4751
# Cancer6_GSE184880         2813      1085
# Cancer7_GSE184880          484      3960
# DTprimary_multi            730      1426
# IM_multi                  3226      8099
# LS_multi                  2454      5210
# N1_normal_GSE181955         44      4605
# Normal1_GSE184880          518      5643
# Normal2_GSE184880         3815      1133
# Normal3_GSE184880         1242      2762
# Normal4_GSE184880         2379      3599
# Normal5_GSE184880          969      2117
# OMT1_CD45+_GSE181955       518      8324
# OMT3_CD45+_GSE181955       139      8989
# PB_multi                  3340      3502
# T1_GSE181955              1199      1231
# T6_GSE181955               689      4740

write.csv(out.query.noresponse.csv, paste(outdir3, 'singleRresults_queryNoResponse_labelsANDscores.csv', row.names = F))

###add the response to the object and re-label the sensitive and resistant in the normal diagnosis to NA
out.query.noresponse.csv <- out.query.noresponse.csv %>% mutate(Predicted.response = ifelse(Diagnosis == "normal", NA, Predicted.response))
##create a table with the same categories from the ovcaH.SandR object 
out.reference.SandR.csv <- data.frame("Cell Name" = colnames(ovcaH.SandR), "Orig.ident" = ovcaH.SandR$orig.ident, "Diagnosis" = ovcaH.SandR$diagnosis, 
                                      "Predicted response" = ovcaH.SandR$predicted.response, "Prune response" = ovcaH.SandR$prune.response, 
                                      "Predicted Score resistant" = ovcaH.SandR$predicted.score.resistant, "Predicted Score sensitive" = ovcaH.SandR$predicted.score.sensitive,
                                      "Delta value" = ovcaH.SandR$delta.value)
##merge the 2 dataframes 
all.query.reference <- rbind(out.query.noresponse.csv, out.reference.SandR.csv)
#also re-label the sensitive and resistant in the prune response to NA if diagnosis is normal 
all.query.reference <- all.query.reference %>% mutate(Prune.response = ifelse(Diagnosis == "normal", NA, Prune.response))

##re-order the Cell.Name on the table to match the colnames from the ovcaH.SandR
all.query.reference.reorder <- all.query.reference %>% arrange(match(Cell.Name, colnames(ovcaH)))
#actually relabel the NA in the Predicted.response and Prune.response to normal
all.query.reference.reorder <- all.query.reference.reorder %>% mutate(Predicted.response = ifelse(Diagnosis == "normal", "normal", Predicted.response))
all.query.reference.reorder <- all.query.reference.reorder %>% mutate(Prune.response = ifelse(Diagnosis == "normal", "normal", Prune.response))
write.csv(all.query.reference.reorder, paste(outdir3,'singleRresults_ALLqueryNoResponseANDSandRordered_labelsANDscores.csv', row.names = F))

##add meta.data variables to the seurat object 
ovcaH$response.2 <- ovcaH$response
ovcaH$response.2[ovcaH$response.2 == 'refractory'] <- 'resistant'

ovcaH@meta.data[["predicted.response"]] <- all.query.reference.reorder$Predicted.response
ovcaH@meta.data[["prune.response"]] <- all.query.reference.reorder$Prune.response
ovcaH@meta.data[["predicted.score.resistant"]] <- all.query.reference.reorder$Predicted.Score.resistant
ovcaH@meta.data[["predicted.score.sensitive"]] <- all.query.reference.reorder$Predicted.Score.sensitive
ovcaH@meta.data[["delta.value"]] <- all.query.reference.reorder$Delta.value

saveRDS(ovcaH, paste(outdir3, "ovca-JoinedLayers_ALL-PostUMAPmodel_withHarmony_withCellTypes ANDPredictedResponse.rds", sep = ""))

##make different useful plots 
response.3 <- data.frame("diagnosis" = ovcaH$diagnosis, "response.2" = ovcaH$response.2)
response.3$response.2[response.3$diagnosis == "normal"] <- "normal" 
ovcaH$response.3 <- response.3$response.2
table(ovcaH$response.3)
# na    normal resistant sensitive 
# 83998     28826     42747      9725 

pdf(paste0(outdir3, 'OvCa-ALL_UMAPbyOGresponse-NormalGrey.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap",group.by = 'response.3', cols = c("normal" = "darkgrey", "sensitive" = "red", "resistant" = "darkgreen", "na" = "darkturquoise"))
dev.off()

pdf(paste0(outdir3, 'OvCa-ALL_UMAPbyPredicted.response_greyANDredANDdarkgreen.pdf'),width = 15, height = 10)
DimPlot(ovcaH, reduction = "umap",group.by = 'predicted.response', cols = c("resistant" = "darkgreen", "sensitive" = "red","normal" = "darkgrey"))
dev.off()

pdf(paste0(outdir3, 'OvCa-ALL_DotPlotByPredicted.Response_MS4A1-CD8A-FCGR3A-CD163-ACTA2-WFDC2-NEW.pdf'),width = 12, height = 4.5)
DotPlot(ovcaH, features = c("MS4A1","CD8A","FCGR3A","CD163","ACTA2","WFDC2"), group.by = "predicted.response", dot.scale = 15)
dev.off()


counts.sensitivity <- as.data.frame(table(ovcaH$predicted.response, ovcaH$predicted.cellType))
write.csv(counts.diagnosis, paste(outdir3, 'Allcounts_perSensitivity_perPredicted.CellType.csv'))

###looked at this table(ovcaH$predicted.response, ovcaH$predicted.cellType, ovcaH$stage)
#to get sensitive vs resistant proportions and then cell type proportions based on sensitive and resistant for each compartment
#saved in an excel sheet

####Make sankey plots with different proportions 
#Code to make the sankey plots but need to be saved manually 
####plot with percentages on treated and untreated samples of sensitive and resistant and from those each cell type 
cell_types <- c("B cell","T cell","NK cell","Macrophages","Fibroblasts",
                "Endothelial","Epithelial")
Iuntreatvstreat_all <- data.frame(
  source = c("Untreated","Untreated",
             rep("Sensitive_u", length(cell_types)),
             rep("Resistant_u", length(cell_types)),
             rep(cell_types, each = 2),
             "Sensitive_t",
             "Resistant_t"),
  target = c("Sensitive_u","Resistant_u",
             cell_types,
             cell_types,
             rep(c("Sensitive_t","Resistant_t"), times = length(cell_types)),
             "Treated","Treated"),
  value  = c(
    0.834889575, 0.165110425,
    c(0.028548211, 0.616940744, 0.152385433, 0.124157988, 0.000748455,
      0.008104698, 0.069114471),
    c(0.003568339, 0.045739619, 0.001621972, 0.227292388, 0.235510381,
     0.056228374, 0.430038927),
    c(0.022996766,0.017694415,
      0.770691101,0.37180634,
      0.054497545,0.033963214,
      0.073062642,0.136831218,
      0.013354893,0.259322413,
      0.008863337,0.016520378,
      0.056533717,0.163862023),
    0.318226864, 0.681773136)
)

nodesI <- data.frame(
 name = c("Untreated",           # layer 1
         "Sensitive_u", "Resistant_u",  # layer 2
        cell_types,                  # layer 3 ← YOUR EXACT ORDER
       "Sensitive_t", "Resistant_t",# layer 4
      "Treated"),
group = c(1,  # Untreated
        2, 2,  # Sensitive_u, Resistant_u
        rep(3, 7),  # cell_types (your exact order!)
       4, 4,  # Sensitive_t, Resistant_t
      5),    # Treated
stringsAsFactors = FALSE
)

# Make group a factor with desired vertical order
nodesI$group <- factor(nodesI$group, levels = 1:5)

# Add ID columns to links (0-based indices)
Iuntreatvstreat_all$IDsource <- match(Iuntreatvstreat_all$source, nodesI$name) - 1
Iuntreatvstreat_all$IDtarget <- match(Iuntreatvstreat_all$target, nodesI$name) - 1

pI <- sankeyNetwork(Links = Iuntreatvstreat_all, Nodes = nodesI,
                 Source = "IDsource", Target = "IDtarget",
                Value = "value", NodeID = "name", NodeGroup = "group", iterations = 0,
               sinksRight=FALSE, fontSize = 0, nodeWidth = 15, nodePadding = 30, margin = 5)
pI

htmlwidgets::saveWidget(pI, "sankey_temp_I.html")

##### make sankey plots per stages ####
###Untreated
#Stage 1
cell_types <- c("B cell","T cell","NK cell","Macrophages","Fibroblasts",
               "Endothelial","Epithelial")
stage1_all <- data.frame(
source = c("Stage 1","Stage 1",
          rep("Sensitive", length(cell_types)),
         rep("Resistant", length(cell_types))),
target = c("Sensitive","Resistant",
         cell_types,
        cell_types),
value  = c(
0.869925214, 0.130074786,
c(0.077167127, 0.680790093, 0.140517859, 0.065909324, 0.000102344,
0.000921093, 0.03459216),
c(0.010951403, 0.06844627, 0.000684463, 0.254620123, 0.045174538,
0.036276523, 0.58384668)
))

#to make a sankey plot instead with teh sankey plot code 
nodesStage1 <- data.frame(
name = c("Stage 1",           # layer 1
        "Sensitive", "Resistant",  # layer 2
       cell_types),
group = c(1,  # Untreated
        2, 2,  # Sensitive_u, Resistant_u
       rep(3, 7)),
stringsAsFactors = FALSE
)

# Make group a factor with desired vertical order
nodesStage1$group <- factor(nodesStage1$group, levels = 1:5)

# Add ID columns to links (0-based indices)
stage1_all$IDsource <- match(stage1_all$source, nodesStage1$name) - 1
stage1_all$IDtarget <- match(stage1_all$target, nodesStage1$name) - 1

pStage1 <- sankeyNetwork(Links = stage1_all, Nodes = nodesStage1,
                  Source = "IDsource", Target = "IDtarget",
                 Value = "value", NodeID = "name", NodeGroup = "group", iterations = 0,
                sinksRight=FALSE, fontSize = 0, nodeWidth = 15, nodePadding = 30, margin = 5)
pStage1

htmlwidgets::saveWidget(pStage1, "sankey_stage1.html")

#Stage 2
cell_types <- c("B cell","T cell","NK cell","Macrophages","Fibroblasts",
              "Endothelial","Epithelial")
stage2_all <- data.frame(
source = c("Stage 2","Stage 2",
          rep("Sensitive", length(cell_types)),
         rep("Resistant", length(cell_types))),
target = c("Sensitive","Resistant",
         cell_types,
        cell_types),
value  = c(
0.81612985, 0.18387015,
c(0.043971411, 0.667339963, 0.12538844, 0.111559975, 0.00077688,
0.019266625, 0.031696706),
c(0.006896552, 0.060689655, 0, 0.269655172, 0.353103448,
0.129655172, 0.18)
))

#to make a sankey plot instead with teh sankey plot code 
nodesStage2 <- data.frame(
name = c("Stage 2",           # layer 1
        "Sensitive", "Resistant",  # layer 2
       cell_types),
group = c(1,  # Untreated
        2, 2,  # Sensitive_u, Resistant_u
       rep(3, 7)),
stringsAsFactors = FALSE
)

# Make group a factor with desired vertical order
nodesStage2$group <- factor(nodesStage2$group, levels = 1:5)

# Add ID columns to links (0-based indices)
stage2_all$IDsource <- match(stage2_all$source, nodesStage2$name) - 1
stage2_all$IDtarget <- match(stage2_all$target, nodesStage2$name) - 1

pStage2 <- sankeyNetwork(Links = stage2_all, Nodes = nodesStage2,
                       Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name", NodeGroup = "group", iterations = 0,
                     sinksRight=FALSE, fontSize = 0, nodeWidth = 15, nodePadding = 30, margin = 5)
pStage2

htmlwidgets::saveWidget(pStage2, "sankey_stage2.html")

#Stage 3
cell_types <- c("B cell","T cell","NK cell","Macrophages","Fibroblasts",
              "Endothelial","Epithelial")
stage3_all <- data.frame(
source = c("Stage 3","Stage 3",
          rep("Sensitive", length(cell_types)),
          rep("Resistant", length(cell_types))),
target = c("Sensitive","Resistant",
         cell_types,
        cell_types),
value  = c(0.278347871, 0.721652129,
c(0.000921659, 0.592626728, 0.098617512, 0.200921659, 0.024884793,
 0.001843318, 0.080184332),
c(0.000355492, 0.039459652, 0.001777462, 0.09953786, 0.540703875,
0.015641664, 0.302523996)
))

#to make a sankey plot instead with teh sankey plot code 
nodesStage3 <- data.frame(
name = c("Stage 3",           # layer 1
       "Sensitive", "Resistant",  # layer 2
      cell_types),
group = c(1,  # Untreated
      2, 2,  # Sensitive_u, Resistant_u
     rep(3, 7)),
stringsAsFactors = FALSE
)

# Make group a factor with desired vertical order
nodesStage3$group <- factor(nodesStage3$group, levels = 1:5)

# Add ID columns to links (0-based indices)
stage3_all$IDsource <- match(stage3_all$source, nodesStage3$name) - 1
stage3_all$IDtarget <- match(stage3_all$target, nodesStage3$name) - 1

pStage3 <- sankeyNetwork(Links = stage3_all, Nodes = nodesStage3,
                       Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name", NodeGroup = "group", iterations = 0,
                     sinksRight=FALSE, fontSize = 0, nodeWidth = 15, nodePadding = 30, margin = 5)
pStage3

htmlwidgets::saveWidget(pStage3, "sankey_stage3.html")

###TREATED
#Stage 3
cell_types <- c("B cell","T cell","NK cell","Macrophages","Fibroblasts",
              "Endothelial","Epithelial")
stage3t_all <- data.frame(
source = c("Stage 3","Stage 3",
          rep("Sensitive", length(cell_types)),
         rep("Resistant", length(cell_types))),
target = c("Sensitive","Resistant",
         cell_types,
        cell_types),
value  = c(0.097463561, 0.902536439,
         c(0, 0.946248004, 0.043108036, 0.01064396, 0,
          0, 0),
       c(0.016264368, 0.437068966, 0.038448276, 0.069310345, 0.390804598,
        0.02091954, 0.027183908)
))

#to make a sankey plot instead with teh sankey plot code 
nodesStage3t <- data.frame(
name = c("Stage 3",           # layer 1
        "Sensitive", "Resistant",  # layer 2
       cell_types),
group = c(1,  # Untreated
        2, 2,  # Sensitive_u, Resistant_u
        rep(3, 7)),
stringsAsFactors = FALSE
)

# Make group a factor with desired vertical order
nodesStage3t$group <- factor(nodesStage3t$group, levels = 1:5)

# Add ID columns to links (0-based indices)
stage3t_all$IDsource <- match(stage3t_all$source, nodesStage3t$name) - 1
stage3t_all$IDtarget <- match(stage3t_all$target, nodesStage3t$name) - 1

pStage3t <- sankeyNetwork(Links = stage3t_all, Nodes = nodesStage3t,
                       Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name", NodeGroup = "group", iterations = 0,
                     sinksRight=FALSE, fontSize = 0, nodeWidth = 15, nodePadding = 30, margin = 5)
pStage3t

htmlwidgets::saveWidget(pStage3t, "sankey_stage3t.html")

#Stage 4
cell_types <- c("B cell","T cell","NK cell","Macrophages","Fibroblasts",
              "Endothelial","Epithelial")
stage4t_all <- data.frame(
source = c("Stage 4","Stage 4",
          rep("Sensitive", length(cell_types)),
         rep("Resistant", length(cell_types))),
target = c("Sensitive","Resistant",
         cell_types,
        cell_types),
value  = c(0.446449553 ,0.553550447,
         c(0.02591268, 0.748431068, 0.055941696, 0.080977124, 0.015048249,
          0.009987179, 0.063702004),
       c(0.019048656, 0.310003265, 0.029715903, 0.200772831, 0.134810058,
        0.012354414, 0.293294873)
))

#to make a sankey plot instead with teh sankey plot code 
nodesStage4t <- data.frame(
name = c("Stage 4",           # layer 1
        "Sensitive", "Resistant",  # layer 2
        cell_types),
group = c(1,  # Untreated
        2, 2,  # Sensitive_u, Resistant_u
       rep(3, 7)),
stringsAsFactors = FALSE
)

# Make group a factor with desired vertical order
nodesStage4t$group <- factor(nodesStage4t$group, levels = 1:5)

# Add ID columns to links (0-based indices)
stage4t_all$IDsource <- match(stage4t_all$source, nodesStage4t$name) - 1
stage4t_all$IDtarget <- match(stage4t_all$target, nodesStage4t$name) - 1

pStage4t <- sankeyNetwork(Links = stage4t_all, Nodes = nodesStage4t,
                        Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name", NodeGroup = "group", iterations = 0,
                      sinksRight=FALSE, fontSize = 0, nodeWidth = 15, nodePadding = 30, margin = 5)
pStage4t

htmlwidgets::saveWidget(pStage4t, "sankey_stage4t.html")


###########################################
###Pathway analysis: resistant vs all
###########################################
outdir4 <- '/results/pathway-analysis_RvsAll/'
dir.create(outdir4)

Idents(ovcaH) <- "predicted.response"
#OR get DEG between resistant and all other cells in dataset 
resistantvsall.markers <- FindMarkers(ovcaH, ident.1 = "resistant")
resistantvsall.markers$gene <- rownames(resistantvsall.markers)
head(resistantvsall.markers) #16460 genes
#            p_val avg_log2FC pct.1 pct.2 p_val_adj     gene
# CDC42SE2     0  -3.109975 0.239 0.712         0 CDC42SE2
# CELF2        0  -2.680689 0.308 0.779         0    CELF2
# SMCHD1       0  -2.506876 0.330 0.798         0   SMCHD1
# FYN          0  -3.126504 0.257 0.721         0      FYN
# CNOT6L       0  -3.514332 0.145 0.592         0   CNOT6L
# TXNIP        0  -2.484243 0.296 0.736         0    TXNIP
write.csv(resistantvsall.markers, paste(outdir4,'differentially_expressed_genes-HigherResistantVSLowerSensitiveANDNormal_scRNAseq.csv'))

sig_resvsall_genes <- resistantvsall.markers %>% filter(p_val_adj <= 0.05) #14988
write.csv(sig_resvsall_genes, paste0(outdir4,'differentially_expressed_genes-SignificantPval<0.05-HigherResistantVSLowerSensitiveANDNormal_scRNAseq.csv'))

####UPREGULATED PATHWAYS 
#we want significant genes with values (padj <= 0.05, log2FC of > or < 0) and background genes (every gene in the differential expression analysis)
sig_resvsall_path_up <- sig_resvsall_genes %>% filter(avg_log2FC > 1) #771 genes upregulated and significant

#need to conver genesymbol to ENTREZID
sig_resvsall_path_up_map <- bitr(geneID = sig_resvsall_path_up$gene,
                                 fromType = "SYMBOL",
                                 toType = "ENTREZID",
                                 OrgDb = "org.Hs.eg.db") #now 683
# 11.41% of input gene IDs are fail to map..., not ideal but will work for now
head(sig_resvsall_path_up_map) #these would be the correct differentially expressed genes
# SYMBOL ENTREZID
# 1   WFDC2    10406
# 2    SLPI     6590
# 3   IGHGP     3505
# 4 RARRES2     5919
# 5    KRT8     3856
# 6   KRT18     3875

#for the background, we would need all the genes from the scRNAseq object which are 33612
background_resvsall_path_up_map <- bitr(geneID = rownames(ovcaH),
                                        fromType = "SYMBOL",
                                        toType = "ENTREZID",
                                        OrgDb = "org.Hs.eg.db")
#34.7% of input gene IDs are fail to map..., now 21950 genes

##pathways from msigdbr
#can pick from the following sets, look at H and C2
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

m_Hpathways <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
head(m_Hpathways)

# m_C2pathways <- msigdbr(species = "Homo sapiens", category = "C2") %>%
#   dplyr::select(gs_name, entrez_gene)
# head(m_C2pathways)

#H pathways
msig_Hpathways_resvsall_upgenes <- enricher(sig_resvsall_path_up_map$ENTREZID, TERM2GENE=m_Hpathways,
                                            universe = background_resvsall_path_up_map$ENTREZID)
head(msig_Hpathways_resvsall_upgenes)
# ID                                Description GeneRatio
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION    44/209
# HALLMARK_COAGULATION                                             HALLMARK_COAGULATION                       HALLMARK_COAGULATION    23/209
# HALLMARK_ESTROGEN_RESPONSE_LATE                       HALLMARK_ESTROGEN_RESPONSE_LATE            HALLMARK_ESTROGEN_RESPONSE_LATE    25/209
# HALLMARK_KRAS_SIGNALING_UP                                 HALLMARK_KRAS_SIGNALING_UP                 HALLMARK_KRAS_SIGNALING_UP    25/209
# HALLMARK_ESTROGEN_RESPONSE_EARLY                     HALLMARK_ESTROGEN_RESPONSE_EARLY           HALLMARK_ESTROGEN_RESPONSE_EARLY    24/209
# HALLMARK_APICAL_JUNCTION                                     HALLMARK_APICAL_JUNCTION                   HALLMARK_APICAL_JUNCTION    22/209       2.146485  5.891280

msig_Hpathways_resvsall_upgenes_df <- msig_Hpathways_resvsall_upgenes@result #41 pathways
write.csv(msig_Hpathways_resvsall_upgenes_df, paste0(outdir4,'mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq.csv'))

pdf(paste0(outdir4,'barPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq.pdf'),width = 12, height = 12)
barplot(msig_Hpathways_resvsall_upgenes, showCategory = 20)
dev.off()

pdf(paste0(outdir4,'barPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq-short.pdf'),width = 12, height = 8)
barplot(msig_Hpathways_resvsall_upgenes)
dev.off()

pdf(paste0(outdir4,'dotPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq.pdf'),width = 12, height = 12)
dotplot(msig_Hpathways_resvsall_upgenes)
dev.off()

####DOWNREGULATED PATHWAYS 
#we want significant genes with values (padj <= 0.05, log2FC of > or < 0) and background genes (every gene in the differential expression analysis)
sig_resvsall_path_down <- sig_resvsall_genes %>% filter(avg_log2FC < -1) #5761 genes upregulated and significant

#need to conver genesymbol to ENTREZID
sig_resvsall_path_down_map <- bitr(geneID = sig_resvsall_path_down$gene,
                                   fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = "org.Hs.eg.db") #now 4746
# 17.62% of input gene IDs are fail to map..., not ideal but will work for now
head(sig_resvsall_path_down_map) #these would be the correct differentially expressed genes
#   SYMBOL ENTREZID
# 1 CDC42SE2    56990
# 2    CELF2    10659
# 3   SMCHD1    23347
# 4      FYN     2534
# 5   CNOT6L   246175
# 6    TXNIP    10628

#for the background, we would need all the genes from the scRNAseq object which are 33612
background_resvsall_path_down_map <- bitr(geneID = rownames(ovcaH),
                                          fromType = "SYMBOL",
                                          toType = "ENTREZID",
                                          OrgDb = "org.Hs.eg.db")
#34.7% of input gene IDs are fail to map..., now 21950 genes

##pathways from msigdbr
#can pick from the following sets, look at H and C2
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

# m_Hpathways <- msigdbr(species = "Homo sapiens", category = "H") %>%
#   dplyr::select(gs_name, entrez_gene)
# head(m_Hpathways)

# m_C2pathways <- msigdbr(species = "Homo sapiens", category = "C2") %>%
#   dplyr::select(gs_name, entrez_gene)
# head(m_C2pathways)

#H pathways
msig_Hpathways_resvsall_downgenes <- enricher(sig_resvsall_path_down_map$ENTREZID, TERM2GENE=m_Hpathways,
                                              universe = background_resvsall_path_down_map$ENTREZID)
head(msig_Hpathways_resvsall_downgenes)
# ID                        Description GeneRatio  BgRatio RichFactor
# HALLMARK_ALLOGRAFT_REJECTION             HALLMARK_ALLOGRAFT_REJECTION       HALLMARK_ALLOGRAFT_REJECTION   80/1020 195/4316  0.4102564
# HALLMARK_IL2_STAT5_SIGNALING             HALLMARK_IL2_STAT5_SIGNALING       HALLMARK_IL2_STAT5_SIGNALING   75/1020 198/4316  0.3787879
# HALLMARK_UV_RESPONSE_DN                       HALLMARK_UV_RESPONSE_DN            HALLMARK_UV_RESPONSE_DN   57/1020 143/4316  0.3986014
# HALLMARK_MITOTIC_SPINDLE                     HALLMARK_MITOTIC_SPINDLE           HALLMARK_MITOTIC_SPINDLE   70/1020 198/4316  0.3535354
# HALLMARK_INTERFERON_GAMMA_RESPONSE HALLMARK_INTERFERON_GAMMA_RESPONSE HALLMARK_INTERFERON_GAMMA_RESPONSE   68/1020 196/4316  0.3469388
# HALLMARK_PI3K_AKT_MTOR_SIGNALING     HALLMARK_PI3K_AKT_MTOR_SIGNALING   HALLMARK_PI3K_AKT_MTOR_SIGNALING   40/1020 105/4316  0.3809524

msig_Hpathways_resvsall_downgenes_df <- msig_Hpathways_resvsall_downgenes@result #50 pathways
write.csv(msig_Hpathways_resvsall_downgenes_df, paste0(figdir1,'mSigHpathways-DEG0.05PvalueANDhigher1FC-DOWNREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq.csv'))

pdf(paste0(outdir4,'barPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-DOWNREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq.pdf'),width = 12, height = 12)
barplot(msig_Hpathways_resvsall_downgenes, showCategory = 20)
dev.off()

pdf(paste0(outdir4,'barPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-DOWNREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq-short.pdf'),width = 12, height = 8)
barplot(msig_Hpathways_resvsall_downgenes)
dev.off()

pdf(paste0(outdir4,'dotPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-DOWNREGULATED_HigherResistantVSLowerSensitiveANDNormal_scRNAseq.pdf'),width = 12, height = 12)
dotplot(msig_Hpathways_resvsall_downgenes)
dev.off()

###########################################
###Extract data for COMET: resistant vs all
###########################################
outdir5 <- '/results/COMETinputs'
dir.created(outdir5)

Idents(ovcaH) <- "predicted.response"
#OR get DEG between resistant and all other cells in dataset 
resistantvsall.markers <- FindMarkers(ovcaH, ident.1 = "resistant")
resistantvsall.markers$gene <- rownames(resistantvsall.markers)

#only use upregulated and significant genes 
sig_resvsall_genes <- resistantvsall.markers %>% filter(p_val_adj <= 0.05) #14988
sig_resvsall_path_up <- sig_resvsall_genes %>% filter(avg_log2FC > 1) #771

#get normalized counts for the entire object
ovcaH_counts_matrix <- GetAssayData(ovcaH, slot = "data", assay = "SCT")

####make the resistant cells and the rest of the cells be the same number
res_cells <- ovcaH$predicted.response == "resistant"
res_ids <- which(res_cells)
#downsample to 65000/2=32500 actually 32499 just in case
res_ids_downsampled <- sample(res_ids, 32499, replace = FALSE) #32499 cells now
#downsample the logical object 
res_cells_downsampled <- res_cells
res_cells_downsampled[setdiff(res_ids, res_ids_downsampled)] <- FALSE

resistant.downsampled.data <- ovcaH_counts_matrix[,res_cells_downsampled] #54772 cells

ovcaH_embeddings <- Embeddings(ovcaH, reduction = "umap")
res.downsampled.umap <- ovcaH_embeddings[res_cells_downsampled,]

#get sensitive and normal data and visualization (umap)
senandnor_cells <- ovcaH$predicted.response != "resistant"
senandnor.data <- ovcaH_counts_matrix[,senandnor_cells] # 110524 cells

senandnor.umap <- ovcaH_embeddings[senandnor_cells,]

#maximum number of cells for COMET is 65K, so the sen and nor have to be downsampled
nSamples.downsampled <- 64998 - dim(resistant.downsampled.data)[2] #32499 samples that are not resistant to add to the dataset
#might be better to have half resistant and half of the rest 

#assign metadata arbitrarily, resistant becomes 1. All others are 0
resistant.downsampled.cluster <- data.frame("X" = rep(1, dim(resistant.downsampled.data)[2]))
senandnor.downsampled.cluster <- data.frame("X" = rep(0, nSamples.downsampled))

#convert to data frames and downsample
resistant.data.df <- data.frame(resistant.downsampled.data[sig_resvsall_path_up$gene,])
senandnor.data.df <- data.frame(senandnor.data[sig_resvsall_path_up$gene,])

cellsKeep.downsampled <- sample(length(senandnor.data.df),nSamples.downsampled)
senandnor.data.downsampled <- senandnor.data.df[,cellsKeep.downsampled]
senandnor.umap.downsampled <- data.frame(senandnor.umap[cellsKeep.downsampled,])

markers.df <- cbind(resistant.data.df,senandnor.data.downsampled) #the colnames got a little changed
umap.df <- rbind(res.downsampled.umap,senandnor.umap.downsampled)
clusters.df <- rbind(resistant.downsampled.cluster,senandnor.downsampled.cluster)
row.names(clusters.df) <- row.names(umap.df) 
colnames(markers.df) <- row.names(umap.df) #match the colnames in markers to the rownames in umap, they don't automatically match
#write outfiles with the correct format for COMET
write.table(markers.df, paste0(outdir5,"ds_markers.txt"), sep="\t",row.names=T, quote = F)
write.table(umap.df, paste0(outdir5,"ds_umap.txt"), sep="\t",row.names=T, col.names = F, quote = F)
write.table(clusters.df, paste0(outdir5,"ds_clusters.txt"), sep="\t",row.names=T,col.names = F, quote = F)
write.table(sig_resvsall_path_up$gene, paste0(outdir5,"ds_genes.txt"), sep="\t",row.names=F, col.names = F, quote = F)

############################################################################################
####Pseudotime analysis of clusters 2 (epithelial cluster) and 3 (fibroblasts cluster) using Monocle
############################################################################################
outdir6 <- '/results/pseudotime/'
dir.create(outdir6)
#Later portions of the monocle code have been commented out because it requires user input and it won't run without the user input

####PSEUDOTIME ON CLUSTER 2 - EPITHELIAL 
clus2.ovcaH <- subset(ovcaH, SCT_snn_res.0.025 %in% "2")
dim(clus2.ovcaH) #[1] 33612 27979

#### Try to combine the sensitive and normal cells to the same name
clus2.ovcaH$predicted.response[clus2.ovcaH$predicted.response == 'normal'] <- 'other'
clus2.ovcaH$predicted.response[clus2.ovcaH$predicted.response == 'sensitive'] <- 'other'
table(clus2.ovcaH$predicted.response)
# other resistant 
# 11177     16802
clus2.resvsall.ovcaH <- clus2.ovcaH
#### CLUSTERING AND CLASSIFYING CELLS
### Step 1: first step is to get the data in a CellDataSet 
#get the different elements for cds
exprs_matrix.clus2.resvsall <- GetAssayData(clus2.resvsall.ovcaH, layer = "data", assay = "SCT") #set to RNA and counts to use raw counts instead
cell_metadata.clus2.resvsall <- clus2.resvsall.ovcaH@meta.data
gene_metadata.clus2.resvsall <- data.frame(gene_short_name = rownames(exprs_matrix.clus2.resvsall))
rownames(gene_metadata.clus2.resvsall) <- rownames(exprs_matrix.clus2.resvsall)

dim(exprs_matrix.clus2.resvsall) #33612 27979
dim(cell_metadata.clus2.resvsall) #27979    46
dim(gene_metadata.clus2.resvsall) #33612     1

#create the cds 
cds.clus2.resvsall <- new_cell_data_set(exprs_matrix.clus2.resvsall,
                               cell_metadata = cell_metadata.clus2.resvsall,
                               gene_metadata = gene_metadata.clus2.resvsall)
# add UMAP coordinates and clusters 
clus2.resvsall_umap_matrix <- Embeddings(clus2.resvsall.ovcaH, reduction = "umap")
reducedDims(cds.clus2.resvsall)$UMAP <- clus2.resvsall_umap_matrix

clus2.resvsall_clusters <- clus2.resvsall.ovcaH$SCT_snn_res.0.025
clus2.resvsall_clusters <- clus2.resvsall_clusters[colnames(cds.clus2.resvsall)] #ensure the order is the same 
colData(cds.clus2.resvsall)$clus2.resvsall_clusters <- as.character(clus2.resvsall_clusters) # add to monocle object 
cds.clus2.resvsall@clusters$UMAP$clusters <- factor(clus2.resvsall_clusters) # add clusters in the monocle inside's structure 

### Step 2: pre-process the data, run Principal component analysis 
cds.clus2.resvsall <- preprocess_cds(cds.clus2.resvsall, num_dim = 100)
pdf(paste(outdir6, "Monocle-PCA-Variance_clus2.resvsall.pdf", sep = ""), width = 28, height = 18)
plot_pc_variance_explained(cds.clus2.resvsall)
dev.off()

### Step 3: reduce dimensionality and visualize 
#Don't think that is necessary since UMAP and clusters were already added 
pdf(paste(outdir6, "Monocle-UMAP_clus2.resvsall_byresponse.pdf", sep = ""), width = 28, height = 18)
plot_cells(cds.clus2.resvsall, color_cells_by = "predicted.response", label_cell_groups = FALSE)
dev.off()

pdf(paste(outdir6, "Monocle-UMAP_clus2.resvsall_bypredicted.cellType.pdf", sep = ""), width = 28, height = 18)
plot_cells(cds.clus2.resvsall, color_cells_by = "predicted.cellType", label_cell_groups = FALSE)
dev.off()

### Step 5: group cells into clusters 
#Don't think it is necessary, already added the previous clusters 
cds.clus2.resvsall <- cluster_cells(cds.clus2.resvsall, resolution = 1e-6) #only 1 cluster
#default parameters should be reduction_method = "UMAP", cluster_method = "leiden"

pdf(paste(outdir6, "Monocle-UMAP_clus2.resvsall_byMonocleCluster.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus2.resvsall, color_cells_by = "partition", label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off()
 

### Step 4: learn trajectory graph 
cds.clus2.resvsall <- learn_graph(cds.clus2.resvsall)
# |==========================================================================================================| 100%
# |==========================================================================================================| 100%
# |==========================================================================================================| 100%

#gray circles are leafs which correspond to a different outcome of the trajectory/ cell fate 
#black circles are nodes in which cells travel to one of the outcomes 
pdf(paste(outdir6, "Monocle-UMAP-trajectories_clus2.resvsall_bypredicted.response.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus2.resvsall, color_cells_by = "predicted.response", label_cell_groups = FALSE, graph_label_size = 5)
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-trajectories_clus2.resvsall_bypredicted.cellType.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus2.resvsall, color_cells_by = "predicted.cellType", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 5)
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-trajectories_clus2.resvsall_byclusters0.025.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus2.resvsall, color_cells_by = "clus2.resvsall_clusters", label_cell_groups = FALSE, graph_label_size = 5)
dev.off()

### Step 5: order cells in pseudotime 
#we would pick where we want the pseudotime to start 
cds.clus2.resvsall <- order_cells(cds.clus2.resvsall)
#a shiny app would open for you to select start nodes 
#plot the cells by pseudotime 
pdf(paste(outdir6, "Monocle-UMAP-pseudotimePlot1startcell_bottomRightCorneronS_clus2.resvsall.#pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus2.resvsall, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = FALSE, graph_label_size = 4)
dev.off()
#gray cells have an infinite pseudotime since they were not reachable from the root nodes picked 
#start point can also be choosen programmatically instead of manually but SKIP for now 

### Plot genes of interest
###try to do 2 separate plots for sensitive and resistant for this pseudotime pseudotimePlot1startcell_bottomRightCorneronS_clus2.resvsall.pdf
#do this for genes SLPI (epithelial resistant) and NNMT (fibroblast resistant) 
SLPI_gene <- "SLPI"

SLPI_lineage_cds <- cds.clus2.resvsall[rowData(cds.clus2.resvsall)$gene_short_name == SLPI_gene, ]
SLPI_lineage_cds <- order_cells(SLPI_lineage_cds)

SLPI_expr_mat <- exprs(SLPI_lineage_cds)
SLPI_pt <- pseudotime(SLPI_lineage_cds)

df_SLPI <- data.frame(
pseudotime = SLPI_pt,
expression = as.numeric(SLPI_expr_mat[1, ]),
group = colData(SLPI_lineage_cds)$predicted.response
)

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_bottomRightCorneronSensitive_clus2.#resvsall_SLPI-gene-withPoints.pdf", sep = ""),width = 8, height = 4)
ggplot(df_SLPI, aes(x = pseudotime, y = expression, color = group)) +
geom_point(size = 0.8, alpha = 0.2) +
geom_smooth(se = TRUE, method = "loess", span = 0.75, linewidth = 2) +
scale_color_manual(values = c("other" = "#EB0F10",
                             "resistant" = "#006502")) +
scale_fill_manual(values = c("other" = "#EB0F10",
                            "resistant" = "#006502")) +
theme_classic() +
ylim(0,7)+
labs(x = "Pseudotime", y = "Expression")
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_bottomRightCorneronSensitive_clus2.#resvsall_SLPI-gene.pdf", sep = ""),width = 8, height = 4)
ggplot(df_SLPI, aes(x = pseudotime, y = expression, color = group)) +
 geom_smooth(se = TRUE, method = "loess", span = 0.75) +
 scale_color_manual(values = c("other" = "#EB0F10",
                               "resistant" = "#006502")) +
 theme_classic() +
 ylim(0,7)+
 labs(x = "Pseudotime", y = "Expression")
dev.off()

NNMT_gene <- "NNMT"

NNMT_lineage_cds <- cds.clus2.resvsall[rowData(cds.clus2.resvsall)$gene_short_name == NNMT_gene, ]
NNMT_lineage_cds <- order_cells(NNMT_lineage_cds)

NNMT_expr_mat <- exprs(NNMT_lineage_cds)
NNMT_pt <- pseudotime(NNMT_lineage_cds)

df_NNMT <- data.frame(
pseudotime = NNMT_pt,
expression = as.numeric(NNMT_expr_mat[1, ]),
group = colData(NNMT_lineage_cds)$predicted.response
)

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_bottomRightCorneronSensitive_clus2.#resvsall_NNMT-gene-withPoints.pdf", sep = ""),width = 8, height = 4)
ggplot(df_NNMT, aes(x = pseudotime, y = expression, color = group)) +
geom_point(size = 0.8, alpha = 0.2) +
geom_smooth(se = TRUE, method = "loess", span = 0.75, linewidth = 2) +
scale_color_manual(values = c("other" = "#EB0F10",
                             "resistant" = "#006502")) +
scale_fill_manual(values = c("other" = "#EB0F10",
                            "resistant" = "#006502")) +
theme_classic() +
ylim(0,7) +
labs(x = "Pseudotime", y = "Expression")
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_bottomRightCorneronSensitive_clus2.#resvsall_NNMT-gene.pdf", sep = ""),width = 8, height = 4)
ggplot(df_NNMT, aes(x = pseudotime, y = expression, color = group)) +
 geom_smooth(se = TRUE, method = "loess", span = 0.75) +
 scale_color_manual(values = c("other" = "#EB0F10",
                              "resistant" = "#006502")) +
theme_classic() +
ylim(0,7) +
labs(x = "Pseudotime", y = "Expression")
dev.off()


#### PSEUDOTIME ON CLUSTER 3 - FIBROBLAST
clus3.ovcaH <- subset(ovcaH, SCT_snn_res.0.025 %in% "3")
dim(clus3.ovcaH) #[1] 33612 26189


pdf(paste0(outdir6, 'OvCa-ALL_UMAPbycellReponse-Resistantvsother_noleyend.pdf'),width = 15, height = 10)
DimPlot(clus3.ovcaH, group.by = "predicted.response",cols = c("other" = "#EB0F10", "resistant" = "#006502"))
dev.off()

#### Try to combine the sensitive and normal cells to the same name
clus3.ovcaH$predicted.response[clus3.ovcaH$predicted.response == 'normal'] <- 'other'
clus3.ovcaH$predicted.response[clus3.ovcaH$predicted.response == 'sensitive'] <- 'other'

clus3.resvsall.ovcaH <- clus3.ovcaH
#### CLUSTERING AND CLASSIFYING CELLS
### Step 1: first step is to get the data in a CellDataSet 
#get the different elements for cds
exprs_matrix.clus3.resvsall <- GetAssayData(clus3.resvsall.ovcaH, layer = "data", assay = "SCT") #set to RNA and counts to use raw counts instead
cell_metadata.clus3.resvsall <- clus3.resvsall.ovcaH@meta.data
gene_metadata.clus3.resvsall <- data.frame(gene_short_name = rownames(exprs_matrix.clus3.resvsall))
rownames(gene_metadata.clus3.resvsall) <- rownames(exprs_matrix.clus3.resvsall)

dim(exprs_matrix.clus3.resvsall) #33612 26189
dim(cell_metadata.clus3.resvsall) #28189    46
dim(gene_metadata.clus3.resvsall) #33612     1

#create the cds 
cds.clus3.resvsall <- new_cell_data_set(exprs_matrix.clus3.resvsall,
                                        cell_metadata = cell_metadata.clus3.resvsall,
                                        gene_metadata = gene_metadata.clus3.resvsall)
#add UMAP coordinates and clusters 
clus3.resvsall_umap_matrix <- Embeddings(clus3.resvsall.ovcaH, reduction = "umap")
reducedDims(cds.clus3.resvsall)$UMAP <- clus3.resvsall_umap_matrix

clus3.resvsall_clusters <- clus3.resvsall.ovcaH$SCT_snn_res.0.025
clus3.resvsall_clusters <- clus3.resvsall_clusters[colnames(cds.clus3.resvsall)] #ensure the order is the same 
colData(cds.clus3.resvsall)$clus3.resvsall_clusters <- as.character(clus3.resvsall_clusters) #add to monocle object 
cds.clus3.resvsall@clusters$UMAP$clusters <- factor(clus3.resvsall_clusters) #add clusters in the monocle inside's structure 

### Step 2: pre-process the data, run Principal component analysis 
cds.clus3.resvsall <- preprocess_cds(cds.clus3.resvsall, num_dim = 100)
pdf(paste(outdir6, "Monocle-PCA-Variance_clus3.resvsall.pdf", sep = ""), width = 28, height = 18)
plot_pc_variance_explained(cds.clus3.resvsall)
dev.off()

### Step 3: reduce dimensionality and visualize 
#Don't think that is necessary since UMAP and clusters were already added 
pdf(paste(outdir6, "Monocle-UMAP_clus3.resvsall_byresponse.pdf", sep = ""), width = 28, height = 18)
plot_cells(cds.clus3.resvsall, color_cells_by = "predicted.response", label_cell_groups = FALSE)
dev.off()

pdf(paste(outdir6, "Monocle-UMAP_clus3.resvsall_byclusters0.025.pdf", sep = ""), width = 28, height = 18)
plot_cells(cds.clus3.resvsall, color_cells_by = "clus3.resvsall_clusters", label_cell_groups = FALSE)
dev.off()

pdf(paste(outdir6, "Monocle-UMAP_clus3.resvsall_bypredicted.cellType.pdf", sep = ""), width = 28, height = 18)
plot_cells(cds.clus3.resvsall, color_cells_by = "predicted.cellType", label_cell_groups = FALSE)
dev.off()

### Step 5: group cells into clusters 
#Don't think it is necessary, already added the previous clusters 
cds.clus3.resvsall <- cluster_cells(cds.clus3.resvsall, resolution = 1e-6) #only 1 cluster
#default parameters should be reduction_method = "UMAP", cluster_method = "leiden"

pdf(paste(outdir6, "Monocle-UMAP_clus3.resvsall_byMonocleCluster.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus3.resvsall, color_cells_by = "partition", label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off()

### Step 4: learn trajectory graph 
cds.clus3.resvsall <- learn_graph(cds.clus3.resvsall)
# |==========================================================================================================| 100%
# |==========================================================================================================| 100%
# |==========================================================================================================| 100%

#gray circles are leafs which correspond to a different outcome of the trajectory/ cell fate 
#black circles are nodes in which cells travel to one of the outcomes 
pdf(paste(outdir6, "Monocle-UMAP-trajectories_clus3.resvsall_bypredicted.response.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus3.resvsall, color_cells_by = "predicted.response", label_cell_groups = FALSE, graph_label_size = 5)
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-trajectories_clus3.resvsall_bypredicted.cellType.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus3.resvsall, color_cells_by = "predicted.cellType", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 5)
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-trajectories_clus3.resvsall_byclusters0.025.pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus3.resvsall, color_cells_by = "clus3.resvsall_clusters", label_cell_groups = FALSE, graph_label_size = 5)
dev.off()

### Step 5: order cells in pseudotime 
#we would pick where we want the pseudotime to start 
cds.clus3.resvsall <- order_cells(cds.clus3.resvsall)
#a shiny app would open for you to select start nodes 
#plot the cells by pseudotime
pdf(paste(outdir6, "Monocle-UMAP-pseudotimePlot1startcell_leftOutonSensitive_clus3.resvsall.#pdf", sep = ""), width = 14, height = 9)
plot_cells(cds.clus3.resvsall, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = FALSE, graph_label_size = 4)
dev.off()

###try to do 2 separate plots for sensitive and resistant for this pseudotime 
#pseudotimePlot1startcell_leftOutonSensitive_clus3.resvsall.pdf
#do this for genes SLPI (epithelial resistant) and NNMT (fibroblast resistant) 
SLPI_gene <- "SLPI"

SLPI_lineage_cds <- cds.clus3.resvsall[rowData(cds.clus3.resvsall)$gene_short_name == SLPI_gene, ]
SLPI_lineage_cds <- order_cells(SLPI_lineage_cds)

SLPI_expr_mat <- exprs(SLPI_lineage_cds)
SLPI_pt <- pseudotime(SLPI_lineage_cds)

df_SLPI <- data.frame(
 pseudotime = SLPI_pt,
 expression = as.numeric(SLPI_expr_mat[1, ]),
 group = colData(SLPI_lineage_cds)$predicted.response
)

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_leftOutonSensitive_clus3.#resvsall_SLPI-gene-withPoints.pdf", sep = ""),width = 8, height = 4)
ggplot(df_SLPI, aes(x = pseudotime, y = expression, color = group)) +
 geom_point(size = 0.8, alpha = 0.2) +
 geom_smooth(se = TRUE, method = "loess", span = 0.75, linewidth = 2) +
 scale_color_manual(values = c("other" = "#EB0F10",
                              "resistant" = "#006502")) +
scale_fill_manual(values = c("other" = "#EB0F10",
                            "resistant" = "#006502")) +
theme_classic() +
ylim(0,6)+
labs(x = "Pseudotime", y = "Expression")
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_leftOutonSensitive_clus3.#resvsall_SLPI-gene.pdf", sep = ""),width = 8, height = 4)
ggplot(df_SLPI, aes(x = pseudotime, y = expression, color = group)) +
geom_smooth(se = TRUE, method = "loess", span = 0.75) +
scale_color_manual(values = c("other" = "#EB0F10",
                             "resistant" = "#006502")) +
theme_classic() +
ylim(0,6)+
labs(x = "Pseudotime", y = "Expression")
dev.off()

NNMT_gene <- "NNMT"

NNMT_lineage_cds <- cds.clus3.resvsall[rowData(cds.clus3.resvsall)$gene_short_name == NNMT_gene, ]
NNMT_lineage_cds <- order_cells(NNMT_lineage_cds)

NNMT_expr_mat <- exprs(NNMT_lineage_cds)
NNMT_pt <- pseudotime(NNMT_lineage_cds)

df_NNMT <- data.frame(
 pseudotime = NNMT_pt,
 expression = as.numeric(NNMT_expr_mat[1, ]),
 group = colData(NNMT_lineage_cds)$predicted.response
)

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_leftOutonSensitive_clus3.#resvsall_NNMT-gene-withPoints.pdf", sep = ""),width = 8, height = 4)
ggplot(df_NNMT, aes(x = pseudotime, y = expression, color = group)) +
 geom_point(size = 0.8, alpha = 0.2) +
geom_smooth(se = TRUE, method = "loess", span = 0.75, linewidth = 2) +
scale_color_manual(values = c("other" = "#EB0F10",
                             "resistant" = "#006502")) +
scale_fill_manual(values = c("other" = "#EB0F10",
                            "resistant" = "#006502")) +
theme_classic() +
ylim(0,6) +
labs(x = "Pseudotime", y = "Expression")
dev.off()

pdf(paste(outdir6, "Monocle-UMAP-expressionAlongPseudotime_1startcell_leftOutonSensitive_clus3.#resvsall_NNMT-gene.pdf", sep = ""),width = 8, height = 4)
ggplot(df_NNMT, aes(x = pseudotime, y = expression, color = group)) +
 geom_smooth(se = TRUE, method = "loess", span = 0.75) +
 scale_color_manual(values = c("other" = "#EB0F10",
                              "resistant" = "#006502")) +
theme_classic() +
ylim(0,6) +
labs(x = "Pseudotime", y = "Expression")
dev.off()


















