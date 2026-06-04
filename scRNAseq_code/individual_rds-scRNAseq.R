#Load libraries
library('ggplot2')
library('Seurat')
library('patchwork')
library('forcats')
library('tidyverse')
library('data.table')
library('future')
library('harmony')
library('clustree')
library('plyr')
library('hdf5r')
library('SingleR')
library(BiocParallel)
library(ArrayExpress)
library('celldex')

#################################################
##  CREATING INDIVIDUAL RDS FILES FOR EACH SAMPLE 
#################################################
outdir <- '/results/scRNAseq_individual/'
dir.create(outdir)

enddir <-"-Unprocessed-RDS_2025-04-14.rds"

###GSE15460#####################################

###T59_bam
T59file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE154600/T59_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
T59.counts <- Read10X(data.dir=T59file)

# create a Seurat object based on the scRNA-seq data
T59 <- CreateSeuratObject(counts=T59.counts, project="T59_GSE154600")
#add metadata
T59@meta.data[["response"]] <- "resistant"
T59@meta.data[["tissue"]] <- "omentum"
T59@meta.data[["stage"]] <- "4"
T59@meta.data[["grade"]] <- "3"
T59@meta.data[["diagnosis"]] <- "HGSOC"
T59@meta.data[["age"]] <- "na"
T59@meta.data[["treatment"]] <- "treated"

#save rds file
saveRDS(T59, paste(outdir,"T59_GSE154600",enddir,sep=""))


###T61_bam
T61file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE154600/T61_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
T61.counts <- Read10X(data.dir=T61file)

# create a Seurat object based on the scRNA-seq data
T61 <- CreateSeuratObject(counts=T61.counts, project="T61_GSE154600")
#add metadata
T61@meta.data[["response"]] <- "refractory"
T61@meta.data[["tissue"]] <- "omentum"
T61@meta.data[["stage"]] <- "3c"
T61@meta.data[["grade"]] <- "3"
T61@meta.data[["diagnosis"]] <- "HGSOC"
T61@meta.data[["age"]] <- "na"
T61@meta.data[["treatment"]] <- "treated"

#save rds file
saveRDS(T61, paste(outdir,"T61_GSE154600",enddir,sep=""))


###T76_bam
T76file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE154600/T76_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
T76.counts <- Read10X(data.dir=T76file)

# create a Seurat object based on the scRNA-seq data
T76 <- CreateSeuratObject(counts=T76.counts, project="T76_GSE154600")
#add metadata
T76@meta.data[["response"]] <- "resistant"
T76@meta.data[["tissue"]] <- "omentum"
T76@meta.data[["stage"]] <- "4b"
T76@meta.data[["grade"]] <- "3"
T76@meta.data[["diagnosis"]] <- "HGSOC"
T76@meta.data[["age"]] <- "na"
T76@meta.data[["treatment"]] <- "treated"
#save rds file
saveRDS(T76, paste(outdir,"T76_GSE154600",enddir,sep=""))


###T89_bam
T89file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE154600/T89_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
T89.counts <- Read10X(data.dir=T89file)

# create a Seurat object based on the scRNA-seq data
T89 <- CreateSeuratObject(counts=T89.counts, project="T89_GSE154600")
#add metadata
T89@meta.data[["response"]] <- "sensitive"
T89@meta.data[["tissue"]] <- "omentum"
T89@meta.data[["stage"]] <- "4"
T89@meta.data[["grade"]] <- "3"
T89@meta.data[["diagnosis"]] <- "HGSOC"
T89@meta.data[["age"]] <- "na"
T89@meta.data[["treatment"]] <- "treated"

#save rds file
saveRDS(T89, paste(outdir,"T89_GSE154600",enddir,sep=""))


###T90_bam
T90file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE154600/T90_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
T90.counts <- Read10X(data.dir=T90file)

# create a Seurat object based on the scRNA-seq data
T90 <- CreateSeuratObject(counts=T90.counts, project="T90_GSE154600")
#add metadata
T90@meta.data[["response"]] <- "sensitive"
T90@meta.data[["tissue"]] <- "omentum"
T90@meta.data[["stage"]] <- "4a"
T90@meta.data[["grade"]] <- "3"
T90@meta.data[["diagnosis"]] <- "HGSOC"
T90@meta.data[["age"]] <- "na"
T90@meta.data[["treatment"]] <- "treated"
#save rds file
saveRDS(T90, paste(outdir,"T90_GSE154600",enddir,sep=""))


###GSE184880#####################################

###SRX12379603_Normal1
Normal1file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379603_Normal1/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Normal1.counts <- Read10X(data.dir=Normal1file)

# create a Seurat object based on the scRNA-seq data
Normal1 <- CreateSeuratObject(counts=Normal1.counts, project="Normal1_GSE184880")
#add metadata
Normal1@meta.data[["response"]] <- "na"
Normal1@meta.data[["tissue"]] <- "normal"
Normal1@meta.data[["stage"]] <- "na"
Normal1@meta.data[["grade"]] <- "na"
Normal1@meta.data[["diagnosis"]] <- "normal"
Normal1@meta.data[["age"]] <- "55"
Normal1@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(Normal1, paste(outdir,"Normal1_GSE184880",enddir,sep=""))


###SRX12379604_Normal2
Normal2file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379604_Normal2/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Normal2.counts <- Read10X(data.dir=Normal2file)

# create a Seurat object based on the scRNA-seq data
Normal2 <- CreateSeuratObject(counts=Normal2.counts, project="Normal2_GSE184880")
#add metadata
Normal2@meta.data[["response"]] <- "na"
Normal2@meta.data[["tissue"]] <- "normal"
Normal2@meta.data[["stage"]] <- "na"
Normal2@meta.data[["grade"]] <- "na"
Normal2@meta.data[["diagnosis"]] <- "normal"
Normal2@meta.data[["age"]] <- "47"
Normal2@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(Normal2, paste(outdir,"Normal2_GSE184880",enddir,sep=""))


###SRX12379605_Normal3
Normal3file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379605_Normal3/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Normal3.counts <- Read10X(data.dir=Normal3file)

# create a Seurat object based on the scRNA-seq data
Normal3 <- CreateSeuratObject(counts=Normal3.counts, project="Normal3_GSE184880")
#add metadata
Normal3@meta.data[["response"]] <- "na"
Normal3@meta.data[["tissue"]] <- "normal"
Normal3@meta.data[["stage"]] <- "na"
Normal3@meta.data[["grade"]] <- "na"
Normal3@meta.data[["diagnosis"]] <- "normal"
Normal3@meta.data[["age"]] <- "46"
Normal3@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(Normal3, paste(outdir,"Normal3_GSE184880",enddir,sep=""))


###SRX12379606_Normal4
Normal4file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379606_Normal4/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Normal4.counts <- Read10X(data.dir=Normal4file)

# create a Seurat object based on the scRNA-seq data
Normal4 <- CreateSeuratObject(counts=Normal4.counts, project="Normal4_GSE184880")
#add metadata
Normal4@meta.data[["response"]] <- "na"
Normal4@meta.data[["tissue"]] <- "normal"
Normal4@meta.data[["stage"]] <- "na"
Normal4@meta.data[["grade"]] <- "na"
Normal4@meta.data[["diagnosis"]] <- "normal"
Normal4@meta.data[["age"]] <- "51"
Normal4@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(Normal4, paste(outdir,"Normal4_GSE184880",enddir,sep=""))


###SRX12379607_Normal5 
Normal5file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379607_Normal5/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Normal5.counts <- Read10X(data.dir=Normal5file)

# create a Seurat object based on the scRNA-seq data
Normal5 <- CreateSeuratObject(counts=Normal5.counts, project="Normal5_GSE184880")
#add metadata
Normal5@meta.data[["response"]] <- "na"
Normal5@meta.data[["tissue"]] <- "normal"
Normal5@meta.data[["stage"]] <- "na"
Normal5@meta.data[["grade"]] <- "na"
Normal5@meta.data[["diagnosis"]] <- "normal"
Normal5@meta.data[["age"]] <- "49"
Normal5@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(Normal5, paste(outdir,"Normal5_GSE184880",enddir,sep=""))


###SRX12379608_Cancer1
Cancer1file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379608_Cancer1/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Cancer1.counts <- Read10X(data.dir=Cancer1file)

# create a Seurat object based on the scRNA-seq data
Cancer1 <- CreateSeuratObject(counts=Cancer1.counts, project="Cancer1_GSE184880")
#add metadata
Cancer1@meta.data[["response"]] <- "na"
Cancer1@meta.data[["tissue"]] <- "primary"
Cancer1@meta.data[["stage"]] <- "3b"
Cancer1@meta.data[["grade"]] <- "na"
Cancer1@meta.data[["diagnosis"]] <- "HGSOC"
Cancer1@meta.data[["age"]] <- "50"
Cancer1@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(Cancer1, paste(outdir,"Cancer1_GSE184880",enddir,sep=""))


###SRX12379609_Cancer2
Cancer2file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379609_Cancer2/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Cancer2.counts <- Read10X(data.dir=Cancer2file)

# create a Seurat object based on the scRNA-seq data
Cancer2 <- CreateSeuratObject(counts=Cancer2.counts, project="Cancer2_GSE184880")
#add metadata
Cancer2@meta.data[["response"]] <- "na"
Cancer2@meta.data[["tissue"]] <- "primary"
Cancer2@meta.data[["stage"]] <- "2b"
Cancer2@meta.data[["grade"]] <- "na"
Cancer2@meta.data[["diagnosis"]] <- "HGSOC"
Cancer2@meta.data[["age"]] <- "51"
Cancer2@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(Cancer2, paste(outdir,"Cancer2_GSE184880",enddir,sep=""))


###SRX12379610_Cancer3
Cancer3file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379610_Cancer3/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Cancer3.counts <- Read10X(data.dir=Cancer3file)

# create a Seurat object based on the scRNA-seq data
Cancer3 <- CreateSeuratObject(counts=Cancer3.counts, project="Cancer3_GSE184880")
#add metadata
Cancer3@meta.data[["response"]] <- "na"
Cancer3@meta.data[["tissue"]] <- "primary"
Cancer3@meta.data[["stage"]] <- "ic2"
Cancer3@meta.data[["grade"]] <- "na"
Cancer3@meta.data[["diagnosis"]] <- "HGSOC"
Cancer3@meta.data[["age"]] <- "41"
Cancer3@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(Cancer3, paste(outdir,"Cancer3_GSE184880",enddir,sep=""))


###SRX12379611_Cancer4
Cancer4file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379611_Cancer4/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Cancer4.counts <- Read10X(data.dir=Cancer4file)

# create a Seurat object based on the scRNA-seq data
Cancer4 <- CreateSeuratObject(counts=Cancer4.counts, project="Cancer4_GSE184880")
#add metadata
Cancer4@meta.data[["response"]] <- "na"
Cancer4@meta.data[["tissue"]] <- "primary"
Cancer4@meta.data[["stage"]] <- "ic2"
Cancer4@meta.data[["grade"]] <- "na"
Cancer4@meta.data[["diagnosis"]] <- "HGSOC"
Cancer4@meta.data[["age"]] <- "47"
Cancer4@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(Cancer4, paste(outdir,"Cancer4_GSE184880",enddir,sep=""))


###SRX12379612_Cancer5
Cancer5file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379612_Cancer5/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Cancer5.counts <- Read10X(data.dir=Cancer5file)

# create a Seurat object based on the scRNA-seq data
Cancer5 <- CreateSeuratObject(counts=Cancer5.counts, project="Cancer5_GSE184880")
#add metadata
Cancer5@meta.data[["response"]] <- "na"
Cancer5@meta.data[["tissue"]] <- "primary"
Cancer5@meta.data[["stage"]] <- "2b"
Cancer5@meta.data[["grade"]] <- "na"
Cancer5@meta.data[["diagnosis"]] <- "HGSOC"
Cancer5@meta.data[["age"]] <- "57"
Cancer5@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(Cancer5, paste(outdir,"Cancer5_GSE184880",enddir,sep=""))


###SRX12379613_Cancer6
Cancer6file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379613_Cancer6/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Cancer6.counts <- Read10X(data.dir=Cancer6file)

# create a Seurat object based on the scRNA-seq data
Cancer6 <- CreateSeuratObject(counts=Cancer6.counts, project="Cancer6_GSE184880")
#add metadata
Cancer6@meta.data[["response"]] <- "na"
Cancer6@meta.data[["tissue"]] <- "primary"
Cancer6@meta.data[["stage"]] <- "2c"
Cancer6@meta.data[["grade"]] <- "na"
Cancer6@meta.data[["diagnosis"]] <- "HGSOC"
Cancer6@meta.data[["age"]] <- "48"
Cancer6@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(Cancer6, paste(outdir,"Cancer6_GSE184880",enddir,sep=""))


###SRX12379614_Cancer7
Cancer7file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE184880/SRX12379614_Cancer7/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Cancer7.counts <- Read10X(data.dir=Cancer7file)

# create a Seurat object based on the scRNA-seq data
Cancer7 <- CreateSeuratObject(counts=Cancer7.counts, project="Cancer7_GSE184880")
#add metadata
Cancer7@meta.data[["response"]] <- "na"
Cancer7@meta.data[["tissue"]] <- "primary"
Cancer7@meta.data[["stage"]] <- "ic2"
Cancer7@meta.data[["grade"]] <- "na"
Cancer7@meta.data[["diagnosis"]] <- "HGSOC"
Cancer7@meta.data[["age"]] <- "53"
Cancer7@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(Cancer7, paste(outdir,"Cancer7_GSE184880",enddir,sep=""))


###GSE181955#####################################
###N1_normal
N1_normalfile <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE181955/N1_normal/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
N1_normal.counts <- Read10X(data.dir=N1_normalfile)

# create a Seurat object based on the scRNA-seq data
N1_normal <- CreateSeuratObject(counts=N1_normal.counts, project="N1_normal_GSE181955")
#add metadata
N1_normal@meta.data[["response"]] <- "na"
N1_normal@meta.data[["tissue"]] <- "normal"
N1_normal@meta.data[["stage"]] <- "na"
N1_normal@meta.data[["grade"]] <- "na"
N1_normal@meta.data[["diagnosis"]] <- "normal"
N1_normal@meta.data[["age"]] <- "na"
N1_normal@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(N1_normal, paste(outdir,"N1_normal_GSE181955",enddir,sep=""))


###OMT-1_CD45+
OMT1_CD45file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE181955/OMT-1_CD45+/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
OMT1_CD45.counts <- Read10X(data.dir=OMT1_CD45file)

# create a Seurat object based on the scRNA-seq data
OMT1_CD45 <- CreateSeuratObject(counts=OMT1_CD45.counts, project="OMT1_CD45+_GSE181955")
#add metadata
OMT1_CD45@meta.data[["response"]] <- "na"
OMT1_CD45@meta.data[["tissue"]] <- "omentum"
OMT1_CD45@meta.data[["stage"]] <- "na"
OMT1_CD45@meta.data[["grade"]] <- "na"
OMT1_CD45@meta.data[["diagnosis"]] <- "HGSOC"
OMT1_CD45@meta.data[["age"]] <- "na"
OMT1_CD45@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(OMT1_CD45, paste(outdir,"OMT1_CD45+_GSE181955",enddir,sep=""))


###OMT-3_CD45+
OMT3_CD45file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE181955/OMT-3_CD45+/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
OMT3_CD45.counts <- Read10X(data.dir=OMT3_CD45file)

# create a Seurat object based on the scRNA-seq data
OMT3_CD45 <- CreateSeuratObject(counts=OMT3_CD45.counts, project="OMT3_CD45+_GSE181955")
#add metadata
OMT3_CD45@meta.data[["response"]] <- "na"
OMT3_CD45@meta.data[["tissue"]] <- "omentum"
OMT3_CD45@meta.data[["stage"]] <- "na"
OMT3_CD45@meta.data[["grade"]] <- "na"
OMT3_CD45@meta.data[["diagnosis"]] <- "HGSOC"
OMT3_CD45@meta.data[["age"]] <- "na"
OMT3_CD45@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(OMT3_CD45, paste(outdir,"OMT3_CD45+_GSE181955",enddir,sep=""))


###T1
T1file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE181955/T1/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
T1.counts <- Read10X(data.dir=T1file)

# create a Seurat object based on the scRNA-seq data
T1 <- CreateSeuratObject(counts=T1.counts, project="T1_GSE181955")
#add metadata
T1@meta.data[["response"]] <- "na"
T1@meta.data[["tissue"]] <- "primary"
T1@meta.data[["stage"]] <- "na"
T1@meta.data[["grade"]] <- "na"
T1@meta.data[["diagnosis"]] <- "HGSOC"
T1@meta.data[["age"]] <- "na"
T1@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(T1, paste(outdir,"T1_GSE181955",enddir,sep=""))


###T6
T6file <- "/blue/ferrallm/00_data/single-cell/OvCa/GSE181955/T6/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
T6.counts <- Read10X(data.dir=T6file)

# create a Seurat object based on the scRNA-seq data
T6 <- CreateSeuratObject(counts=T6.counts, project="T6_GSE181955")
#add metadata
T6@meta.data[["response"]] <- "na"
T6@meta.data[["tissue"]] <- "primary"
T6@meta.data[["stage"]] <- "na"
T6@meta.data[["grade"]] <- "na"
T6@meta.data[["diagnosis"]] <- "HGSOC"
T6@meta.data[["age"]] <- "na"
T6@meta.data[["treatment"]] <- "untreated"
#save rds file
saveRDS(T6, paste(outdir,"T6_GSE181955",enddir,sep=""))


###CORE samples#####################################
###UF_1
UF_1file <- "/blue/ferrallm/00_data/single-cell/OvCa/GE7897-MFerral/UF_1/cellranger_multi_tutorial/outs/per_sample_outs/cellranger_multi_tutorial/count/sample_filtered_feature_bc_matrix"
UF_1.counts <- Read10X(data.dir=UF_1file)

# create a Seurat object based on the scRNA-seq data
UF_1 <- CreateSeuratObject(counts=UF_1.counts, project="UF_1_multi")
#add metadata
UF_1@meta.data[["response"]] <- "na"
UF_1@meta.data[["tissue"]] <- "primary"
UF_1@meta.data[["stage"]] <- "na"
UF_1@meta.data[["grade"]] <- "na"
UF_1@meta.data[["diagnosis"]] <- "HGSOC" #potentially, leave it like this until we know the diagnosis
UF_1@meta.data[["age"]] <- "na"
UF_1@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(UF_1, paste(outdir,"UF_1_multi",enddir,sep=""))


###UF_2
UF_2file <- "/blue/ferrallm/00_data/single-cell/OvCa/GE7897-MFerral/UF_2/cellranger_multi_tutorial/outs/per_sample_outs/cellranger_multi_tutorial/count/sample_filtered_feature_bc_matrix"
UF_2.counts <- Read10X(data.dir=UF_2file)

# create a Seurat object based on the scRNA-seq data
UF_2 <- CreateSeuratObject(counts=UF_2.counts, project="UF_2_multi")
#add metadata
UF_2@meta.data[["response"]] <- "na"
UF_2@meta.data[["tissue"]] <- "primary"
UF_2@meta.data[["stage"]] <- "na"
UF_2@meta.data[["grade"]] <- "na"
UF_2@meta.data[["diagnosis"]] <- "HGSOC" #potentially, leave it like this until we know the diagnosis
UF_2@meta.data[["age"]] <- "na"
UF_2@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(UF_2, paste(outdir,"UF_2_multi",enddir,sep=""))


###UF_3
UF_3file <- "/blue/ferrallm/00_data/single-cell/OvCa/GE7897-MFerral/UF_3/cellranger_multi_tutorial/outs/per_sample_outs/cellranger_multi_tutorial/count/sample_filtered_feature_bc_matrix"
UF_3.counts <- Read10X(data.dir=UF_3file)

# create a Seurat object based on the scRNA-seq data
UF_3 <- CreateSeuratObject(counts=UF_3.counts, project="UF_3_multi")
#add metadata
UF_3@meta.data[["response"]] <- "na"
UF_3@meta.data[["tissue"]] <- "primary"
UF_3@meta.data[["stage"]] <- "na"
UF_3@meta.data[["grade"]] <- "na"
UF_3@meta.data[["diagnosis"]] <- "HGSOC" #potentially, leave it like this until we know the diagnosis
UF_3@meta.data[["age"]] <- "na"
UF_3@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(UF_3, paste(outdir,"UF_3_multi",enddir,sep=""))


###UF_4
UF_4file <- "/blue/ferrallm/00_data/single-cell/OvCa/GE7897-MFerral/UF_4/cellranger_multi_tutorial/outs/per_sample_outs/cellranger_multi_tutorial/count/sample_filtered_feature_bc_matrix"
UF_4.counts <- Read10X(data.dir=UF_4file)

# create a Seurat object based on the scRNA-seq data
UF_4 <- CreateSeuratObject(counts=UF_4.counts, project="UF_4_multi")
#add metadata
UF_4@meta.data[["response"]] <- "na"
UF_4@meta.data[["tissue"]] <- "primary"
UF_4@meta.data[["stage"]] <- "na"
UF_4@meta.data[["grade"]] <- "na"
UF_4@meta.data[["diagnosis"]] <- "HGSOC" #potentially, leave it like this until we know the diagnosis
UF_4@meta.data[["age"]] <- "na"
UF_4@meta.data[["treatment"]] <- "na"
#save rds file
saveRDS(UF_4, paste(outdir,"UF_4_multi",enddir,sep=""))


##################################################################################
###CREATE RDS FILE FOR THE REFERENCE WITH THE CELL TYPES
##################################################################################
####Need to reformat the reference first so it is normalized using SCT
###Reference won't be scaled since it contains the right number of cells that we have the metadata for (13646 cells in total)
################GSE173682
####download  Patient_8_3BAE2L_RNA AND	Patient_9_3E5CFL_RNA, in total they should have 13647 cells 
pt8file <- "/blue/ferrallm/00_data/single-cell/OvCa/reference_for_singleR/GSE173682/Patient_8_3BAE2L"
pt8.counts <- Read10X(data.dir=pt8file)
pt8 <- CreateSeuratObject(counts=pt8.counts, project="pt8_GSE173682")

pt9file <- "/blue/ferrallm/00_data/single-cell/OvCa/reference_for_singleR/GSE173682/Patient_9_3E5CFL"
pt9.counts <- Read10X(data.dir=pt9file)
pt9 <- CreateSeuratObject(counts=pt9.counts, project="pt9_GSE173682")
#merge both patients 
reference_GSE173682 <- merge(x=pt8, y=list(pt9), merge.data=TRUE, project="reference_GSE173682")
reference_GSE173682 <- JoinLayers(reference_GSE173682, assay="RNA", new.layer="counts")
saveRDS(reference_GSE173682, paste(outdir,"OvCa-scRNAseq-Pilot-reference-GSE173682-JoinedLayers-unprocessed_2025-04-21.rds"))

