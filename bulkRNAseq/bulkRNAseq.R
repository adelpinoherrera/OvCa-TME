#This script shows the preprocessing, quality control, normalization and downstream analysis for the bulkRNAseq cohort used in the study 

library(edgeR)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)

#set up output directories
outdir <- '/results/bulkRNAseq/'
dir.create(outdir)
dir <- '/data/bulkRNAseq/'

################
###Preprocessing
################
outdir1 <- '/results/bulkRNAseq/preprocessing/'
dir.create(outdir1)

counts_6142 <- read.table(file.path(dir, "6142_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
#skip = 4 deletes the first 4 rows so it gets rid off the N mapped data 
rownames(counts_6142) <- counts_6142[,1] #make the gene names the row names 
counts_6142 <- counts_6142[,2, drop = FALSE]#get the second colum (which is the first count column which contains 
#the unstranded reads which combines the forward-stranded and the reverse-stranded reads)

counts_6159 <- read.table(file.path(dir, "6159_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_6159) <- counts_6159[,1]
counts_6159 <- counts_6159[,2, drop = FALSE] 

counts_6161 <- read.table(file.path(dir, "6161_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_6161) <- counts_6161[,1]
counts_6161 <- counts_6161[,2, drop = FALSE] 

counts_6242 <- read.table(file.path(dir, "6242_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_6242) <- counts_6242[,1]
counts_6242 <- counts_6242[,2, drop = FALSE] 

counts_6273 <- read.table(file.path(dir, "6273_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_6273) <- counts_6273[,1]
counts_6273 <- counts_6273[,2, drop = FALSE] 

counts_6589 <- read.table(file.path(dir, "6589_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_6589) <- counts_6589[,1]
counts_6589 <- counts_6589[,2, drop = FALSE] 

counts_6823 <- read.table(file.path(dir, "6823_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_6823) <- counts_6823[,1]
counts_6823 <- counts_6823[,2, drop = FALSE] 

counts_7472 <- read.table(file.path(dir, "7472_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_7472) <- counts_7472[,1]
counts_7472 <- counts_7472[,2, drop = FALSE] 

counts_12799 <- read.table(file.path(dir, "12799_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_12799) <- counts_12799[,1]
counts_12799 <- counts_12799[,2, drop = FALSE] 

counts_12838 <- read.table(file.path(dir, "12838_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_12838) <- counts_12838[,1]
counts_12838 <- counts_12838[,2, drop = FALSE] 

counts_12970_1AN2 <- read.table(file.path(dir, "12970-1AN2_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_12970_1AN2) <- counts_12970_1AN2[,1]
counts_12970_1AN2 <- counts_12970_1AN2[,2, drop = FALSE] 

counts_14257 <- read.table(file.path(dir, "14257_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_14257) <- counts_14257[,1]
counts_14257 <- counts_14257[,2, drop = FALSE] 

counts_14379_1AN2 <- read.table(file.path(dir, "14379-1AN2_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_14379_1AN2) <- counts_14379_1AN2[,1]
counts_14379_1AN2 <- counts_14379_1AN2[,2, drop = FALSE] 

counts_14379_1BC1 <- read.table(file.path(dir, "14379-1BC1_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_14379_1BC1) <- counts_14379_1BC1[,1]
counts_14379_1BC1 <- counts_14379_1BC1[,2, drop = FALSE] 

counts_38716 <- read.table(file.path(dir, "38716_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_38716) <- counts_38716[,1]
counts_38716 <- counts_38716[,2, drop = FALSE] 

counts_39170 <- read.table(file.path(dir, "39170_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_39170) <- counts_39170[,1]
counts_39170 <- counts_39170[,2, drop = FALSE] 

counts_39998 <- read.table(file.path(dir, "39998_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_39998) <- counts_39998[,1]
counts_39998 <- counts_39998[,2, drop = FALSE] 

counts_40960 <- read.table(file.path(dir, "40960_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_40960) <- counts_40960[,1]
counts_40960 <- counts_40960[,2, drop = FALSE] 

counts_41427 <- read.table(file.path(dir, "41427_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_41427) <- counts_41427[,1]
counts_41427 <- counts_41427[,2, drop = FALSE] 

counts_46828_1AC1 <- read.table(file.path(dir, "46828-1AC1_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_46828_1AC1) <- counts_46828_1AC1[,1]
counts_46828_1AC1 <- counts_46828_1AC1[,2, drop = FALSE] 

counts_46828_1AN2 <- read.table(file.path(dir, "46828-1AN2_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_46828_1AN2) <- counts_46828_1AN2[,1]
counts_46828_1AN2 <- counts_46828_1AN2[,2, drop = FALSE] 

counts_48420 <- read.table(file.path(dir, "48420_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_48420) <- counts_48420[,1]
counts_48420 <- counts_48420[,2, drop = FALSE] 

counts_56302 <- read.table(file.path(dir, "56302_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_56302) <- counts_56302[,1]
counts_56302 <- counts_56302[,2, drop = FALSE] 

counts_62941 <- read.table(file.path(dir, "62941_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_62941) <- counts_62941[,1]
counts_62941 <- counts_62941[,2, drop = FALSE] 

counts_66062 <- read.table(file.path(dir, "66062_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_66062) <- counts_66062[,1]
counts_66062 <- counts_66062[,2, drop = FALSE] 

counts_71201 <- read.table(file.path(dir, "71201_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_71201) <- counts_71201[,1]
counts_71201 <- counts_71201[,2, drop = FALSE] 

counts_99855 <- read.table(file.path(dir, "99855_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_99855) <- counts_99855[,1]
counts_99855 <- counts_99855[,2, drop = FALSE] 

counts_UF_5 <- read.table(file.path(dir, "UF_5_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_UF_5) <- counts_UF_5[,1]
counts_UF_5 <- counts_UF_5[,2, drop = FALSE] 

counts_UF_4 <- read.table(file.path(dir, "UF_4_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_UF_4) <- counts_UF_4[,1]
counts_UF_4 <- counts_UF_4[,2, drop = FALSE]

counts_UF_3 <- read.table(file.path(dir, "UF_3_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_UF_3) <- counts_UF_3[,1]
counts_UF_3 <- counts_UF_3[,2, drop = FALSE] 

counts_UF_2 <- read.table(file.path(dir, "UF_2_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_UF_2) <- counts_UF_2[,1]
counts_UF_2 <- counts_UF_2[,2, drop = FALSE] 

counts_UF_1 <- read.table(file.path(dir, "UF_1_ReadsPerGene.out.tab"), sep = '\t', header = F, skip = 4, quote = '')
rownames(counts_UF_1) <- counts_UF_1[,1]
counts_UF_1 <- counts_UF_1[,2, drop = FALSE] 

#join all the objects
all_counts <- cbind(counts_6142, counts_6159, counts_6161, counts_6242, counts_6273, counts_6589, counts_6823, 
                    counts_7472, counts_12799, counts_12838, counts_12970_1AN2, counts_14257, counts_14379_1AN2,
                    counts_14379_1BC1, counts_38716, counts_39170, counts_39998, counts_40960, counts_41427,
                    counts_46828_1AC1, counts_46828_1AN2, counts_48420, counts_56302, counts_62941, counts_66062,
                    counts_71201, counts_99855, counts_UF_5, counts_UF_4, counts_UF_3, counts_UF_2, counts_UF_1)

#change colnames in all_counts
colnames(all_counts) <- c("6142", "6159","6161","6242","6273","6589","6823",
                          "7472","12799","12838","12970.1AN2","14257","14379.1AN2",
                          "14379.1BC1","38716","39170","39998","40960","41427",
                          "46828.1AC1","46828.1AN2","48420","56302","62941","66062",
                          "71201","99855","UF_5","UF_4","UF_3","UF_2","UF_1")

#save the data.frame as a csv file 
write.csv(all_counts, paste(outdir1, "/results/bulkRNAseq/preprocessing/unprocessed_counts_readspergene.csv"))

####convert ensembl symbols to gene symbols
#add gene column to the data
all_counts$Gene <- rownames(all_counts)
all_counts$Gene <- mapIds(org.Hs.eg.db, keys = row.names(all_counts), keytype = "ENSEMBL", 
                          column = "SYMBOL", multiVals = "first")

table(is.na(all_counts$Gene))
#FALSE  TRUE 
#24555 12046 

#delete the rows which gene wasnt mapped and add gene symbols to the first column 
all_counts <- all_counts[!is.na(all_counts$Gene), ]
dim(all_counts)
#24555    33

all_counts <- all_counts[, c("Gene", setdiff(names(all_counts), "Gene"))]
write.csv(all_counts, "/results/bulkRNAseq/preprocessing/unprocessed_countsWITHgenesymbol_readspergene.csv")

#delete the rownames and add new ones 
rownames(all_counts) <- NULL

#remove duplicates 
table(duplicated(all_counts$Gene))
# FALSE  TRUE 
# 24492    63 

# Add mean expression column
all_counts$mean_expr <- rowMeans(all_counts[, -which(names(all_counts) == "Gene")], na.rm = TRUE)

# Keep row with highest expression for each gene
all_counts_unique <- all_counts %>%
  group_by(Gene) %>%
  slice_max(order_by = mean_expr, n = 1) %>%
  ungroup()
dim(all_counts_unique)
# 24494    34

#delete expression column and add gene names as rownames 
all_counts_unique$mean_expr <- NULL
all_counts_unique <- as.data.frame(all_counts_unique)
#these 2 genes are also duplicated with 0 expression each so delete both of them ‘DEFB130A’, ‘TSPY3’ 
all_counts_unique <- all_counts_unique[!all_counts_unique$Gene %in% c("DEFB130A","TSPY3"),]
dim(all_counts_unique)
#24490    33

#set gene names as rownames
rownames(all_counts_unique) <- all_counts_unique[,1]
all_counts_unique$Gene <- NULL
write.csv(all_counts_unique, "/results/bulkRNAseq/preprocessing/unprocessed_countsWITHgenesymbolNoDuplicates_readspergene.csv")

#################################################
###Normalization and differential gene expression
#################################################
outdir2 <- '/results/bulkRNAseq/normalization/'
dir.create(outdir2)

####prepare data for normalization
#calculate standard deviation across all samples, no need to exclude the first column because the gene symbols are the rownames
all_counts_unique_sd <- all_counts_unique[order(-apply(all_counts_unique, 1, sd)), ]

#convert data (excluding the first column) to a numeric matrix 
all_counts_matrix <- as.matrix(all_counts_unique_sd)

#filtering the count matrix, based on a threshold for CPM
minCPM <- 0.5 
nLibraries <- 1
all_counts_matrix_filtered <- all_counts_matrix[ which( apply( cpm(DGEList(counts = all_counts_matrix)),  1,
                                                               function(y) sum(y >= minCPM) ) >= nLibraries ), ] 
dim(all_counts_matrix_filtered) #21487    32 after filtering 

###build DESeq2, to compare data in the treatment group, if wanting to compare the sensitive group then change "design"
##DESeq2 needs colData to be able to run but it will ignore it if design = ~1 so create a dummy that it can run with 
metadata_all_counts <- read.csv("/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/OVCA_allsamples/metadata_unprocessed_readspergene_ptinrows.csv")

#make a metadata factor based on tfi_group
rownames(metadata_all_counts) <- metadata_all_counts$X
#add another column with non-cancerous vs cancerous 
metadata_all_counts$Diagnosis2 <- metadata_all_counts$Diagnosis
metadata_all_counts$Diagnosis2[metadata_all_counts$Diagnosis2 == 'Normal'] <- "noncancerous"
metadata_all_counts$Diagnosis2[metadata_all_counts$Diagnosis2 == 'Benign'] <- "noncancerous"
metadata_all_counts$Diagnosis2[metadata_all_counts$Diagnosis2 == 'Benign '] <- "noncancerous"
metadata_all_counts$Diagnosis2[metadata_all_counts$Diagnosis2 == 'Primary'] <- "cancerous"
all(colnames(all_counts_matrix_filtered_corrected) %in% rownames(metadata_all_counts)) #true 
all(rownames(metadata_all_counts) == colnames(all_counts_matrix_filtered_corrected)) #true, in the same order 

#run the dds on the matrix first
dds <- DESeqDataSetFromMatrix(countData = all_counts_matrix_filtered, colData = metadata_all_counts, design = ~ Diagnosis2) #use design 1 if not comparing any groups and eliminate colData
dds <- DESeq(dds)
#write the test first and then the control: primary, normal then positive FC would be high in cancer and low in normal
#cancerous vs non-cancerous
res.cvsn <- results(dds, contrast = c("Diagnosis2", "cancerous","noncancerous"))
#create table with gene names, log2FC, pval and padj 
res.cvsn$GeneSymbol <- rownames(res.cvsn)
res_df_cvsn <- res.cvsn[,c("GeneSymbol", "log2FoldChange", "pvalue", "padj")]
res_df_cvsn_table <- as.data.frame(res_df_cvsn@listData)
write.csv(res_df_cvsn_table,"/results/bulkRNAseq/normalization/differetially_expressed_genes-HigherCancerousVSLowerNon-Cancerous_OVCAallsamples32.csv" )

#normalize the data
c <- 4 #pseudocount is 4
all_counts_matrix_log_norm <- log2(counts(dds, normalized=TRUE) + c)
write.csv(all_counts_matrix_log_norm, "/results/bulkRNAseq/normalization/processed_countsWITHgenesymbolNoDuplicates_logNorm.csv")

################################
###Prepare matrix for CIBERSORTx
################################
outdir3 <- '/results/bulkRNAseq/CIBERSORT/'
dir.create(outdir3)
###Use raw count matrix 
#convert to counts per million 
CPM_all_counts_matrix <- apply(all_counts_unique, 2, function(x) x / sum(x)*1e6)
dim(CPM_all_counts_matrix) #24490    32

#add gene names as a column
Gene <- rownames(CPM_all_counts_matrix)
CPM_all_counts_matrix <- cbind(Gene, CPM_all_counts_matrix)

write.table(CPM_all_counts_matrix, "/results/bulkRNAseq/CIBERSORT/CPM-count-matrix_OVCAallsamples32_forCIBERSORTx.txt", sep = "\t", na="", row.names = FALSE, quote=FALSE)


###Prepare reference for CIBERSORTx
###################################
#read reference after single-cell normalization
subset_reference_GSE17.norm <- readRDS("/results/scRNAseq/singleR_for_cellType_classification/OvCa-scRNAseq-Pilot-subset-reference-GSE173683-JoinedLayers-PostPCA-withCellTypes.rds")

count_matrix_cell_type <- GetAssayData(subset_reference_GSE17.norm, assay = "RNA",layer = "counts") #picked RNA and counts because those are the raw data, the data slot is log-normalized

#calculate the CPM for CIBERSORTx
counts_per_million_count_matrix <- apply(count_matrix_cell_type, 2, function(x) {
  x/sum(x)*1e6
})

#change column names to cell types 
colnames(counts_per_million_count_matrix) <- subset_reference_GSE17.norm$cell.type

#change column names so repeated columns don't contain .1, .2, ect
colnames(counts_per_million_count_matrix) <- sub("\\.\\d+$", "", colnames(counts_per_million_count_matrix))

#add gene names to a new columns and export the table 
counts_per_million_count_matrix_df <- as.data.frame(counts_per_million_count_matrix)
counts_per_million_count_matrix_df <- cbind(GeneSymbol = rownames(counts_per_million_count_matrix_df), counts_per_million_count_matrix_df)
#checking I have the right matrix 
colnames(counts_per_million_count_matrix_df)
rownames(counts_per_million_count_matrix_df)
#save as a txt file 
write.table(counts_per_million_count_matrix_df, "/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim1-scRNAseq/results_rds/referenceGSE173682_individualCells_postPCA_forCIBERSORTx.txt", sep = "\t", na="", row.names = FALSE, col.names = TRUE, quote=FALSE)

#originally included 33538 genes 
dim(counts_per_million_count_matrix_df) #33538 13647

#calculate the standard deviation for each gene 
reference_numeric <- counts_per_million_count_matrix_df[ , sapply(counts_per_million_count_matrix_df, is.numeric) ]

#calculate the standard deviation 
gene_sd <- apply(reference_numeric, 1, sd)

#keep the top 3500 most variable genes 
top_n <- 3500
top_gene_indices <- order(gene_sd, decreasing = TRUE)[1:top_n]
reference_filtered <- reference[top_gene_indices, ] 
dim(reference_filtered) #3500 13647

#check size
object.size(reference_filtered) / 1024^2 #366.1 bytes 

#delete the number after the column names for the cell types 
colnames(reference_filtered) <- sub("\\.\\d+$", "", colnames(reference_filtered))

#create a vector of names for the celltypes
celltypes <- colnames(reference_filtered)[-1]
table(celltypes)
# B.cell Empty.Epithelial.cell      Endothelial.cell       Epithelial.cell            Fibroblast 
# 192                   308                   334                  5743                  3597 
# Macrophage               NK.cell                T.cell 
# 1947                   418                  1107 
n_cells_per_type <- 192 #to keep all the B cells 

#create a dataframe of col indices and cell types 
col_df <- data.frame(index = 2:ncol(reference_filtered), cellType = celltypes)

#set a seed for reproducibility 
set.seed(123)

selected_indices <- col_df %>%
  group_by(cellType) %>%
  group_modify(~ {
    n <- min(n_cells_per_type, nrow(.x))
    .x[sample(nrow(.x), n), ]
  }) %>%
  pull(index)

#subset to keep gene symbol column 
reference_downsampled <- reference_filtered[, c(1, selected_indices)]
dim(reference_downsampled) #3500  1537

#check size
object.size(reference_downsampled) / 1024^2 #41.4 bytes, might be small enough 

#delete the numbers and dots on the colnames
colnames(reference_downsampled) <- sub("\\.\\d+$", "", colnames(reference_downsampled))


#export table 
write.table(reference_downsampled, "/results/bulkRNAseq/CIBERSORT/referenceGSE173682_individualCells_3500mostVariableGenes-192cellperCellType_postPCA_forCIBERSORTx.txt", sep = "\t", na="", row.names = FALSE, col.names = TRUE, quote=FALSE)
