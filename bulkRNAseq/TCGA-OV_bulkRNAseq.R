#This script shows the preprocessing, quality control, normalization and downstream analysis for the TCGA-OV bulkRNAseq cohort used in the study 
#load libraries
library('TCGAbiolinks') #only for downloading data 
library('SummarizedExperiment') #only for downloading data 
library('dplyr')
library('readr')
library(ggplot2)
library(org.Hs.eg.db)
#install.packages("janitor")
library(janitor)
library(DESeq2)
library(tidyverse)
library(tximport)
library(BiocManager)
library(GSVA)

#set up output directories
outdir <- '/results/TCGA-OV_bulkRNAseq/'
dir.create(outdir)

#BiocParallel::register(BiocParallel::MulticoreParam(workers = 8)) 
#options(future.globals.maxSize = Inf) #if not SCTransform won't run, 50 * 1024^3

#############################
###Download and preprocessing
#############################
outdir1 <- '/results/TCGA-OV_bulkRNAseq/preprocessing/'
dir.create(outdir1)

OVCA_query <- GDCquery(project = "TCGA-OV",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

GDCdownload(OVCA_query)

OVCA_res = getResults(OVCA_query)
colnames(OVCA_res)
# [1] "id"                        "data_format"               "cases"                     "access"                   
# [5] "file_name"                 "submitter_id"              "data_category"             "type"                     
# [9] "platform"                  "file_size"                 "created_datetime"          "md5sum"                   
# [13] "updated_datetime"          "file_id"                   "data_type"                 "state"                    
# [17] "experimental_strategy"     "version"                   "data_release"              "project"                  
# [21] "analysis_id"               "analysis_state"            "analysis_submitter_id"     "analysis_workflow_link"   
# [25] "analysis_workflow_type"    "analysis_workflow_version" "sample_type"               "is_ffpe"                  
# [29] "cases.submitter_id"        "sample.submitter_id"      
table(OVCA_res$sample_type)

#to make a summarizedExperiment object 
OVCA_data <- GDCprepare(OVCA_query)

#get the metadata only first 
barcodes <- OVCA_data@colData$barcode
days_to_diagnosis <- OVCA_data@colData@listData[["days_to_diagnosis"]]
days_to_last_follow_up <- OVCA_data@colData@listData[["days_to_last_follow_up"]]
age_at_diagnosis <- OVCA_data@colData@listData[["age_at_diagnosis"]]
year_of_diagnosis <- OVCA_data@colData@listData[["year_of_diagnosis"]]
vital_status <- OVCA_data@colData@listData[["vital_status"]]
age_at_index <- OVCA_data@colData@listData[["age_at_index"]]
days_to_birth <- OVCA_data@colData@listData[["days_to_birth"]]
days_to_death <- OVCA_data@colData@listData[["days_to_death"]]

metadata <- cbind(barcodes,days_to_diagnosis,days_to_last_follow_up,age_at_diagnosis,year_of_diagnosis,vital_status,
                  age_at_index,days_to_birth,days_to_death)
metadata_df <- data.frame(metadata)
rownames(metadata_df) <- metadata_df$barcodes

write.csv(metadata_df, '/results/TCGA-OV_bulkRNAseq/preprocessing/metadata_all_TCGA-OVpatients.csv')

#format a table to get the therapeutic agents and when they were delivered
col_data <- colData(OVCA_data)
barcodes <- col_data$barcode
days_to_diagnosis <- col_data$days_to_diagnosis
treatment_list <- col_data$treatments
prior_treatment <- col_data$prior_treatment #originally 278 NO and 149 YES 
days_to_last_follow_up <- col_data$days_to_last_follow_up


# Now unpack treatments
all_treatments <- lapply(seq_along(treatment_list), function(i) {
  treatments <- treatment_list[[i]]
  
  # Skip if no treatment data
  if (is.null(treatments) || nrow(treatments) == 0) return(NULL)
  
  n <- nrow(treatments)
  
  # Extract and clean days_to_treatment_start
  days_start <- treatments$days_to_treatment_start
  if (is.null(days_start) || length(days_start) != n) {
    days_start <- rep(NA, n)
  }
  
  # --- Handle treatment end ---
  days_end <- treatments$days_to_treatment_end
  if (is.null(days_end) || length(days_end) != n) {
    days_end <- rep(NA, n)
  }

  # add number of cycles
  number_of_cycles <- if (!is.null(treatments$number_of_cycles) && length(treatments$number_of_cycles) == n)
    treatments$number_of_cycles else rep(NA, n)
  
  # Safely handle therapeutic_agents
  agents <- treatments$therapeutic_agents
  if (is.null(agents) || length(agents) != n) {
    agents_collapsed <- rep(NA_character_, n)
  } else if (is.list(agents)) {
    agents_collapsed <- sapply(agents, function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      paste(x, collapse = ", ")
    })
  } else {
    agents_collapsed <- as.character(agents)
  }
  
  # Create the data frame with matching row counts
  data.frame(
    barcode = rep(barcodes[i], n),
    prior_treatment = rep(prior_treatment[i], n),
    days_to_diagnosis = rep(days_to_diagnosis[i], n),
    days_to_treatment_start = days_start,
    days_to_treatment_end = days_end,
    days_to_last_follow_up = rep(days_to_last_follow_up[i], n),
    #treatment_or_therapy = treatment_therapy,
    therapeutic_agents = agents_collapsed,
    number_of_cycles = number_of_cycles,
    #initial_disease_status = initial_disease_status,
    stringsAsFactors = FALSE
  )
})


treatment_df <- do.call(rbind, all_treatments)#table with treatment information 

#save table with the previous information: barcodes, prior treatment, days to diagnosis, days to treatment start, 
#days to treatment end, days to last follow up, therapeutic agents, number of cycles 
write.csv(treatment_df, '/results/TCGA-OV_bulkRNAseq/preprocessing/all_treatmentsANDcyclesANDdaystolastfollowup_per_patient.csv')

#from that table calculate treatment free intervals on excel and delete patients with missing information

################
###Normalization 
################
outdir2 <- '/results/TCGA-OV_bulkRNAseq/normalization/'
dir.create(outdir2)
###Upload table with 102 patients and TFI after formatting it in excel
patient102 <- read.csv(dir,"allTreatments_102patients_noNA-cis_withTFI.csv")
table(patient102$platin_second)
# no  no (nonresponder) yes (responder) 
# 1  54                 47 

###add 180 groups 
patient102 <- patient102 %>%
  mutate(tfi_180group = ifelse(tfi > 180, "above_180", "below_180"))
table(patient102$tfi_180group)
# above_180 (sensitive) below_180 (resistant)
# 43 (0.42)             59 (0.578)

###Extract RNA data from the OVCA object 
OVCA_count_matrix <- assay(OVCA_data)
OVCA_count_matrix_df <- as.data.frame(OVCA_count_matrix)

OVCA_genes <- rownames(OVCA_count_matrix) %>%
  tibble::enframe() %>%
  mutate(ENSEMBL = stringr::str_replace(value, "\\.[0-9]+",""))

row.names(OVCA_count_matrix) <- OVCA_genes$ENSEMBL
OVCA_count_matrix_df <- data.frame(OVCA_count_matrix)

#change from ENSEMBL to symbol 
clusterProfiler::bitr(OVCA_genes$ENSEMBL,
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db) %>%
  janitor::get_dupes(SYMBOL) %>%
  head()

#keep only one gene out of the duplicates
OVCA_genes_map <- clusterProfiler::bitr(OVCA_genes$ENSEMBL,
                                        fromType = "ENSEMBL",
                                        toType = "SYMBOL",
                                        OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE)

keep_ensembl <- intersect(OVCA_genes_map$ENSEMBL, rownames(OVCA_count_matrix))
OVCA_count_matrix <- OVCA_count_matrix[keep_ensembl, ]

symbol_map <- OVCA_genes_map$SYMBOL[match(rownames(OVCA_count_matrix), OVCA_genes_map$ENSEMBL)]
rownames(OVCA_count_matrix) <- symbol_map

OVCA_count_matrix_df <- data.frame(OVCA_count_matrix)
write.csv(OVCA_count_matrix_df,'/results/TCGA-OV_bulkRNAseq/normalization/TCGA-OVCA_count_matrix.csv',row.names=TRUE)

#change the colnames to have - between the numbers and not .
colnames(OVCA_count_matrix_df) <- gsub("\\.", "-", colnames(OVCA_count_matrix_df))

#collect barcodes for the 102 patients 
barcodes_102patients <- patient102$barcode

count_matrix_reduced <- OVCA_count_matrix_df[,colnames(OVCA_count_matrix_df) %in% barcodes_102patients] #102 columns 
dim(count_matrix_reduced) #36940   102
write.csv(count_matrix_reduced, '/results/TCGA-OV_bulkRNAseq/normalization/raw-count-matrix_102patients_TFI-from-firstPlatinumTOnextTreat.csv', row.names= TRUE)

#normalize the count matrix 
#calculate standard deviation across all samples, no need to exclude the first column because the gene symbols are the rownames
count_matrix_sd <- count_matrix_reduced[order(-apply(count_matrix_reduced, 1, sd)), ]
dim(count_matrix_sd) #36940   102

#remove duplicated genes, in this case all of this was done prior 
count_matrix_sd <- count_matrix_sd[!duplicated(rownames(count_matrix_sd)) ,] #36940 leftover genes, no duplicated genes
#remove genes that were not mapped 
count_matrix_sd <- count_matrix_sd[!is.na(rownames(count_matrix_sd)),] #36940 leftover genes, all genes were mapped

#set gene names as row names, these are already the rownames 
#rownames(count_matrix_sd) <- count_matrix_sd[,1]

#convert data (excluding the first column) to a numeric matrix 
count_matrix_sd <- as.matrix(count_matrix_sd)

#filtering the count matrix, based on a threshold for CPM
minCPM <- 0.5 
nLibraries <- 1
count_matrix_sd_matrix_filtered <- count_matrix_sd[ which( apply( cpm(DGEList(counts = count_matrix_sd)),  1,
                                                                  function(y) sum(y >= minCPM) ) >= nLibraries ), ] 
dim(count_matrix_sd_matrix_filtered) #22597 102 after filtering 

###build DESeq2, to compare data in the treatment group, if wanting to compare the sensitive group then change "design"
##DESeq2 needs colData to be able to run but it will ignore it if design = ~1 so create a dummy that it can run with 

#make a metadata factor based on tfi_group
patient102 <- data.frame(lapply(patient102, trimws, which = "left"), stringsAsFactors = FALSE)
rownames(patient102) <- patient102$barcode
all(colnames(count_matrix_sd_matrix_filtered) %in% rownames(patient102)) #true 
all(rownames(patient102) == colnames(count_matrix_sd_matrix_filtered)) #false, not in the same order 
#re order columns in count matrix 
patient102 <- patient102[colnames(count_matrix_sd_matrix_filtered), ]
all(rownames(patient102) == colnames(count_matrix_sd_matrix_filtered)) #true, in the same order 

#run dds based on response
dds <- DESeqDataSetFromMatrix(countData = count_matrix_sd_matrix_filtered, colData = patient102, design = ~ tfi_180group) 
dds <- DESeq(dds)
res <- results(dds, contrast = c("tfi_180group", "above_180","below_180")) #upregulated in above_180, downregulated in below_180
#create table with gene names, log2FC, pval and padj 
res$GeneSymbol <- rownames(res)
#reorder columns: 
res_df <- res[,c("GeneSymbol", "log2FoldChange", "pvalue", "padj")]
res_df_table <- as.data.frame(res_df@listData)
write.csv(res_df_table,"/results/TCGA-OV_bulkRNAseq/normalization/differetially_expressed_genes-HigherAbove180VSLowerBelow180-tfi_102patients.csv" )

#lastly do log transformation 
c <- 4 #pseudocount is 4
count_matrix_log_norm <- log2(counts(dds, normalized=TRUE) + c)
write.csv(count_matrix_log_norm,"/results/TCGA-OV_bulkRNAseq/normalization/processed_data_102patients_logNorm.csv")


#### Attempt to add metadata to the log normalized count matrix table ####
#transpose matrix
count_matrix_log_norm_t <- t(count_matrix_log_norm)

#add a barcode column on the transpose count matrix 
count_matrix_log_norm_t <- data.frame(barcode = rownames(count_matrix_log_norm_t), count_matrix_log_norm_t)

#delete non-relevant information from 90 patient dataframe
patient102$days_to_treatment_start <- NULL
patient102$days_to_treatment_end <- NULL
patient102$therapeutic_agents <- NULL
patient102$prior_treatment <- NULL
patient102$platin_first <- NULL
patient102$platin_type <- NULL

#merge count matrix with metadata table
count_matrix_with_102patients <- merge(patient102, count_matrix_log_norm_t, by = "barcode", all.x = TRUE)
write_csv(count_matrix_with_102patients, '/results/TCGA-OV_bulkRNAseq/normalization/lognorm_102patients_count_matrix_WITH_tfi_firstPlatinumTOnexTreat.csv')

###################################
###Scoring based on COMET signature 
###################################
outdir3 <- '/results/TCGA-OV_bulkRNAseq/COMET_scoring/'
dir.create(outdir3)

#normalized TCGA OV table to matrix format
count_matrix_log_norm_mat <- as.matrix(count_matrix_log_norm)

#read COMET genes 
COMET_genes <- read.csv('/data/bulkRNAseq/cluster_1_singleton_all_ranked.csv')
COMET_genes$gene_1
gene_set <- list(COMET = COMET_genes$gene_1) #makes a list of 1 where all the genes are part of 1 list, get only 1 score per sample
#gene_set <- as.list(COMET_genes$gene_1) the list has a slot per gene so get a score per gene per sample
#run gsva
res_sig <- gsva(gsvaParam(TCGAOV_lognorm_mat, gene_set))

res_sig_df <- as.data.frame(res_sig)

res_sig_df_t <- t(res_sig_df)

#save the matrix
write.csv(res_sig_df_t,paste0(outdir3, 'TCGA-OV-cis_res-signaturescoresfromCOMET.csv'))



