#This script shows the pathway analysis done between cancerous and normal samples in both bulk RNA- and scRNA-seq used in the study 
library('ggplot2')
library('Seurat') 
library('patchwork')
library('tidyverse') 
library('data.table')
#BiocManager::install('glmGamPoi')
library('glmGamPoi') 
library('hdf5r')
# only 1 time BiocManager::install("SingleR")
library('SingleR')
library(BiocParallel)
#only 1 time BiocManager::install("ArrayExpress")
library(ArrayExpress)
library('plyr')
library('celldex')
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(msigdbr)
library(ggrepel)

#set up output directories
outdir <- '/results/bulkandscRNAseq/'
dir.create(outdir)

###############################################
#Differential gene expression, cancer vs normal
###############################################

###scRNAseq genes
ovcaH <- readRDS("/results/singleR_for_sensitivity_classification/OvCa-scRNAseq-ovca-JoinedLayers_ALL-PostUMAPmodel_withHarmony_withPredictedResponseANDcellTypes.rds")
Idents(ovcaH) <- "diagnosis"

#don't use FindAllMarkers if comparing 2 specific groups
cancervsnormal.markers <- FindMarkers(ovcaH, ident.1 = "HGSOC", ident.2 = "normal") #17713 genes
cancervsnormal.markers$gene <- rownames(cancervsnormal.markers)
write.csv(cancervsnormal.markers, paste0(outdir,'differentially_expressed_genes-HigherHGSOCVSLowerNormal_scRNAseq.csv'))

sig_sc_cvsn <- cancervsnormal.markers %>% filter(p_val_adj <= 0.05) #16896
write.csv(sig_sc_cvsn, paste0(outdir,'differentially_expressed_genes-SignificantPval<0.05-HigherHGSOCVSLowerNormal_scRNAseq.csv'))

###bulk genes 
res_df_cvsn_table <- read_csv('/results/bulkRNAseq/normalization/differetially_expressed_genes-HigherCancerousVSLowerNon-Cancerous_OVCAallsamples32.csv')
sig_bulk_cvsn <- res_df_cvsn_table %>% filter(padj <= 0.05) #10403 before 21487
write.csv(sig_bulk_cvsn, paste0(outdir,'differentially_expressed_genes-SignificantPval<0.05-HigherHGSOCVSLowerNormal_bulkRNAseq.csv'))


#create a new table with FC and adj pval for all the matched genes 
intersect_cvsn_all <- intersect(sig_sc_cvsn$gene, sig_bulk_cvsn$GeneSymbol) #7220

intersect_cvsn_sc_match <- match(intersect_cvsn_all, sig_sc_cvsn$gene)

sig_cvsn_match <- data.frame(gene = intersect_cvsn_all, sc_avg_log2FC = sig_sc_cvsn$avg_log2FC[intersect_cvsn_sc_match], 
                             sc_p_val_adj = sig_sc_cvsn$p_val_adj[intersect_cvsn_sc_match], stringsAsFactors = FALSE)

instersect_cvsn_bulk_match <- match(intersect_cvsn_all, sig_bulk_cvsn$GeneSymbol)

sig_cvsn_bulk_match <- data.frame(gene = intersect_cvsn_all, bulk_avg_log2FC = sig_bulk_cvsn$log2FoldChange[instersect_cvsn_bulk_match], 
                                  bulk_p_val_adj = sig_bulk_cvsn$padj[instersect_cvsn_bulk_match], stringsAsFactors = FALSE)

intersect_cvsn_all_data <- cbind(sig_cvsn_match, sig_cvsn_bulk_match)
write.csv(intersect_cvsn_all_data, paste0(outdir,'significant-DEG-HigherCancerousVSLowerNon-Cancerous_matchingBulkRNAseqANDscRNAseq_ALL-withDATA.csv'))

####pathway analysis sc and bulk cancer vs normal UPREGULATED log2FC > 1 ####
#we want significant genes with values (padj <= 0.05, log2FC of > or < 0) and background genes (every gene in the differential expression analysis) 
sig_sc_cvsn_path_up <- sig_sc_cvsn %>% filter(avg_log2FC > 1) #1106 genes upregulated and significant
sig_bulk_cvsn_path_up <- sig_bulk_cvsn %>% filter(log2FoldChange > 1) #4845

###sc
#need to conver genesymbol to ENTREZID
sig_sc_cvsn_path_up_map <- bitr(geneID = sig_sc_cvsn_path_up$gene,
                                fromType = "SYMBOL",
                                toType = "ENTREZID",
                                OrgDb = "org.Hs.eg.db") #now 937
# 15.28% of input gene IDs are fail to map..., not ideal but will work for now
head(sig_sc_cvsn_path_up_map) #these would be the correct differentially expressed genes
# SYMBOL ENTREZID
# 1    IGKC     3514
# 2    CD74      972
# 3   IGLC2     3538
# 4 HLA-DRA     3122
# 5    RGS1     5996
# 6   WFDC2    10406

#for the background, we would need all the genes from the scRNAseq object which are 33612
background_sc_cvsn_path_map <- bitr(geneID = rownames(ovcaH),
                                    fromType = "SYMBOL",
                                    toType = "ENTREZID",
                                    OrgDb = "org.Hs.eg.db")
#34.7% of input gene IDs are fail to map..., now 21950 genes  

###bulk
sig_bulk_cvsn_path_up_map <- bitr(geneID = sig_bulk_cvsn_path_up$GeneSymbol,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = "org.Hs.eg.db") #now 4830
# 0.33% of input gene IDs are fail to map...
head(sig_bulk_cvsn_path_up_map) #these would be the correct differentially expressed genes
# SYMBOL ENTREZID
# 2  NEAT1   283131
# 3  MUC16    94025
# 4    CLU     1191
# 5   CD74      972
# 6  MUC5B   727897
# 7   VCAN     1462

#for the background, we would need all the genes from the bulk DEG list object which are 21487
background_bulk_genes_cvsn_path_map <- bitr(geneID = res_df_cvsn_table$GeneSymbol,
                                            fromType = "SYMBOL",
                                            toType = "ENTREZID",
                                            OrgDb = "org.Hs.eg.db")
#0.41% of input gene IDs are fail to map..., now 21401
background_all <- intersect(background_sc_cvsn_path_map$ENTREZID, background_bulk_genes_cvsn_path_map$ENTREZID)

###sc
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
msig_Hpathways_sc_upgenes <- enricher(sig_sc_cvsn_path_up_map$ENTREZID, TERM2GENE=m_Hpathways, 
                                      universe = background_sc_cvsn_path_map$ENTREZID)
head(msig_Hpathways_sc_upgenes)
# ID                      Description GeneRatio  BgRatio RichFactor FoldEnrichment   zScore
# HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB    36/277 199/4316  0.1809045       2.818715 6.878504
# HALLMARK_INFLAMMATORY_RESPONSE     HALLMARK_INFLAMMATORY_RESPONSE   HALLMARK_INFLAMMATORY_RESPONSE    36/277 200/4316  0.1800000       2.804621 6.843160
# HALLMARK_ALLOGRAFT_REJECTION         HALLMARK_ALLOGRAFT_REJECTION     HALLMARK_ALLOGRAFT_REJECTION    29/277 195/4316  0.1487179       2.317208 4.929056
# HALLMARK_ESTROGEN_RESPONSE_LATE   HALLMARK_ESTROGEN_RESPONSE_LATE  HALLMARK_ESTROGEN_RESPONSE_LATE    29/277 197/4316  0.1472081       2.293683 4.866968
# HALLMARK_IL6_JAK_STAT3_SIGNALING HALLMARK_IL6_JAK_STAT3_SIGNALING HALLMARK_IL6_JAK_STAT3_SIGNALING    17/277  87/4316  0.1954023       3.044608 5.044802
# HALLMARK_ESTROGEN_RESPONSE_EARLY HALLMARK_ESTROGEN_RESPONSE_EARLY HALLMARK_ESTROGEN_RESPONSE_EARLY    28/277 198/4316  0.1414141       2.203406 4.539363

msig_Hpathways_sc_upgenes_df <- msig_Hpathways_sc_upgenes@result #45 pathways
write.csv(msig_Hpathways_sc_upgenes_df, paste0(outdir,'mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerVSLowerNormal_scRNAseq.csv'))

pdf(paste0(outdir,'barPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerVSLowerNormal_scRNAseq-short.pdf'),width = 12, height = 8)
barplot(msig_Hpathways_sc_upgenes)
dev.off()

pdf(paste0(outdir,'dotPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerVSLowerNormal_scRNAseq.pdf'),width = 12, height = 12)
dotplot(msig_Hpathways_sc_upgenes)
dev.off()

###bulk 
#H pathways
msig_Hpathways_bulk_upgenes <- enricher(sig_bulk_cvsn_path_up_map$ENTREZID, TERM2GENE=m_Hpathways, 
                                        universe = background_bulk_genes_cvsn_path_map$ENTREZID)
head(msig_Hpathways_bulk_upgenes)
# ID                        Description GeneRatio  BgRatio RichFactor
# HALLMARK_ALLOGRAFT_REJECTION             HALLMARK_ALLOGRAFT_REJECTION       HALLMARK_ALLOGRAFT_REJECTION  128/1250 196/4310  0.6530612
# HALLMARK_INFLAMMATORY_RESPONSE         HALLMARK_INFLAMMATORY_RESPONSE     HALLMARK_INFLAMMATORY_RESPONSE  128/1250 199/4310  0.6432161
# HALLMARK_TNFA_SIGNALING_VIA_NFKB     HALLMARK_TNFA_SIGNALING_VIA_NFKB   HALLMARK_TNFA_SIGNALING_VIA_NFKB  125/1250 198/4310  0.6313131
# HALLMARK_INTERFERON_GAMMA_RESPONSE HALLMARK_INTERFERON_GAMMA_RESPONSE HALLMARK_INTERFERON_GAMMA_RESPONSE  105/1250 199/4310  0.5276382
# HALLMARK_ESTROGEN_RESPONSE_LATE       HALLMARK_ESTROGEN_RESPONSE_LATE    HALLMARK_ESTROGEN_RESPONSE_LATE   96/1250 200/4310  0.4800000
# HALLMARK_KRAS_SIGNALING_UP                 HALLMARK_KRAS_SIGNALING_UP         HALLMARK_KRAS_SIGNALING_UP   95/1250 198/4310  0.4797980

msig_Hpathways_bulk_upgenes_df <- msig_Hpathways_bulk_upgenes@result #50 pathways
write.csv(msig_Hpathways_bulk_upgenes_df, paste0(outdir,'mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerousVSLowerNon-Cancerous_bulkRNAseq.csv'))

pdf(paste0(outdir, 'barPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerousVSLowerNon-Cancerous_bulkRNAseq-short.pdf'),width = 12, height = 8)
barplot(msig_Hpathways_bulk_upgenes)
dev.off()

pdf(paste0(outdir,'dotPlot-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerousVSLowerNon-Cancerous_bulkRNAseq.pdf'),width = 12, height = 12)
dotplot(msig_Hpathways_bulk_upgenes)
dev.off()

#intersection of pathways 
sig_Hpath_sc_upgenes <- msig_Hpathways_sc_upgenes_df %>% filter(p.adjust <= 0.05) #12
sig_Hpath_bulk_upgenes <- msig_Hpathways_bulk_upgenes_df %>% filter(p.adjust <= 0.05) #17
intersect_Hpath_up_scandbulk <- intersect(sig_Hpath_sc_upgenes$Description, sig_Hpath_bulk_upgenes$Description) #11 shared pathways
#sc 11/12 = 92% , bulk 11/17 = 65% of pathways for each are shared 
write.csv(intersect_Hpath_up_scandbulk, paste0(outdir,'Significant-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerousVSLowerNon-Cancerous_intersectionscANDbulk.csv'))

#create a new table with FC and adj pval for all the matched pathways 
intersect_Hpath_up_cvsn_sc_match <- match(intersect_Hpath_up_scandbulk, sig_Hpath_sc_upgenes$Description)

data_Hpath_up_cvsn_sc_match <- data.frame(sc_Hpath = intersect_Hpath_up_scandbulk, sc_FoldEnrichment = sig_Hpath_sc_upgenes$FoldEnrichment[intersect_Hpath_up_cvsn_sc_match], 
                                          sc_p_val_adj = sig_Hpath_sc_upgenes$p.adjust[intersect_Hpath_up_cvsn_sc_match], stringsAsFactors = FALSE)

instersect_Hpath_up_cvsn_bulk_match <- match(intersect_Hpath_up_scandbulk, sig_Hpath_bulk_upgenes$Description)

data_Hpath_up_cvsn_bulk_match <- data.frame(bulk_Hpath = intersect_Hpath_up_scandbulk, bulk_FoldEnrichment = sig_Hpath_bulk_upgenes$FoldEnrichment[instersect_Hpath_up_cvsn_bulk_match], 
                                            bulk_p_val_adj = sig_Hpath_bulk_upgenes$p.adjust[instersect_Hpath_up_cvsn_bulk_match], stringsAsFactors = FALSE)

Hpath_up_intersect_cvsn_all_data <- cbind(data_Hpath_up_cvsn_sc_match, data_Hpath_up_cvsn_bulk_match)
write.csv(Hpath_up_intersect_cvsn_all_data, paste0(outdir, 'Significant-mSigHpathways-DEG0.05PvalueANDhigher1FC-UPREGULATED_HigherCancerousVSLowerNon-Cancerous_intersectionscANDbulk-withDATA.csv'))

####pathway analysis sc and bulk, cancer vs normal DOWNREGULATED log2FC < -1 #### 
#we want significant genes with values (padj <= 0.05, log2FC of > or < 0) and background genes (every gene in the differential expression analysis) 
sig_sc_cvsn_path_down <- sig_sc_cvsn %>% filter(avg_log2FC <= -1) #8575 genes downregulated and significant
sig_bulk_cvsn_path_down <- sig_bulk_cvsn %>% filter(log2FoldChange <= -1) #1641

###sc
#need to conver genesymbol to ENTREZID
sig_sc_cvsn_path_down_map <- bitr(geneID = sig_sc_cvsn_path_down$gene,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = "org.Hs.eg.db") #now 7199
# 16.05% of input gene IDs are fail to map..., not ideal but will work for now
head(sig_sc_cvsn_path_down_map) #these would be the correct differentially expressed genes
# SYMBOL ENTREZID
# 1    IGKC     3514
# 2    CD74      972
# 3   IGLC2     3538
# 4 HLA-DRA     3122
# 5    RGS1     5996
# 6   WFDC2    10406

#for the background, we would need all the genes from the scRNAseq object which are 33612
background_sc_cvsn_path_map <- bitr(geneID = rownames(ovcaH),
                                    fromType = "SYMBOL",
                                    toType = "ENTREZID",
                                    OrgDb = "org.Hs.eg.db")
#34.7% of input gene IDs are fail to map..., now 21950 genes  

###bulk
sig_bulk_cvsn_path_down_map <- bitr(geneID = sig_bulk_cvsn_path_down$GeneSymbol,
                                    fromType = "SYMBOL",
                                    toType = "ENTREZID",
                                    OrgDb = "org.Hs.eg.db") #now 1632
# 0.55% of input gene IDs are fail to map...
head(sig_bulk_cvsn_path_down_map) #these would be the correct differentially expressed genes
# SYMBOL ENTREZID
# 2  NEAT1   283131
# 3  MUC16    94025
# 4    CLU     1191
# 5   CD74      972
# 6  MUC5B   727897
# 7   VCAN     1462

#for the background, we would need all the genes from the bulk DEG list object which are 21487
background_bulk_genes_cvsn_path_map <- bitr(geneID = res_df_cvsn_table$GeneSymbol,
                                            fromType = "SYMBOL",
                                            toType = "ENTREZID",
                                            OrgDb = "org.Hs.eg.db")
#0.41% of input gene IDs are fail to map..., now 21401

##sc
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
msig_Hpathways_sc_downgenes <- enricher(sig_sc_cvsn_path_down_map$ENTREZID, TERM2GENE=m_Hpathways, 
                                        universe = background_sc_cvsn_path_map$ENTREZID)
head(msig_Hpathways_sc_downgenes)
# ID                                Description GeneRatio  BgRatio RichFactor
# HALLMARK_UV_RESPONSE_DN                                       HALLMARK_UV_RESPONSE_DN                    HALLMARK_UV_RESPONSE_DN   93/1578 143/4316  0.6503497
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION  118/1578 200/4316  0.5900000
# HALLMARK_ADIPOGENESIS                                           HALLMARK_ADIPOGENESIS                      HALLMARK_ADIPOGENESIS  106/1578 199/4316  0.5326633
# HALLMARK_FATTY_ACID_METABOLISM                         HALLMARK_FATTY_ACID_METABOLISM             HALLMARK_FATTY_ACID_METABOLISM   81/1578 156/4316  0.5192308
# HALLMARK_OXIDATIVE_PHOSPHORYLATION                 HALLMARK_OXIDATIVE_PHOSPHORYLATION         HALLMARK_OXIDATIVE_PHOSPHORYLATION   99/1578 200/4316  0.4950000
# HALLMARK_MYOGENESIS                                               HALLMARK_MYOGENESIS                        HALLMARK_MYOGENESIS   96/1578 197/4316  0.4873096

msig_Hpathways_sc_downgenes_df <- msig_Hpathways_sc_downgenes@result #50 pathways
write.csv(msig_Hpathways_sc_downgenes_df, paste0(outdir,'mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerVSLowerNormal_scRNAseq.csv'))

pdf(paste0(outdir,'barPlot-mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerVSLowerNormal_scRNAseq-short.pdf'),width = 12, height = 8)
barplot(msig_Hpathways_sc_downgenes)
dev.off()

pdf(paste0(outdir,'dotPlot-mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerVSLowerNormal_scRNAseq.pdf'),width = 12, height = 12)
dotplot(msig_Hpathways_sc_downgenes)
dev.off()

###bulk 
#H pathways
msig_Hpathways_bulk_downgenes <- enricher(sig_bulk_cvsn_path_down_map$ENTREZID, TERM2GENE=m_Hpathways, 
                                          universe = background_bulk_genes_cvsn_path_map$ENTREZID)
head(msig_Hpathways_bulk_downgenes)
# ID                                Description GeneRatio  BgRatio RichFactor
# HALLMARK_UV_RESPONSE_DN                                       HALLMARK_UV_RESPONSE_DN                    HALLMARK_UV_RESPONSE_DN    27/281 143/4310  0.1888112
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION    30/281 200/4310  0.1500000
# HALLMARK_APICAL_JUNCTION                                     HALLMARK_APICAL_JUNCTION                   HALLMARK_APICAL_JUNCTION    29/281 199/4310  0.1457286
# HALLMARK_MYOGENESIS                                               HALLMARK_MYOGENESIS                        HALLMARK_MYOGENESIS    26/281 192/4310  0.1354167
# HALLMARK_COAGULATION                                             HALLMARK_COAGULATION                       HALLMARK_COAGULATION    19/281 132/4310  0.1439394
# HALLMARK_BILE_ACID_METABOLISM                           HALLMARK_BILE_ACID_METABOLISM              HALLMARK_BILE_ACID_METABOLISM    15/281 108/4310  0.1388889

msig_Hpathways_bulk_downgenes_df <- msig_Hpathways_bulk_downgenes@result #50 pathways
write.csv(msig_Hpathways_bulk_downgenes_df, paste0(outdir,'mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerousVSLowerNon-Cancerous_bulkRNAseq.csv'))

pdf(paste0(outdir, 'barPlot-mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerousVSLowerNon-Cancerous_bulkRNAseq-short.pdf'),width = 12, height = 8)
barplot(msig_Hpathways_bulk_downgenes)
dev.off()

pdf(paste0(outdir,'dotPlot-mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerousVSLowerNon-Cancerous_bulkRNAseq.pdf'),width = 12, height = 12)
dotplot(msig_Hpathways_bulk_downgenes)
dev.off()

#intersection of pathways 
sig_Hpath_sc_downgenes <- msig_Hpathways_sc_downgenes_df %>% filter(p.adjust <= 0.05) #8
sig_Hpath_bulk_downgenes <- msig_Hpathways_bulk_downgenes_df %>% filter(p.adjust <= 0.05) #6
intersect_Hpath_down_scandbulk <- intersect(sig_Hpath_sc_downgenes$Description, sig_Hpath_bulk_downgenes$Description) #4 shared pathways
#sc 4/8 = 80% , bulk 4/6 = 66% of pathways for each are shared 
write.csv(intersect_Hpath_down_scandbulk, paste0(outdir,'Significant-mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerousVSLowerNon-Cancerous_intersectionscANDbulk.csv'))

#create a new table with FC and adj pval for all the matched pathways 
intersect_Hpath_down_cvsn_sc_match <- match(intersect_Hpath_down_scandbulk, sig_Hpath_sc_downgenes$Description)

data_Hpath_down_cvsn_sc_match <- data.frame(sc_Hpath = intersect_Hpath_down_scandbulk, sc_FoldEnrichment = sig_Hpath_sc_downgenes$FoldEnrichment[intersect_Hpath_down_cvsn_sc_match], 
                                            sc_p_val_adj = sig_Hpath_sc_downgenes$p.adjust[intersect_Hpath_down_cvsn_sc_match], stringsAsFactors = FALSE)

instersect_Hpath_down_cvsn_bulk_match <- match(intersect_Hpath_down_scandbulk, sig_Hpath_bulk_downgenes$Description)

data_Hpath_down_cvsn_bulk_match <- data.frame(bulk_Hpath = intersect_Hpath_down_scandbulk, bulk_FoldEnrichment = sig_Hpath_bulk_downgenes$FoldEnrichment[instersect_Hpath_down_cvsn_bulk_match], 
                                              bulk_p_val_adj = sig_Hpath_bulk_downgenes$p.adjust[instersect_Hpath_down_cvsn_bulk_match], stringsAsFactors = FALSE)

Hpath_down_intersect_cvsn_all_data <- cbind(data_Hpath_down_cvsn_sc_match, data_Hpath_down_cvsn_bulk_match)
write.csv(Hpath_down_intersect_cvsn_all_data, paste0(outdir, 'Significant-mSigHpathways-DEG0.05PvalueANDlower1FC-DOWNREGULATED_HigherCancerousVSLowerNon-Cancerous_intersectionscANDbulk-withDATA.csv'))
