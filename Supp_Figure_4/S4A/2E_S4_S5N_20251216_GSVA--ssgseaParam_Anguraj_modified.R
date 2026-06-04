################################################################################
# Figure Panels : 2E, Supp Figure 4, Supp Figure 5N
# Description  : GSVA ssGSEA parameter run across multiple cohorts. Earlier
#                version (Dec 2025) of the multi-cohort integration script.
#                Kept for reproducibility.
# Input data   : Cohort-specific expression matrices; kegg_legacy gmt
# Output       : Cohort-specific PDFs (scatter + boxplot)
# R packages   : GSVA, GSEABase, ggplot2, dplyr, cowplot
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : December 2025
################################################################################

pacman::p_load(magrittr, dplyr, tidyr, GSVA, GSEABase, ggplot2, ggthemes, ggprism, ggsci, reshape2)

# Colossus PDX
Colossus_PDX = read.table("GSE204805_vsd_LMX_BASALE_all_RNAseq.txt",
                          sep='\t', header=TRUE, row.names=1, check.names=FALSE)
# COLOSUSS processed Proxy
COLOSUSS_processed_Proxy = read.table("Supplementary data/clinical_proxy/20250521_COLOSUSS_processed_Proxy_1_2_trial_before_AFTER_correction.txt",
                                      sep='\t', header=TRUE, row.names=1, check.names=FALSE)
# COLOSUS RNAseq
COLOSUS_EXP = read.table("20200830_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_sd0_DESeq2_rlog_filtered_5_ALL_Combat_6_batches_corrected_Validation_sd0.txt",
                         sep='\t', check.names=FALSE, header=TRUE, row.names=1)
# COLOSUS proteomics
COLOSUS_protein = read.table("20200928_metabolomics_median_replicates_data_80_Non_zero_sample_median_centred_imputed_MARS_combat_lysis_batch_protein_names_only_sd0.txt",
                             sep='\t', header=TRUE, check.names=FALSE,, row.names=1) #%>%
 # mutate(Gene=gsub('.*_', '', proteins)) %>%
 # separate_rows(Gene, sep=';') %>%
 #dplyr::select(-proteins) %>%
  #aggregate(.~Gene, ., mean) %>%
 # tibble::column_to_rownames('Gene')
# TCGA
#TCGA_CRC = read.table("Figure 2 modification+validation/DATA/HiSeqV2.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE)

#MT_barcode = read.table("Figure 2 modification+validation/DATA/tcga_coadread_pancanatlas_RAS_mutant_MSS_TCGA_barcodes.txt",
 #                       sep='\t', header=TRUE, check.names=FALSE)
TCGA_MT = TCGA_CRC %>% dplyr::select(intersect(MT_barcode$TCGA_RAS_mutant_samples, colnames(.)))

#WT_barcode = read.table("Figure 2 modification+validation/DATA/tcga_coadread_pancanatlas_RAS_wildtype_MSS_TCGA_barcodes.txt",
 #                       sep='\t', header=TRUE, check.names=FALSE)
TCGA_WT = TCGA_CRC %>% dplyr::select(intersect(WT_barcode$TCGA_RAS_wt_samples, colnames(.)))

# Public metabolic data - Anguraj added
PUBLIC_RNA = read.table("Supplementary data/GSE89076/GSE89076.Agilent8x60K.Matched_expression_tumours_matched.csv",
                             sep=',', header=TRUE, check.names=FALSE,, row.names=1) 

# Public metabolic data - Anguraj added
CELLLINE_mut_RNA = read.table("20251220_OmicsExpressionExpectedCountHumanProteinCodingGenes_CRC_RAS_mutant_expectedCount_transposed_genesXsamples_deseq2_rlog_data.txt",
                        sep='\t', header=TRUE, check.names=FALSE,, row.names=1) 

# Public metabolic data - Anguraj added
CELLLINE_wt_RNA = read.table("20251220_OmicsExpressionExpectedCountHumanProteinCodingGenes_CRC_RAS_wildtype_expectedCount_transposed_genesXsamples_deseq2_rlog_data.txt",
                              sep='\t', header=TRUE, check.names=FALSE,, row.names=1) 

#Tumouroids RAS mut
Colossus_tumoroids_RAS_mut_all <- read.table("20251230_GSE204805_vsd_LMO_organoids_RAS_mutant_RNAseq_data.txt",sep='\t', header=TRUE, check.names=FALSE,, row.names=1)

#Tumouroids RAS wildtype
Colossus_tumoroids_RAS_wt_all <- read.table("20260322_GSE204805_vsd_LMO_organoids_RAS_wt_RNAseq_data.txt",sep='\t', header=TRUE, check.names=FALSE,, row.names=1)

#9 RAS mt tumoroids
tumoroids_9_RAS_mt <- read.table("S4I_20260322_LMO_9_organoids_for_metabolism_duplicate_mean_gene_exp.txt",sep='\t', header=TRUE, check.names=FALSE,row.names=1)

#9 RAS mt tumoroids
tumoroids_3_RAS_wt <- read.table("20251224_GSE204805_vsd_1331_322_399_LMO_organoids_Basale_genes.txt",sep='\t', header=TRUE, check.names=FALSE,row.names=1)
# Extract gene symbols (remove everything up to and including first underscore)
rownames(tumoroids_3_RAS_wt) <- sub("^.*?_", "", rownames(tumoroids_3_RAS_wt))

## MSigDB KEGG ####
ALANINE_ASPARTATE_AND_GLUTAMATE = readGMT("Figure 2 modification+validation/DATA/KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM.v2025.1.Hs.gmt")
  
PURINE_METABOLISM = readGMT("Figure 2 modification+validation/DATA/KEGG_PURINE_METABOLISM.v2025.1.Hs.gmt")

ARGININE_AND_PROLINE = readGMT("Figure 2 modification+validation/DATA/KEGG_ARGININE_AND_PROLINE_METABOLISM.v2025.1.Hs.gmt")

BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS = readGMT("Figure 2 modification+validation/DATA/KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS.v2025.1.Hs.gmt")

BUTANOATE = readGMT("Figure 2 modification+validation/DATA/KEGG_BUTANOATE_METABOLISM.v2025.1.Hs.gmt")

CITRATE_CYCLE_TCA_CYCLE = readGMT("Figure 2 modification+validation/DATA/KEGG_CITRATE_CYCLE_TCA_CYCLE.v2025.1.Hs.gmt")

GLYCINE_SERINE_AND_THREONINE = readGMT("Figure 2 modification+validation/DATA/KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM.v2025.1.Hs.gmt")

GLYOXYLATE_AND_DICARBOXYLATE = readGMT("Figure 2 modification+validation/DATA/KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM.v2025.1.Hs.gmt")

PYRUVATE = readGMT("Figure 2 modification+validation/DATA/KEGG_PYRUVATE_METABOLISM.v2025.1.Hs.gmt")

# GSEA analysis -----
# combine
MSigDB_KEGG = data.frame(pathway='ALANINE_ASPARTATE_AND_GLUTAMATE', gene=ALANINE_ASPARTATE_AND_GLUTAMATE$KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM) %>%
  rbind(data.frame(pathway='PURINE_METABOLISM', gene=PURINE_METABOLISM$KEGG_PURINE_METABOLISM)) %>%
  rbind(data.frame(pathway='ARGININE_AND_PROLINE', gene=ARGININE_AND_PROLINE$KEGG_ARGININE_AND_PROLINE_METABOLISM)) %>%
  rbind(data.frame(pathway='BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS', gene=BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS$KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS)) %>%
  rbind(data.frame(pathway='BUTANOATE', gene=BUTANOATE$KEGG_BUTANOATE_METABOLISM)) %>%
  rbind(data.frame(pathway='CITRATE_CYCLE_TCA_CYCLE', gene=CITRATE_CYCLE_TCA_CYCLE$KEGG_CITRATE_CYCLE_TCA_CYCLE)) %>%
  rbind(data.frame(pathway='GLYCINE_SERINE_AND_THREONINE', gene=GLYCINE_SERINE_AND_THREONINE$KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM)) %>%
  rbind(data.frame(pathway='GLYOXYLATE_AND_DICARBOXYLATE', gene=GLYOXYLATE_AND_DICARBOXYLATE$KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM)) %>%
  rbind(data.frame(pathway='PYRUVATE', gene=PYRUVATE$KEGG_PYRUVATE_METABOLISM)) %>%
  separate_rows(gene, sep=',') %>%
  as.data.frame()

pathway_order = c('ALANINE_ASPARTATE_AND_GLUTAMATE', 'GLYOXYLATE_AND_DICARBOXYLATE', 'BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS',
                  'BUTANOATE', 'CITRATE_CYCLE_TCA_CYCLE', 'ARGININE_AND_PROLINE', 'PYRUVATE', 'GLYCINE_SERINE_AND_THREONINE',
                  'PURINE_METABOLISM')

## Colossus tumoroids ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(Colossus_tumoroids), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_Colossus_tumoroids = gsva(gsvapar)
# export
write.table(gsva_Colossus_tumoroids, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/Colossus tumoroids-replicate_mean_GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)

## Colossus tumoroids all RAS mutant samples####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(Colossus_tumoroids_RAS_mut_all), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_Colossus_tumoroids_RAS_mut_all = gsva(gsvapar)
# export
#write.table(gsva_Colossus_tumoroids_RAS_mut_all, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/Colossus tumoroids-RAS_mutant_all_samples_replicate_mean_GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)

## Colossus tumoroids all RAS WILDTYPE samples####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(Colossus_tumoroids_RAS_wt_all), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_Colossus_tumoroids_RAS_wt_all = gsva(gsvapar)
# export
#write.table(gsva_Colossus_tumoroids_RAS_wt_all, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/Colossus tumoroids-RAS_wildtype_all_samples_replicate_mean_GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)

## Colossus tumoroids 9 selected RAS mt samples####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(tumoroids_9_RAS_mt), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_tumoroids_9_RAS_mt = gsva(gsvapar)
# export
#write.table(gsva_tumoroids_9_RAS_mt, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/Colossus tumoroids-9_selected_RAS_mutant_samples_replicate_mean_GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)

## Colossus tumoroids 3 selected RAS WILDTYPE samples####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(tumoroids_3_RAS_wt), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_tumoroids_e_RAS_wt = gsva(gsvapar)
# export
write.table(gsva_tumoroids_e_RAS_wt, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/Colossus tumoroids-3_selected_RAS_wildtype_samples_replicate_mean_GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)



## Colossus PDX ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(Colossus_PDX), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_Colossus_PDX = gsva(gsvapar)
# export
#write.table(gsva_Colossus_PDX, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/Colossus PDX-LMX_BASALE_ALL_GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)

## COLOSUSS processed Proxy ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(COLOSUSS_processed_Proxy), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
  gsva_COLOSUSS_processed_Proxy = gsva(gsvapar)
# export
#write.table(gsva_COLOSUSS_processed_Proxy, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/COLOSUSS processed Proxy-GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)

## COLOSUS RNAseq ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(COLOSUS_EXP), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_COLOSUS_EXP = gsva(gsvapar)
# export
#write.table(gsva_COLOSUS_EXP, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/COLOSUS RNAseq-GSEA_Anguraj_modified.txt", sep="\t", quote=FALSE, row.names=TRUE)

## COLOSUS proteomics ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(COLOSUS_protein), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_COLOSUS_protein = gsva(gsvapar)
# export
#write.table(gsva_COLOSUS_protein, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/COLOSUS proteomics-GSEA_Anguraj_modified.txt", sep="\t", quote=FALSE, row.names=TRUE)

## TCGA MT ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
gene.sets <- gene.sets[!names(gene.sets) %in% "BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS"]
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(TCGA_MT), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_TCGA_MT = gsva(gsvapar)
# export
#write.table(gsva_TCGA_MT, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/TCGA MT-GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)

## TCGA WT ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(TCGA_WT), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_TCGA_WT = gsva(gsvapar)
# export
#write.table(gsva_TCGA_WT, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/TCGA WT-GSEA.txt", sep="\t", quote=FALSE, row.names=TRUE)


## PUBLIC RNAseq ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(PUBLIC_RNA), geneSets=gene.sets)
# estimate GSVA enrichment scores for the background genesets
gsva_PUBLIC_RNA = gsva(gsvapar)
# export
#write.table(gsva_PUBLIC_RNA, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/PUBLIC RNAseq-GSEA_Anguraj_modified.txt", sep="\t", quote=FALSE, row.names=TRUE)


## CELL LINES RAS mutant RNAseq ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])

gene.sets1 <- gene.sets
#NOTICE CITRATE CYCLE - removed one gene - "SDHCP5"
gene.sets1$CITRATE_CYCLE_TCA_CYCLE <- gene.sets$CITRATE_CYCLE_TCA_CYCLE[-27]
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(CELLLINE_mut_RNA), geneSets=gene.sets1)
# estimate GSVA enrichment scores for the background genesets
gsva_CELLLINE_mut_RNA = gsva(gsvapar)
# export
write.table(gsva_CELLLINE_mut_RNA, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/CELLLINE_RAS_mut_RNAseq-GSEA_Anguraj_modified.txt", sep="\t", quote=FALSE, row.names=TRUE)

## CELL LINES RAS wildtype RNAseq ####
# background gene set
gene.sets = split(MSigDB_KEGG[, 2], MSigDB_KEGG[, 1])
gene.sets1 <- gene.sets
#NOTICE CITRATE CYCLE - removed one gene - "SDHCP5"
gene.sets1$CITRATE_CYCLE_TCA_CYCLE <- gene.sets$CITRATE_CYCLE_TCA_CYCLE[-27]
# build GSVA parameter object
gsvapar = ssgseaParam(exprData=as.matrix(CELLLINE_wt_RNA), geneSets=gene.sets1)
# estimate GSVA enrichment scores for the background genesets
gsva_CELLLINE_wt_RNA = gsva(gsvapar)
# export
write.table(gsva_CELLLINE_wt_RNA, "Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/CELLLINE_RAS_wt_RNAseq-GSEA_Anguraj_modified.txt", sep="\t", quote=FALSE, row.names=TRUE)



#------------------------------------------- boxplot -------------------------------------------####
## Colossus tumoroids -----
# data prepare
boxplot_data = gsva_Colossus_tumoroids %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p1 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p1

# Save output figure.
ggsave(filename="Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/boxplot--Colossus tumoroids_replicate_mean.pdf", plot=p1, width=6, height=15)


## Colossus tumoroids RAS mutant all samples -----
# data prepare
boxplot_data = gsva_Colossus_tumoroids_RAS_mut_all %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p1.1 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p1.1

write.table(boxplot_data,"S4G_tumoroids_RAS_mutant_all_replicate_mean_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4G_boxplot_tumoroids_RAS_mutant_all_replicate_mean.pdf", plot=p1.1, width=6, height=15)



## Colossus tumoroids RAS wildtype all samples -----
# data prepare
boxplot_data = gsva_Colossus_tumoroids_RAS_wt_all %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p1.2 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p1.2

write.table(boxplot_data,"S4H_tumoroids_RAS_wt_all_replicate_mean_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4H_boxplot_tumoroids_RAS_wt_all_replicate_mean.pdf", plot=p1.2, width=6, height=15)


## Colossus tumoroids 9 selected RAS mt all samples -----
# data prepare
boxplot_data = gsva_tumoroids_9_RAS_mt %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p1.3 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p1.3

write.table(boxplot_data,"S4I_20260322_LMO_9_organoids_for_metabolism_duplicate_mean_gene_exp_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4I_boxplot_9_tumoroids_RAS_mt_selected_replicate_mean.pdf", plot=p1.3, width=6, height=15)


## Colossus tumoroids 3 selected RAS wt all samples -----
# data prepare
boxplot_data = gsva_tumoroids_e_RAS_wt %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p1.4 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p1.4

write.table(boxplot_data,"S4J_20260322_LMO_3_RAS_wt_organoids_for_metabolism_duplicate_mean_gene_exp_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4J_boxplot_3_tumoroids_RAS_wt_selected_replicate_mean.pdf", plot=p1.4, width=6, height=15)



## Colossus PDX -----
# data prepare
boxplot_data = gsva_Colossus_PDX %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p2 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p2

write.table(boxplot_data,"S4F_PDX_LMX_BASALE_All_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4F_boxplot_PDX_LMX_BASALE_All.pdf", plot=p2, width=6, height=15)

## COLOSUSS processed Proxy -----
# data prepare
boxplot_data = gsva_COLOSUSS_processed_Proxy %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p3 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p3

write.table(boxplot_data,"S4E_COLOSUSS_processed_Proxy_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4E_boxplot_COLOSUSS_processed_Proxy_GSVA.pdf", plot=p3, width=6, height=15)

## COLOSUS RNAseq -----
# data prepare
boxplot_data = gsva_COLOSUS_EXP %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p4 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p4

ggsave(filename="2E_boxplot_ssGSEA_COLOSUS RNAseq.pdf", plot=p4, width=6, height=15)
write.table(boxplot_data,"2E_20260322_COLOSSUS_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


## COLOSUS proteomics -----
# data prepare
boxplot_data = gsva_COLOSUS_protein %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p5 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p5

write.table(boxplot_data,"2F_20260322_COLOSSUS_Proteomics_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="2F_boxplot-GSVA_COLOSUS proteomics.pdf", plot=p5, width=6, height=15)

## TCGA MT -----
# data prepare
boxplot_data = gsva_TCGA_MT %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p6 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p6

write.table(boxplot_data,"S4C_20260322_TCGA MT_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4C_boxplot_TCGA MT_GSVA.pdf", plot=p6, width=6, height=15)

## TCGA WT -----
# data prepare
boxplot_data = gsva_TCGA_WT %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)

# plot
p7 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p7

write.table(boxplot_data,"S4D_20260322_TCGA WT_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4D_boxplot_TCGA WT_GSVA.pdf", plot=p7, width=6, height=15)


## PUBLIC RNAseq -----
# data prepare
boxplot_data = gsva_PUBLIC_RNA %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)


# plot
p8 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p8

write.table(boxplot_data,"S4B_20260322_Public_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores_Rachel.txt", sep = "\t", quote = FALSE, row.names = FALSE)


ggsave(filename="S4B_boxplot--PUBLIC RNAseq.pdf", plot=p8, width=6, height=15)

## CELL LINE RAS Mutant RNAseq -----
# data prepare
boxplot_data = gsva_CELLLINE_mut_RNA %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)


# plot
p9 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p9

ggsave(filename="Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/boxplot--PUBLIC RNAseq.pdf", plot=p4, width=6, height=15)

## CELL LINE RAS Wildtype RNAseq -----
# data prepare
boxplot_data = gsva_CELLLINE_wt_RNA %>%
  melt() %>%
  set_colnames(c('Pathway', 'Sample', 'Score'))

boxplot_data$Pathway = factor(boxplot_data$Pathway, levels=pathway_order)


# plot
p10 = ggplot(data=boxplot_data, mapping=aes(x=Pathway, y=Score, fill=Pathway)) +
  geom_boxplot(outlier.shape=21, color='black') +
  # scale_fill_manual(values=c('#6CA6CD', '#CD5C5C', 'navyblue')) +
  theme_bw() +
  labs(x=NULL, y='GSVA Gaussian score') +
  scale_color_npg() + scale_fill_npg() +
  theme(legend.position='none',
        axis.text.x=element_text(angle=90, hjust=1),
        axis.text=element_text(color='black', size=12))

p10

ggsave(filename="Figure 2 modification+validation/GSVA Gaussian+ssGSEA/output/ssgseaParam/data/Anguraj/boxplot--PUBLIC RNAseq.pdf", plot=p4, width=6, height=15)
