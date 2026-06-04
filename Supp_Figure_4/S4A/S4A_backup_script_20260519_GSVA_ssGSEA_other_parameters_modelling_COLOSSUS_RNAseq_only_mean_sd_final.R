################################################################################
# Figure Panel : Supplementary Figure 4A (GSVA backup)
# Description  : Alternative scoring using GSVA ssGSEA parameter for comparison
#                against singscore. Kept for reproducibility audit.
# Input data   : COLOSSUS RNA-seq DESeq2 rlog matrix; kegg_legacy gmt
# Output       : S4A_..._ssGSEA_violin_plot.pdf
# R packages   : GSVA, GSEABase, ggplot2, dplyr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : Primary figure uses singscore (S4A_new_script). This is the
#                GSVA backup for methods comparison.
################################################################################

## install if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GSVA")
# BiocManager::install("GSEABase")

# GSVA backup script for S4A — kept for methods comparison with singscore.
# Uses GSVA ssGSEA parameter on COLOSSUS RNA-seq.
# Primary figure uses singscore (S4A_new_script_singscore_analysis.R).
library(GSVA)
library(GSEABase)

# ── Input file paths ─────────────────────────────────────────────────────────
expr_file <- "20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt"

gmt_file  <- "~/Downloads/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

out_gct   <- "20201014_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_sd0_DESeq2_rlog_filtered_5_ALL_Combat_6_batches_corrected_Validation_sd0_GSVA_ssGSEA_results_KEGG_Legacy.gct"

expr_df <- read.delim(expr_file,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

## If first column is gene names and the rest are samples:
#rownames(expr_df) <- expr_df[, 1]
#expr_df <- expr_df[, -1, drop = FALSE]

expr_mat <- as.matrix(expr_df)

#Read KEGG gene sets
gmt <- getGmt(gmt_file)         # GeneSetCollection
gs_list <- geneIds(gmt)         # named list: pathway → genes

#Run GSVA in ssGSEA mode (to mimic your previous ssGSEA)
gs_list <- geneIds(gmt) 

# ssGSEA parameter set — same configuration as Fig 2E for consistency.
ssgsea_par <- ssgseaParam(
  exprData = expr_mat,
  geneSets = gs_list,
  minSize  = 1,
  maxSize  = Inf,
  alpha    = 0.25,
  normalize = FALSE
)

#ssgsea_par <- gsvaParam(
 # exprData = expr_mat,
  #geneSets = gs_list,
  #minSize  = 1,
  #maxSize  = Inf,
#)

#ssgsea_par <- plageParam(
 # exprData = expr_mat,
  #geneSets = gs_list,
  #minSize  = 1,
  #maxSize  = Inf,
#)

#ssgsea_par <- zscoreParam(
 # exprData = expr_mat,
  #geneSets = gs_list,
  #minSize  = 1,
  #maxSize  = Inf,
#)


# Run ssGSEA on COLOSSUS RNA-seq.
gsva_ssgsea <- gsva(ssgsea_par, verbose = FALSE)

# This is your ssGSEA-like matrix
GSVA_RNAseq_m <- gsva_ssgsea
dim(GSVA_RNAseq_m)

GSVA_RNAseq_sig <- c( GSVA_RNAseq_m[5,], GSVA_RNAseq_m[76,],GSVA_RNAseq_m[26,],GSVA_RNAseq_m[28,],GSVA_RNAseq_m[37,],GSVA_RNAseq_m[16,],GSVA_RNAseq_m[145,],GSVA_RNAseq_m[66,],GSVA_RNAseq_m[143,])
names <- c("Alanine, aspartate and glutamate metabolism","Glyoxylate and dicarboxylate metabolism", "Biosynthesis of unsaturated fatty acids","Butanoate metabolism","Citrate cycle & TCA cycle","Arginine and proline metabolism","Pyruvate metabolism", "Glycine, serine and threonine metabolism","Purine metabolism")
length(GSVA_RNAseq_sig)

cl <- rep(names, each = 134)
length(cl)

boxplot(GSVA_RNAseq_sig~factor(cl))

df <- data.frame(cbind(GSVA_RNAseq_sig),cl)
colnames(df) <- c("GSVA_RNAseq_sig","class")

# Explicitly define the order of levels in 'class'
df$class <- factor(df$class, levels = c(
  "Alanine, aspartate and glutamate metabolism",
  "Glyoxylate and dicarboxylate metabolism",
  "Biosynthesis of unsaturated fatty acids",
  "Butanoate metabolism",
  "Citrate cycle & TCA cycle",
  "Arginine and proline metabolism",
  "Pyruvate metabolism",
  "Glycine, serine and threonine metabolism",
  "Purine metabolism")
)

library(ggplot2)

#pdf("20241014_COLOSSUS_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot.pdf")
ggplot(df, aes(x=class, y=GSVA_RNAseq_sig)) + geom_violin(aes(fill = class), trim = FALSE) + 
  geom_boxplot(width = 0.2)
#dev.off()

#mean of mean and sd across the selected pathways
## ---- prerequisites --------------------------------------------------------
## Assumes you already have:
## expr_mat : numeric matrix (genes x samples), rownames = gene symbols
## gs_list  : named list of KEGG pathways (e.g. KEGG_PYRUVATE_METABOLISM -> genes)

kegg_selected <- c(
  "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM",
  "KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM",
  "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
  "KEGG_BUTANOATE_METABOLISM",
  "KEGG_CITRATE_CYCLE_TCA_CYCLE",
  "KEGG_ARGININE_AND_PROLINE_METABOLISM",
  "KEGG_PYRUVATE_METABOLISM",
  "KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM",
  "KEGG_PURINE_METABOLISM"
)

## Keep only those KEGG IDs that exist in gs_list
kegg_selected <- kegg_selected[kegg_selected %in% names(gs_list)]
kegg_selected

################################################################################
# Define Custom Color Palette (matching your reference figure)
################################################################################

# Pathway-specific colors matching your figure
pathway_colors <- c(
  "Alanine, aspartate and\nglutamate metabolism" = "#E76F51",      # Red/coral
  "Glyoxylate and\ndicarboxylate metabolism" = "#F4A261",          # Orange
  "Biosynthesis of\nunsaturated fatty acids" = "#A9A965",          # Olive/yellow-green
  "Butanoate metabolism" = "#2A9D8F",                               # Teal-green
  "Citrate cycle\n(TCA cycle)" = "#4ECDC4",                         # Cyan/teal
  "Arginine and\nproline metabolism" = "#1982C4",                   # Blue
  "Pyruvate metabolism" = "#6A9FB5",                                # Light blue
  "Glycine, serine and\nthreonine metabolism" = "gray50",         # Pink/magenta                           # Yellow (for new pathway)
  "Purine metabolism" = "gray50"                                   # Purple (for new pathway)
)

## ---------------------- Pathway Summary Stats --------------------------

#Global thresholds: what is “low” and “high” expression?

# Mean expression per gene across all samples
all_gene_means <- rowMeans(expr_mat)

# Global thresholds (you can tweak these)
low_thresh  <- quantile(all_gene_means, 0.25)  # bottom 25% = "low"
high_thresh <- quantile(all_gene_means, 0.90)  # top 10%   = "very high"

#For each pathway, summarise the gene mean distribution and other metrics
summarise_pathway_expression <- function(pathway_id, expr_mat, gs_list,
                                         low_thresh, high_thresh) {
  genes <- intersect(gs_list[[pathway_id]], rownames(expr_mat))
  if (length(genes) == 0) return(NULL)
  
  submat <- expr_mat[genes, , drop = FALSE]
  gene_means <- rowMeans(submat)
  
  data.frame(
    pathway            = pathway_id,
    n_genes            = length(genes),
    mean_of_means      = mean(gene_means),
    median_of_means    = median(gene_means),
    prop_low_genes     = mean(gene_means <= low_thresh),
    prop_high_genes    = mean(gene_means >= high_thresh),
    min_mean           = min(gene_means),
    max_mean           = max(gene_means),
    range_mean         = max(gene_means) - min(gene_means)
  )
}

#Apply this to your selected KEGG pathways:
pathway_expr_summary_list <- lapply(kegg_selected, summarise_pathway_expression,
                                    expr_mat = expr_mat,
                                    gs_list = gs_list,
                                    low_thresh = low_thresh,
                                    high_thresh = high_thresh)

pathway_expr_summary <- do.call(rbind, pathway_expr_summary_list)
pathway_expr_summary

#Visualise distributions per pathway (systematic, ggplot)
library(dplyr)
library(ggplot2)

gene_means_df <- lapply(kegg_selected, function(pw) {
  genes <- intersect(gs_list[[pw]], rownames(expr_mat))
  if (length(genes) == 0) return(NULL)
  submat <- expr_mat[genes, , drop = FALSE]
  gene_means <- rowMeans(submat)
  
  data.frame(
    pathway   = pw,
    gene      = genes,
    gene_mean = gene_means
  )
}) %>%
  bind_rows()

#A. Violin/boxplot of gene means per pathway
# Order pathways by median gene mean
gene_means_df$pathway <- factor(
  gene_means_df$pathway,
  levels = gene_means_df %>%
    group_by(pathway) %>%
    summarise(median_mean = median(gene_mean)) %>%
    arrange(desc(median_mean)) %>%
    pull(pathway)
)


ggplot(gene_means_df, aes(x = pathway, y = gene_mean)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  #coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Distribution of gene mean expression per KEGG pathway",
    x = "KEGG Pathway",
    y = "Gene mean expression across samples"
  )

#Overlay global thresholds (to “show” low vs high genes)

pdf("20260104_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_mean_of_genes_across_pathways_violin_plot.pdf")
ggplot(gene_means_df, aes(x = pathway, y = gene_mean)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  geom_hline(yintercept = low_thresh, linetype = "dashed") +
  geom_hline(yintercept = high_thresh, linetype = "dotted") +
#  coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Gene mean expression per pathway\nDashed = low threshold, dotted = high threshold",
    x = "KEGG Pathway",
    y = "Gene mean expression"
  )
dev.off()

# Median analysis

# Median expression per gene across all samples
all_gene_median <- apply(expr_mat,1, median)

# Global thresholds (you can tweak these)
low_thresh_median  <- quantile(all_gene_median, 0.25)  # bottom 25% = "low"
high_thresh_median <- quantile(all_gene_median, 0.90)  # top 10%   = "very high"


#Visualise distributions per pathway (systematic, ggplot)
library(dplyr)
library(ggplot2)

gene_median_df <- lapply(kegg_selected, function(pw) {
  genes <- intersect(gs_list[[pw]], rownames(expr_mat))
  if (length(genes) == 0) return(NULL)
  submat <- expr_mat[genes, , drop = FALSE]
  gene_median <- apply(submat, 1, median, na.rm = TRUE)
  
  data.frame(
    pathway     = pw,
    gene        = genes,
    gene_median = gene_median
  )
}) %>%
  bind_rows()

#A. Violin/boxplot of gene medians per pathway
# Order pathways by "median of gene medians" — HIGHEST FIRST
gene_median_df$pathway <- factor(
  gene_median_df$pathway,
  levels = gene_median_df %>%
    group_by(pathway) %>%
    summarise(median_of_medians = median(gene_median, na.rm = TRUE)) %>%
    arrange(desc(median_of_medians)) %>%   # highest median-of-medians first
    pull(pathway)
)

ggplot(gene_median_df, aes(x = pathway, y = gene_median)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Distribution of gene median expression per KEGG pathway",
    x = "KEGG Pathway (highest median-of-medians first)",
    y = "Gene median expression across samples"
  )


#Overlay global thresholds (to “show” low vs high genes)

pdf("20260104_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_median_of_genes_across_pathways_violin_plot.pdf")


ggplot(gene_median_df, aes(x = pathway, y = gene_median)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  geom_hline(yintercept = low_thresh_median, linetype = "dashed") +
  geom_hline(yintercept = high_thresh_median, linetype = "dotted") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Gene median expression per pathway\nDashed = low threshold, dotted = high threshold",
    x = "KEGG Pathway (highest median-of-medians first)",
    y = "Gene median expression"
  )

dev.off()

# -------------------------------------------------------------------------
# SD of genes across samples per selected KEGG pathway (gene_sd_df)
# -------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

gene_sd_df <- lapply(kegg_selected, function(pw) {
  genes <- intersect(gs_list[[pw]], rownames(expr_mat))
  if (length(genes) == 0) return(NULL)
  submat <- expr_mat[genes, , drop = FALSE]
  
  # SD per gene across all samples
  gene_sds <- apply(submat, 1, sd, na.rm = TRUE)
  
  data.frame(
    pathway  = pw,
    gene     = genes,
    gene_sd  = gene_sds
  )
}) %>%
  bind_rows()

# Order pathways by median gene SD (same style as gene_means_df ordering)
gene_sd_df$pathway <- factor(
  gene_sd_df$pathway,
  levels = gene_sd_df %>%
    group_by(pathway) %>%
    summarise(median_sd = median(gene_sd)) %>%
    arrange(median_sd) %>%
    pull(pathway)
)

# A. Violin/boxplot of gene SDs per pathway (same style as mean plot)
ggplot(gene_sd_df, aes(x = pathway, y = gene_sd)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Distribution of gene SD (variability) per KEGG pathway",
    x = "KEGG Pathway",
    y = "Gene SD across samples"
  )

## ---------------------------
## 4. Rank genes (singscore)
## ---------------------------
library(singscore)
# rankGenes expects genes in rows, samples in columns
rank_mat <- rankGenes(expr_mat)

## ---------------------------
## 5. Compute singscores for each KEGG pathway
## ---------------------------

# For KEGG, we typically just use them as "up" gene sets (no downSet)
scores_list <- lapply(gs_list, function(genes) {
  common_genes <- intersect(rownames(rank_mat), genes)
  
  if (length(common_genes) == 0) {
    # no overlap → NAs
    return(rep(NA_real_, ncol(rank_mat)))
  }
  
  sc <- simpleScore(rank_mat,
                    upSet       = common_genes,
                    centerScore = TRUE)  # centered total score
  
  # simpleScore returns a list; we take the per-sample total scores
  sc$TotalScore
})

# Combine into a matrix: pathways x samples
singscore_mat <- do.call(rbind, scores_list)
rownames(singscore_mat) <- names(gs_list)
colnames(singscore_mat) <- colnames(expr_mat)

# This is the direct analogue of your ssGSEA_RNAseq_m / singscore_RNAseq_m
singscore_RNAseq_m <- singscore_mat

singscore_RNAseq_sig <- c( singscore_RNAseq_m[5,], singscore_RNAseq_m[76,],singscore_RNAseq_m[26,],singscore_RNAseq_m[28,],singscore_RNAseq_m[37,],singscore_RNAseq_m[16,],singscore_RNAseq_m[145,],singscore_RNAseq_m[66,],singscore_RNAseq_m[143,])
names <- c("Alanine, aspartate and glutamate metabolism","Glyoxylate and dicarboxylate metabolism", "Biosynthesis of unsaturated fatty acids","Butanoate metabolism","Citrate cycle & TCA cycle","Arginine and proline metabolism","Pyruvate metabolism","Glcine, serine and threonine metabolism", "Purine metabolism")
length(singscore_RNAseq_sig)

cl <- rep(names, each = 134)
length(cl)

boxplot(singscore_RNAseq_sig~factor(cl))

df <- data.frame(cbind(singscore_RNAseq_sig),cl)
colnames(df) <- c("singscore_RNAseq_sig","class")

# Explicitly define the order of levels in 'class'
df$class <- factor(df$class, levels = c(
  "Alanine, aspartate and glutamate metabolism",
  "Glyoxylate and dicarboxylate metabolism",
  "Biosynthesis of unsaturated fatty acids",
  "Butanoate metabolism",
  "Citrate cycle & TCA cycle",
  "Arginine and proline metabolism",
  "Pyruvate metabolism",
  "Glcine, serine and threonine metabolism",
  "Purine metabolism"
))

library(ggplot2)

ggplot(df, aes(x=class, y=singscore_RNAseq_sig)) + geom_violin(aes(fill = class), trim = FALSE) + 
  geom_boxplot(width = 0.2)



# Create box/violin plot
# Create box/violin plot with pathway colors
panel_b <- ggplot(df, aes(x = class, y = singscore_RNAseq_sig)) +
  geom_violin(aes(fill = class), trim = FALSE, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = 16, outlier.size = 1, 
               alpha = 0.8, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_manual(
    values = pathway_colors,
    name = "Pathway"
  ) +
  labs(
    x = NULL,
    y = "sigScore Enrichment Score",
    title = "sigScore Distribution Across Metabolically Dysregulated Pathways",
    subtitle = sprintf("MSS RAS Mutant CRC (n=%d samples)", ncol(singscore_RNAseq_sig))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    legend.position = "none",  # Hide legend as colors already on x-axis
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

library(dplyr)

df <- df %>%
  mutate(
    class = factor(
      class,
      levels = df %>%
        group_by(class) %>%
        summarise(med = median(singscore_RNAseq_sig, na.rm = TRUE)) %>%
        arrange(desc(med)) %>%   # ⬅ reverse order
        pull(class)
    )
  )

panel_b_ordered <- ggplot(df, aes(x = class, y = singscore_RNAseq_sig)) +
  geom_violin(aes(fill = class), trim = FALSE, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = 16, outlier.size = 1, 
               alpha = 0.8, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_manual(
    values = pathway_colors,
    name = "Pathway"
  ) +
  labs(
    x = NULL,
    y = "sigScore Enrichment Score",
    title = "sigScore Distribution Across Metabolically Dysregulated Pathways",
    subtitle = sprintf("MSS RAS Mutant CRC (n=%d samples)", ncol(singscore_RNAseq_sig))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

pdf("S4A_20260519_COLOSSUS_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_sigScore.pdf")
panel_b_ordered
dev.off()


# plot like mean enrichment plot
pdf("20260104_COLOSSUS_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_like_mean_plot_sigScore.pdf")
ggplot(df, aes(x = class, y = singscore_RNAseq_sig)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  #  coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Gene mean expression per pathway\nDashed = low threshold, dotted = high threshold",
    x = "KEGG Pathway",
    y = "Gene mean expression"
  )
dev.off()


#----------------
# reorder mean enrichment plot
#--------------

pdf("20260104_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_mean_of_genes_across_pathways_violin_plot.pdf")
ggplot(gene_means_df, aes(x = pathway, y = gene_mean)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  geom_hline(yintercept = low_thresh, linetype = "dashed") +
  geom_hline(yintercept = high_thresh, linetype = "dotted") +
  #  coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Gene mean expression per pathway\nDashed = low threshold, dotted = high threshold",
    x = "KEGG Pathway",
    y = "Gene mean expression"
  )
dev.off()

## ---------------------------
## 6. Write out scores
## ---------------------------

# If you want a simple matrix file: pathway x samples (tab-delimited)
write.table(singscore_mat,
            file      = out_gct,
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

# Optional: quick sanity checks in R
dim(singscore_RNAseq_m)
head(rownames(singscore_RNAseq_m))  # KEGG pathways
head(colnames(singscore_RNAseq_m))  # samples
summary(as.numeric(singscore_RNAseq_m))
