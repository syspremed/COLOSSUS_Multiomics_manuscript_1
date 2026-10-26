################################################################################
# Figure Panel : Supp S10E, S10I, S10P, S11I
# Description  : Tumoroid bulk RNA-seq (RAS-mt and RAS-wt, n = 131): mean MCM3/4/7 expression vs
#                replication stress (S10E), AURKB (S10I), mismatch repair (S10P) and Wong
#                proteasome genes (S11I).
# Input data   : 20260523_GSE204805_vsd_organoids_RAS_mutant_wildtype_combined_131_samples_replicate_mean_expr_data.tsv
#                2024_Rep_stress_genes_Takashi_Cancer_Res_Comm.txt
#                WONG_PROTEASOME_GENE_MODULE.v2024.1.Hs.gmt
# Output       : S9E_20260217_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_MCM_proteins_MCM_RepStress_scatter_plot.pdf
#                S9P_220260427_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_MCM_proteins_MCM_MMR_genes_scatter_plot.pdf
#                S10H_20260427_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_MCM_proteins_MCM_Wong_Proteosome_genes_scatter_plot.pdf
#                S9I_20260427_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_mcm_vs_AURKB_with_pvalues_barplot.pdf
# R packages   : ggplot2, reshape2, dplyr, ggpubr, GSEABase
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : Converted from a chunk-less .Rmd to .R. Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Load expression data
expr <- read.delim("20260523_GSE204805_vsd_organoids_RAS_mutant_wildtype_combined_131_samples_replicate_mean_expr_data.tsv", stringsAsFactors = FALSE)
cat("Expression matrix dimensions:", dim(expr), "\n")

# Keep only matched samples in expression data
expr_filtered <- expr
dim(expr_filtered)

# Load replication stress genes
repstress <- read.delim("2024_Rep_stress_genes_Takashi_Cancer_Res_Comm.txt", stringsAsFactors = FALSE)

# Define MCM genes of interest
mcm_genes <- c("MCM3", "MCM4", "MCM7")

# Filter expression data for MCM and rep-stress genes
mcm_expr <- expr_filtered[rownames(expr_filtered) %in% mcm_genes, , drop = FALSE]
repstress_genes <- unique(repstress$Gene)  # adjust if needed
repstress_expr <- expr_filtered[rownames(expr_filtered) %in% repstress_genes, , drop = FALSE]

# Calculate mean expression per sample
mean_mcm <- colMeans(mcm_expr, na.rm = TRUE)
mean_repstress <- colMeans(repstress_expr, na.rm = TRUE)

# Create data frame for plotting
df_plot <- data.frame(
  Sample = names(mean_mcm),
  Mean_MCM = mean_mcm,
  Mean_RepStress = mean_repstress
)

# Melt for ggplot
df_melt <- melt(df_plot, id.vars = "Sample", variable.name = "Gene_Set", value.name = "Mean_Expression")

library(ggpubr)
# Compute Pearson correlation
cor_test <- cor.test(df_plot$Mean_MCM, df_plot$Mean_RepStress, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 3)
label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot with manual annotation
pdf("S9E_20260217_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_MCM_proteins_MCM_RepStress_scatter_plot.pdf")
ggplot(df_plot, aes(x = Mean_MCM, y = Mean_RepStress)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, label = label_text, hjust = 1.1, vjust = 2,
           size = 5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Expression: MCM vs Replication Stress Genes",
    x = "Mean MCM Expression (MCM3–7)",
    y = "Mean Replication Stress Gene Expression"
  )
  dev.off()
  
  # MMR
  #MMR signature - KEGG MsigDB

mmr_genes <- c("EXO1",	"LIG1",	"MLH1",	"MLH3",	"MSH2",	"MSH3",	"MSH6",	"PCNA",	"PMS2",	"POLD1",	"POLD2",	"POLD3",	"POLD4",	"RFC1",	"RFC2",	"RFC3",	"RFC4",	"RFC5",	"RPA1",	"RPA2",	"RPA3",	"RPA4",	"SSBP1")

# Filter expression data
mcm_expr <- expr_filtered[rownames(expr_filtered) %in% mcm_genes, , drop = FALSE]
mmr_expr <- expr_filtered[rownames(expr_filtered) %in% mmr_genes, , drop = FALSE]

# Compute mean expression per sample
mean_mcm <- colMeans(mcm_expr, na.rm = TRUE)
mean_mmr <- colMeans(mmr_expr, na.rm = TRUE)

# Combine into one dataframe
common_samples <- intersect(names(mean_mcm), names(mean_mmr))
df_mmr <- data.frame(
  Sample = common_samples,
  Mean_MCM = mean_mcm[common_samples],
  Mean_MMR = mean_mmr[common_samples]
)

# 4. Pearson correlation
cor_result <- cor.test(df_mmr$Mean_MCM, df_mmr$Mean_MMR, method = "pearson")
label_text <- sprintf("R = %.2f\nP = %.2g", cor_result$estimate, cor_result$p.value)

# Scatter plot

# Plot with manual annotation
pdf("S9P_220260427_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_MCM_proteins_MCM_MMR_genes_scatter_plot.pdf")
ggplot(df_mmr, aes(x = Mean_MCM, y = Mean_MMR)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, label = label_text, hjust = 1.1, vjust = 2,
           size = 5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Expression: MCM vs MMR Gene Expression",
    x = "Mean MCM Gene Expression (MCM3–7)",
    y = "Mean MMR Gene Expression"
  )
  dev.off()
  
  
  # Wong proteosome
  # Load necessary libraries
library(ggplot2)
library(dplyr)
library(GSEABase)  # for read.gmt()

# 1. Load NFKB gene module (Proteasome gene set in GMT format)
gmt_data <- getGmt("WONG_PROTEASOME_GENE_MODULE.v2024.1.Hs.gmt")
nfkb_genes <- geneIds(gmt_data[[1]])  # assuming single gene set in file

# 2. Filter for matching genes in expression data
nfkb_expr <- expr_filtered[rownames(expr_filtered) %in% nfkb_genes, , drop = FALSE]

# 3. Match samples between MCM and NFKB expression matrices
common_samples <- intersect(colnames(mcm_expr), colnames(nfkb_expr))
mcm_expr_aligned <- mcm_expr[, common_samples]
nfkb_expr_aligned <- nfkb_expr[, common_samples]

# 4. Calculate mean expression
mean_mcm <- colMeans(mcm_expr_aligned, na.rm = TRUE)
mean_nfkb <- colMeans(nfkb_expr_aligned, na.rm = TRUE)

# 5. Create dataframe
df_nfkb <- data.frame(
  Sample = common_samples,
  Mean_MCM = mean_mcm,
  Mean_NFKB = mean_nfkb
)

# 6. Correlation
cor_result <- cor.test(df_nfkb$Mean_MCM, df_nfkb$Mean_NFKB, method = "pearson")
label_text <- sprintf("R = %.2f\nP = %.2g", cor_result$estimate, cor_result$p.value)

# 7. Plot
pdf("S10H_20260427_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_MCM_proteins_MCM_Wong_Proteosome_genes_scatter_plot.pdf")
ggplot(df_nfkb, aes(x = Mean_MCM, y = Mean_NFKB)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, label = label_text, hjust = 1.1, vjust = 2,
           size = 5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Expression: MCM vs Proteasome/NFKB Gene Expression",
    x = "Mean MCM Gene Expression (MCM3–7)",
    y = "Mean Proteasome Gene Expression (WONG NFKB Module)"
  )
  dev.off()
  
  
  
# AURKB vs MCMs
# Define the genes of interest
genes_of_interest <- c("AURKB", "MCM3", "MCM4", "MCM7")

# Subset the expression matrix to include only those genes
expr_selected <- expr_filtered[rownames(expr_filtered) %in% genes_of_interest, ]

# Optional: Check if any genes were not found
missing_genes <- setdiff(genes_of_interest, rownames(expr_filtered))
if (length(missing_genes) > 0) {
  message("The following genes were not found in expression data: ", paste(missing_genes, collapse = ", "))
}

# Transpose for correlation analysis (genes become columns)
expr_selected_t <- t(expr_selected)

# Ensure required genes are present
genes_of_interest <- c("AURKB", "MCM3", "MCM4", "MCM7")
expr_selected <- expr_filtered[genes_of_interest, ]
expr_selected_t <- t(expr_selected)

# Compute Pearson correlations — FIXED with simplify = FALSE
cor_results <- sapply(
  genes_of_interest[-1],  # exclude AURKB
  function(gene) cor.test(expr_selected_t[, "AURKB"], expr_selected_t[, gene], method = "pearson"),
  simplify = FALSE
)

# Build a data frame of correlations and adjusted p-values
cor_df <- data.frame(
  Gene = names(cor_results),
  Correlation = sapply(cor_results, function(x) x$estimate),
  p.value = sapply(cor_results, function(x) x$p.value)
)
cor_df$adj.p.value <- p.adjust(cor_df$p.value, method = "BH")
cor_df$p.label <- paste0("Adj p=", signif(cor_df$adj.p.value, 2))

# Plotting
library(ggplot2)

pdf("S9I_20260427_20251230_GSE204805_vsd_LMO_organoids_all_samples_RNAseq_data_mcm_vs_AURKB_with_pvalues_barplot.pdf")
library(ggplot2)

ggplot(cor_df, aes(x = reorder(Gene, Correlation), y = Correlation, fill = Correlation)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0("Adj p=", signif(adj.p.value, 2))), vjust = -0.5) +
    labs(title = paste("Correlation with AURKB"),
         y = "Pearson Correlation",
         x = "MCM Genes") +
    theme_minimal() +
    coord_flip()
  dev.off()

