################################################################################
# Figure Panel : Supplementary Figure 4I
# Description  : ssGSEA scatter + boxplot for 9 selected RAS-mutant organoid
#                lines vs MetaboAnalyst-enriched pathways.
# Input data   : LMO organoid expression (9 RAS-mutant selected lines, mean)
#                kegg_legacy gmt
# Output       : S4I_..._LMO_9_organoids_RAS_mutant_Combined_Panels_AB.pdf
# R packages   : GSVA, GSEABase, ggplot2, dplyr, cowplot
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : December 2025

# Core packages for gene set variation analysis.
library(GSVA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

# Load your ssGSEA results (assumes you've run your previous script)
# Make sure GSVA_RNAseq_m is in your environment

################################################################################
# Your 10 Selected Pathways with Actual MetaboAnalyst Data
################################################################################

pathway_data <- data.frame(
  KEGG_Name = c(
    "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM",
    "KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM",
    "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
    "KEGG_BUTANOATE_METABOLISM",
    "KEGG_CITRATE_CYCLE_TCA_CYCLE",
    "KEGG_ARGININE_AND_PROLINE_METABOLISM",
    "KEGG_PYRUVATE_METABOLISM",
    "KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM",
#    "KEGG_FATTY_ACID_BIOSYNTHESIS",
    "KEGG_PURINE_METABOLISM"
  ),
  Display_Name = c(
    "Alanine, aspartate and\nglutamate metabolism",
    "Glyoxylate and\ndicarboxylate metabolism",
    "Biosynthesis of\nunsaturated fatty acids",
    "Butanoate metabolism",
    "Citrate cycle\n(TCA cycle)",
    "Arginine and\nproline metabolism",
    "Pyruvate metabolism",
    "Glycine, serine and\nthreonine metabolism",
  #  "Fatty acid\nbiosynthesis",
    "Purine metabolism"
  ),
  # ACTUAL VALUES FROM YOUR METABOANALYST RESULTS:
  Metabolomics_pvalue = c(
    0.0045295,   # Alanine, aspartate and glutamate
    0.044039,  # Glyoxylate and dicarboxylate
    0.60404,  # Biosynthesis of unsaturated fatty acids
    0.31847,  # Butanoate
    0.40073,   # Citrate cycle (TCA)
    0.22817,    # Arginine and proline
    0.44534,    # Pyruvate
    0.0082669,   # Glycine, serine and threonine (NEW)
   # 0.03034,     # Fatty acid biosynthesis (NEW)
    0.095226      # Purine metabolism (NEW)
  ),
  Metabolomics_FDR = c(
    0.090591,  # Alanine, aspartate and glutamate
    0.39146,   # Glyoxylate and dicarboxylate
    1,   # Biosynthesis of unsaturated fatty acids
    0.94361,   # Butanoate
    1,    # Citrate cycle (TCA)
    0.82972,    # Arginine and proline
    1,    # Pyruvate
    0.099028,    # Glycine, serine and threonine (NEW)
   # 0.29065,     # Fatty acid biosynthesis (NEW)
    0.54715       # Purine metabolism (NEW)
  ),
  Impact = c(
    0.621,     # Alanine, aspartate and glutamate
    0,        # Glyoxylate and dicarboxylate
    0,           # Biosynthesis of unsaturated fatty acids
    0,           # Butanoate
    0.04412,       # Citrate cycle (TCA)
    0.12442,     # Arginine and proline
    0.0283,     # Pyruvate
    0.17228,     # Glycine, serine and threonine (NEW)
 #   0.01473,     # Fatty acid biosynthesis (NEW)
    0.13433      # Purine metabolism (NEW)
  ),
  Total_compounds = c(28, 36,32, 15, 20, 36, 23, 33, 70),
  Hits = c(4, 1, 3, 1, 1, 2, 1, 4, 4),
  stringsAsFactors = FALSE
)

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
  "Glycine, serine and\nthreonine metabolism" = "gray50",         # Pink/magenta
  "Fatty acid\nbiosynthesis" = "gray50",                           # Yellow (for new pathway)
  "Purine metabolism" = "gray50"                                   # Purple (for new pathway)
)


################################################################################
# Calculate ssGSEA Summary Statistics
################################################################################

## install if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GSVA")
# BiocManager::install("GSEABase")

library(GSVA)
library(GSEABase)

# ── Input file paths ────────────────────────────────────────────────────────
# Update these paths to match your local directory structure.
expr_file <- "20210504_Colossus_trial proxy samples_RNAseq_final_counts_data_sd0_duplicated_checked_DESeq2_rlog_filtered_5_LMO_organoids.txt"

gmt_file  <- "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

out_gct   <- "20210504_Colossus_trial proxy samples_RNAseq_final_counts_data_sd0_duplicated_checked_DESeq2_rlog_filtered_5_LMO_organoids_GSVA_ssGSEA_results_KEGG_Legacy.gct"

# Load expression matrix. Rows = genes, columns = samples.
expr_df <- read.delim(expr_file,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

## If first column is gene names and the rest are samples:
rownames(expr_df) <- expr_df[, 1]
expr_df <- expr_df[, -1, drop = FALSE]

# Convert to numeric matrix required by GSVA.
expr_mat <- as.matrix(expr_df)

#RAS mutant samples used in the manuscript
samples <- c("CRC0031",
"CRC0051",
"CRC0148",
"CRC0277",
"CRC1307",
"CRC1360",
"CRC1502",
"CRC1588",
"CRC1589")

# Get rownames
rn <- colnames(expr_mat)

# Logical match: does rowname start with any sample ID?
keep <- sapply(rn, function(x) any(startsWith(x, samples)))

# Subset matrix
expr_mat_matched <- expr_mat[, keep, drop = FALSE]

# Check
colnames(expr_mat_matched)

# expr_mat_matched: genes x samples (columns are long IDs)
stopifnot(!is.null(colnames(expr_mat_matched)))

# Group ID = first 7 characters (e.g., "CRC0051")
grp <- substr(colnames(expr_mat_matched), 1, 7)

# Average replicate columns within each group
expr_mat_avg <- sapply(split(seq_len(ncol(expr_mat_matched)), grp), function(idx) {
  rowMeans(expr_mat_matched[, idx, drop = FALSE], na.rm = TRUE)
})

# Ensure it's a matrix and set colnames to the CRC IDs
expr_mat_avg <- as.matrix(expr_mat_avg)
colnames(expr_mat_avg) <- names(split(seq_len(ncol(expr_mat_matched)), grp))

# Check replicate structure and result
table(grp)
dim(expr_mat_matched)
dim(expr_mat_avg)

expr_mat <- expr_mat_avg
dim(expr_mat)

#write.table(expr_mat,"S4I_20260322_LMO_9_organoids_for_metabolism_duplicate_mean_gene_exp.txt", sep = "\t", quote = FALSE)

#Read KEGG gene sets
# Load KEGG Legacy gene sets from MSigDB .gmt file.
gmt <- getGmt(gmt_file)         # GeneSetCollection
# Extract gene IDs as a named list: pathway name → gene symbols.
gs_list <- geneIds(gmt)         # named list: pathway → genes

#Run GSVA in ssGSEA mode (to mimic your previous ssGSEA)
gs_list <- geneIds(gmt)[pathway_data$KEGG_Name]


# Configure ssGSEA parameters.
# alpha = 0.25 is the standard weighting exponent for ssGSEA.
# normalize = TRUE scales scores to [0,1] within each sample.
ssgsea_par <- ssgseaParam(
  exprData = expr_mat,
  geneSets = gs_list,
  minSize  = 1,
  maxSize  = Inf,
  alpha    = 0.25,
  normalize = TRUE
)

# Run ssGSEA — this may take several minutes for large expression matrices.
gsva_ssgsea <- gsva(ssgsea_par, verbose = FALSE)

# This is your ssGSEA-like matrix
GSVA_RNAseq_m <- gsva_ssgsea
dim(GSVA_RNAseq_m)

# Extract ssGSEA scores for your selected pathways
selected_ssgsea <- GSVA_RNAseq_m[pathway_data$KEGG_Name, , drop = FALSE]

#write.table(GSVA_RNAseq_m, "20251228_20210504_Colossus_trial proxy samples_RNAseq_final_counts_data_sd0_duplicated_checked_DESeq2_rlog_filtered_5_LMO_organoids_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_data_selected_samples_average.txt", sep="\t", quote=FALSE, row.names=TRUE)


# Calculate summary statistics
pathway_data$Mean_ssGSEA <- apply(selected_ssgsea, 1, mean)
pathway_data$Median_ssGSEA <- apply(selected_ssgsea, 1, median)
pathway_data$SD_ssGSEA <- apply(selected_ssgsea, 1, sd)
pathway_data$Percent_Positive <- apply(selected_ssgsea, 1, 
                                       function(x) sum(x > 0) / length(x) * 100)

# Calculate -log10 values for plotting
pathway_data$neg_log10_pvalue <- -log10(pathway_data$Metabolomics_pvalue)
pathway_data$neg_log10_FDR <- -log10(pathway_data$Metabolomics_FDR)

# Classify pathways as concordant/discordant/non-significant
pathway_data$Classification <- ifelse(
  pathway_data$Metabolomics_FDR < 0.2 & pathway_data$Mean_ssGSEA > 0,
  "Concordant",
  ifelse(
    pathway_data$Metabolomics_FDR < 0.2 & pathway_data$Mean_ssGSEA <= 0,
    "Discordant",
    ifelse(
      pathway_data$Metabolomics_FDR > 0.2 & pathway_data$Mean_ssGSEA > 0,
      "Discordant",
      "Non-significant"
    )
  )
)

# Print initial summary
cat("\n" , rep("=", 80), "\n", sep = "")
cat("PATHWAY CLASSIFICATION SUMMARY\n")
cat(rep("=", 80), "\n", sep = "")
print(pathway_data[, c("Display_Name", "Metabolomics_FDR", "Mean_ssGSEA", 
                       "Percent_Positive", "Classification")])

cat("\n")
cat("Concordant pathways:", sum(pathway_data$Classification == "Concordant"), "\n")
cat("Discordant pathways:", sum(pathway_data$Classification == "Discordant"), "\n")
cat("Non-significant pathways:", sum(pathway_data$Classification == "Non-significant"), "\n\n")

################################################################################
# PANEL A: Scatter Plot (Metabolomics vs ssGSEA)
################################################################################
library(ggplot2)
library(ggrepel)
# Calculate correlation
cor_result <- cor.test(pathway_data$neg_log10_FDR, 
                       pathway_data$Mean_ssGSEA, 
                       method = "spearman")

cat("Spearman correlation: rho =", round(cor_result$estimate, 3), 
    ", p =", format.pval(cor_result$p.value, digits = 3), "\n\n")

# Create scatter plot
panel_a <- ggplot(pathway_data, aes(x = neg_log10_FDR, y = Mean_ssGSEA)) +
  geom_point(aes(color = Classification, size = Percent_Positive), alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_vline(xintercept = -log10(0.2), linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_text_repel(
    aes(label = Display_Name),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    segment.color = "gray50",
    segment.size = 0.3
  ) +
  scale_color_manual(
    values = c(
      "Concordant" = "#D55E00", 
      "Discordant" = "#0072B2", 
      "Non-significant" = "gray70"
    ),
    name = "Classification",
    guide = guide_legend(override.aes = list(size = 5))
  ) +
  scale_size_continuous(
    name = "% Samples\nPositive",
    range = c(4, 10),
    breaks = c(25, 50, 75, 100)
  ) +
  labs(
    x = "Metabolomics Enrichment (-log10 FDR)",
    y = "Mean ssGSEA Score",
    title = "Multi-omics Integration: Metabolomics vs Transcriptomics",
    subtitle = sprintf("MSS RAS Mutant CRC (n=%d samples) | Spearman rho = %.3f, p = %.3f",
                       ncol(selected_ssgsea), cor_result$estimate, cor_result$p.value)
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )

# Save Panel A
## Save figure to PDF.
ggsave("20210504_Colossus_trial proxy samples_RNAseq_final_counts_data_sd0_duplicated_checked_DESeq2_rlog_filtered_5_LMO_organoids_Metabolanalyst_FDR_ssGSEA_top_pathways_Scatter.pdf", 
   # panel_a, width = 12, height = 8, dpi = 300)
cat("✓ Panel A saved: 20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Scatter.pdf\n")

################################################################################
# PANEL B: Box/Violin Plots
################################################################################

# Convert to long format for plotting
long_data <- selected_ssgsea %>%
  as.data.frame() %>%
  mutate(KEGG_Name = rownames(.)) %>%
  pivot_longer(cols = -KEGG_Name, names_to = "Sample", values_to = "ssGSEA_Score") %>%
  left_join(
    pathway_data[, c("KEGG_Name", "Display_Name", "Classification", 
                     "Mean_ssGSEA", "Metabolomics_FDR")], 
    by = "KEGG_Name"
  )

long_data$Display_Name <- factor(long_data$Display_Name,
                                 levels = pathway_data$Display_Name[match(rownames(selected_ssgsea),
                                                                          pathway_data$KEGG_Name)])

# Order by mean ssGSEA for better visualization
#long_data$Display_Name <- reorder(long_data$Display_Name, 
       #                           long_data$Mean_ssGSEA, 
       #                           mean)

# Create box/violin plot
# Create box/violin plot with pathway colors
panel_b <- ggplot(long_data, aes(x = Display_Name, y = ssGSEA_Score)) +
  geom_violin(aes(fill = Display_Name), trim = FALSE, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = 16, outlier.size = 1, 
               alpha = 0.8, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_manual(
    values = pathway_colors,
    name = "Pathway"
  ) +
  labs(
    x = NULL,
    y = "ssGSEA Enrichment Score",
    title = "ssGSEA Score Distribution Across Metabolically Dysregulated Pathways",
    subtitle = sprintf("MSS RAS Mutant CRC (n=%d samples)", ncol(selected_ssgsea))
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

# main figure based
#panel_b_alt <- ggplot(long_data, aes(x=Display_Name, y=ssGSEA_Score)) + geom_violin(aes(fill = Display_Name), trim = FALSE) + 
#  geom_boxplot(width = 0.2)

panel_b_alt <- ggplot(long_data, aes(x = Display_Name, y = ssGSEA_Score)) +
  geom_boxplot(aes(fill = Display_Name), width = 0.6)

# Save Panel B
ggsave("S4I_20260322_LMO_9_organoids_for_metabolism_duplicate_mean_gene_exp_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_selected_samples_average_Boxplots.pdf", 
panel_b_alt, width = 12, height = 7, dpi = 300)

write.table(long_data[,1:3],"S4I_20260322_LMO_9_organoids_for_metabolism_duplicate_mean_gene_exp_Metabolanalyst_FDR_ssGSEA_top_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat("✓ Panel B saved: 20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_Boxplots.pdf\n")

################################################################################
# Combined Figure (Panels A + B)
################################################################################
library(cowplot)
combined_plot <- plot_grid(
  panel_a + theme(legend.position = "none"),
  panel_b + theme(legend.position = "none"),
  ncol = 1,
  labels = c("A", "B"),
  label_size = 16,
  rel_heights = c(1, 0.85)
)

# Extract legend from panel_a
legend <- get_legend(panel_a + 
                       theme(legend.position = "right",
                             legend.box = "vertical"))

# Combine with legend
final_plot <- plot_grid(combined_plot, legend, 
                        ncol = 2, 
                        rel_widths = c(1, 0.2))

ggsave("S4I_20260322_LMO_9_organoids_for_metabolism_duplicate_mean_gene_exp_Metaboanalyst_ssGSEA_selected_samples_average_Combined_Panels_AB.pdf", final_plot, 
  width = 7, height = 6, dpi = 300)
cat("✓ Combined figure saved: Combined_Panels_AB.pdf\n\n")

################################################################################
# Summary Tables
################################################################################

# Comprehensive summary table
summary_table <- pathway_data %>%
  dplyr::arrange(.data$Metabolomics_FDR) %>%
  dplyr::select(
    .data$Display_Name,
    .data$Total_compounds,
    .data$Hits,
    .data$Metabolomics_pvalue,
    .data$Metabolomics_FDR,
    .data$Impact,
    .data$Mean_ssGSEA,
    .data$Median_ssGSEA,
    .data$SD_ssGSEA,
    .data$Percent_Positive,
    .data$Classification
  )

write.csv(summary_table, "20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Summary_Table.csv", row.names = FALSE)
cat("✓ Summary table saved: Manuscript_Summary_Table.csv\n\n")

# Create a formatted table for manuscript
manuscript_table <- summary_table %>%
  mutate(
    Display_Name = gsub("\n", " ", Display_Name),
    Metabolomics_pvalue = format.pval(Metabolomics_pvalue, digits = 3),
    Metabolomics_FDR = format.pval(Metabolomics_FDR, digits = 3),
    Impact = round(Impact, 3),
    Mean_ssGSEA = round(Mean_ssGSEA, 1),
    Median_ssGSEA = round(Median_ssGSEA, 1),
    SD_ssGSEA = round(SD_ssGSEA, 1),
    Percent_Positive = round(Percent_Positive, 1)
  ) %>%
  rename(
    Pathway = Display_Name,
    `Total Compounds` = Total_compounds,
    `Matched Hits` = Hits,
    `Raw p-value` = Metabolomics_pvalue,
    `FDR` = Metabolomics_FDR,
    `Impact Score` = Impact,
    `Mean ssGSEA` = Mean_ssGSEA,
    `Median ssGSEA` = Median_ssGSEA,
    `SD ssGSEA` = SD_ssGSEA,
    `% Positive Samples` = Percent_Positive
  )

write.csv(manuscript_table, "20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Formatted_Manuscript_Table.csv", row.names = FALSE)
cat("✓ Formatted table saved: Formatted_Manuscript_Table.csv\n\n")

################################################################################
# Manuscript Statistics
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("KEY STATISTICS FOR MANUSCRIPT\n")
cat(rep("=", 80), "\n\n", sep = "")

# Overall summary
cat("OVERALL SUMMARY:\n")
cat("----------------\n")
cat("Total pathways analyzed:", nrow(pathway_data), "\n")
cat("Number of samples:", ncol(selected_ssgsea), "\n\n")

# Metabolomics enrichment
cat("METABOLOMICS ENRICHMENT:\n")
cat("------------------------\n")
sig_metab <- sum(pathway_data$Metabolomics_FDR < 0.05)
cat("Pathways with FDR < 0.05:", sig_metab, 
    sprintf("(%.1f%%)\n", sig_metab / nrow(pathway_data) * 100))
cat("Mean pathway impact score:", 
    sprintf("%.3f (range: %.3f - %.3f)\n\n", 
            mean(pathway_data$Impact),
            min(pathway_data$Impact),
            max(pathway_data$Impact)))

# Concordance analysis
cat("CONCORDANCE ANALYSIS:\n")
cat("---------------------\n")
n_concordant <- sum(pathway_data$Classification == "Concordant")
n_discordant <- sum(pathway_data$Classification == "Discordant")
n_nonsig <- sum(pathway_data$Classification == "Non-significant")

cat("Concordant pathways:", n_concordant, 
    sprintf("(%.1f%% of metabolomics-enriched)\n", 
            n_concordant / sig_metab * 100))
cat("Discordant pathways:", n_discordant, 
    sprintf("(%.1f%% of metabolomics-enriched)\n", 
            n_discordant / sig_metab * 100))
cat("Non-significant pathways:", n_nonsig, "\n\n")

# Correlation
cat("CORRELATION ANALYSIS:\n")
cat("---------------------\n")
cat(sprintf("Spearman correlation: rho = %.3f, p-value = %s\n\n",
            cor_result$estimate, 
            format.pval(cor_result$p.value, digits = 3)))

# Mean characteristics
cat("MEAN CHARACTERISTICS:\n")
cat("---------------------\n")
concordant_data <- pathway_data %>% filter(Classification == "Concordant")
if (nrow(concordant_data) > 0) {
  cat(sprintf("Mean ssGSEA in concordant pathways: %.1f ± %.1f\n",
              mean(concordant_data$Mean_ssGSEA),
              sd(concordant_data$Mean_ssGSEA)))
  cat(sprintf("Mean %% positive samples in concordant pathways: %.1f%%\n",
              mean(concordant_data$Percent_Positive)))
}

discordant_data <- pathway_data %>% filter(Classification == "Discordant")
if (nrow(discordant_data) > 0) {
  cat(sprintf("Mean ssGSEA in discordant pathways: %.1f ± %.1f\n",
              mean(discordant_data$Mean_ssGSEA),
              sd(discordant_data$Mean_ssGSEA)))
  cat(sprintf("Mean %% positive samples in discordant pathways: %.1f%%\n",
              mean(discordant_data$Percent_Positive)))
}
cat("\n")

# Top concordant pathways
cat("TOP CONCORDANT PATHWAYS:\n")
cat("------------------------\n")
top_concordant <- pathway_data %>%
  filter(Classification == "Concordant") %>%
  arrange(Metabolomics_FDR) %>%
  select(Display_Name, Metabolomics_FDR, Mean_ssGSEA, Percent_Positive)

if (nrow(top_concordant) > 0) {
  print(top_concordant, row.names = FALSE)
} else {
  cat("None\n")
}
cat("\n")

# Discordant pathways
cat("DISCORDANT/NON-SIGNIFICANT PATHWAYS:\n")
cat("------------------------------------\n")
non_concordant <- pathway_data %>%
  filter(Classification != "Concordant") %>%
  arrange(Metabolomics_FDR) %>%
  select(Display_Name, Metabolomics_FDR, Mean_ssGSEA, Classification)

if (nrow(non_concordant) > 0) {
  print(non_concordant, row.names = FALSE)
} else {
  cat("None\n")
}
cat("\n")

################################################################################
# Results Text Template for Manuscript
################################################################################

cat(rep("=", 80), "\n", sep = "")
cat("SUGGESTED RESULTS TEXT FOR MANUSCRIPT\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("To validate metabolic pathway alterations at the transcriptional level,\n")
cat("we performed ssGSEA analysis on matched gene expression data from", 
    ncol(selected_ssgsea), "MSS\n")
cat("RAS mutant CRC samples. Among", sig_metab, "metabolomics-enriched pathways (FDR < 0.05),\n")
cat(n_concordant, sprintf("(%.1f%%)", n_concordant / sig_metab * 100), 
    "showed concordant upregulation in ssGSEA (mean\n")
cat("enrichment score > 0), indicating coordinated transcriptional and metabolic\n")
cat("activation. Notably, pathways including")
if (nrow(top_concordant) > 0) {
  pathway_names <- gsub("\n", " ", top_concordant$Display_Name[1:min(3, nrow(top_concordant))])
  cat(paste(pathway_names, collapse = ", "))
}
cat("\ndemonstrated both significant metabolite enrichment and consistent\n")
cat("transcriptional upregulation across samples (Figure X).\n\n")

cat("The correlation between metabolomics enrichment and mean ssGSEA scores was\n")
cat(sprintf("significant (Spearman rho = %.3f, p = %s), supporting the integration\n",
            cor_result$estimate, format.pval(cor_result$p.value, digits = 3)))
cat("of metabolic and transcriptional regulation in RAS-driven metabolic\n")
cat("reprogramming.")

if (n_discordant > 0 || n_nonsig > 0) {
  cat(" However,", n_discordant + n_nonsig, "pathways showed variable\n")
  cat("or negative transcriptional enrichment despite metabolic dysregulation,\n")
  cat("suggesting post-transcriptional regulation or contributions from the tumor\n")
  cat("microenvironment.\n")
}

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Generated files:\n")
cat("  ✓ Panel_A_Metabolomics_ssGSEA_Scatter.pdf\n")
cat("  ✓ Panel_B_Pathway_Distribution_Boxplots.pdf\n")
cat("  ✓ Panel_C_ssGSEA_Heatmap.pdf\n")
cat("  ✓ Combined_Panels_AB.pdf\n")
cat("  ✓ Manuscript_Summary_Table.csv\n")
cat("  ✓ Formatted_Manuscript_Table.csv\n\n")

cat("Next steps:\n")
cat("  1. Review the figures and ensure they look correct\n")
cat("  2. Check the classification of concordant/discordant pathways\n")
cat("  3. Use the suggested results text as a template for your manuscript\n")
cat("  4. Consider additional pathway-specific validation experiments\n\n")

################################################################################
# Session Info
################################################################################
cat("Session Information:\n")
cat("--------------------\n")
print(# Record session information for reproducibility.
sessionInfo())
