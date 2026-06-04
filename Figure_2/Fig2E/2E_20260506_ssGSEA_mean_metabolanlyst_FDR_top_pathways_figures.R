################################################################################
# Figure Panel : 2E
# Description  : Multi-omics integration – metabolomics (MetaboAnalyst) vs
#                transcriptomics (ssGSEA) for COLOSSUS MSS RAS-mutant CRC.
#                Generates: (A) scatter plot of concordance, (B) boxplot of
#                ssGSEA scores across selected KEGG pathways.
# Input data   : (1) COLOSSUS RNA-seq VST matrix (Combat-corrected)
#                (2) c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt
#                (3) MetaboAnalyst pathway_results.csv (hardcoded values)
# Output       : 2E_..._ssGSEA_top_pathways_Pathway_Distribution_Boxplots.pdf
#                2E_..._Combined_Panels_AB.pdf
# R packages   : GSVA, GSEABase, ggplot2, dplyr, tidyr, pheatmap,
#                RColorBrewer, ggrepel, cowplot
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : Update expr_file, gmt_file, and out_gct paths before running.
################################################################################

################################################################################
# updated from Multi-omics Integration: MSS RAS Mutant CRC
# Metabolomics (MetaboAnalyst) vs Transcriptomics (ssGSEA)
# Author:Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date: December 2024
################################################################################

# Core packages: GSVA for gene set scoring, tidyverse tools for data wrangling,
# ggrepel for non-overlapping labels, cowplot for multi-panel assembly.
library(GSVA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

# NOTE: GSVA_RNAseq_m will be computed from scratch below (see 'Run GSVA' section).
# This comment is kept for historical context only.

################################################################################
# Selected KEGG pathways — MetaboAnalyst enrichment results
# These pathways were identified as significantly enriched (FDR < 0.2) in the
# COLOSSUS RAS-mutant patient metabolomics data using MetaboAnalyst pathway
# analysis. Values (p-values, FDR, impact scores) are hardcoded from the
# MetaboAnalyst output CSV to ensure reproducibility without re-running the
# web tool.
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
    3.5198e-06,   # Alanine, aspartate and glutamate
    0.00011063,  # Glyoxylate and dicarboxylate
    0.00022024,  # Biosynthesis of unsaturated fatty acids
    0.00042163,  # Butanoate
    0.0013627,   # Citrate cycle (TCA)
    0.001847,    # Arginine and proline
    0.0023545,    # Pyruvate
    0.0090275,   # Glycine, serine and threonine (NEW)
   # 0.03034,     # Fatty acid biosynthesis (NEW)
   0.10216      # Purine metabolism (NEW)
  ),
  Metabolomics_FDR = c(
    0.00028159,  # Alanine, aspartate and glutamate
    0.0029501,   # Glyoxylate and dicarboxylate
    0.0044048,   # Biosynthesis of unsaturated fatty acids
    0.0067461,   # Butanoate
    0.018169,    # Citrate cycle (TCA)
    0.021109,    # Arginine and proline
    0.023545,    # Pyruvate
    0.072751,    # Glycine, serine and threonine (NEW)
   # 0.29065,     # Fatty acid biosynthesis (NEW)
   0.5108       # Purine metabolism (NEW)
  ),
  Impact = c(
    0.39984,     # Alanine, aspartate and glutamate
    0.35,        # Glyoxylate and dicarboxylate
    0,           # Biosynthesis of unsaturated fatty acids
    0,           # Butanoate
    0.153,       # Citrate cycle (TCA)
    0.14884,     # Arginine and proline
    0.21967,     # Pyruvate
    0.25981,     # Glycine, serine and threonine (NEW)
 #   0.01473,     # Fatty acid biosynthesis (NEW)
    0.06569      # Purine metabolism (NEW)
  ),
  Total_compounds = c(28, 32, 36, 15, 20, 36, 23, 33, 70),
  Hits = c(7, 6, 6, 4, 4, 5, 4, 4, 4),
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
# Update these paths to point to your local copies of the data files.
expr_file <- "20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt"

gmt_file  <- "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

out_gct   <- "20201014_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_sd0_DESeq2_rlog_filtered_5_ALL_Combat_6_batches_corrected_Validation_sd0_GSVA_ssGSEA_results_KEGG_Legacy.gct"

# Load the Combat-corrected RNA-seq VST expression matrix.
# Rows = genes (gene symbols), columns = patient samples.
expr_df <- read.delim(expr_file,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

## If first column is gene names and the rest are samples:
#rownames(expr_df) <- expr_df[, 1]
#expr_df <- expr_df[, -1, drop = FALSE]

# Convert to numeric matrix — required format for GSVA.
expr_mat <- as.matrix(expr_df)

# Read KEGG Legacy gene sets from MSigDB .gmt file.
gmt <- getGmt(gmt_file)         # GeneSetCollection
gs_list <- geneIds(gmt)  # named list: pathway name → character vector of gene symbols

# Run GSVA in ssGSEA mode.
# ssGSEA (single-sample GSEA) scores each sample independently against each
# gene set, producing a per-sample enrichment score.
# Subset gene sets to only the 9 MetaboAnalyst-enriched pathways.
gs_list <- geneIds(gmt)[pathway_data$KEGG_Name]


# ssGSEA parameters:
#   alpha = 0.25  — standard rank weighting exponent
#   normalize = TRUE — scale scores within each sample to improve comparability
ssgsea_par <- ssgseaParam(
  exprData = expr_mat,
  geneSets = gs_list,
  minSize  = 1,
  maxSize  = Inf,
  alpha    = 0.25,
  normalize = TRUE
)

# Run ssGSEA. This may take a few minutes for large cohorts.
gsva_ssgsea <- gsva(ssgsea_par, verbose = FALSE)

# This is your ssGSEA-like matrix
GSVA_RNAseq_m <- gsva_ssgsea
dim(GSVA_RNAseq_m)

# Extract ssGSEA scores for the 9 selected pathways (rows) across all samples (columns).
selected_ssgsea <- GSVA_RNAseq_m[pathway_data$KEGG_Name, , drop = FALSE]

#write.table(GSVA_RNAseq_m, "20251216_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_data.txt", sep="\t", quote=FALSE, row.names=TRUE)


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
# PANEL A: Scatter Plot — Metabolomics enrichment vs mean ssGSEA score
# Each point = one KEGG pathway.
# x-axis = metabolomics enrichment (-log10 FDR from MetaboAnalyst)
# y-axis = mean ssGSEA score across all COLOSSUS patients
# Point size = % of samples with positive ssGSEA score
# Concordant = metabolomically enriched AND transcriptionally upregulated
################################################################################

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
#ggsave("20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Scatter.pdf", 
    #  panel_a, width = 12, height = 8, dpi = 300)
cat("✓ Panel A saved: 20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Scatter.pdf\n")

################################################################################
# PANEL B: Distribution of ssGSEA scores per pathway
# Violin + boxplot showing the spread of per-sample ssGSEA scores for each
# of the 9 MetaboAnalyst-enriched pathways across all COLOSSUS patients.
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
ggsave("2E_20260506_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_Boxplots.pdf", 
  panel_b_alt, width = 12, height = 7, dpi = 300)
cat("✓ Panel B saved: 2E_20260506_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_Boxplots.pdf\n")


write.table(long_data[,1:3],"2E_20260506_COLOSSUS_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)


################################################################################
# Assemble final figure: Panel A (scatter) above Panel B (boxplots)
# cowplot::plot_grid() arranges panels and extracts a shared legend.
################################################################################

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

ggsave("2E_20260506_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_scatter_Boxplots_ssGSEA_Combined_Panels_AB.pdf", final_plot, 
       width = 14, height = 13, dpi = 300)
cat("✓ Combined figure saved: Combined_Panels_AB.pdf\n\n")

################################################################################
# Summary tables for manuscript supplementary materials
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
# Print key statistics to console — copy into manuscript Results section.
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
# Record R session details for reproducibility reporting.
################################################################################
cat("Session Information:\n")
cat("--------------------\n")
print(sessionInfo())
