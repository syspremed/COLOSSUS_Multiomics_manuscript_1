################################################################################
# Figure Panel : 2F
# Description  : Multi-omics integration – metabolomics vs proteomics (ssGSEA).
#                Mirror of Fig2E but using COLOSSUS RPPA/proteomics data instead
#                of RNA-seq. Generates scatter + boxplot panels.
# Input data   : COLOSSUS proteomics expression matrix; kegg_legacy gmt;
#                MetaboAnalyst pathway_results (hardcoded values)
# Output       : 2F_..._ssGSEA_PROTEIN_top_pathways_Scatter_violin_plot_
#                Combined_Panels_AB.pdf
# R packages   : GSVA, GSEABase, ggplot2, dplyr, tidyr, pheatmap,
#                RColorBrewer, ggrepel, cowplot
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : December 2025
# Note         : Update expr_file to point to the proteomics matrix.

# Same pipeline as Fig 2E but applied to COLOSSUS RPPA proteomics data
# instead of RNA-seq. This tests whether metabolic pathway enrichment is
# also reflected at the protein level.
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
    2.543e-06,   # Alanine, aspartate and glutamate
    0.00011063,  # Glyoxylate and dicarboxylate
    0.00022024,  # Biosynthesis of unsaturated fatty acids
    0.00042163,  # Butanoate
    0.0013627,   # Citrate cycle (TCA)
    0.001847,    # Arginine and proline
    0.002345,    # Pyruvate
    0.0090275,   # Glycine, serine and threonine (NEW)
   # 0.03034,     # Fatty acid biosynthesis (NEW)
    0.10216      # Purine metabolism (NEW)
  ),
  Metabolomics_FDR = c(
    0.00020159,  # Alanine, aspartate and glutamate
    0.0029501,   # Glyoxylate and dicarboxylate
    0.0044048,   # Biosynthesis of unsaturated fatty acids
    0.0067461,   # Butanoate
    0.018169,    # Citrate cycle (TCA)
    0.021109,    # Arginine and proline
    0.023545,    # Pyruvate
    0.072751,    # Glycine, serine and threonine (NEW)
   # 0.29065,     # Fatty acid biosynthesis (NEW)
    0.5103       # Purine metabolism (NEW)
  ),
  Impact = c(
    0.39994,     # Alanine, aspartate and glutamate
    0.35,        # Glyoxylate and dicarboxylate
    0,           # Biosynthesis of unsaturated fatty acids
    0,           # Butanoate
    0.153,       # Citrate cycle (TCA)
    0.14884,     # Arginine and proline
    0.21667,     # Pyruvate
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
  "Fatty acid biosynthesis" = "gray50",                           # Yellow (for new pathway)
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


## ── Input file paths ─────────────────────────────────────────────────────────
# Update to point to your local proteomics matrix and GMT file.
expr_file <- "20200928_metabolomics_median_replicates_data_80_Non_zero_sample_median_centred_imputed_MARS_combat_lysis_batch_protein_names_only.txt"

gmt_file  <- "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

out_gct   <- "20201014_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_sd0_DESeq2_rlog_filtered_5_ALL_Combat_6_batches_corrected_Validation_sd0_GSVA_ssGSEA_results_KEGG_Legacy.gct"

#screenExpr(expr_file,0)


#name of the protein file was incorrectly written.

# here is a test for it 

d1 <- read.delim("20200928_metabolomics_median_replicates_data_80_Non_zero_sample_median_centred_imputed_MARS_combat_lysis_batch_protein_names_only.txt")
dim(d1)

d2 <- read.delim("20200607_colossus_retro_protein_data_log2_median_normalized_imputed_cohort_lysis_batch_corrected.txt")
dim(d2)

cor(c(as.numeric(unlist(d1[,-1]))),c(as.numeric(unlist(d2[,-1]))))


# check the protein file name, which is correct, but incorrectly labelled as metabolomics
expr_file <- "20200928_metabolomics_median_replicates_data_80_Non_zero_sample_median_centred_imputed_MARS_combat_lysis_batch_protein_names_only_sd0.txt"


# Load RPPA protein expression matrix.
# Rows = protein antibody probes, columns = patient samples.
expr_df <- read.delim(expr_file,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

## If first column is gene names and the rest are samples:
rownames(expr_df) <- expr_df[,1]
expr_df <- expr_df[, -1, drop = FALSE]

expr_mat <- as.matrix(expr_df)

#Read KEGG gene sets
gmt <- getGmt(gmt_file)         # GeneSetCollection
gs_list <- geneIds(gmt)         # named list: pathway → genes

#Run GSVA in ssGSEA mode (to mimic your previous ssGSEA)
# Subset gene sets to MetaboAnalyst-enriched pathways only.
# Note: RPPA covers a limited antibody panel; not all pathway genes will be
# present. ssGSEA handles this gracefully by using available gene overlap.
gs_list <- geneIds(gmt)[pathway_data$KEGG_Name] 

# ssGSEA parameters (identical to RNA-seq script for comparability).
ssgsea_par <- ssgseaParam(
  exprData = expr_mat,
  geneSets = gs_list,
  minSize  = 1,
  maxSize  = Inf,
  alpha    = 0.25,
  normalize = TRUE
)

# Run ssGSEA on proteomics data.
gsva_ssgsea <- gsva(ssgsea_par, verbose = FALSE)

# This is your ssGSEA-like matrix
GSVA_RNAseq_m <- gsva_ssgsea
dim(GSVA_RNAseq_m)

# Extract ssGSEA scores for your selected pathways
selected_ssgsea <- GSVA_RNAseq_m[pathway_data$KEGG_Name, , drop = FALSE]

#write.table(GSVA_RNAseq_m, "20251216_COLOSSUS_Metabolanalyst_FDR_ssGSEA_Protein_top_pathways_Pathway_Distribution_data.txt", sep="\t", quote=FALSE, row.names=TRUE)


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

# Calculate correlation
cor_result <- cor.test(pathway_data$neg_log10_FDR, 
                       pathway_data$Mean_ssGSEA, 
                       method = "spearman")

cat("Spearman correlation: rho =", round(cor_result$estimate, 3), 
    ", p =", format.pval(cor_result$p.value, digits = 3), "\n\n")

size_vals <- sort(unique(pathway_data$Percent_Positive))

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
    name   = "% Samples\nPositive",
    range  = c(4, 10),
    breaks = size_vals,                  # exact data values
    labels = round(size_vals)            # pretty legend labels (e.g. 90, 99, 100)
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
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )

# Save Panel A
#ggsave("20251216_COLOSSUS_Metabolanalyst_FDR_ssGSEA_PROTEIN_top_pathways_Scatter.pdf", 
 #  panel_a, width = 12, height = 8, dpi = 300)
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
ggsave("2F_20260322_COLOSSUS_Metabolanalyst_FDR_ssGSEA_PROTEIN_top_pathways_Pathway_Distribution_Boxplots.pdf", 
   panel_b_alt, width = 12, height = 7, dpi = 300)
cat("✓ Panel B saved: 20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_Boxplots.pdf\n")

write.table(long_data[,1:3],"2F_20260322_COLOSSUS_Protein_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# Combined Figure (Panels A + B)
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

ggsave("2F_20260322_COLOSSUS_Metabolanalyst_FDR_ssGSEA_PROTEIN_top_pathways_Scatter_violin_plot_Combined_Panels_AB.pdf", final_plot, 
       width = 14, height = 13, dpi = 300)
cat("✓ Combined figure saved: Combined_Panels_AB.pdf\n\n")

################################################################################
# Summary Tables
################################################################################

# Comprehensive summary table
summary_table <- pathway_data %>%
  arrange(Metabolomics_FDR) %>%
  select(
    Display_Name,
    Total_compounds,
    Hits,
    Metabolomics_pvalue,
    Metabolomics_FDR,
    Impact,
    Mean_ssGSEA,
    Median_ssGSEA,
    SD_ssGSEA,
    Percent_Positive,
    Classification
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
# Session Info
################################################################################
cat("Session Information:\n")
cat("--------------------\n")
print(sessionInfo())
