################################################################################
# Figure Panel : Supplementary Figure 4E
# Description  : ssGSEA validation in the COLOSSUS clinical proxy cohort
#                (RAS-mutant, all samples, batch-corrected).
# Input data   : COLOSSUS clinical proxy RNA-seq (Combat-corrected)
#                kegg_legacy gmt
# Output       : S4E_..._clinical_Proxy_RAS_mutant_ssGSEA_boxplot.pdf
# R packages   : GSVA, GSEABase, ggplot2, dplyr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026

# Core packages for gene set variation analysis.
library(GSVA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

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
  )
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
expr_file <- "20260426_Colossus_trial proxy samples_RNAseq_final_counts_data_names_corrected_DESeq2_VST_filtered_5_combat_2_library_prep_sequencing_lane_batches.txt"

gmt_file  <- "~/Downloads/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

out_gct   <- "20250521_COLOSUSS_processed_Proxy_1_2_trial_before_AFTER_correction_GSVA_ssGSEA_data.txt"

# Load expression matrix. Rows = genes, columns = samples.
expr_df <- read.delim(expr_file,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

## If first column is gene names and the rest are samples:
#rownames(expr_df) <- expr_df[, 1]
#expr_df <- expr_df[, -1, drop = FALSE]

# Convert to numeric matrix required by GSVA.
expr_mat <- as.matrix(expr_df)

# Check
head(rownames(expr_mat))

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

#write.table(GSVA_RNAseq_m, "20250521_COLOSUSS_processed_Proxy_1_2_trial_before_AFTER_correction_GSVA_ssGSEA_top_pathways_Pathway_Distribution_data.txt", sep="\t", quote=FALSE, row.names=TRUE)




################################################################################
# PANEL B: Box/Violin Plots
################################################################################

# Convert selected_ssgsea (pathways x samples) into long format for ggplot
long_data <- selected_ssgsea %>%
  as.data.frame() %>%
  tibble::rownames_to_column("KEGG_Name") %>%
  tidyr::pivot_longer(
    cols = -KEGG_Name,
    names_to = "Sample",
    values_to = "ssGSEA_Score"
  ) %>%
  dplyr::left_join(pathway_data, by = "KEGG_Name") %>%
  dplyr::mutate(
    Display_Name = factor(Display_Name, levels = pathway_data$Display_Name),
    KEGG_Name    = factor(KEGG_Name, levels = pathway_data$KEGG_Name)
  )

# sanity check
stopifnot(nrow(long_data) == nrow(pathway_data) * ncol(selected_ssgsea))

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
    subtitle = sprintf("MSS RAS Wildtype CRC (n=%d samples)", ncol(selected_ssgsea))
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

panel_b

# Save Panel B
# Save figure to PDF.
ggsave("S4E_20250521_COLOSUSS_processed_Proxy_1_2_trial_before_AFTER_correction_ssGSEA_top_pathways_Pathway_Distribution_data_violin_plot.pdf", 
panel_b, width = 12, height = 7, dpi = 300)

# main figure based
#panel_b_alt <- ggplot(long_data, aes(x=Display_Name, y=ssGSEA_Score)) + geom_violin(aes(fill = Display_Name), trim = FALSE) + 
#  geom_boxplot(width = 0.2)

panel_b_alt <- ggplot(long_data, aes(x = Display_Name, y = ssGSEA_Score)) +
  geom_boxplot(aes(fill = Display_Name), width = 0.6)

# Save Panel B
ggsave("S4E_20250521_COLOSUSS_processed_Proxy_1_2_trial_before_AFTER_correction_ssGSEA_top_pathways_Pathway_Distribution_data_boxplot.pdf", 
panel_b_alt, width = 12, height = 7, dpi = 300)
cat("✓ Panel B saved: 20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_Boxplots.pdf\n")

write.table(long_data[,1:3],"S4E_20250521_COLOSUSS_processed_Proxy_1_2_trial_before_AFTER_correction_ssGSEA_top_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)
