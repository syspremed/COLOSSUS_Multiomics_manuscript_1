################################################################################
# Figure Panel : Supplementary Figure 4A
# Description  : singscore-based single-sample gene set scoring for COLOSSUS
#                RNA-seq. Scores each sample against selected KEGG metabolic
#                pathways and generates violin plots of signature scores.
# Input data   : COLOSSUS RNA-seq DESeq2 rlog-normalised matrix
#                c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt
# Output       : S4A_..._violin_plot_singscore.pdf
# R packages   : GSEABase, singscore, ggplot2, dplyr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : singscore is the primary method for S4A (see backup script
#                for GSVA/ssGSEA alternative).
################################################################################

## =========================================================
## COLOSSUS RNA-seq KEGG pathway  analysis
## Rewritten for robustness and clarity
## =========================================================

## ---- Install if needed ----
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("GSEABase"))
# install.packages(c("singscore", "ggplot2", "dplyr"))

# ── Required libraries ────────────────────────────────────────────────────────
# GSEABase: read GMT gene set files
# singscore: single-sample gene set scoring (alternative to ssGSEA)
# ggplot2/dplyr: visualisation and data manipulation
## ---- Load libraries ----
library(GSEABase)
library(singscore)
library(ggplot2)
library(dplyr)

# ── Input / output file paths ────────────────────────────────────────────────
# Update these paths to match your local directory structure.
## ---- Input / output files ----
expr_file <- "20200830_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_sd0_DESeq2_rlog_filtered_5_ALL_Combat_6_batches_corrected_Validation_sd0.txt"

gmt_file  <- "~/Downloads/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

pdf_file  <- "S4A_20260404_COLOSSUS_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_singscore.pdf"

# ── Step 1: Load RNA-seq expression matrix ───────────────────────────────────
# DESeq2 rlog-normalised, Combat batch-corrected.
# Rows = genes (symbols), columns = patient samples.
## ---- Read expression data ----
expr_df <- read.delim(
  expr_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## Assume first column contains gene symbols
rownames(expr_df) <- expr_df[[1]]
expr_df <- expr_df[, -1, drop = FALSE]

## Convert to numeric matrix safely
expr_df[] <- lapply(expr_df, as.numeric)
expr_mat <- as.matrix(expr_df)
mode(expr_mat) <- "numeric"

## Remove rows with missing/blank gene names
expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & rownames(expr_mat) != "", , drop = FALSE]

## If duplicated gene symbols exist, keep first occurrence
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]

cat("Expression matrix dimensions:", dim(expr_mat), "\n")
cat("Number of samples:", ncol(expr_mat), "\n")
cat("Number of genes:", nrow(expr_mat), "\n")

## ---- Read KEGG gene sets ----
gmt <- getGmt(gmt_file)
gs_list <- geneIds(gmt)

cat("Number of KEGG pathways in GMT:", length(gs_list), "\n")

## ---- Rank genes for singscore ----
## rankGenes expects genes in rows, samples in columns
rank_mat <- rankGenes(expr_mat)

## ---- Compute singscore for each KEGG pathway ----
singscore_list <- lapply(gs_list, function(genes) {
  common_genes <- intersect(rownames(rank_mat), genes)
  
  if (length(common_genes) == 0) {
    return(rep(NA_real_, ncol(rank_mat)))
  }
  
  sc <- simpleScore(
    rankData    = rank_mat,
    upSet       = common_genes,
    centerScore = TRUE
  )
  
  sc$TotalScore
})

singscore_mat <- do.call(rbind, singscore_list)
rownames(singscore_mat) <- names(gs_list)
colnames(singscore_mat) <- colnames(expr_mat)

cat("singscore matrix dimensions:", dim(singscore_mat), "\n")

## ---- Define pathways of interest by exact GMT names ----
## Check your GMT file to ensure these names match exactly.
selected_pathways <- c(
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

missing_pathways <- setdiff(selected_pathways, rownames(singscore_mat))
if (length(missing_pathways) > 0) {
  stop(
    "These selected pathways were not found in the singscore matrix:\n",
    paste(missing_pathways, collapse = "\n")
  )
}

singscore_selected <- singscore_mat[selected_pathways, , drop = FALSE]

## ---- Display labels for plotting ----
display_labels <- c(
  "Alanine, aspartate and glutamate metabolism",
  "Glyoxylate and dicarboxylate metabolism",
  "Biosynthesis of unsaturated fatty acids",
  "Butanoate metabolism",
  "Citrate cycle (TCA cycle)",
  "Arginine and proline metabolism",
  "Pyruvate metabolism",
  "Glycine, serine and threonine metabolism",
  "Purine metabolism"
)

names(display_labels) <- selected_pathways

## ---- Define pathway colors with labels matching plot labels exactly ----
pathway_colors <- c(
  "Alanine, aspartate and glutamate metabolism"      = "#E76F51",
  "Glyoxylate and dicarboxylate metabolism"          = "#F4A261",
  "Biosynthesis of unsaturated fatty acids"          = "#A9A965",
  "Butanoate metabolism"                             = "#2A9D8F",
  "Citrate cycle (TCA cycle)"                        = "#4ECDC4",
  "Arginine and proline metabolism"                  = "#1982C4",
  "Pyruvate metabolism"                              = "#6A9FB5",
  "Glycine, serine and threonine metabolism"         = "gray50",
  "Purine metabolism"                                = "gray30"
)

## ---- Convert selected score matrix to long-format data frame ----
df <- data.frame(
  pathway_id = rep(rownames(singscore_selected), each = ncol(singscore_selected)),
  sample     = rep(colnames(singscore_selected), times = nrow(singscore_selected)),
  score      = as.vector(t(singscore_selected)),
  stringsAsFactors = FALSE
)

df$pathway <- display_labels[df$pathway_id]

## ---- Order pathways by median score ----
pathway_order <- df %>%
  group_by(pathway) %>%
  summarise(median_score = median(score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_score)) %>%
  pull(pathway)

df$pathway <- factor(df$pathway, levels = pathway_order)

## ---- Optional: inspect medians ----
median_table <- df %>%
  group_by(pathway) %>%
  summarise(
    median_score = median(score, na.rm = TRUE),
    mean_score   = mean(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_score))

print(median_table)

## ---- Plot ----
panel_b_ordered <- ggplot(df, aes(x = pathway, y = score)) +
  geom_violin(aes(fill = pathway), trim = FALSE, alpha = 0.6, scale = "width") +
  geom_boxplot(
    width = 0.2,
    outlier.shape = 16,
    outlier.size = 1,
    alpha = 0.8,
    fill = "white"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = pathway_colors) +
  labs(
    x = NULL,
    y = "singscore enrichment score",
    title = "singscore distribution across selected KEGG metabolic pathways",
    subtitle = sprintf("MSS RAS-mutant CRC (n = %d samples)", ncol(expr_mat))
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

## ---- Save PDF ----
pdf(pdf_file, width = 10, height = 6)
print(panel_b_ordered)
dev.off()

cat("Plot saved to:\n", pdf_file, "\n")
