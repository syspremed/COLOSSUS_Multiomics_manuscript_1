################################################################################
# Figure Panel : Supplementary Figure 4A (Claude-assisted version)
# Description  : Updated singscore analysis script (May 2026 revision) with
#                improved robustness and clarity. Functionally equivalent to
#                S4A_new_script_singscore_analysis.R.
# Input data   : COLOSSUS RNA-seq DESeq2 rlog matrix; kegg_legacy gmt
# Output       : S4A_..._singscore_v2.pdf
# R packages   : GSEABase, singscore, ggplot2, dplyr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
################################################################################

## =========================================================
##  COLOSSUS RNA-seq — KEGG pathway singscore analysis
##  Independent implementation for cross-validation
##  Reference: Foroutan et al., BMC Bioinformatics 2018
## =========================================================

library(GSEABase)
library(singscore)
library(ggplot2)
library(dplyr)
library(tidyr)

## ---- Paths -----------------------------------------------
# ── Input file paths — update to your local paths ───────────────────────────
expr_file <- "20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt"

gmt_file  <- "~/Downloads/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"

pdf_file  <- "S4A_20260519_COLOSSUS_RNAseq_ssGSEA_analysis_KEGG_Legacy_pathways_selected_metabolites_enriched_violin_plot_singscore_v2.pdf"

## ---- Step 1: Load expression matrix ----------------------
## Read with base R read.table for transparency
expr_raw <- read.table(
  expr_file,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  check.names      = FALSE,
  comment.char     = ""
)

## Store gene symbols from column 1, then drop it
#gene_symbols <- as.character(expr_raw[, 1])
#expr_raw     <- expr_raw[, -1, drop = FALSE]

## Convert all columns to numeric explicitly
#expr_raw <- as.data.frame(lapply(expr_raw, function(x) suppressWarnings(as.numeric(x))))
expr_mat <- as.matrix(expr_raw)
#rownames(expr_mat) <- gene_symbols

## Quality filters
keep <- !is.na(rownames(expr_mat)) &          # no NA gene names
  nchar(trimws(rownames(expr_mat))) > 0 & # no blank gene names
  rowSums(is.na(expr_mat)) < ncol(expr_mat) # at least one non-NA value
expr_mat <- expr_mat[keep, , drop = FALSE]

## Remove duplicated gene symbols — keep the one with highest mean expression
## (different from your script which keeps the first occurrence)
dup_genes <- duplicated(rownames(expr_mat))
if (any(dup_genes)) {
  gene_means <- rowMeans(expr_mat, na.rm = TRUE)
  # For each gene name, keep the row with the highest mean
  expr_mat <- expr_mat[order(rownames(expr_mat), -gene_means), , drop = FALSE]
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  cat("Duplicated genes removed (kept highest-mean row):", sum(dup_genes), "\n")
}

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

## ---- Step 2: Rank all genes per sample -------------------
## Full transcriptome ranking — this is the reference universe for singscore.
## Do NOT subset to pathway genes before ranking.
rank_mat <- rankGenes(expr_mat)

cat("Rank matrix dimensions:", dim(rank_mat), "\n")

## ---- Step 3: Load GMT and extract selected pathways ------
gmt     <- getGmt(gmt_file)
gs_list <- geneIds(gmt)   # named list: pathway_name -> character vector of genes

cat("Total KEGG pathways in GMT:", length(gs_list), "\n")

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

## Hard stop if any pathway is missing from GMT
not_found <- setdiff(selected_pathways, names(gs_list))
if (length(not_found) > 0) {
  stop("Pathway(s) not found in GMT:\n", paste(not_found, collapse = "\n"))
}

## ---- Step 4: Score only the 9 selected pathways ----------
## Ranking universe = full expr_mat (all genes).
## Scoring loop restricted to 9 pathways only — efficient and correct.
## centerScore = TRUE (default) centres TotalScore around 0.
## knownDirection not set — uses default (FALSE, unidirectional).

score_pathway <- function(pathway_name) {
  pathway_genes <- gs_list[[pathway_name]]
  overlap       <- intersect(rownames(rank_mat), pathway_genes)
  n_overlap     <- length(overlap)
  
  cat(sprintf("  %-55s %d / %d genes overlap\n",
              pathway_name, n_overlap, length(pathway_genes)))
  
  if (n_overlap == 0) return(rep(NA_real_, ncol(rank_mat)))
  
  ## Pass gene names as character vector — singscore handles subsetting internally
  result <- simpleScore(rankData = rank_mat, upSet = overlap)
  result$TotalScore
}

scores_list <- lapply(selected_pathways, score_pathway)
names(scores_list) <- selected_pathways

## Build sample x pathway data frame (transposed relative to your script's matrix)
## — then pivot to long format directly, avoiding an intermediate wide matrix
scores_df <- as.data.frame(scores_list)   # samples x pathways
scores_df$SampleID <- colnames(expr_mat)

cat("\nSingscore data frame:", nrow(scores_df), "samples x",
    ncol(scores_df) - 1, "pathways\n")

## ---- Step 5: Pivot to long format ------------------------
plot_df <- scores_df %>%
  pivot_longer(
    cols      = all_of(selected_pathways),
    names_to  = "pathway_id",
    values_to = "score"
  ) %>%
  filter(!is.na(score))

## ---- Step 6: Display labels ------------------------------
display_labels <- c(
  "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM" = "Alanine, aspartate and glutamate metabolism",
  "KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM"    = "Glyoxylate and dicarboxylate metabolism",
  "KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS"    = "Biosynthesis of unsaturated fatty acids",
  "KEGG_BUTANOATE_METABOLISM"                        = "Butanoate metabolism",
  "KEGG_CITRATE_CYCLE_TCA_CYCLE"                    = "Citrate cycle (TCA cycle)",
  "KEGG_ARGININE_AND_PROLINE_METABOLISM"             = "Arginine and proline metabolism",
  "KEGG_PYRUVATE_METABOLISM"                         = "Pyruvate metabolism",
  "KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM"    = "Glycine, serine and threonine metabolism",
  "KEGG_PURINE_METABOLISM"                           = "Purine metabolism"
)

plot_df$pathway <- display_labels[plot_df$pathway_id]

## Order by ascending median score (lowest enrichment on left)
## — deliberate difference from your script's descending order,
##   so you can compare which direction better suits the figure
pathway_order <- plot_df %>%
  group_by(pathway) %>%
  summarise(med = median(score, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(pathway)

plot_df$pathway <- factor(plot_df$pathway, levels = pathway_order)

## ---- Step 7: Inspect medians and score range -------------
summary_tbl <- plot_df %>%
  group_by(pathway) %>%
  summarise(
    n        = n(),
    median   = round(median(score, na.rm = TRUE), 4),
    mean     = round(mean(score,   na.rm = TRUE), 4),
    sd       = round(sd(score,     na.rm = TRUE), 4),
    min      = round(min(score,    na.rm = TRUE), 4),
    max      = round(max(score,    na.rm = TRUE), 4),
    .groups  = "drop"
  ) %>%
  arrange(desc(median))

print(summary_tbl, n = Inf)

## ---- Step 8: Pathway colours (same palette as reference) --
pathway_colors <- c(
  "Alanine, aspartate and glutamate metabolism" = "#E76F51",
  "Glyoxylate and dicarboxylate metabolism"     = "#F4A261",
  "Biosynthesis of unsaturated fatty acids"     = "#A9A965",
  "Butanoate metabolism"                         = "#2A9D8F",
  "Citrate cycle (TCA cycle)"                   = "#4ECDC4",
  "Arginine and proline metabolism"              = "#1982C4",
  "Pyruvate metabolism"                          = "#6A9FB5",
  "Glycine, serine and threonine metabolism"    = "gray50",
  "Purine metabolism"                            = "gray30"
)

## ---- Step 9: Violin plot ---------------------------------
p <- ggplot(plot_df, aes(x = pathway, y = score)) +
  geom_violin(
    aes(fill = pathway),
    trim  = FALSE,
    alpha = 0.6,
    scale = "width"
  ) +
  geom_boxplot(
    width         = 0.2,
    outlier.shape = 16,
    outlier.size  = 1,
    alpha         = 0.8,
    fill          = "white"
  ) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "gray40",
    linewidth  = 0.8
  ) +
  scale_fill_manual(values = pathway_colors) +
  labs(
    x        = NULL,
    y        = "singscore enrichment score",
    title    = "singscore distribution across selected KEGG metabolic pathways",
    subtitle = sprintf("MSS RAS-mutant CRC (n = %d samples)", ncol(expr_mat))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

## ---- Step 10: Save ---------------------------------------
pdf(pdf_file, width = 10, height = 6)
print(p)
dev.off()

cat("Plot saved to:\n", pdf_file, "\n")
