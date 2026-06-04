################################################################################
# Figure Panel : Supplementary Figure 4F
# Description  : ssGSEA validation in LMX/BASALE PDX cohort (GSE204805).
#                All samples; generates boxplot/violin of pathway scores.
# Input data   : GSE204805 VST expression matrix (LMX_BASALE)
#                kegg_legacy gmt
# Output       : S4F_..._LMX_BASALE_ALL_GSVA_ssGSEA_top_pathways.pdf
# R packages   : GSVA, GSEABase, ggplot2, dplyr
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
expr_file <- "GSE204805_vsd.tsv"

gmt_file  <- "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"


out_gct   <- "GSE204805_vsd_LMX_RAS_wildtype_GSVA_ssGSEA_results_KEGG_Legacy.gct"

# Load expression matrix. Rows = genes, columns = samples.
expr_df <- read.delim(expr_file,
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

## If first column is gene names and the rest are samples:
rownames(expr_df) <- expr_df[, 1]
expr_df <- expr_df[, -1, drop = FALSE]

# Convert to numeric matrix required by GSVA.
expr_mat <- as.matrix(expr_df)

select <- read.delim("20251230_LMX_BASALE_all_IDs_meta_data.txt", stringsAsFactors = FALSE)

#RAS mutant samples used in the manuscript
samples <- substr(select$sample_id_R,1,10)

# Get rownames
rn <- colnames(expr_mat)

# Logical match: does rowname start with any sample ID?
keep <- sapply(rn, function(x) any(startsWith(x, samples)))

# Subset matrix
expr_mat_matched <- expr_mat[, keep, drop = FALSE]

# Check
colnames(expr_mat_matched)

library(dplyr)
library(readr)

# metadata file (tab-delimited). Change path as needed:
meta <- read_tsv("GSE204805_selected_metadata_annot_final_nolinfo_nooutlier_ctx.tsv", show_col_types = FALSE)

## 1) Exact-match expression columns to metadata
matched_tbl <- tibble(sample_id_R = colnames(expr_mat_matched)) %>%
  inner_join(meta, by = "sample_id_R")

## 2) Keep ONLY LMX_BASALE and LMX_BASALE.1
basale_tbl <- matched_tbl %>%
  filter(type %in% c("LMX_BASALE", "LMX_BASALE.1")) %>%
  arrange(sample_id_R)

basale_tbl

#Subset the expression matrix (this is usually what you want)
expr_mat_LMX_BASALE <- expr_mat_matched[
  ,
  basale_tbl$sample_id_R,
  drop = FALSE
]

#Sanity checks (recommended)
# How many samples retained?
ncol(expr_mat_LMX_BASALE)

# Confirm types
table(basale_tbl$type)

# Ensure no accidental mismatches
stopifnot(all(colnames(expr_mat_LMX_BASALE) == basale_tbl$sample_id_R))


# expr_mat_matched: genes x samples (columns are long IDs)
stopifnot(!is.null(colnames(expr_mat_matched)))

expr_mat_matched <- expr_mat_LMX_BASALE

# Group ID = first 7 characters (e.g., "CRC0051")
grp <- substr(colnames(expr_mat_matched), 1, 7)

# Average replicate columns within each group
expr_mat_avg <- sapply(split(seq_len(ncol(expr_mat_matched)), grp), function(idx) {
  rowMeans(expr_mat_matched[, idx, drop = FALSE], na.rm = TRUE)
})

#find the duplicate samples 
# grp is your grouping vector (same length as ncol(expr_mat_matched))
sp <- split(seq_len(ncol(expr_mat_matched)), grp)

# Groups with more than one entry
multi <- sp[lengths(sp) > 1]

# (A) show the group names
names(multi)

# (B) show counts per group (how many columns each has)
sort(lengths(multi), decreasing = TRUE)

# (C) show the actual column names for those groups
multi_colnames <- lapply(multi, function(idx) colnames(expr_mat_matched)[idx])
multi_colnames

# Ensure it's a matrix and set colnames to the CRC IDs
expr_mat_avg <- as.matrix(expr_mat_avg)
colnames(expr_mat_avg) <- names(split(seq_len(ncol(expr_mat_matched)), grp))

# Check replicate structure and result
table(grp)
dim(expr_mat_matched)
dim(expr_mat_avg)

expr_mat <- expr_mat_avg
dim(expr_mat)

# Extract current row names
rn <- rownames(expr_mat)

# Remove leading "H_"
rn_clean <- sub("^H_", "", rn)

# Assign back
rownames(expr_mat) <- rn_clean

# Check
head(rownames(expr_mat))

dim(expr_mat)

#write.table(expr_mat,"GSE204805_vsd_LMX_BASALE_all_RNAseq.txt", sep = "\t", quote = FALSE)

#ssGSEA analysis
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
# 10 Selected Pathways with Actual MetaboAnalyst Data
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
  ) )

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
  "Purine metabolism" = "gray50"                                   # Purple (for new pathway)
)

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

#write.table(GSVA_RNAseq_m, "GSE204805_vsd_LMX_BASALE_all_GSVA_ssGSEA_top_pathways_Pathway_Distribution_data_selected_samples_average.txt", sep="\t", quote=FALSE, row.names=TRUE)




################################################################################
# PANEL B: Box/Violin Plots
################################################################################

# Convert to long format for plotting
long_data <- selected_ssgsea %>%
  as.data.frame() %>%
  mutate(KEGG_Name = rownames(.)) %>%
  pivot_longer(cols = -KEGG_Name, names_to = "Sample", values_to = "ssGSEA_Score") %>%
  left_join(
    pathway_data[, c("KEGG_Name", "Display_Name")], 
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
    subtitle = sprintf("ALL LMX BASALE CRC (n=%d samples)", ncol(selected_ssgsea))
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
ggsave("S4F_GSE204805_vsd_LMX_BASALE_ALL_GSVA_ssGSEA_top_pathways_Pathway_Distribution_data_selected_samples_average_violin_plot.pdf", 
panel_b, width = 5, height = 4, dpi = 300)

# main figure based
#panel_b_alt <- ggplot(long_data, aes(x=Display_Name, y=ssGSEA_Score)) + geom_violin(aes(fill = Display_Name), trim = FALSE) + 
#  geom_boxplot(width = 0.2)

panel_b_alt <- ggplot(long_data, aes(x = Display_Name, y = ssGSEA_Score)) +
  geom_boxplot(aes(fill = Display_Name), width = 0.6)

# Save Panel B
ggsave("S4F_GSE204805_vsd_LMX_BASALE_ALL_GSVA_ssGSEA_top_pathways_Pathway_Distribution_data_selected_samples_average_boxplot.pdf", 
panel_b_alt, width = 12, height = 7, dpi = 300)
cat("✓ Panel B saved: 20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_Boxplots.pdf\n")

write.table(long_data[,1:3],"S4F_GSE204805_vsd_LMX_BASALE_ALL_GSVA_ssGSEA_top_pathways_selected_metabolites_enriched_violin_plot_GSVA_ssGSEA_ordered_scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)


################################################################################
# PANEL C: Heatmap
################################################################################

# Prepare data for heatmap
heatmap_data <- selected_ssgsea
rownames(heatmap_data) <- gsub("\n", " ", pathway_data$Display_Name)

# Row annotation
row_annotation <- data.frame(
  Classification = pathway_data$Classification,
  FDR = pathway_data$Metabolomics_FDR,
  row.names = gsub("\n", " ", pathway_data$Display_Name)
)

# Define annotation colors
ann_colors <- list(
  Classification = c(
    "Concordant" = "#D55E00",
    "Discordant" = "#0072B2",
    "Non-significant" = "gray70"
  ),
  FDR = colorRampPalette(c("darkred", "white"))(100)
)
library(pheatmap)
library(RColorBrewer)
# Create heatmap
#pdf("GSE204805_vsd_LMX_BASALE_ALL_GSVA_ssGSEA_top_pathways_Pathway_Distribution_data_selected_samples_average_Heatmap.pdf", width = 14, height = 8)
pheatmap(
  heatmap_data,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-max(abs(heatmap_data)), max(abs(heatmap_data)), length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 9,
  main = "ssGSEA Enrichment of Metabolically Dysregulated Pathways\nMSS RAS Mutant CRC",
  border_color = NA,
  cellheight = 20
)
#dev.off()
cat("✓ Panel C saved: 20251210_COLOSSUS_Metabolanalyst_FDR_ssGSEA_top_pathways_Pathway_Distribution_Boxplots_ssGSEA_Heatmap.pdf\n\n")
