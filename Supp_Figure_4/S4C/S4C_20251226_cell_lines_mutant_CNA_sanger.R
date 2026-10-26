################################################################################
# Figure Panel : Supp S4C
# Description  : Categorical copy-number heatmap (gain/loss/neutral) of selected genes in CRC
#                RAS-mutant cell lines (Sanger Cell Model Passports, PureCN).
# Input data   : WES_pureCN_CNV_genes_cn_category_20250207.csv
# Output       : CRC_CN_Sanger_categorical.csv
#                20260301_CRC_CELL_LINES_CNA_Sanger_categorical_heatmap.pdf
#                CRC_CN_Sanger_categorical_heatmap.png
#                CRC_CN_Sanger_oncoprint.pdf
#                CRC_CN_Sanger_summary.csv
#                CRC_CN_Sanger_cellline_summary.csv
# R packages   : tidyverse, ComplexHeatmap, circlize, grid, dplyr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : December 2025
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# CRC Copy Number Analysis - Sanger Cell Model Passports Data
# ===========================================================================
# Data format: Genes in rows, Cell lines in columns
# This is ALREADY CATEGORICAL (Gain/Loss/Neutral)

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Define genes of interest
genes_of_interest <- c("EGFR","PTK2","APPBP2","JPH2","TOX2",
                       "PTEN","B2M","MAF","WWOX","CDH2","SMAD4","BCL2","ZNF516")

cat("Reading Sanger Cell Model Passports data...\n")
# Read the Sanger data
sanger_cn <- read_csv("WES_pureCN_CNV_genes_cn_category_20250207.csv")

cat("\nData dimensions:", nrow(sanger_cn), "genes x", ncol(sanger_cn)-3, "cell lines\n")

# The file structure is:
# Row 1: model_name | model_id | source | [cell line names...]
# Row 2+: gene_symbol | values... 

# Check the structure
cat("\nFirst few column names:\n")
print(head(names(sanger_cn), 10))

# The first few columns are metadata (model_name, model_id, source likely contain gene info)
# Then the rest are cell line names

# Find where the cell line columns start
# Usually after columns like "symbol" or "gene"
metadata_cols <- c("model_name", "model_id", "source", "symbol", "gene", "gene_symbol")
metadata_present <- names(sanger_cn)[names(sanger_cn) %in% metadata_cols]

cat("\nMetadata columns found:", paste(metadata_present, collapse=", "), "\n")

# Assume the gene symbol is in the first column (adjust if needed)
gene_col <- names(sanger_cn)[1]
cat("Gene column:", gene_col, "\n")

# Get all cell line columns (everything except metadata)
cell_line_cols <- setdiff(names(sanger_cn), metadata_cols)
cat("\nNumber of cell lines in dataset:", length(cell_line_cols), "\n")

# Show first few cell line names
cat("First 10 cell lines:", paste(head(cell_line_cols, 10), collapse=", "), "\n")

# OPTION 1: Filter for specific CRC cell lines by name
crc_cell_lines <- c(
  "LS-513", "CL-11","LS1034","COLO-678","T84","SKCO1","NCI-H747","SW837","HCC56","SW1463","SW1116",
  "LS-123","SNU-61","RCM-1", "SW620", "SW948","CL-40","SNU-81","HCC2998","SW626", "SNU-1033", "SW403")

length(crc_cell_lines)

#cell_line_cols[1099] <- "SKCO1"
#colnames(sanger_cn)[1099] <- "SKCO1"
#cell_line_cols[130] <- "CL-11"
#colnames(sanger_cn)[130] <- "CL-11"
#cell_line_cols[1177] <- "SNU-61"
#colnames(sanger_cn)[1177] <- "SNU-61"

cell_line_cols[1100] <- "SKCO1"
colnames(sanger_cn)[1100] <- "SKCO1"
cell_line_cols[131] <- "CL-11"
colnames(sanger_cn)[131] <- "CL-11"
cell_line_cols[1178] <- "SNU-61"
colnames(sanger_cn)[1178] <- "SNU-61"


# Find which of the selected cell lines are in the dataset
available_cell_lines <- intersect(crc_cell_lines, cell_line_cols)
missing_cell_lines <- setdiff(crc_cell_lines, cell_line_cols)

cat("\nFound", length(available_cell_lines), "of", length(crc_cell_lines), "requested CRC cell lines\n")
cat("Available:", paste(head(available_cell_lines, 10), collapse=", "), "\n")
if(length(missing_cell_lines) > 0) {
  cat("Missing:", paste(head(missing_cell_lines, 5), collapse=", "), "...\n")
}

# Filter for genes of interest
available_genes <- sanger_cn %>%
  filter(.data[[gene_col]] %in% genes_of_interest) %>%
  pull(gene_col)

missing_genes <- setdiff(genes_of_interest, available_genes)

cat("\nGenes available:", length(available_genes), "of", length(genes_of_interest), "\n")
if(length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse=", "), "\n")
}

# Create the matrix: genes in rows, cell lines in columns
# Then transpose to get cell lines in rows, genes in columns
cn_matrix <- sanger_cn %>%
  filter(.data[[gene_col]] %in% genes_of_interest) %>%
  dplyr::select(all_of(gene_col), all_of(available_cell_lines)) %>%
  column_to_rownames(gene_col)

# Transpose so cell lines are rows, genes are columns
cn_categorical <- as.data.frame(t(cn_matrix))

cat("\nMatrix dimensions:", nrow(cn_categorical), "cell lines x", ncol(cn_categorical), "genes\n")
cat("Cell lines:", paste(head(rownames(cn_categorical), 5), collapse=", "), "...\n")
cat("Genes:", paste(head(colnames(cn_categorical), 5), collapse=", "), "...\n")

# Standardize category names: Gain → Amplification, Loss → Deletion
#cn_categorical <- cn_categorical %>%
 # mutate(across(everything(), ~case_when(
  #  . == "Gain" ~ "Gain",
   # . == "Loss" ~ "Loss",
  #  . == "Neutral" ~ "Neutral",
   # is.na(.) ~ "Neutral",
  #  TRUE ~ as.character(.)
  #)))

# Create summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
alteration_summary <- cn_categorical %>%
  rownames_to_column("Cell_Line") %>%
  pivot_longer(-Cell_Line, names_to = "Gene", values_to = "Status") %>%
  count(Gene, Status) %>%
  complete(Gene, Status = c("Amplification", "Deletion", "Neutral"), 
           fill = list(n = 0)) %>%
  pivot_wider(names_from = Status, values_from = n) %>%
  mutate(
    Total_Altered = Amplification + Deletion,
    Pct_Amp = round(Amplification / nrow(cn_categorical) * 100, 1),
    Pct_Del = round(Deletion / nrow(cn_categorical) * 100, 1),
    Pct_Altered = round(Total_Altered / nrow(cn_categorical) * 100, 1)
  ) %>%
  arrange(desc(Pct_Altered))

print(alteration_summary)

# Sanity check
cat("\n=== SANITY CHECK ===\n")
cat("Tumor suppressors (expect DELETIONS):\n")
for(gene in c("SMAD4", "PTEN", "B2M")) {
  if(gene %in% names(cn_categorical)) {
    dels <- sum(cn_categorical[[gene]] == "Deletion", na.rm = TRUE)
    amps <- sum(cn_categorical[[gene]] == "Amplification", na.rm = TRUE)
    cat(sprintf("  %s: %d deletions (%.1f%%), %d amplifications\n", 
                gene, dels, dels/nrow(cn_categorical)*100, amps))
  }
}

cat("\nPotential oncogenes (might show AMPLIFICATIONS):\n")
for(gene in c("EGFR", "AKT3")) {
  if(gene %in% names(cn_categorical)) {
    dels <- sum(cn_categorical[[gene]] == "Deletion", na.rm = TRUE)
    amps <- sum(cn_categorical[[gene]] == "Amplification", na.rm = TRUE)
    cat(sprintf("  %s: %d deletions, %d amplifications (%.1f%%)\n", 
                gene, dels, amps, amps/nrow(cn_categorical)*100))
  }
}

# ===================================================================
# CREATE CATEGORICAL HEATMAP
# ===================================================================

cat("\nGenerating categorical heatmap...\n")

# Define colors
cat_colors <- c("Amplification" = "#d88900", 
                 "Gain" = "#f7c97f",
                "Deletion" = "#5ec2bd", 
                "Loss" = "#a4dcd9", 
                "Neutral" = "white")

cn_categorical_ord_t <- t(cn_categorical[,genes_of_interest])

# Create heatmap
ht_categorical <- Heatmap(
  as.matrix(cn_categorical_ord_t),
  name = "CN Status",
  col = cat_colors,
  
  # Clustering
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  # Styling
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 2),
  
  # Legend
  heatmap_legend_param = list(
    title = "Copy Number",
    at = c("Amplification", "Gain", "Neutral", "Loss", "Deletion"),
    labels = c("Amplification", "Gain", "Neutral", "Loss", "Deletion"),
    legend_height = unit(4, "cm")
  ),
  
  # Titles
  column_title = "Copy Number Alterations in CRC Cell Lines (Sanger Data)",
  row_title = "Cell Lines"
)

# Display
draw(ht_categorical, heatmap_legend_side = "right")

# Save as PDF
#pdf("20260301_CRC_CELL_LINES_CNA_Sanger_categorical_heatmap.pdf", width = 12, height = 14)
draw(ht_categorical, heatmap_legend_side = "right")
#dev.off()

# Save as PNG
#png("CRC_CN_Sanger_categorical_heatmap.png", 
 #   width = 12, height = 14, units = "in", res = 300)
#draw(ht_categorical, heatmap_legend_side = "right")
#dev.off()

cat("\nHeatmap saved!\n")

# ===================================================================
# CREATE ONCOPRINT
# ===================================================================

cat("\nGenerating OncoPrint...\n")

alter_fun <- list(
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = "#d88900", col = NA))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = "#f7c97f", col = NA))
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = "#a4dcd9", col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = "#5ec2bd", col = NA))
  },
  Neutral = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, 
              gp = gpar(fill = "white", col = NA))
  }
)

#pdf("CRC_CN_Sanger_oncoprint.pdf", width = 12, height = 14)
oncoPrint(
  as.matrix(cn_categorical_ord_t),
  alter_fun = alter_fun,
  col = cat_colors,
  column_title = "Copy Number Alterations (Sanger Data)",
  heatmap_legend_param = list(
    title = "Alterations",
    at = c("Amplification", "Gain", "Loss", "Deletion", "Neutral")
  ),
  top_annotation = HeatmapAnnotation(
    cbar = anno_oncoprint_barplot()
  ),
  right_annotation = rowAnnotation(
    rbar = anno_oncoprint_barplot()
  )
)
#dev.off()

cat("OncoPrint saved!\n")

###############################################################################
# CLEAN ORDERED COPY NUMBER ONCOPRINT (SANGER DATA)
# Neutral retained
# Cell line names displayed
###############################################################################

library(ComplexHeatmap)
library(grid)

###############################################################################
# 1️⃣  Prepare matrix
###############################################################################

mat <- as.matrix(cn_categorical_ord_t)
n_samples <- ncol(mat)

# Logical matrix for burden (exclude Neutral)
altered_mat <- mat != "Neutral"

###############################################################################
# 2️⃣  Compute metrics
###############################################################################

row_count <- rowSums(altered_mat)
row_pct   <- round(row_count / n_samples * 100, 1)
col_count <- colSums(altered_mat)

###############################################################################
# 3️⃣  Order rows and columns
###############################################################################

row_order <- order(row_pct, decreasing = TRUE)
col_order <- order(col_count, decreasing = TRUE)

mat        <- mat[row_order, col_order]
row_count  <- row_count[row_order]
row_pct    <- row_pct[row_order]
col_count  <- col_count[col_order]

###############################################################################
# 4️⃣  Define colors
###############################################################################

cat_colors <- c(
  Amplification = "#d88900",
  Gain          = "#f7c97f",
  Loss          = "#a4dcd9",
  Deletion      = "#5ec2bd",
  Neutral       = "white"
)

###############################################################################
# 5️⃣  Alteration drawing
###############################################################################

alter_fun <- list(
  background = function(x, y, w, h)
    grid.rect(x, y, w, h, gp = gpar(fill = "white", col = NA)),
  
  Amplification = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cat_colors["Amplification"], col = NA)),
  
  Gain = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cat_colors["Gain"], col = NA)),
  
  Loss = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cat_colors["Loss"], col = NA)),
  
  Deletion = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cat_colors["Deletion"], col = NA)),
  
  Neutral = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "white", col = NA))
)

###############################################################################
# 6️⃣  Top annotation (sample burden)
###############################################################################

top_anno <- HeatmapAnnotation(
  Alterations = anno_barplot(
    col_count,
    border = FALSE,
    gp = gpar(fill = "grey40"),
    height = unit(3, "cm"),
    axis_param = list(side = "right")
  )
)

###############################################################################
# 7️⃣  Right annotation (gene burden + percentage)
###############################################################################

right_anno <- rowAnnotation(
  
  Count = anno_barplot(
    row_count,
    border = FALSE,
    gp = gpar(fill = "grey40"),
    width = unit(3, "cm"),
    axis_param = list(side = "top")
  ),
  
  Percent = anno_text(
    paste0(row_pct, "%"),
    just = "left",
    gp = gpar(fontsize = 9)
  ),
  
  annotation_name_side = "top"
)

###############################################################################
# 8️⃣  Generate Ordered OncoPrint
###############################################################################

oncoPrint(
  mat,
  alter_fun = alter_fun,
  alter_fun_is_vectorized = FALSE,
  col = cat_colors,
  column_title = "Copy Number Alterations (Sanger Data)",
  show_pct = FALSE,
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  top_annotation = top_anno,
  right_annotation = right_anno,
  heatmap_legend_param = list(
    title = "Copy Number"
  )
)

# ===================================================================
# EXPORT DATA FILES
# ===================================================================

cat("\nExporting data files...\n")

# 1. Categorical data (cn_categorical)
cn_categorical %>%
  rownames_to_column("Cell_Line") %>%
  write_csv("CRC_CN_Sanger_categorical.csv")

# 2. Summary statistics
#write_csv(alteration_summary, "CRC_CN_Sanger_summary.csv")

# 3. Per cell line summary
cellline_summary <- cn_categorical %>%
  rownames_to_column("Cell_Line") %>%
  pivot_longer(-Cell_Line, names_to = "Gene", values_to = "Status") %>%
  group_by(Cell_Line) %>%
  summarise(
    Total_Amps = sum(Status == "Amplification", na.rm = TRUE),
    Total_Dels = sum(Status == "Deletion", na.rm = TRUE),
    Total_Altered = Total_Amps + Total_Dels,
    Genes_Amplified = paste(Gene[Status == "Amplification"], collapse = ", "),
    Genes_Deleted = paste(Gene[Status == "Deletion"], collapse = ", ")
  ) %>%
  arrange(desc(Total_Altered))

#write_csv(cellline_summary, "CRC_CN_Sanger_cellline_summary.csv")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files saved:\n")
cat("  1. CRC_CN_Sanger_categorical_heatmap.pdf\n")
cat("  2. CRC_CN_Sanger_categorical_heatmap.png\n")
cat("  3. CRC_CN_Sanger_oncoprint.pdf\n")
cat("  4. CRC_CN_Sanger_categorical.csv (categorical copy-number matrix)\n")
cat("  5. CRC_CN_Sanger_summary.csv\n")
cat("  6. CRC_CN_Sanger_cellline_summary.csv\n")

cat("\nDone!\n")
