################################################################################
# Figure Panels : Supplementary Figures 5M & 5N
# Description  : Validation of glutamine–glutamate axis in the public matched
#                metabolomic–transcriptomic cohort GSE89076. Generates:
#                (5M) scatter plot of glutamine vs glutamate metabolites;
#                (5N) scatter/violin plots of GLS, GOT2, GLUL gene expression
#                     vs glutamate/glutamine ratio.
# Input data   : GSE89076.Agilent8x60K.log2_transformed.gene_symbol.csv
#                PreprocessedData_COAD_tumours_only.txt (metabolomics)
#                Sample annotation file
# Output       : S5M_..._scatter_plot.pdf
#                S5N_..._violin_plot.pdf / linear_regression_plot.pdf
# R packages   : ggplot2, dplyr, ggpubr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : April 2026
################################################################################



#Publicly available COAD metabolic and expression data preparation

## ===================================================
## 0. File names
## ===================================================

sample_file <- "Supplementary data/GSE89076/20251210_GSE89076_sample_GSE_ID.txt"
expr_file   <- "Supplementary data/GSE89076/GSE89076.Agilent8x60K.log2_transformed.gene_symbol.csv"
metab_file  <- "Supplementary data/PreprocessedData_COAD_tumours_only.txt"

expr_tumor_out  <- "Supplementary data/GSE89076/GSE89076.Agilent8x60K.log2_transformed.gene_symbol_Tumour_only.csv"
expr_match_out  <- "Supplementary data/GSE89076/GSE89076.Agilent8x60K.Matched_expression_tumours_matched.csv"
metab_match_out <- "Supplementary data/GSE89076/GSE89076.Agilent8x60K.Matched_metabolites_tumours_matched.txt"


## ===================================================
## 1. Read sample annotation (GSE ID / Sample ID)
##    Columns: "GSE ID"    "Sample ID"
## ===================================================

sample_info <- read.table(
  sample_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Define column names explicitly
gse_col    <- "GSE ID"
sample_col <- "Sample ID"

# Check they exist
stopifnot(gse_col %in% colnames(sample_info),
          sample_col %in% colnames(sample_info))

# Tumour samples = Sample ID ending with "T"
tumour_info <- sample_info[grep("T$", sample_info[[sample_col]]), ]

cat("Tumour samples in sample_info:", nrow(tumour_info), "\n")


## ===================================================
## 2. Read expression matrix and keep tumour GSMs
##    Format: gene_symbol, GSM2358437, GSM2358438, ...
## ===================================================

expr <- read.csv(
  expr_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_col <- "gene_symbol"
stopifnot(gene_col %in% colnames(expr))

# GSM IDs for tumours from tumour_info
tumour_gsms <- tumour_info[[gse_col]]

# Keep only those GSMs that are columns in expr
tumour_gsms_in_expr <- intersect(tumour_gsms, colnames(expr))

cat("Tumour GSMs in expression matrix:", length(tumour_gsms_in_expr), "\n")

if (length(tumour_gsms_in_expr) == 0) {
  stop("No tumour GSM IDs from sample_info are present in expression matrix columns.")
}

# Subset expression to tumour-only data
expr_tumor <- expr[, c(gene_col, tumour_gsms_in_expr), drop = FALSE]

# Optional: save tumour-only expression matrix
#write.csv(expr_tumor, expr_tumor_out, row.names = FALSE)
cat("Tumour-only expression written to:", expr_tumor_out, "\n")


## ===================================================
## 3. Read metabolomics tumour data
##    Format: metabolites   6T  1T  2T  ...
## ===================================================

metab <- read.table(
  metab_file,
  header       = TRUE,
  sep          = "\t",
  stringsAsFactors = FALSE,
  check.names  = FALSE,
  quote        = "",      # <– ignore all quote characters
  comment.char = "",      # <– don't treat # or others as comments
  fill         = TRUE     # <– pad rows with fewer fields using NA
)

metab_feature_col <- colnames(metab)[1]   # usually "metabolites"
metab_samples     <- colnames(metab)[-1]  # e.g. 6T, 1T, 2T, ...

cat("Metabolite tumour samples:", length(metab_samples), "\n")


## ===================================================
## 4. Map expression GSM IDs to Sample IDs (e.g. 6T, 1T)
## ===================================================

# GSM columns actually present in tumour expression matrix
expr_gsms <- colnames(expr_tumor)[-1]

# Map GSM -> Sample ID using sample_info (not only tumour_info,
# in case you want to reuse the mapping later)
expr_sample_ids <- sample_info[[sample_col]][match(expr_gsms, sample_info[[gse_col]])]

# Name the vector by GSM ID
names(expr_sample_ids) <- expr_gsms

# Check mapping
cat("Mapped GSM -> Sample ID for expression columns:\n")
print(head(data.frame(GSM = expr_gsms, SampleID = expr_sample_ids), 5))


## ===================================================
## 5. Find common tumour Sample IDs between expression & metabolomics
## ===================================================

common_sample_ids <- intersect(expr_sample_ids, metab_samples)

cat("Common tumour Sample IDs in BOTH expression & metabolomics:",
    length(common_sample_ids), "\n")

if (length(common_sample_ids) == 0) {
  stop("No overlapping tumour Sample IDs between expression and metabolomics data.")
}

# GSMs corresponding to these common Sample IDs
matched_gsms <- names(expr_sample_ids)[expr_sample_ids %in% common_sample_ids]


## ===================================================
## 6. Subset and align matrices
## ===================================================

## 6a. Expression: gene_symbol + matched GSMs
expr_matched <- expr_tumor[, c(gene_col, matched_gsms), drop = FALSE]

# Rename GSM columns to Sample IDs (e.g. GSM2358438 -> 6T)
colnames(expr_matched)[-1] <- expr_sample_ids[matched_gsms]

cat("Expression samples after matching:\n")
print(colnames(expr_matched))


## 6b. Metabolomics: metabolites + same Sample IDs, same order
metab_matched <- metab[, c(metab_feature_col, colnames(expr_matched)[-1]), drop = FALSE]

cat("Metabolomics samples after matching (should match expression):\n")
print(colnames(metab_matched))


## ===================================================
## 7. Save matched matrices
## ===================================================

#write.csv(
#  expr_matched,
#  expr_match_out,
#  row.names = FALSE
#)
cat("Matched tumour expression matrix written to:", expr_match_out, "\n")

#write.table(
#  metab_matched,
#  metab_match_out,
#  sep = "\t",
#  quote = FALSE,
#  row.names = FALSE
#)
cat("Matched tumour metabolite matrix written to:", metab_match_out, "\n")


#sd for metabolites
metab_known <- metab
dim(metab_known)

sd_metab <- apply(metab_known[,-1],1,sd)
w_sd_metab <- which(sd_metab > 0.8)

length(w_sd_metab)

metab_known_sd <- metab_known[w_sd_metab,]
dim(metab_known_sd)

#added 31st Dec 2025
## ===================================================
## 8. Correlations: GLUL / GOT2 / GLS vs Glutamine metabolism and metabolites
##    1) GLUL vs Glutamine metabolism score
##    2) GOT2 vs Glutamine
##    3) GOT2 vs Glutamic acid
##    4) GLS  vs Glutamic acid
##    5) GLS  vs Glutamine
## ===================================================

# Helper: extract a gene vector from expr_matched (handles multiple probes by median)
get_gene_vec <- function(expr_df, gene_symbol) {
  stopifnot(gene_col %in% colnames(expr_df))
  idx <- which(toupper(expr_df[[gene_col]]) == toupper(gene_symbol))
  if (length(idx) == 0) {
    stop(paste0("Gene not found in expression matrix: ", gene_symbol))
  }
  mat <- as.matrix(expr_df[idx, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  if (nrow(mat) == 1) {
    vec <- as.numeric(mat[1, ])
  } else {
    vec <- apply(mat, 2, median, na.rm = TRUE)
  }
  names(vec) <- colnames(expr_df)[-1]
  return(vec)
}

# Helper: extract a metabolite vector from metab_matched (flexible matching; handles duplicates by median)
get_metab_vec <- function(metab_df, pattern) {
  feature_names <- metab_df[[metab_feature_col]]
  idx <- grep(pattern, feature_names, ignore.case = TRUE)
  if (length(idx) == 0) {
    stop(paste0("Metabolite not found (pattern): ", pattern))
  }
  mat <- as.matrix(metab_df[idx, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  if (nrow(mat) == 1) {
    vec <- as.numeric(mat[1, ])
  } else {
    vec <- apply(mat, 2, median, na.rm = TRUE)
  }
  names(vec) <- colnames(metab_df)[-1]
  return(vec)
}

# Helper: correlation runner (pairwise complete obs)
run_cor <- function(x, y, label_x, label_y) {
  common_ids <- intersect(names(x), names(y))
  x2 <- x[common_ids]
  y2 <- y[common_ids]
  ok <- is.finite(x2) & is.finite(y2)
  x2 <- x2[ok]
  y2 <- y2[ok]
  cat("\n--------------------------------------------------\n")
  cat("Correlation:", label_x, "vs", label_y, "\n")
  cat("N:", length(x2), "\n")
  if (length(x2) < 3) {
    cat("Too few samples to run cor.test.\n")
    return(invisible(NULL))
  }
  ct <- suppressWarnings(cor.test(x2, y2, method = "pearson", exact = FALSE))
  cat("Pearson rho:", unname(ct$estimate), "\n")
  cat("P-value:", ct$p.value, "\n")
  invisible(ct)
}

# 8a) Extract gene expression vectors
glul_vec <- get_gene_vec(expr_matched, "GLUL")
got2_vec <- get_gene_vec(expr_matched, "GOT2")
gls_vec  <- get_gene_vec(expr_matched, "GLS")

# 8b) Define a simple "Glutamine metabolism" expression score (edit gene list if needed)
# Score = mean of listed genes per sample (log2 expression already in expr_matched input)
glutamine_metabolism_genes <- c(
  "GLS","GLUL","GOT1","GOT2","ASNS","GLUD1","GLUD2",
  "SLC1A5","SLC7A5","SLC3A2","CAD","GFPT1","GFPT2"
)

# Build score using genes present; multiple probes handled by median in get_gene_vec-like logic
present_genes <- glutamine_metabolism_genes[
  toupper(glutamine_metabolism_genes) %in% toupper(expr_matched[[gene_col]])
]

cat("\nGlutamine metabolism genes requested:", length(glutamine_metabolism_genes), "\n")
cat("Glutamine metabolism genes present:", length(present_genes), "\n")
if (length(present_genes) < 2) {
  stop("Too few glutamine metabolism genes found in expression matrix to compute a score.")
}

score_mat <- sapply(present_genes, function(g) get_gene_vec(expr_matched, g))
# sapply can return a vector if only 1 gene; enforce matrix
if (is.null(dim(score_mat))) {
  score_mat <- matrix(score_mat, ncol = 1)
  colnames(score_mat) <- present_genes
  rownames(score_mat) <- names(glul_vec)
}

glutamine_metabolism_score <- rowMeans(score_mat, na.rm = TRUE)
names(glutamine_metabolism_score) <- rownames(score_mat)

# 8c) Extract metabolite vectors (patterns are flexible; adjust if your feature names differ)
glutamine_vec    <- get_metab_vec(metab_matched, "^glutamine$|glutamine")
glutamicacid_vec <- get_metab_vec(metab_matched, "^glutamate|glutamic")

# 8d) Run requested correlations
# 1) GLUL vs Glutamine metabolism
run_cor(glul_vec, glutamine_vec, "GLUL (expression)", "Glutamine (metabolite)")
summary(lm(glutamine_vec~glul_vec))

# 2) GOT2 vs Glutamine
run_cor(got2_vec, glutamine_vec, "GOT2 (expression)", "Glutamine (metabolite)")
summary(lm(glutamine_vec~got2_vec))

# 3) GOT2 vs Glutamic acid
run_cor(got2_vec, glutamicacid_vec, "GOT2 (expression)", "Glutamic acid / Glutamate (metabolite)")
summary(lm(glutamicacid_vec~got2_vec))

# 4) GLS vs Glutamic acid
run_cor(gls_vec, glutamicacid_vec, "GLS (expression)", "Glutamic acid / Glutamate (metabolite)")
summary(lm(glutamicacid_vec~gls_vec))

# 5) GLS vs Glutamine
run_cor(gls_vec, glutamine_vec, "GLS (expression)", "Glutamine (metabolite)")
summary(lm(glutamine_vec~gls_vec))

# 6) glutamate vs Glutamine
run_cor(glutamicacid_vec, glutamine_vec, "Glutamic acid / Glutamate (metabolite)", "Glutamine (metabolite)")
summary(lm(glutamine_vec~glutamicacid_vec))


# 6) glutamate vs Glutamine
run_cor(glutamicacid_vec, glutamine_vec, "Glutamic acid / Glutamate (metabolite)", "Glutamine (metabolite)")

plot(glutamicacid_vec, glutamine_vec)

## ===================================================
## 9b. Violin plot (styled): Glutamine vs Glutamate
##     – matching reference figure colours & layout
## ===================================================

library(ggplot2)

# Prepare paired data
common_ids <- intersect(names(glutamicacid_vec), names(glutamine_vec))

glu <- glutamicacid_vec[common_ids]
gln <- glutamine_vec[common_ids]

ok <- is.finite(glu) & is.finite(gln)
glu <- glu[ok]
gln <- gln[ok]

violin_df <- data.frame(
  Metabolite = factor(
    c(rep("Glutamine", length(gln)),
      rep("Glutamate", length(glu))),
    levels = c("Glutamine", "Glutamate")
  ),
  Abundance = c(gln, glu)
)

# Wilcoxon paired test
wilcox_p <- wilcox.test(gln, glu, paired = TRUE)$p.value
p_label  <- paste0("Wilcoxon, p = ", formatC(wilcox_p, format = "e", digits = 2))

# Colours matched to example
metab_cols <- c(
  "Glutamine" = "#66C2A5",   # teal
  "Glutamate" = "#FC8D62"   # orange
)

# Plot
p_violin <- ggplot(violin_df, aes(x = Metabolite, y = Abundance, fill = Metabolite)) +
  geom_violin(
    color = "black",
    linewidth = 1,
    trim = FALSE,
    alpha = 0.9
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    color = "black",
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.08,
    size = 1.5,
    alpha = 0.6,
    color = "grey30"
  ) +
  scale_fill_manual(values = metab_cols) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Abundance of glutamine\nand glutamate",
    y = "Abundance"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(violin_df$Abundance, na.rm = TRUE) * 1.03,
    label = p_label,
    size = 5
  )

#pdf("20251231_Public_metabolic_cohort_glutamine_vs_glutamate_violin_plot.pdf")
print(p_violin)
#dev.off()

## ===================================================
## 9c. Scatter plot (styled): Glutamine vs Glutamate
##     – matching reference figure (points + red lm line + r/p label)
## ===================================================

library(ggplot2)

# Prepare paired data
common_ids <- intersect(names(glutamicacid_vec), names(glutamine_vec))

glu <- glutamicacid_vec[common_ids]
gln <- glutamine_vec[common_ids]

ok <- is.finite(glu) & is.finite(gln)
glu <- glu[ok]
gln <- gln[ok]

scatter_df <- data.frame(
  Glutamine = as.numeric(gln),
  Glutamate = as.numeric(glu)
)

# Pearson correlation (to match "R" style in reference)
ct <- cor.test(scatter_df$Glutamine, scatter_df$Glutamate, method = "pearson")
r_val <- unname(ct$estimate)
p_val <- ct$p.value

p_txt <- if (p_val < 2.2e-16) {
  "p < 2.2e-16"
} else {
  paste0("p = ", formatC(p_val, format = "e", digits = 2))
}

label_txt <- paste0("r = R = ", sprintf("%.2f", r_val), " ,  ", p_txt)

# Plot
p_scatter <- ggplot(scatter_df, aes(x = Glutamine, y = Glutamate)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 2.2, color = "red") +
  theme_classic(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text    = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1.2)
  ) +
  labs(
    x = "Glutamine",
    y = "Glutamate"
  ) +
  annotate(
    "text",
    x = min(scatter_df$Glutamine, na.rm = TRUE) + 0.1,
    y = max(scatter_df$Glutamate, na.rm = TRUE) - 0.1,
    label = label_txt,
    hjust = 0,
    vjust = 1,
    size = 4.5
  )

#pdf("20251231_Public_metabolic_cohort_glutamine_vs_glutamate_scatter_plot.pdf")
print(p_scatter)
#dev.off()
