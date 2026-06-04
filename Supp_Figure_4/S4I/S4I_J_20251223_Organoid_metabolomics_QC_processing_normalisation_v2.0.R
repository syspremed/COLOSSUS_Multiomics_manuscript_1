################################################################################
# Figure Panels : Supplementary Figures 4I & 4J (QC/pre-processing step)
# Description  : Quality control, processing, and normalisation of organoid
#                metabolomics data (v2.0). Run this before S4I and S4J scripts.
# Input data   : Raw organoid metabolomics data
# Output       : Processed/normalised organoid metabolomics matrix
# R packages   : ggplot2, dplyr, tidyr (standard QC packages)
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : December 2025
# Note         : Must be run first; outputs feed into S4I and S4J analyses.
################################################################################

############################################################
# Metabolomics pipeline (samples + Mock + QC)
# - Imports a wide table (samples in rows, metabolites in columns)
# - Separates Sample / Mock / QC by row name patterns
# - Basic QC: missingness, Mock background filter, QC CV filter
# - Log2 transform with pseudocount
# - QC-anchored normalization using sweep(..., MARGIN=1, FUN="-")
# - Diagnostics: density/boxplot, PCA, QC CV plot
# - Exports normalized matrix + QC metrics
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

### 0) USER SETTINGS ###########################################################
# Input file: your pasted table saved as TSV (tab-delimited)
# It must contain:
#   - rows: samples (e.g., MCF002329_CC01_1, ... Mock..., ... QC...)
#   - columns: metadata rows ("Formula", "RT") may exist at top (we will remove)
#   - metabolite columns are numeric
#infile  <- "MCF002329_Colossus_results_organoids.txt"
infile  <- "MCF002329_Colossus_results_organoids__metabolite_HMDB_merged.tsv"

# Output prefix
outpref <- "MCF002329_metabolomics"

# Row name column (first column in your pasted data is sample ID)
sample_id_col <- "Sample"

# Patterns used to detect sample types from sample IDs
pat_qc   <- "QC"
pat_mock <- "Mock"

# Pseudocount for log transform (helps when zeros exist)
# Use 1 for peak areas; for very small intensities you can use 0.5 or 10.
pseudocount <- 1

# Filters (typical starting points)
max_missing_frac_per_feature <- 0.20   # drop metabolites missing in >20% of biological samples
qc_cv_threshold              <- 0.30   # keep metabolites with QC CV <= 30% (on raw or log scale; here raw)
mock_to_sample_ratio_thresh  <- 0.10   # drop metabolites where mock median >= 10% of sample median

do_pca <- TRUE
do_impute <- TRUE
################################################################################

### 1) LOAD DATA ###############################################################
raw <- read.delim(infile, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

if (!sample_id_col %in% colnames(raw)) {
  stop("Cannot find sample ID column: ", sample_id_col,
       "\nColumns are: ", paste(colnames(raw), collapse = ", "))
}

# Remove non-sample metadata rows (like "Formula", "RT", blanks) if present
meta_rows <- raw[[sample_id_col]] %in% c("Formula", "RT", "")
raw2 <- raw[!meta_rows, , drop = FALSE]

# Set rownames = sample IDs; drop the ID column from matrix
rownames(raw2) <- raw2[[sample_id_col]]
raw2[[sample_id_col]] <- NULL

# Coerce metabolite columns to numeric
mdat_raw <- as.data.frame(lapply(raw2[,-1], function(x) as.numeric(as.character(x))))
mdat_raw <- as.matrix(mdat_raw)
rownames(mdat_raw) <- rownames(raw2)
cat("Loaded matrix dim (samples x metabolites): ", dim(mdat_raw), "\n")

setwd("metabolomics")

### 2) DEFINE GROUPS (ROWS = samples) #########################################
sids <- rownames(mdat_raw)

is_qc   <- grepl(pat_qc,   sids, ignore.case = TRUE)
is_mock <- grepl(pat_mock, sids, ignore.case = TRUE)
is_bio  <- !(is_qc | is_mock)

cat("Counts: bio =", sum(is_bio), " QC =", sum(is_qc), " Mock =", sum(is_mock), "\n")

### 3) CLEANING ################################################################
mdat_raw[mdat_raw < 0] <- NA

# Missingness per feature (COLUMN) evaluated only on biological samples (ROWS)
feat_missing_bio <- colMeans(is.na(mdat_raw[is_bio, , drop = FALSE]))
keep_missing <- feat_missing_bio <= max_missing_frac_per_feature

cat("Features kept after missingness filter: ", sum(keep_missing), "/", ncol(mdat_raw), "\n")
mdat1 <- mdat_raw[, keep_missing, drop = FALSE]

#Density plot (RAW scale, sample-wise)
# identify sample groups
grp <- ifelse(grepl("QC", rownames(mdat_raw), ignore.case=TRUE), "QC",
              ifelse(grepl("Mock", rownames(mdat_raw), ignore.case=TRUE), "Mock", "Bio"))

cols <- c(Bio="grey40", QC="black", Mock="dodgerblue3")

# select samples to plot (e.g. first 12)
sel <- 1:min(12, nrow(mdat_raw))

par(mfrow=c(1,2))

# Density of first sample
plot(density(mdat_raw[sel[1], ], na.rm=TRUE),
     main="RAW intensity density (samples)",
     xlab="Raw peak intensity",
     lwd=3,
     col=cols[grp[sel[1]]])

# Overlay remaining samples
for (i in sel) {
  lines(density(mdat_raw[i, ], na.rm=TRUE),
        lwd=2,
        col=cols[grp[i]])
}

legend("topright", legend=names(cols), col=cols, lwd=3, bty="n")

#Boxplot (RAW scale, sample-wise)
boxplot(
  t(mdat_raw[sel, ]),
  outline=FALSE,
  las=2,
  col=cols[grp[sel]],
  border=cols[grp[sel]],
  main="RAW intensity boxplot (samples)",
  ylab="Raw peak intensity"
)

legend("topright", legend=names(cols), fill=cols, bty="n")

### 4) MOCK BACKGROUND FILTER (COLUMN-wise features) ###########################
if (sum(is_mock) >= 1 && sum(is_bio) >= 1) {
  bio_med  <- apply(mdat1[is_bio,  , drop = FALSE], 2, median, na.rm = TRUE)  # per COLUMN
  mock_med <- apply(mdat1[is_mock, , drop = FALSE], 2, median, na.rm = TRUE)  # per COLUMN
  
  ratio <- mock_med / pmax(bio_med, 1)
  keep_mock <- ratio < mock_to_sample_ratio_thresh
  
  cat("Features kept after Mock filter: ", sum(keep_mock), "/", ncol(mdat1), "\n")
  mdat2 <- mdat1[, keep_mock, drop = FALSE]
} else {
  warning("No Mock or no Bio samples detected; skipping Mock background filter.")
  mdat2 <- mdat1
}

### 5) QC CV FILTER (per feature COLUMN, using QC rows) ########################
if (sum(is_qc) >= 2) {
  qc_mat  <- mdat2[is_qc, , drop = FALSE]                 # QC rows
  qc_mean <- apply(qc_mat, 2, mean, na.rm = TRUE)         # per COLUMN
  qc_sd   <- apply(qc_mat, 2, sd,   na.rm = TRUE)         # per COLUMN
  qc_cv   <- qc_sd / qc_mean
  
  keep_qc_cv <- is.finite(qc_cv) & !is.na(qc_cv) & (qc_cv <= qc_cv_threshold)
  
  cat("Features kept after QC CV filter: ", sum(keep_qc_cv), "/", ncol(mdat2), "\n")
  mdat3 <- mdat2[, keep_qc_cv, drop = FALSE]
  qc_cv <- qc_cv[keep_qc_cv]
} else {
  warning("Not enough QCs for QC CV filter; skipping.")
  mdat3 <- mdat2
  qc_cv <- rep(NA_real_, ncol(mdat3)); names(qc_cv) <- colnames(mdat3)
}

### 6) LOG2 TRANSFORM ##########################################################
mdat_log <- log2(mdat3 + pseudocount)

# Log2 density + boxplot for first 12 samples (ROWS)
#pdf(paste0(outpref, "_01_log2_density_boxplot.pdf"), width = 9, height = 6)
par(mfrow=c(1,2))

plot(density(mdat_log[1, ], na.rm=TRUE),
     main="Log2 density (samples 1–12)",
     xlab="log2(intensity)", lwd=2, xlim=c(7,35))
for (i in 1:min(12, nrow(mdat_log))) {
  lines(density(mdat_log[i, ], na.rm=TRUE), lwd=1, col=i+1)
}

boxplot(t(mdat_log[1:min(12, nrow(mdat_log)), ]),
        outline=FALSE, las=2,
        main="Log2 boxplot (samples 1–12)",
        ylab="log2(intensity)")
#dev.off()

### 7) QC-ANCHORED NORMALIZATION (IMPORTANT: MARGIN=2) #########################
# Because: rows = samples, columns = metabolites
# So we subtract one value PER COLUMN (metabolite) from all rows (samples).
if (sum(is_qc) >= 1) {
  qc_median <- apply(mdat_log[is_qc, , drop = FALSE], 2, median, na.rm = TRUE) # per COLUMN
  mdat_norm <- sweep(mdat_log, 2, qc_median, FUN = "-")  # MARGIN=2 = columns
} else {
  warning("No QCs detected; skipping QC-anchored normalization.")
  mdat_norm <- mdat_log
}

#plot
library(reshape2)
library(ggplot2)

# If you created an imputed matrix, use the NON-imputed version for plotting if possible:
# mdat_plot <- mdat_norm_noimp   # ideal
# If you only have mdat_norm (imputed), we’ll filter out the imputation floor:

mdat_plot <- mdat_norm

# Identify imputation floor per metabolite (if you used min-1)
# Values <= (min observed - 0.5) are likely imputed floor; tune threshold if needed
col_floor <- apply(mdat_plot, 2, function(x) min(x, na.rm=TRUE))
floor_mask <- sweep(mdat_plot, 2, col_floor + 0.01, FUN="<=")  # near-min values
mdat_plot[floor_mask] <- NA  # remove floor values for density plotting

# Long format
df <- melt(mdat_plot, varnames=c("Sample","Metabolite"), value.name="log2_fc")
df <- df[is.finite(df$log2_fc), ]

grp <- ifelse(grepl("QC", rownames(mdat_plot), ignore.case=TRUE), "QC",
              ifelse(grepl("Mock", rownames(mdat_plot), ignore.case=TRUE), "Mock", "Bio"))
df$Group <- grp[df$Sample]

# Group-level density
p <- ggplot(df, aes(x=log2_fc, color=Group, fill=Group)) +
  geom_density(alpha=0.15, linewidth=1.2, adjust=1.1) +
  geom_vline(xintercept=0, linetype="dashed") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_bw(base_size=13) +
  labs(title="After QC normalisation (QC-anchored): group densities",
       x="log2 fold-change vs QC", y="Density")

p

# Order samples: Bio first, then QC, then Mock (or any order you prefer)
sample_levels <- rownames(mdat_plot)
df$Sample <- factor(df$Sample, levels=sample_levels)

p2 <- ggplot(df, aes(x=Sample, y=log2_fc, fill=Group)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.5) +
  geom_boxplot(width=0.15, outlier.shape=NA, alpha=0.9) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(title="After QC normalisation: per-sample distributions",
       x=NULL, y="log2 fold-change vs QC")

p2

#density per sample (but readable): facet by group
p3 <- ggplot(df, aes(x=log2_fc, group=Sample)) +
  geom_density(alpha=0.15, linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed") +
  coord_cartesian(xlim=c(-6,6)) +
  facet_wrap(~Group, ncol=1, scales="free_y") +
  theme_bw(base_size=13) +
  labs(title="After QC normalisation: sample densities by group",
       x="log2 fold-change vs QC", y="Density")

p3

### 8) IMPUTATION (optional) ###################################################
if (do_impute) {
  feat_min <- apply(mdat_norm, 2, function(x) min(x, na.rm = TRUE))  # per COLUMN
  for (j in seq_len(ncol(mdat_norm))) {
    mdat_norm[is.na(mdat_norm[, j]), j] <- feat_min[j] - 1
  }
}

### 9) DIAGNOSTICS #############################################################
# QC CV plot (raw scale CV computed earlier)
if (sum(is_qc) >= 2 && all(!is.na(qc_cv))) {
  qc_df <- data.frame(feature = names(qc_cv), qc_cv = qc_cv) |>
    arrange(desc(qc_cv))
  
  p <- ggplot(qc_df, aes(x = reorder(feature, qc_cv), y = qc_cv)) +
    geom_point() +
    coord_flip() +
    labs(title = "QC CV per metabolite (raw scale)", x = "Metabolite", y = "CV") +
    theme_bw(base_size = 10)
  
  # Save figure.
ggsave(paste0(outpref, "_02_qc_cv_plot.pdf"), p, width = 8, height = 10)
}

# PCA on normalized data (rows=samples)
if (do_pca) {
  pca <- prcomp(mdat_norm, center = TRUE, scale. = TRUE)
  
  grp <- ifelse(is_qc, "QC", ifelse(is_mock, "Mock", "Bio"))
  pca_df <- data.frame(sample = rownames(mdat_norm),
                       PC1 = pca$x[,1], PC2 = pca$x[,2],
                       type = grp)
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = type)) +
    geom_point(size = 3) +
    theme_bw(base_size = 12) +
    labs(title="PCA: QC-anchored normalized (log2 fold vs QC)",
         x="PC1", y="PC2")
  
  ggsave(paste0(outpref, "_03_pca.pdf"), p, width = 7, height = 5)
}

### 10) EXPORT #################################################################
write.table(
  cbind(Sample = rownames(mdat_norm), as.data.frame(mdat_norm)),
  file = paste0(outpref, "_normalized_log2_QCanchored.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

feat_report <- data.frame(
  feature = colnames(mdat_norm),
  missing_frac_bio = colMeans(is.na(mdat_raw[is_bio, colnames(mdat_norm), drop=FALSE])),
  qc_cv_raw = if (sum(is_qc) >= 2) qc_cv[colnames(mdat_norm)] else NA_real_
)

write.table(
  feat_report,
  file = paste0(outpref, "_feature_QC_report.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nDONE.\n",
    "Outputs:\n",
    " - ", outpref, "_normalized_log2_QCanchored.tsv\n",
    " - ", outpref, "_feature_QC_report.tsv\n",
    " - ", outpref, "_01_log2_density_boxplot.pdf\n",
    " - ", outpref, "_02_qc_cv_plot.pdf (if QCs>=2)\n",
    " - ", outpref, "_03_pca.pdf\n", sep="")

#sd

library(dplyr)
library(stringr)

# ------------ INPUTS ------------
met_mat   <- mdat_norm  # or mdat_log; samples x metabolites
hmdb_file <- "HMDB_IDS_from_Bart.txt"

sd_thresh <- 0.5

# choose which samples to compute SD over:
# usually Bio only (recommended)
is_qc   <- grepl("QC", rownames(met_mat),   ignore.case=TRUE)
is_mock <- grepl("Mock", rownames(met_mat), ignore.case=TRUE)
is_bio  <- !(is_qc | is_mock)

mat_for_sd <- met_mat[is_bio, , drop=FALSE]   # SD across Bio samples

# ------------ CANONICAL NAME FUNCTION ------------
canon_name <- function(x) {
  x %>%
    as.character() %>%
    str_squish() %>%
    str_replace_all("\\.", " ") %>%
    str_replace_all("[_\\-]+", " ") %>%
    str_replace_all("[^[:alnum:] ]+", "") %>%
    str_to_lower()
}

# If your column names contain suffixes like _RT_mz, uncomment this stripper:
# strip_suffix <- function(x) sub("_[0-9.]+_[0-9.]+$", "", x)
strip_suffix <- function(x) x

# ------------ 1) SD FILTER ------------
met_sd <- apply(mat_for_sd, 2, sd, na.rm=TRUE)  # per metabolite (column)
keep <- met_sd > sd_thresh

sel_names <- names(met_sd)[keep]
sel_sd    <- met_sd[keep]

cat("Selected metabolites (SD >", sd_thresh, "): ", length(sel_names), "\n")

# ------------ 2) LOAD + PREP HMDB ------------
hmdb <- read.delim(hmdb_file, sep="\t", check.names=FALSE, stringsAsFactors=FALSE) %>%
  mutate(
    name_key = canon_name(name),
    dbId_clean = str_squish(dbId),
    ppm = suppressWarnings(as.numeric(`m/z error (ppm)`)),
    msigma = suppressWarnings(as.numeric(mSigma))
  )

# pick best HMDB hit per canonical name (lowest |ppm| then lowest mSigma)
hmdb_best <- hmdb %>%
  filter(!is.na(name_key), name_key != "", !is.na(dbId_clean), dbId_clean != "") %>%
  group_by(name_key) %>%
  arrange(abs(ppm), msigma, .by_group=TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(name_key, HMDB_ID = dbId_clean, formula, ion, mz, ppm, msigma)

# ------------ 3) MAP SELECTED METABOLITES ------------
res <- tibble(metabolite = sel_names) %>%
  mutate(
    metabolite_sd = as.numeric(sel_sd[metabolite]),
    name_key = canon_name(strip_suffix(metabolite))
  ) %>%
  left_join(hmdb_best, by="name_key") %>%
  arrange(desc(metabolite_sd)) %>%
  select(metabolite, metabolite_sd, HMDB_ID, formula, ion, mz, ppm, msigma)

cat("Mapped HMDB IDs: ", sum(!is.na(res$HMDB_ID)), "/", nrow(res), "\n")

# ------------ 4) WRITE OUTPUT ------------
write.table(res,
            file = paste0("metabolites_SDgt", sd_thresh, "_with_HMDB.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

res

# additional column with HMDB Ids

library(dplyr)
library(stringr)

# res is your data.frame / tibble shown in the screenshot
# columns include: metabolite, metabolite_sd, HMDB_ID, ...

res2 <- res %>%
  mutate(
    # extract HMDB id from metabolite string
    HMDB_from_name = str_extract(metabolite, "HMDB\\d{5}"),
    
    # make a clean metabolite name by removing the HMDB suffix/pattern
    metabolite_name = metabolite %>%
      str_replace("_?HMDB\\d{5}", "") %>%  # remove "_HMDB12345" or "HMDB12345"
      str_replace_all("_+$", "")           # drop trailing underscores if any
  ) %>%
  # If HMDB_ID column exists but is NA, fill it using extracted HMDB
  mutate(
    HMDB_ID = if ("HMDB_ID" %in% names(.)) coalesce(HMDB_ID, HMDB_from_name) else HMDB_from_name
  ) %>%
  select(metabolite_name, HMDB_ID, everything(), -HMDB_from_name)

# quick check
head(res2[, c("metabolite", "metabolite_name", "HMDB_ID")], 20)

# write out
write.table(res2, "metabolites_SDgt0.5_with_HMDB_with_extracted_HMDB.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



############################################################
# Select metabolites with SD > 0.5 across 3 RAS-WT samples
# + map to HMDB IDs (handles '.' vs space name issues)
#
# Assumptions:
#  - met_mat: samples in rows, metabolites in columns (use mdat_norm ideally)
#  - rownames(met_mat) are sample IDs
#  - HMDB file has columns: name, dbId, mSigma, m/z error (ppm), etc.
############################################################

library(dplyr)
library(stringr)

# -------------------- INPUTS --------------------
met_mat   <- mdat_norm  # <- change if needed (samples x metabolites)
sd_thresh <- 0.5

hmdb_file <- "HMDB_IDS_from_Bart.txt"

# Put your 3 RAS wild-type sample IDs here (must match rownames(met_mat) exactly)
ras_wt_samples <- c(
  "MCF002329_CC01_1",
  "MCF002329_CC02_1",
  "MCF002329_CC03_1"
)

# -------------------- HELPERS --------------------
canon_name <- function(x) {
  x %>%
    as.character() %>%
    str_squish() %>%
    str_replace_all("\\.", " ") %>%          # Penicilloic.acid -> Penicilloic acid
    str_replace_all("[_\\-]+", " ") %>%      # underscores/hyphens -> space
    str_replace_all("[^[:alnum:] ]+", "") %>%# drop punctuation
    str_to_lower()
}

# -------------------- CHECKS --------------------
missing_ids <- setdiff(ras_wt_samples, rownames(met_mat))
if (length(missing_ids) > 0) {
  stop("These RAS-WT sample IDs are not found in rownames(met_mat):\n  ",
       paste(missing_ids, collapse = "\n  "))
}

# -------------------- 1) SUBSET TO RAS-WT --------------------
mat_wt <- met_mat[ras_wt_samples, , drop = FALSE]

# -------------------- 2) SD PER METABOLITE (ACROSS WT SAMPLES) ----------------
wt_sd <- apply(mat_wt, 2, sd, na.rm = TRUE)

# -------------------- 3) FILTER SD > THRESHOLD ----------------
keep <- wt_sd > sd_thresh & is.finite(wt_sd)
sel_names <- names(wt_sd)[keep]
sel_sd    <- wt_sd[keep]

cat("Selected metabolites (RAS-WT SD >", sd_thresh, "): ",
    length(sel_names), " / ", ncol(met_mat), "\n", sep="")

# -------------------- 4) LOAD HMDB + CHOOSE BEST HIT PER NAME -----------------
hmdb <- read.delim(hmdb_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    name_key   = canon_name(name),
    HMDB_ID    = str_squish(dbId),
    ppm        = suppressWarnings(as.numeric(`m/z error (ppm)`)),
    msigma     = suppressWarnings(as.numeric(mSigma)),
    mz_numeric = suppressWarnings(as.numeric(mz))
  )

# Best hit per canonical name: lowest |ppm| then lowest mSigma
hmdb_best <- hmdb %>%
  filter(!is.na(name_key), name_key != "", !is.na(HMDB_ID), HMDB_ID != "") %>%
  group_by(name_key) %>%
  arrange(abs(ppm), msigma, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(name_key, HMDB_ID, formula, ion, mz = mz_numeric, ppm, msigma)

# -------------------- 5) MAP SELECTED METABOLITES TO HMDB ---------------------
res <- tibble(
  metabolite = sel_names,
  SD_RAS_WT  = as.numeric(sel_sd[sel_names]),
  name_key   = canon_name(sel_names)
) %>%
  left_join(hmdb_best, by = "name_key") %>%
  arrange(desc(SD_RAS_WT)) %>%
  select(metabolite, SD_RAS_WT, HMDB_ID, formula, ion, mz, ppm, msigma)

cat("Mapped HMDB IDs: ", sum(!is.na(res$HMDB_ID)), " / ", nrow(res), "\n", sep="")

# -------------------- 6) WRITE OUTPUTS -------------------------
out1 <- paste0("RAS_WT_metabolites_SDgt", sd_thresh, "_with_HMDB.tsv")
write.table(res, file = out1, sep = "\t", quote = FALSE, row.names = FALSE)

# additional column with HMDB Ids

library(dplyr)
library(stringr)

# res is your data.frame / tibble shown in the screenshot
# columns include: metabolite, metabolite_sd, HMDB_ID, ...

res2 <- res %>%
  mutate(
    # extract HMDB id from metabolite string
    HMDB_from_name = str_extract(metabolite, "HMDB\\d{5}"),
    
    # make a clean metabolite name by removing the HMDB suffix/pattern
    metabolite_name = metabolite %>%
      str_replace("_?HMDB\\d{5}", "") %>%  # remove "_HMDB12345" or "HMDB12345"
      str_replace_all("_+$", "")           # drop trailing underscores if any
  ) %>%
  # If HMDB_ID column exists but is NA, fill it using extracted HMDB
  mutate(
    HMDB_ID = if ("HMDB_ID" %in% names(.)) coalesce(HMDB_ID, HMDB_from_name) else HMDB_from_name
  ) %>%
  select(metabolite_name, HMDB_ID, everything(), -HMDB_from_name)

# quick check
head(res2[, c("metabolite", "metabolite_name", "HMDB_ID")], 20)

# write out
write.table(res2, "RAS_WT_metabolites_SDgt0.5_with_HMDB_with_extracted_HMDB.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# Optional: save the filtered matrix (WT samples x selected metabolites)
mat_wt_sel <- mat_wt[, sel_names, drop = FALSE]
out2 <- paste0("RAS_WT_matrix_SDgt", sd_thresh, ".tsv")
write.table(cbind(Sample = rownames(mat_wt_sel), as.data.frame(mat_wt_sel)),
            file = out2, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote:\n  ", out1, "\n  ", out2, "\n", sep="")

############################################################
# Differential abundance analysis: Mutant vs Wildtype
# - QC and Mock samples are excluded
# - WT samples are explicitly defined
# - ALL remaining biological samples are Mutant
#
# Data format:
#   rows    = samples
#   columns = metabolites
#
# Input should be QC-normalised log2 data (mdat_norm)
############################################################

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
})

# --------------------- INPUTS ------------------------------
met_mat <- mdat_norm   # samples x metabolites (QC-normalised log2)

# Define WT sample IDs (EXACT rownames)
wt_samples <- c(
  "MCF002329_CC01_1",
  "MCF002329_CC02_1",
  "MCF002329_CC03_1"
)

# Analysis thresholds
min_non_na_per_group <- 2
fdr_cutoff  <- 0.2
logfc_cutoff <- 0.5

# ------------------- STEP 1: REMOVE QC & MOCK --------------
is_qc   <- grepl("QC",   rownames(met_mat), ignore.case = TRUE)
is_mock <- grepl("Mock", rownames(met_mat), ignore.case = TRUE)

bio_mat <- met_mat[!(is_qc | is_mock), , drop = FALSE]

cat("Biological samples retained:", nrow(bio_mat), "\n")

# ------------------- STEP 2: DEFINE GROUPS -----------------
missing_wt <- setdiff(wt_samples, rownames(bio_mat))
if (length(missing_wt) > 0) {
  stop("WT sample IDs not found:\n", paste(missing_wt, collapse = "\n"))
}

mut_samples <- setdiff(rownames(bio_mat), wt_samples)

cat("WT samples :", length(wt_samples), "\n")
cat("Mut samples:", length(mut_samples), "\n")

# Subset to analysis samples
use_samples <- c(wt_samples, mut_samples)
X <- bio_mat[use_samples, , drop = FALSE]

group <- factor(
  c(rep("WT", length(wt_samples)),
    rep("Mut", length(mut_samples))),
  levels = c("WT", "Mut")
)

# ---------------- STEP 3: FEATURE FILTERING ----------------
nonNA_WT  <- colSums(!is.na(X[group == "WT",  , drop = FALSE]))
nonNA_Mut <- colSums(!is.na(X[group == "Mut", , drop = FALSE]))

keep_feat <- (nonNA_WT >= min_non_na_per_group) &
  (nonNA_Mut >= min_non_na_per_group)

Xf <- X[, keep_feat, drop = FALSE]

cat("Metabolites tested:", ncol(Xf), "/", ncol(X), "\n")

# ---------------- STEP 4: IMPUTATION (limma) ---------------
feat_min <- apply(Xf, 2, function(v) min(v, na.rm = TRUE))
for (j in seq_len(ncol(Xf))) {
  Xf[is.na(Xf[, j]), j] <- feat_min[j] - 1
}

# ---------------- STEP 5: LIMMA ----------------------------
Y <- t(Xf)  # metabolites x samples

design <- model.matrix(~ group)
fit <- lmFit(Y, design)
fit <- eBayes(fit)

tt <- topTable(fit, coef = "groupMut", number = Inf, sort.by = "P")

res <- tt %>%
  tibble::rownames_to_column("metabolite") %>%
  rename(log2FC = logFC) %>%
  mutate(
    direction = case_when(
      adj.P.Val < fdr_cutoff & log2FC >=  logfc_cutoff ~ "Up_in_Mut",
      adj.P.Val < fdr_cutoff & log2FC <= -logfc_cutoff ~ "Down_in_Mut",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(adj.P.Val)

# ---------------- STEP 6: OUTPUT TABLE ---------------------
write.table(
  res,
  file = "Differential_metabolites_WT_vs_Mut.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Significant (FDR <", fdr_cutoff, "): ",
    sum(res$adj.P.Val < fdr_cutoff), "\n")

# ---------------- STEP 7: VOLCANO PLOT ---------------------
p_vol <- ggplot(res, aes(x = log2FC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    title = "Differential metabolites: Mutant vs WT",
    x = "log2 fold-change (Mut − WT)",
    y = "-log10(FDR)"
  )

ggsave("Volcano_WT_vs_Mut.pdf", p_vol, width = 7, height = 5)

# ---------------- STEP 8: TOP 50 HEATMAP -------------------
topN <- 50
top_feats <- head(res$metabolite, min(topN, nrow(res)))

Z <- t(Xf[, top_feats, drop = FALSE])
Z <- t(scale(t(Z)))   # z-score per metabolite

pdf("Top50_heatmap_WT_vs_Mut.pdf", width = 10, height = 7)
heatmap(
  Z,
  Colv = NA,
  scale = "none",
  labCol = use_samples,
  labRow = top_feats,
  main = "Top differential metabolites (WT vs Mut)"
)
dev.off()

cat("DONE.\n",
    "Outputs:\n",
    " - Differential_metabolites_WT_vs_Mut.tsv\n",
    " - Volcano_WT_vs_Mut.pdf\n",
    " - Top50_heatmap_WT_vs_Mut.pdf\n")


############################################################
# Select metabolites with SD > 0.5 across MUTANT samples
# + map to HMDB IDs
#
# Assumptions:
#  - met_mat: samples in rows, metabolites in columns
#  - Use QC-normalised log2 data (mdat_norm)
#  - WT samples are explicitly listed
#  - QC and Mock samples are excluded
############################################################

library(dplyr)
library(stringr)

# -------------------- INPUTS --------------------
met_mat <- mdat_norm   # samples x metabolites (QC-normalised log2)

# Define WT sample IDs (EXACT rownames)
wt_samples <- c(
  "MCF002329_CC01_1",
  "MCF002329_CC02_1",
  "MCF002329_CC03_1"
)

sd_thresh <- 0.5

hmdb_file <- "HMDB_IDS_from_Bart.txt"

# -------------------- HELPERS --------------------
canon_name <- function(x) {
  x %>%
    as.character() %>%
    str_squish() %>%
    str_replace_all("\\.", " ") %>%          # Penicilloic.acid -> Penicilloic acid
    str_replace_all("[_\\-]+", " ") %>%
    str_replace_all("[^[:alnum:] ]+", "") %>%
    str_to_lower()
}

# -------------------- 1) REMOVE QC & MOCK -----------------
is_qc   <- grepl("QC",   rownames(met_mat), ignore.case = TRUE)
is_mock <- grepl("Mock", rownames(met_mat), ignore.case = TRUE)

bio_mat <- met_mat[!(is_qc | is_mock), , drop = FALSE]

# -------------------- 2) DEFINE MUTANT SAMPLES ------------
missing_wt <- setdiff(wt_samples, rownames(bio_mat))
if (length(missing_wt) > 0) {
  stop("WT sample IDs not found:\n", paste(missing_wt, collapse = "\n"))
}

mut_samples <- setdiff(rownames(bio_mat), wt_samples)

cat("Mutant samples:", length(mut_samples), "\n")

mat_mut <- bio_mat[mut_samples, , drop = FALSE]

# -------------------- 3) SD PER METABOLITE ----------------
mut_sd <- apply(mat_mut, 2, sd, na.rm = TRUE)

# -------------------- 4) FILTER SD > 0.5 ------------------
keep <- mut_sd > sd_thresh & is.finite(mut_sd)

sel_names <- names(mut_sd)[keep]
sel_sd    <- mut_sd[keep]

cat("Selected metabolites (Mut SD >", sd_thresh, "): ",
    length(sel_names), "\n", sep="")

# -------------------- 5) LOAD HMDB + BEST HIT --------------
hmdb <- read.delim(hmdb_file, sep = "\t",
                   check.names = FALSE,
                   stringsAsFactors = FALSE) %>%
  mutate(
    name_key = canon_name(name),
    HMDB_ID  = str_squish(dbId),
    ppm      = suppressWarnings(as.numeric(`m/z error (ppm)`)),
    msigma   = suppressWarnings(as.numeric(mSigma)),
    mz       = suppressWarnings(as.numeric(mz))
  ) %>%
  filter(!is.na(name_key), name_key != "",
         !is.na(HMDB_ID), HMDB_ID != "")

# Best HMDB hit per name: lowest |ppm| → lowest mSigma
hmdb_best <- hmdb %>%
  group_by(name_key) %>%
  arrange(abs(ppm), msigma, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(name_key, HMDB_ID, formula, ion, mz, ppm, msigma)

# -------------------- 6) MAP TO HMDB ----------------------
res <- tibble(
  metabolite = sel_names,
  SD_Mutant  = as.numeric(sel_sd[sel_names]),
  name_key   = canon_name(sel_names)
) %>%
  left_join(hmdb_best, by = "name_key") %>%
  arrange(desc(SD_Mutant)) %>%
  select(metabolite, SD_Mutant, HMDB_ID, formula, ion, mz, ppm, msigma)

cat("Mapped to HMDB IDs:", sum(!is.na(res$HMDB_ID)),
    "/", nrow(res), "\n")

# -------------------- 7) SAVE OUTPUT ----------------------
out_file <- paste0("Mutant_metabolites_SDgt", sd_thresh, "_with_HMDB.tsv")

write.table(res,
            file = out_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Written:", out_file, "\n")

# additional column with HMDB Ids

library(dplyr)
library(stringr)

# res is your data.frame / tibble shown in the screenshot
# columns include: metabolite, metabolite_sd, HMDB_ID, ...

res2 <- res %>%
  mutate(
    # extract HMDB id from metabolite string
    HMDB_from_name = str_extract(metabolite, "HMDB\\d{5}"),
    
    # make a clean metabolite name by removing the HMDB suffix/pattern
    metabolite_name = metabolite %>%
      str_replace("_?HMDB\\d{5}", "") %>%  # remove "_HMDB12345" or "HMDB12345"
      str_replace_all("_+$", "")           # drop trailing underscores if any
  ) %>%
  # If HMDB_ID column exists but is NA, fill it using extracted HMDB
  mutate(
    HMDB_ID = if ("HMDB_ID" %in% names(.)) coalesce(HMDB_ID, HMDB_from_name) else HMDB_from_name
  ) %>%
  select(metabolite_name, HMDB_ID, everything(), -HMDB_from_name)

# quick check
head(res2[, c("metabolite", "metabolite_name", "HMDB_ID")], 20)

# write out
write.table(res2, "Mutant_metabolites_SDgt_0.5_with_HMDB_with_extracted_HMDB.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
