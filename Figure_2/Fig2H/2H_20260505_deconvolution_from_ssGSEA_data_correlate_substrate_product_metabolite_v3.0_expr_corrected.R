################################################################################
# Figure Panel : 2H
# Description  : Cell-type deconvolution (MuSiC) of COLOSSUS bulk RNA-seq using
#                scRNA-seq reference, then linear regression of cell-type
#                fractions against substrate/product metabolite concentrations
#                (glutamine–glutamate axis). Generates forest-plot style
#                coefficient summary.
# Input data   : (1) Seurat object: query_seurat_mss_kras_with_scATOMIC_epi_
#                    malignancy.Rds
#                (2) COLOSSUS RNA-seq VST matrix (Combat-corrected)
#                (3) Metabolomics batch-corrected data
# Output       : 2H_..._coeff_summary_glutamine_glutamate_regression.pdf
# R packages   : Seurat, MuSiC, SingleCellExperiment, SummarizedExperiment,
#                dplyr, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : Seurat .Rds file and bulk expression matrix must be accessible.
#                MuSiC requires matching gene names between bulk and single-cell.
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2H — Cell-type deconvolution + metabolite regression
# ══════════════════════════════════════════════════════════════════════════════
# Strategy:
#   1. Use MuSiC to estimate cell-type composition of COLOSSUS bulk RNA-seq,
#      using a matched scRNA-seq dataset as reference.
#   2. Correlate estimated cell-type fractions with metabolite levels
#      (glutamine, glutamate) using linear regression.
#   3. Visualise regression coefficients as a forest plot.
#

# Seurat: single-cell RNA-seq processing; MuSiC: deconvolution algorithm.
# SingleCellExperiment/SummarizedExperiment: Bioconductor data containers.
library(Seurat)
library(MuSiC)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)


#Prepare the single-cell reference from the Seurat object

# ── Step 1: Load scRNA-seq reference ─────────────────────────────────────────
# Seurat object containing annotated COLOSSUS MSS KRAS-mutant single cells.
# Cell-type labels are stored in the 'main_annot3' metadata column.
seurat_object1<-readRDS("query_seurat_mss_kras_with_scATOMIC_epi_malignancy.Rds") #seurat object

seu <- seurat_object1
# Use 'main_annot3' as the cell-type label for MuSiC deconvolution.
seu$cell_type <- seu$main_annot3

# ── Step 2: Convert Seurat object to SingleCellExperiment ────────────────────
# MuSiC requires a SingleCellExperiment object as its reference.
# Convert Seurat -> SingleCellExperiment
sce <- as.SingleCellExperiment(seu)

# Add MuSiC-required metadata
# clusters: cell-type label
# samples:  sample/patient ID
# MuSiC uses 'clusters' for cell-type labels and 'samples' for donor IDs.
colData(sce)$clusters <- sce$cell_type
colData(sce)$samples  <- sce$orig.ident

# ── Step 3: Load bulk RNA-seq expression matrix ───────────────────────────────
# Combat-corrected VST matrix. Genes in rows, samples in columns.
# MuSiC requires overlapping gene names between bulk and scRNA-seq.
#Prepare bulk expression as SummarizedExperiment
expr_file <- "20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt"
#expr_file <- "20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_0.5M_reads_percent_3_5_Combat_Centre_Lane_0.5M_reads_batches_corrected_no_reads_involved.txt"

expr_df <- read.delim(expr_file,
                      stringsAsFactors = FALSE,
                      check.names = FALSE, row.names = 1)
dim(expr_df)

## If first column is gene names and the rest are samples:
#rownames(expr_df) <- expr_df[, 1]
#expr_df <- expr_df[, -1, drop = FALSE]

# Convert to plain numeric matrix for MuSiC.
bulk_mat <- as.matrix(expr_df)

min_val <- min(bulk_mat, na.rm = TRUE)
# Shift all values to be non-negative (MuSiC requires non-negative input).
bulk_shifted <- bulk_mat - min_val   # now minimum is 0

any(bulk_shifted < 0)  # should be FALSE


# ── Step 4: Align genes between bulk and single-cell ─────────────────────────
# Only genes present in both datasets can be used for deconvolution.
# Keep only genes that overlap
common_genes <- intersect(rownames(bulk_shifted), rownames(sce))
length(common_genes)

bulk_mtx <- as.matrix(bulk_shifted[common_genes, ])
sce_use  <- sce[common_genes, ]

# ------------------------------- #
# 3. Run MuSiC
# ------------------------------- #

# ── Step 5: Run MuSiC deconvolution ─────────────────────────────────────────
# Estimates cell-type proportions for each bulk sample.
# Output: Est.prop.weighted = weighted cell-type fraction matrix
#         (samples × cell types)
music_res <- music_prop(
  bulk.mtx = bulk_mtx,        # bulk matrix
  sc.sce   = sce_use,         # SingleCellExperiment
  clusters = "clusters",      # cell-type column
  samples  = "samples",       # patient/sample ID column
  verbose  = TRUE
)

# Output: cell-type fractions
# Extract the weighted proportion estimates (primary MuSiC output).
frac_est <- music_res$Est.prop.weighted
head(frac_est)

cell_fractions <- as.data.frame(frac_est)
#write.table(cell_fractions, "2H_20260103_Music_cell_fractions.txt", sep = "\t", quote = FALSE)
#cell_fractions$sample <- rownames(cell_fractions)

# Alternative: load pre-computed MuSiC fractions from Hari's analysis
# for comparison/validation purposes.
#Hari comparison
cell_fractions <- read.csv("MuSiC_celltype_proportions_134_rnaseq_samples_latest.csv")

  # ── Step 6: Load metabolomics data ───────────────────────────────────────
  # Batch-corrected metabolomics matrix — same file used in Fig 2B/2K/2L.
  #metabolomics
  metab <- read.delim("20200902_ColossusResults_Metabolomics_2_tissue_weight_Centre_run_batch_corrected_data.txt", stringsAsFactors = FALSE)
  dim(metab)
  
  unk <- grep("Unknown",metab[,1])
  
  metab_known <- metab[-unk,]
  dim(metab_known)
  
  sd_metab <- apply(metab_known[,-1],1,sd)
  w_sd_metab <- which(sd_metab > 1)
  
  length(w_sd_metab)
  
  metab_known_sd <- metab_known[w_sd_metab,]
  dim(metab_known_sd)

# Align metabolite samples with bulk
  mm <- match(colnames(bulk_mat),colnames(metab_known_sd))
  ww <- which(!is.na(mm))
metab_use <- metab_known_sd[,mm[ww]]
rownames(metab_use) <- metab_known_sd[,1]
dim(metab_use)

bulk_mat_1 <- bulk_mat[,ww]
dim(bulk_mat_1)
#correlation with gene expression

gee <- which(rownames(bulk_mat_1) == "GLUL")

## Core step: correlate vector with each column of the matrix of cell fractions
correlations <- cor(bulk_mat[gee,], data.matrix(cell_fractions), use = "complete.obs", method = "pearson")

correlations

summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat_1[gee,])))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Myeloid cells` ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$Fibroblasts ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[49,]) ~ cell_fractions$`Epithelial (Malignant)` ))


summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat_1[gee,])))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Myeloid cells` ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$Fibroblasts ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[51,]) ~ cell_fractions$`Epithelial (Malignant)` ))


gee <- which(rownames(bulk_mat_1) == "GLS")
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat_1[gee,])))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Myeloid cells` ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$Fibroblasts ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[49,]) ~ cell_fractions$`Epithelial (Malignant)` ))


summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat_1[gee,])))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Myeloid cells` ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$Fibroblasts ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[51,]) ~ cell_fractions$`Epithelial (Malignant)` ))


gee <- which(rownames(bulk_mat_1) == "GOT2")
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat_1[gee,])))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Myeloid cells` ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$Fibroblasts ))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[49,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[49,]) ~ cell_fractions$`Epithelial (Malignant)` ))


summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat_1[gee,])))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Myeloid cells` ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$Fibroblasts ))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[51,]) ~ as.numeric(bulk_mat[gee,]) + cell_fractions$`Epithelial (Malignant)` + cell_fractions$`Epithelial (Normal)`))
summary(lm(as.numeric(metab_use[51,]) ~ cell_fractions$`Epithelial (Malignant)` ))


## ============================================================
## Coefficient (estimate) summary plots for glutamine & glutamate
## Genes: GLUL, GLS, GOT2
## - Main figure: significant gene effects only (p < 0.05)
## - Supplement: all gene effects (incl. non-significant)
## ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

## -----------------------------
## 0) Inputs (assumed in env)
## -----------------------------
# bulk_mat: genes x samples matrix
# metab_use: metabolites x samples matrix
# (optional) cell_fractions: data.frame with sample rows aligned to columns

## -----------------------------
## 1) Define outcomes & genes
## -----------------------------
metab_map <- c(
  "Glutamine" = 49,
  "Glutamate" = 51
)

genes_use <- c("GLUL", "GLS", "GOT2")

## -----------------------------
## 2) Ensure sample alignment
## -----------------------------
# Require same samples/order across bulk_mat and metab_use
stopifnot(all(colnames(bulk_mat_1) %in% colnames(metab_use)))
stopifnot(all(colnames(metab_use) %in% colnames(bulk_mat_1)))

# Force the same order (use bulk_mat as reference)
metab_use_aligned <- metab_use[, colnames(bulk_mat_1), drop = FALSE]

## If you plan to adjust for cell fractions, also align them:
has_cf <- exists("cell_fractions")
if (has_cf) {
  # rownames(cell_fractions) should be sample IDs matching colnames(bulk_mat)
  if (!is.null(rownames(cell_fractions))) {
    cell_fractions_aligned <- cell_fractions[colnames(bulk_mat), , drop = FALSE]
  } else {
    # If no rownames, assume already aligned to bulk_mat columns
    cell_fractions_aligned <- cell_fractions
  }
}

## -----------------------------
## 3) Helper to fit & extract gene coefficient
## -----------------------------
fit_extract <- function(y, x_gene, gene_name,
                        adjust_epi = FALSE,
                        epi_mal = NULL, epi_nor = NULL) {
  
  df <- data.frame(
    y = as.numeric(y),
    gene = as.numeric(x_gene)
  )
  
  if (adjust_epi) {
    df$epi_mal <- as.numeric(epi_mal)
    df$epi_nor <- as.numeric(epi_nor)
    m <- lm(y ~ gene + epi_mal + epi_nor, data = df)
  } else {
    m <- lm(y ~ gene, data = df)
  }
  
  # Extract only the gene term
  tt <- broom::tidy(m, conf.int = TRUE) %>%
    dplyr::filter(term == "gene") %>%
    dplyr::mutate(
      gene = gene_name,
      model = ifelse(adjust_epi, "Adjusted: gene + epi_mal + epi_nor", "Unadjusted: gene only")
    ) %>%
    dplyr::select(model, gene, estimate, conf.low, conf.high, std.error, statistic, p.value)
  
  tt
}

## -----------------------------
## 4) Build coefficient tables
## -----------------------------
# Toggle this:
# - FALSE: gene-only models (cleanest coefficient summary)
# - TRUE: gene + epithelial fractions (composition-adjusted)
adjust_epi <- FALSE

if (adjust_epi) {
  if (!has_cf) stop("adjust_epi=TRUE but cell_fractions not found in environment.")
  
  epi_mal_vec <- cell_fractions_aligned$`Epithelial (Malignant)`
  epi_nor_vec <- cell_fractions_aligned$`Epithelial (Normal)`
}

coef_df <- lapply(names(metab_map), function(mname) {
  midx <- metab_map[[mname]]
  y <- metab_use_aligned[midx, ]
  
  do.call(rbind, lapply(genes_use, function(g) {
    gidx <- which(rownames(bulk_mat_1) == g)
    if (length(gidx) != 1) stop(paste("Gene not found uniquely in bulk_mat:", g))
    
    x <- bulk_mat_1[gidx, ]
    
    if (adjust_epi) {
      fit_extract(y, x, g, adjust_epi = TRUE,
                  epi_mal = epi_mal_vec, epi_nor = epi_nor_vec) %>%
        dplyr::mutate(metabolite = mname)
    } else {
      fit_extract(y, x, g, adjust_epi = FALSE) %>%
        dplyr::mutate(metabolite = mname)
    }
  }))
}) %>% dplyr::bind_rows()


## Standardize effect sizes? (optional)
## If you prefer comparable betas across genes, uncomment this and re-fit models with scale():
# (Keeping as raw coefficients by default, matching your summaries)

## -----------------------------
## 5) Create MAIN and SUPP datasets
## -----------------------------
coef_df <- coef_df %>%
  dplyr::mutate(
    sig = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE            ~ "ns"
    ),
    gene = factor(gene, levels = genes_use),
    metabolite = factor(metabolite, levels = c("Glutamine", "Glutamate"))
  )

coef_df <- coef_df %>%
  dplyr::mutate(
    p_label = ifelse(
      p.value < 0.001,
      "p < 0.001",
      paste0("p = ", formatC(p.value, format = "f", digits = 3))
    )
  )

coef_main <- coef_df %>% dplyr::filter(p.value < 0.05)   # significant only
coef_supp <- coef_df       

# all

## If you want to *force* positives only in main:
# coef_main <- coef_df %>% filter(p.value < 0.05, estimate > 0)

## -----------------------------
## 6) Plot functions
## -----------------------------
plot_coef <- function(df, title_text) {
  
  # Compute offset so p-values sit consistently outside CI
  x_range <- range(c(df$conf.low, df$conf.high), na.rm = TRUE)
  x_offset <- diff(x_range) * 0.12
  
  ggplot(df, aes(x = gene, y = estimate, colour = gene)) +
    
    # Zero line
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    
    # Point + CI
    geom_pointrange(
      aes(ymin = conf.low, ymax = conf.high),
      size = 0.8
    ) +
    
    # Significance stars near point
    geom_text(
      aes(label = sig),
      nudge_y = 0.02 * diff(x_range),
      size = 4,
      show.legend = FALSE
    ) +
    
    # p-values to the right of CI
    geom_text(
      aes(
        y = conf.high + x_offset,
        label = p_label
      ),
      hjust = 0,
      size = 3.4,
      show.legend = FALSE
    ) +
    
    coord_flip(clip = "off") +
    
    facet_wrap(~ metabolite, scales = "free_x") +
    
    scale_colour_manual(
      values = c(
        "GLUL" = "#009E73",  # green (synthesis)
        "GLS"  = "#D55E00",  # orange/red (consumption)
        "GOT2" = "#0072B2"   # blue (flux)
      )
    ) +
    
    labs(
      x = NULL,
      y = "Regression coefficient (estimate ± 95% CI)",
      title = title_text,
      subtitle = unique(df$model),
      colour = "Gene"
    ) +
    
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 40, 5.5, 5.5),
      legend.position = "top"
    )
}
p_main <- plot_coef(
  coef_main,
  title_text = "Gene–metabolite regression coefficients (significant terms only)"
)

p_supp <- plot_coef(
  coef_supp,
  title_text = "Gene–metabolite regression coefficients (all terms)"
)

## -----------------------------
## 7) Save figures
## -----------------------------
# Main
ggsave("2H_20260103_Fig_MAIN_coeff_summary_glutamine_glutamate_regression.pdf", p_main, width = 8, height = 4.5, useDingbats = FALSE)

# Supplement
ggsave("2H_20260103_Fig_SUPP_coeff_summary_glutamine_glutamate_regression_all.pdf", p_supp, width = 8, height = 4.5, useDingbats = FALSE)

## -----------------------------
## 8) Optional: export the coefficient table
## -----------------------------
write.csv(coef_df, "2H_20260103_CoeffTable_GLUL_GLS_GOT2_glutamine_glutamate.csv", row.names = FALSE)




## 2) Glutamate ~ gene + covariate   (gene = GLS/GLUL/GOT2)
##
## Output:
##   - MAIN: positive & significant (p < 0.05) gene coefficients only
##   - SUPP: all gene coefficients (all covariates, all genes)
##
## Plot style:
##   Forest plot (estimate ± 95% CI) with exact p-values annotated.
##
## Assumes in environment:
##   bulk_mat (genes x samples), metab_use (metabolites x samples),
##   cell_fractions (samples x cell-types) with rownames = sample IDs
## ============================================================

library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)

## -----------------------------
## 0) Settings you can edit
## -----------------------------
metab_map <- c("Glutamine" = 49, "Glutamate" = 51)
genes_use <- c("GLUL", "GLS", "GOT2")

# Covariates to test (must match colnames(cell_fractions))
covariates_use <- c(
  "Epithelial (Normal)",
  "Epithelial (Malignant)",
  "Myeloid cells",
  "Fibroblasts"
  # add/remove as needed
)

# MAIN filter
main_alpha <- 0.05
main_positive_only <- TRUE

# Output filenames
out_main_pdf <- "20260505_Fig_MAIN_oneGene_oneCov_forest.pdf"
out_supp_pdf <- "20260505_Fig_SUPP_oneGene_oneCov_forest.pdf"

## -----------------------------
## 1) Align samples
## -----------------------------
stopifnot(all(colnames(bulk_mat_1) %in% colnames(metab_use)))
metab_use_aligned <- metab_use[, colnames(bulk_mat_1), drop = FALSE]

if (!exists("cell_fractions")) stop("cell_fractions not found in environment.")
if (!is.null(rownames(cell_fractions))) {
  cell_frac_aligned <- cell_fractions[colnames(bulk_mat_1), , drop = FALSE]
} else {
  cell_frac_aligned <- cell_fractions
}

missing_cov <- setdiff(covariates_use, colnames(cell_frac_aligned))
if (length(missing_cov) > 0) {
  stop(paste("Missing covariates in cell_fractions:", paste(missing_cov, collapse = ", ")))
}

## -----------------------------
## 2) Fit all one-gene + one-cov models and extract gene term
## -----------------------------
fit_one_gene_one_cov <- function(y, x_gene, x_cov, gene_name, cov_name) {
  
  df <- data.frame(
    y = as.numeric(y),
    gene = as.numeric(x_gene),
    cov = as.numeric(x_cov)
  )
  
  m <- lm(y ~ gene + cov, data = df)
  
  broom::tidy(m, conf.int = TRUE) %>%
    dplyr::filter(term == "gene") %>%
    dplyr::mutate(
      gene = gene_name,
      covariate = cov_name,
      sig = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE ~ "ns"
      ),
      p_label = ifelse(
        p.value < 0.001,
        "p < 0.001",
        paste0("p = ", formatC(p.value, format = "f", digits = 3))
      )
    ) %>%
    dplyr::select(gene, covariate, estimate, conf.low, conf.high, p.value, sig, p_label)
}


all_coef <- lapply(names(metab_map), function(mname) {
  
  y <- as.numeric(metab_use_aligned[metab_map[[mname]], ])
  
  do.call(rbind, lapply(genes_use, function(g) {
    
    gi <- which(rownames(bulk_mat_1) == g)
    if (length(gi) != 1) stop(paste("Gene not found uniquely in bulk_mat:", g))
    x_gene <- as.numeric(bulk_mat_1[gi, ])
    
    do.call(rbind, lapply(covariates_use, function(cv) {
      x_cov <- as.numeric(cell_frac_aligned[[cv]])
      fit_one_gene_one_cov(y, x_gene, x_cov, gene_name = g, cov_name = cv) %>%
        dplyr::mutate(metabolite = mname)
    }))
  }))
}) %>% dplyr::bind_rows()

all_coef <- all_coef %>%
  dplyr::mutate(
    metabolite = factor(metabolite, levels = c("Glutamine", "Glutamate")),
    gene = factor(gene, levels = genes_use),
    covariate = factor(covariate, levels = covariates_use),
    label = paste0(gene, " + ", covariate)
  )

## -----------------------------
## 3) MAIN vs SUPP tables
## -----------------------------
coef_supp <- all_coef

coef_main <- all_coef %>%
  dplyr::filter(p.value < main_alpha) %>%
  { if (main_positive_only) dplyr::filter(., estimate > 0) else . }

if (nrow(coef_main) == 0) {
  warning("No terms passed MAIN filter. Consider main_positive_only=FALSE or relaxing main_alpha.")
}

## -----------------------------
## 4) Forest plot function (with p-values)
## -----------------------------
plot_one_gene_one_cov <- function(df, title_text) {
  
  # For consistent label placement, compute facet-wise offsets
  df <- df %>%
    group_by(metabolite) %>%
    dplyr::mutate(
      rng = max(conf.high, na.rm = TRUE) - min(conf.low, na.rm = TRUE),
      offset = ifelse(rng == 0, 0.1, rng * 0.18),
      y_p = conf.high + offset
    ) %>%
    ungroup()
  
  ggplot(df, aes(x = label, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.6) +
    geom_text(aes(label = sig), nudge_y = 0.02, size = 3.8) +
    geom_text(aes(y = y_p, label = p_label), hjust = 0, size = 3.1) +
    coord_flip(clip = "off") +
    facet_wrap(~ metabolite, scales = "free_x") +
    labs(
      x = NULL,
      y = "Gene coefficient (estimate ± 95% CI)",
      title = title_text,
      subtitle = "Each point: coefficient of the gene term from lm(metabolite ~ gene + covariate)"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 70, 5.5, 5.5)
    )
}

## -----------------------------
## 5) Make plots
## -----------------------------
p_main <- plot_one_gene_one_cov(
  coef_main,
  title_text = "One-gene + one-covariate models (positive & significant gene effects)"
)

p_supp <- plot_one_gene_one_cov(
  coef_supp,
  title_text = "One-gene + one-covariate models (all gene effects)"
)

print(p_main)
print(p_supp)

## -----------------------------
## 6) Save figures and tables
## -----------------------------
ggsave(out_main_pdf, p_main, width = 10, height = 6, useDingbats = FALSE)
ggsave(out_supp_pdf, p_supp, width = 10, height = 6, useDingbats = FALSE)

write.csv(all_coef, "20260505_CoeffTable_oneGene_oneCov_GLN_GLU.csv", row.names = FALSE)
write.csv(coef_main, "20260505_CoeffTable_MAIN_oneGene_oneCov_positive_sig.csv", row.names = FALSE)
