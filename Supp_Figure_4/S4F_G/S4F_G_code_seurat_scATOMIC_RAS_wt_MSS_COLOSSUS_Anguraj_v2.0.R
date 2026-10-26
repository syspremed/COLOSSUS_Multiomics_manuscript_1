################################################################################
# Figure Panel : Supp S4F, S4G
# Description  : scATOMIC classification of epithelial cells into cancer vs normal (per-sample
#                run) in COLOSSUS RAS-wildtype MSI/MSS snRNA-seq (4 samples). Adds scATOMIC_pred_epi
#                and epi_malignancy to the Seurat metadata and saves the updated object.
# Input data   : seurat_merged_harmony_main_annot_kras_wt_4_samples_20260302.rds
# Output       : query_seurat_WILDTYPE_MSI_MSS_with_scATOMIC_epi_malignancy_20260302.Rds
# R packages   : reticulate, Seurat, Matrix, scATOMIC
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Requires a Python environment where `import magic` works; set RETICULATE_PYTHON.
################################################################################

suppressPackageStartupMessages({
  library(reticulate)
  library(Seurat)
  library(Matrix)
  library(scATOMIC)
})


# ----------------------------
# 0) (Optional) Confirm python + magic module
#     Run before scATOMIC in a fresh session
# ----------------------------

# Set to the Python binary of an environment where `import magic` works
Sys.setenv(RETICULATE_PYTHON = "path/to/python")
library(reticulate)
py_config()
stopifnot(py_module_available("magic"))
cat("magic available\n")

cat("\n--- Python config ---\n")
print(reticulate::py_config())
stopifnot(reticulate::py_module_available("magic"))
cat("Python module 'magic' is available.\n")

# ----------------------------
# 0) Read original Seurat object
# ----------------------------

#seu_obj <- readRDS("seurat_merged_harmony_main_annot.Rds")


#length(table(seu_obj$orig.ident))
#table(seu_obj$orig.ident)

#samples_keep <- c("snrMPP302_03",
#                  "snrMPP302_15",
 #                 "snrMPP302_19",
  #                "snrMPP302_20")

#seu_subset <- subset(
 # seu_obj,
#  subset = orig.ident %in% samples_keep
#)

#table(seu_subset$orig.ident)

#saveRDS(seu_subset, file = "seurat_merged_harmony_main_annot_kras_wt_4_samples_20260302.rds")

# ----------------------------
# 1) Read Seurat object
# ----------------------------
seurat_object1 <- readRDS(
  "seurat_merged_harmony_main_annot_kras_wt_4_samples_20260302.rds"
)

cat("\nCell-type distribution (main_annot2):\n")
print(table(seurat_object1$main_annot2))

# ----------------------------
# 2) Subset epithelial cells
# ----------------------------
epi <- subset(seurat_object1, subset = main_annot2 == "Epithelial cells")
DefaultAssay(epi) <- "RNA"

if (!"RNA" %in% names(epi@assays)) stop("RNA assay not found in epithelial subset.")
if (ncol(epi) == 0) stop("No epithelial cells found (main_annot2 == 'Epithelial cells').")

cat("\nEpithelial cells:", ncol(epi), "\n")

# ----------------------------
# 3) Choose sample column for per-sample run
# ----------------------------
split_col <- "orig.ident"
if (!split_col %in% colnames(epi@meta.data)) {
  stop(paste0("Metadata column '", split_col, "' not found. Available columns:\n",
              paste(colnames(epi@meta.data), collapse = ", ")))
}

samples <- unique(as.character(epi@meta.data[[split_col]]))
samples <- samples[!is.na(samples)]
cat("Number of samples (", split_col, "): ", length(samples), "\n", sep = "")

# ----------------------------
# 4) Find scATOMIC runner function
# ----------------------------
get_runner <- function() {
  ns <- asNamespace("scATOMIC")
  if (exists("run_scATOMIC", where = ns, inherits = FALSE)) {
    return(get("run_scATOMIC", envir = ns))
  }
  if (exists("scATOMIC_pipeline", where = ns, inherits = FALSE)) {
    return(get("scATOMIC_pipeline", envir = ns))
  }
  if (exists("run_scATOMIC_pipeline", where = ns, inherits = FALSE)) {
    return(get("run_scATOMIC_pipeline", envir = ns))
  }
  stop("No scATOMIC runner found in namespace.")
}
runner_fun <- get_runner()

# ----------------------------
# 5) Parameters (edit if needed)
# ----------------------------
qc_max_pct_mt      <- 25      # drop cells with >= this % MT
qc_min_features    <- 500     # drop cells with <= this many detected genes
min_cells_per_gene <- 10      # keep genes expressed in >= this many cells
max_cells_per_run  <- 15000   # downsample within sample above this (scATOMIC guidance)
use_CNVs           <- FALSE   # set TRUE for CNV refinement (slower; extra deps)
mc_cores           <- 1       # keep 1 to avoid parallel python surprises
min_prop           <- 0.5     # summary threshold

# ----------------------------
# 6) Helper: QC + filtering on sparse counts
# ----------------------------
filter_counts <- function(x) {
  stopifnot(inherits(x, "dgCMatrix"))
  
  total_umi <- Matrix::colSums(x)
  mt_idx <- grep("^MT-", rownames(x))
  
  pct_mt <- rep(0, ncol(x))
  names(pct_mt) <- colnames(x)
  if (length(mt_idx) > 0) {
    pct_mt <- Matrix::colSums(x[mt_idx, , drop = FALSE]) / total_umi * 100
  }
  
  nFeature <- Matrix::colSums(x > 0)
  
  keep_cells <- names(which(pct_mt < qc_max_pct_mt & nFeature > qc_min_features))
  x <- x[, keep_cells, drop = FALSE]
  
  if (ncol(x) == 0) return(x)
  
  keep_genes <- Matrix::rowSums(x > 0) >= min_cells_per_gene
  x <- x[keep_genes, , drop = FALSE]
  
  x
}

# ----------------------------
# 7) Output vectors (for epithelial cells only)
# ----------------------------
epi_cells_all <- colnames(epi)
pred_scATOMIC <- setNames(rep(NA_character_, length(epi_cells_all)), epi_cells_all)
malignancy    <- setNames(rep(NA_character_, length(epi_cells_all)), epi_cells_all)

# ----------------------------
# 8) Run scATOMIC per sample
# ----------------------------
for (s in samples) {
  cat("\n============================================================\n")
  cat("Running scATOMIC for ", split_col, " = ", s, "\n", sep = "")
  
  cells_s <- rownames(epi@meta.data)[as.character(epi@meta.data[[split_col]]) == s]
  if (length(cells_s) == 0) next
  
  epi_s <- subset(epi, cells = cells_s)
  x <- epi_s@assays$RNA@counts

  
  cat("Raw: ", nrow(x), " genes x ", ncol(x), " cells\n", sep = "")
  
  x <- filter_counts(x)
  cat("After QC/filter: ", nrow(x), " genes x ", ncol(x), " cells\n", sep = "")
  
  if (ncol(x) < 50) {
    cat("Skipping (too few cells after QC)\n")
    next
  }
  if (nrow(x) < 500) {
    cat("Skipping (too few genes after filter)\n")
    next
  }
  
  if (ncol(x) > max_cells_per_run) {
    set.seed(1)
    keep <- sample(colnames(x), max_cells_per_run)
    x <- x[, keep, drop = FALSE]
    cat("Downsampled to: ", ncol(x), " cells\n", sep = "")
  }
  
  # Run scATOMIC
  pred_list <- try(runner_fun(x), silent = TRUE)
  if (inherits(pred_list, "try-error")) {
    cat("scATOMIC runner failed for ", s, ":\n", sep = "")
    cat(as.character(pred_list), "\n")
    next
  }
  
  # Summarize predictions
  res_sum <- try(
    create_summary_matrix(
      prediction_list = pred_list,
      raw_counts      = x,
      use_CNVs        = use_CNVs,
      modify_results  = TRUE,
      min_prop        = min_prop,
      mc.cores        = mc_cores
    ),
    silent = TRUE
  )
  
  if (inherits(res_sum, "try-error")) {
    cat("create_summary_matrix failed for ", s, ":\n", sep = "")
    cat(as.character(res_sum), "\n")
    next
  }
  
  if (!all(c("cell_names", "scATOMIC_pred") %in% colnames(res_sum))) {
    cat("Unexpected columns in scATOMIC summary for ", s, ":\n", sep = "")
    print(colnames(res_sum))
    next
  }
  
  # Map epithelial subtype prediction
  cc <- intersect(res_sum$cell_names, names(pred_scATOMIC))
  pred_scATOMIC[cc] <- as.character(res_sum$scATOMIC_pred[match(cc, res_sum$cell_names)])
  
  # Preferred malignancy label if available
  if ("pan_cancer_cluster" %in% colnames(res_sum)) {
    malignancy[cc] <- as.character(res_sum$pan_cancer_cluster[match(cc, res_sum$cell_names)])
  } else {
    # Fallback: infer from prediction label text
    labs <- pred_scATOMIC[cc]
    malignancy[cc] <- ifelse(grepl("Cancer", labs, ignore.case = TRUE), "Cancer", "Normal")
  }
  
  cat("Done: ", s, "\n", sep = "")
  cat("  Pred assigned: ", sum(!is.na(pred_scATOMIC[cc])), "/", length(cc), "\n", sep = "")
  cat("  Malig assigned: ", sum(!is.na(malignancy[cc])), "/", length(cc), "\n", sep = "")
}

# ----------------------------
# 9) Attach results back to FULL Seurat object
# ----------------------------
seurat_object1$scATOMIC_pred_epi <- NA_character_
seurat_object1$epi_malignancy    <- NA_character_

seurat_object1$scATOMIC_pred_epi[epi_cells_all] <- pred_scATOMIC[epi_cells_all]
seurat_object1$epi_malignancy[epi_cells_all]    <- malignancy[epi_cells_all]

# ----------------------------
# 10) Summaries
# ----------------------------
cat("\nFinal epi_malignancy table:\n")
print(table(seurat_object1$epi_malignancy, useNA = "ifany"))

cat("\nTop scATOMIC_pred_epi levels:\n")
print(head(sort(table(seurat_object1$scATOMIC_pred_epi, useNA = "ifany"), decreasing = TRUE), 25))

# Optional: plot if UMAP exists
# DimPlot(seurat_object1, group.by = "epi_malignancy", reduction = "umap", label = TRUE)


# =========================================================
# 11) Create main_annot3: same as main_annot2, but epithelial split
# =========================================================

# Start from main_annot2
seurat_object1$main_annot3 <- as.character(seurat_object1$main_annot2)

is_epi <- seurat_object1$main_annot2 == "Epithelial cells"

# Prefer Script 1 label if available
if ("epi_malignancy" %in% colnames(seurat_object1@meta.data)) {
  # epi_malignancy should be "Cancer"/"Normal" for epithelial cells
  seurat_object1$main_annot3[is_epi & seurat_object1$epi_malignancy == "Cancer"] <- "Epithelial (Malignant)"
  seurat_object1$main_annot3[is_epi & seurat_object1$epi_malignancy == "Normal"] <- "Epithelial (Normal)"
  
  # if any epithelial are NA, label as uncertain
  seurat_object1$main_annot3[is_epi & is.na(seurat_object1$epi_malignancy)] <- "Epithelial (Uncertain)"
  
} else if ("epi_state" %in% colnames(seurat_object1@meta.data)) {
  # Script 2 label: "Cancer_epithelial"/"Normal_epithelial"/"Uncertain_epithelial"
  seurat_object1$main_annot3[is_epi & seurat_object1$epi_state == "Cancer_epithelial"] <- "Epithelial (Malignant)"
  seurat_object1$main_annot3[is_epi & seurat_object1$epi_state == "Normal_epithelial"] <- "Epithelial (Normal)"
  seurat_object1$main_annot3[is_epi & seurat_object1$epi_state == "Uncertain_epithelial"] <- "Epithelial (Uncertain)"
  
  # if any epithelial are NA, label as uncertain
  seurat_object1$main_annot3[is_epi & is.na(seurat_object1$epi_state)] <- "Epithelial (Uncertain)"
  
} else if ("scATOMIC_epi_label" %in% colnames(seurat_object1@meta.data)) {
  # Alternative label column (used if epi_state is not present)
  seurat_object1$main_annot3[is_epi & seurat_object1$scATOMIC_epi_label == "Cancer_epithelial"] <- "Epithelial (Malignant)"
  seurat_object1$main_annot3[is_epi & seurat_object1$scATOMIC_epi_label == "Normal_epithelial"] <- "Epithelial (Normal)"
  seurat_object1$main_annot3[is_epi & is.na(seurat_object1$scATOMIC_epi_label)] <- "Epithelial (Uncertain)"
  
} else {
  stop("No epithelial malignancy labels found. Need epi_malignancy (Script 1) or epi_state/scATOMIC_epi_label (Script 2).")
}

# Make it a factor with a nice order (optional)
lvl <- c("Epithelial (Malignant)", "Epithelial (Normal)", "Epithelial (Uncertain)",
         "Myeloid cells", "Endothelial cells", "Fibroblasts", "Bcells", "Tcells")
seurat_object1$main_annot3 <- factor(seurat_object1$main_annot3, levels = intersect(lvl, unique(seurat_object1$main_annot3)))

# Check counts
cat("\nmain_annot2:\n")
print(table(seurat_object1$main_annot2))

cat("\nmain_annot3:\n")
print(table(seurat_object1$main_annot3, useNA = "ifany"))

# ----------------------------
# 12) Save updated object
# ----------------------------
out_rds <- "query_seurat_WILDTYPE_MSI_MSS_with_scATOMIC_epi_malignancy_20260302.Rds"

saveRDS(seurat_object1, file = out_rds)
cat("\nSaved updated Seurat object:\n", out_rds, "\n", sep = "")
################################################################################
