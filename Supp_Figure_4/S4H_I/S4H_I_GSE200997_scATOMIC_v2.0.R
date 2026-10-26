################################################################################
# Figure Panel : Supp S4H, S4I (step 2 of 2)
# Description  : scATOMIC classification of epithelial cells into cancer vs normal (per-sample
#                run) in GSE200997; UMAP of cell types and EGFR/PTK2 dot plot.
# Input data   : GSE200997_integrated_publication_style_scDblFinder_SingleR_only.rds
# Output       : 20260226_GSE200997_tumor_MSS_only_scRNASeq_scATOMIC_EGFR_PTK2_dotplot.pdf
#                GSE200997_integrated_publication_style_scDblFinder_scATOMIC_epi_malignancy.Rds
#                20260302_GSE200997_tumor_MSS_only_scRNASeq_scATOMIC_UMAP_cell_types_plot.pdf
# R packages   : reticulate, Seurat, Matrix, scATOMIC, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Uses the RDS written by the processing script. Requires a Python
#                environment where `import magic` works; set RETICULATE_PYTHON.
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
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# Set to the Python binary of an environment where `import magic` works
Sys.setenv(RETICULATE_PYTHON = 
             "path/to/python")

library(reticulate)
py_config()
stopifnot(reticulate::py_module_available("magic"))
cat("Python module 'magic' is available.\n")

# ----------------------------
# 1) Read Seurat object
# ----------------------------
seurat_object1 <- readRDS(
  "GSE200997_integrated_publication_style_scDblFinder_SingleR_only.rds"
)

cat("\nCell-type distribution (SingleR_6):\n")
print(table(seurat_object1$SingleR_6))

# ----------------------------
# 2) Subset epithelial cells
# ----------------------------
epi <- subset(seurat_object1, subset = SingleR_6 == "Epithelial")
DefaultAssay(epi) <- "RNA"

if (!"RNA" %in% names(epi@assays)) stop("RNA assay not found in epithelial subset.")
if (ncol(epi) == 0) stop("No epithelial cells found (main_annot2 == 'Epithelial cells').")

cat("\nEpithelial cells:", ncol(epi), "\n")

# ----------------------------
# 3) Choose sample column for per-sample run
# ----------------------------
split_col <- "samples"
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

get_counts_safe <- function(seu, assay = "RNA") {
  DefaultAssay(seu) <- assay
  
  # Seurat v5 layers
  x <- tryCatch(LayerData(seu, assay = assay, layer = "counts"),
                error = function(e) NULL)
  
  # fallback: older slots / Seurat compat
  if (is.null(x)) {
    x <- tryCatch(GetAssayData(seu, assay = assay, slot = "counts"),
                  error = function(e) NULL)
  }
  
  # if still null, fallback to data (not ideal but prevents crash)
  if (is.null(x)) {
    x <- tryCatch(GetAssayData(seu, assay = assay, slot = "data"),
                  error = function(e) NULL)
  }
  
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
  cat("Cells in epi_s: ", ncol(epi_s), "\n", sep = "")
  cat("Assays: ", paste(names(epi_s@assays), collapse=", "), "\n", sep = "")
  cat("Layers in RNA: ",
      paste(tryCatch(names(epi_s@assays$RNA@layers), error=function(e) "NA"), collapse=", "),
      "\n", sep = "")
  
 # x <- LayerData(epi_s, assay = "RNA", layer = "counts")
  x <- get_counts_safe(epi_s, assay = "RNA")
  
  # Guardrails: must be a 2D sparse matrix
  if (is.null(x)) {
    cat("Skipping: could not retrieve counts/data for sample ", s, "\n", sep="")
    next
  }
  if (!inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
  if (length(dim(x)) != 2) {
    cat("Skipping: counts is not 2D for sample ", s, "\n", sep="")
    next
  }
  if (ncol(x) < 2 || nrow(x) < 2) {
    cat("Skipping: too small matrix for sample ", s, "\n", sep="")
    next
  }
  
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
seurat_object1$main_annot3 <- as.character(seurat_object1$SingleR_6)

is_epi <- seurat_object1$SingleR_6 == "Epithelial"

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
         "Myeloid", "Endothelial", "Fibroblast", "B_cells", "T_cells")
seurat_object1$main_annot3 <- factor(seurat_object1$main_annot3, levels = intersect(lvl, unique(seurat_object1$main_annot3)))

# Check counts
cat("\nmain_annot2:\n")
print(table(seurat_object1$SingleR_6))

cat("\nmain_annot3:\n")
print(table(seurat_object1$main_annot3, useNA = "ifany"))

DefaultAssay(seurat_object1) <- "RNA"

# Convert to character first (important)
seurat_object1$main_annot3 <- as.character(seurat_object1$main_annot3)

# Replace NA with "Uncertain"
seurat_object1$main_annot3[is.na(seurat_object1$main_annot3)] <- "Uncertain"

# Convert back to factor and drop empty levels
seurat_object1$main_annot3 <- droplevels(as.factor(seurat_object1$main_annot3))

table(seurat_object1$main_annot3, useNA = "ifany")

##############
# Figures for the paper
#################

DotPlot(
  seurat_object1,
  features = c( "EGFR","PTK2"),
  group.by = "main_annot3"
) +
  RotatedAxis() +
  theme_classic() +
  labs(title = "PTK2 and EGFR Expression Across Major Cell Types")

library(Seurat)
library(ggplot2)

DotPlot(
  seurat_object1,
  features = c("EGFR", "PTK2"),
  group.by = "main_annot3"
) +
  scale_colour_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  RotatedAxis() +
  theme_classic() +
  labs(title = "PTK2 and EGFR Expression Across Major Cell Types")


library(Seurat)
library(ggplot2)

# Ensure correct ordering
seurat_object1$main_annot3 <- factor(
  seurat_object1$main_annot3,
  levels = c(
    "Epithelial (Malignant)",
    "Epithelial (Normal)",
    "Epithelial (Uncertain)",
    "Myeloid",
    "Endothelial",
    "Fibroblast",
    "B_cells",
    "T_cells"
  )
)

p <- DotPlot(
  seurat_object1,
  features = c("EGFR", "PTK2"),
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2, 2)
  ) +
  scale_size(range = c(4, 12)) +
  RotatedAxis() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "PTK2 and EGFR Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )

pdf("20260226_GSE200997_tumor_MSS_only_scRNASeq_scATOMIC_EGFR_PTK2_dotplot.pdf")

p

dev.off()

#pdf("20260302_GSE200997_tumor_MSS_only_scRNASeq_scATOMIC_UMAP_cell_types_plot.pdf")
DimPlot(
  seurat_object1,
  group.by = "main_annot3",
  reduction = "umap",
  pt.size = 0.3
) +
  theme_classic() +
  ggtitle("Malignant vs Normal Epithelial Cells")
#dev.off()

# ----------------------------
# 12) Save updated object
# ----------------------------
out_rds <- "GSE200997_integrated_publication_style_scDblFinder_scATOMIC_epi_malignancy.Rds"

saveRDS(seurat_object1, file = out_rds)
cat("\nSaved updated Seurat object:\n", out_rds, "\n", sep = "")
################################################################################
