#!/usr/bin/env Rscript
################################################################################
# Figure Panel : Supp S4J, S4K (step 1 of 2)
# Description  : GSE178318 (CRC liver metastasis) scRNA-seq QC following the original paper:
#                nGene >= 500, per-batch +/-3 SD filters on nUMI and nGene, mitochondrial <= 15%,
#                normalisation, clustering and SingleR annotation.
# Input data   : GSE178318_matrix.mtx.gz
#                GSE178318_genes.tsv.gz
#                GSE178318_barcodes.tsv.gz
# Output       : GSE178318_Seurat_QC_filtered.rds
#                GSE178318_QC_Seurat_Modern_CRC_LM_3_3_samples_SingleR_only.rds
# R packages   : Matrix, Seurat, biomaRt, dplyr, tidyr, ggplot2, SingleR, SingleCellExperiment, celldex
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Place the GSE178318 matrix, genes and barcodes files in data/.
################################################################################

# ============================================================
# GSE178318 (Che et al., Cell 2021) QC pipeline (Seurat v2.3.4)
# Input: GSE178318_matrix.mtx.gz, GSE178318_genes.tsv.gz, GSE178318_barcodes.tsv.gz
#
# Filters (as described in paper):
#  (1) nGene >= 500
#  (2) per batch: log10(nUMI) within mean ± 3 SD
#  (3) per batch: nGene within mean ± 3 SD
#  (4) percent.mt <= 15
#
# Optional: keep protein-coding genes only (biomaRt)
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)   # v2.3.4
})

# -----------------------------
# 0) Paths
# -----------------------------
data_dir <- "data"
mtx_gz   <- file.path(data_dir, "GSE178318_matrix.mtx.gz")
genes_gz <- file.path(data_dir, "GSE178318_genes.tsv.gz")
bc_gz    <- file.path(data_dir, "GSE178318_barcodes.tsv.gz")

out_rds  <- file.path(data_dir, "GSE178318_Seurat_QC_filtered.rds")

# -----------------------------
# 1) Read MTX + features + barcodes
# -----------------------------
cat("Reading matrix...\n")
counts <- readMM(gzfile(mtx_gz))

cat("Reading genes/features...\n")
genes_df <- read.table(gzfile(genes_gz), header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")
# Typically 10x v2 format: V1=gene_id, V2=gene_name
# Column 1 = gene ID, column 2 = gene name
gene_ids   <- genes_df[[1]]
gene_names <- genes_df[[2]]

cat("Reading barcodes...\n")
barcodes <- read.table(gzfile(bc_gz), header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")[[1]]

# Apply dimnames
rownames(counts) <- make.unique(gene_names)
colnames(counts) <- barcodes

#rm(genes_df, gene_ids, gene_names, barcodes)
#gc()

# -----------------------------
# 2) Create Seurat object (v2)
# -----------------------------

# Keep only cells ending with "_CRC"
crc_cells <- grepl("_CRC$", colnames(counts))

# Keep only cells ending with "_LM"
LM_cells <- grepl("_LM$", colnames(counts))

cat("Total cells:", ncol(counts), "\n")
cat("CRC cells:", sum(crc_cells), "\n")

# Subset matrix
counts_crc1 <- counts[, crc_cells]
counts_crc2 <- counts[, LM_cells]
counts_crc3 <- cbind(counts_crc1,counts_crc2)

# Keep only cells treatment naive samples - COL07, COL12, COL16
col07 <- grepl("_COL07", colnames(counts_crc3))
col12 <- grepl("_COL12", colnames(counts_crc3))
col16 <- grepl("_COL16", colnames(counts_crc3))

counts_crc4 <- counts_crc3[,col07]
counts_crc5 <- counts_crc3[,col12]
counts_crc6 <- counts_crc3[,col16]

counts_crc <- cbind(counts_crc4, counts_crc5, counts_crc6)
dim(counts_crc)
  
library(Matrix)

counts_crc <- as(counts_crc, "dgCMatrix")
gc()

# Optional: free memory
#rm(counts)
#gc()

cat("Creating Seurat object...\n")
so <- CreateSeuratObject(
  counts = counts_crc,
  min.cells = 0,
  min.features = 0
)
#rm(counts)
#gc()

# -----------------------------
# 4. ADD MITOCHONDRIAL %
# -----------------------------
so[["percent.mt"]] <- PercentageFeatureSet(
  so,
  pattern = "^MT-"
)

# -----------------------------
# 5. STEP 1 — Remove <500 genes
# -----------------------------
n_before <- ncol(so)
so <- subset(so, subset = nFeature_RNA >= 500)
cat("Step 1 removed:", n_before - ncol(so), "cells\n")

# -----------------------------
# 6. DEFINE BATCH (CRITICAL)
# -----------------------------
# Barcode format: AAACCTGAGAAACCTA_COL07_CRC
so$batch <- sub(".*_(COL[0-9]+)_.*", "\\1", colnames(so))

cat("Batch distribution:\n")
print(table(so$batch))

# -----------------------------
# 7. STEP 2 & 3 — Per-batch ±3 SD
# -----------------------------
meta <- so@meta.data
keep_cells <- rep(TRUE, nrow(meta))

for (b in unique(meta$batch)) {
  idx <- which(meta$batch == b)
  
  # Step 2: log10(nUMI)
  logUMI <- log10(meta$nCount_RNA[idx] + 1)
  keep_log <- logUMI >= (mean(logUMI) - 3*sd(logUMI)) &
    logUMI <= (mean(logUMI) + 3*sd(logUMI))
  
  # Step 3: nFeature
  gene <- meta$nFeature_RNA[idx]
  keep_gene <- gene >= (mean(gene) - 3*sd(gene)) &
    gene <= (mean(gene) + 3*sd(gene))
  
  keep_cells[idx] <- keep_log & keep_gene
}

n_before23 <- ncol(so)
so <- so[, keep_cells]
cat("Step 2–3 removed:", n_before23 - ncol(so), "cells\n")

# -----------------------------
# 8. STEP 4 — Remove >15% mito
# -----------------------------
n_before4 <- ncol(so)
so <- subset(so, subset = percent.mt <= 15)
cat("Step 4 removed:", n_before4 - ncol(so), "cells\n")

# -----------------------------
# 8b. FEATURE QC — Remove low-quality genes
# -----------------------------
# The paper reports "protein-coding genes" remaining after QC.
# Before protein-coding filtering, genes that are not detected
# (or detected in very few cells) are dropped.

MIN_CELLS_PER_GENE <- 3   # minimum number of cells expressing a gene

counts_mat <- GetAssayData(so, slot = "counts")
gene_detected_cells <- Matrix::rowSums(counts_mat > 0)

n_genes_before <- nrow(so)
keep_genes_basic <- gene_detected_cells >= MIN_CELLS_PER_GENE

so <- subset(so, features = names(keep_genes_basic)[keep_genes_basic])

cat("Feature QC removed:", n_genes_before - nrow(so),
    "genes (kept genes detected in >=", MIN_CELLS_PER_GENE, "cells)\n")

#rm(counts_mat, gene_detected_cells, keep_genes_basic)
gc()

# -----------------------------
# 9. OPTIONAL — Keep Protein-Coding Genes Only
# -----------------------------
DO_PROTEIN_CODING <- TRUE  # Set TRUE if desired

if (DO_PROTEIN_CODING) {
  suppressPackageStartupMessages(library(biomaRt))
  
  cat("Filtering protein-coding genes...\n")
  
  gene_symbols <- rownames(so)
  
  mart <- useEnsembl(biomart = "genes",
                     dataset = "hsapiens_gene_ensembl")
  
  ann <- getBM(
    attributes = c("hgnc_symbol", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = mart
  )
  
  protein_coding <- unique(
    ann$hgnc_symbol[ann$gene_biotype == "protein_coding"]
  )
  
  so <- subset(so, features = protein_coding)
}

# -----------------------------
# 10. FINAL SUMMARY
# -----------------------------
cat("\nFinal dimensions:\n")
cat("Genes:", nrow(so), "\n")
cat("Cells:", ncol(so), "\n")


###############################################################################
# 11) NORMALISATION + DIM REDUCTION + CLUSTERING (needed for markers/UMAP)
###############################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(ggplot2)
})

DefaultAssay(so) <- "RNA"

# ---- Normalise
so <- NormalizeData(
  so,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ---- HVGs (use standard 2000)
so <- FindVariableFeatures(
  so,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# ---- Scale (include percent.mt if present)
vars_to_regress <- intersect(c("percent.mt"), colnames(so@meta.data))
so <- ScaleData(
  so,
  features = VariableFeatures(so),
#  vars.to.regress = vars_to_regress,
  verbose = FALSE
)

# ---- PCA
so <- RunPCA(
  so,
  features = VariableFeatures(so),
  npcs = 50,
  verbose = FALSE
)

# ---- Neighbours + clustering
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
so <- FindClusters(so, resolution = 0.6, verbose = FALSE)

# ---- UMAP + tSNE (both used in plots)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)
so <- RunTSNE(so, dims = 1:30, verbose = FALSE)


###############################################################################
# 12) SingleR
###############################################################################
#SingleR annotation (orthogonal evidence)
suppressPackageStartupMessages({
  library(SingleR)
  library(SingleCellExperiment)
  library(celldex)
})

# Convert to SCE for SingleR
sce <- as.SingleCellExperiment(so)

# Reference (pick one; HumanPrimaryCellAtlas is broad, Monaco is immune-focused)
ref <- HumanPrimaryCellAtlasData()

# Use log-normalized data if present; otherwise SingleR will use assay data
pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

# Add SingleR labels to Seurat
so$SingleR_label_main <- pred$labels
so$SingleR_pruned     <- pred$pruned.labels

table(so$SingleR_label_main)

DimPlot(so, reduction = "umap", group.by = "SingleR_label_main")


#Final forced assignment into the 6 major types (cell-level)

# Map SingleR broad labels into six major cell classes
map_singleR_to_6 <- function(x) {
  # if (is.na(x)) return("Uncertain")
  if (grepl("T cells|T_cell|T-cell|NK_cell", x, ignore.case = TRUE)) return("T_cells")
  if (grepl("B cells|B_cell|B-cell|Plasma| Pre-B_cell_CD34-|Pro-B_cell_CD34+", x, ignore.case = TRUE)) return("B_cells")
  if (grepl("Monocyte|Myelocyte|Macrophage|DC|Neutrophils|Pro-Myelocyte|CMP|GMP|MEP|Erythroblast|BM & Prog.|HSC_-G-CSF|HSC_CD34+", x, ignore.case = TRUE)) return("Myeloid")
  if (grepl("Endothelial|Smooth_muscle_cells", x, ignore.case = TRUE)) return("Endothelial")
  if (grepl("Fibroblast|Stromal|MSC|Osteoblasts|Chondrocytes", x, ignore.case = TRUE)) return("Fibroblast")
  if (grepl("Epithelial|Tissue_stem_cells|Keratinocytes|Neuroepithelial_cell|Stem|Hepatocytes|iPS_cells|Neurons|Astrocyte|Embryonic_stem_cells|Gametocytes", x, ignore.case = TRUE)) return("Epithelial")
  return("Uncertain")
}


so$SingleR_6 <- vapply(so$SingleR_label_main, map_singleR_to_6, character(1))
table(so$SingleR_6, useNA = "ifany")

DimPlot(so, reduction = "umap", group.by = "SingleR_6")


###############################################################################
# 13) MALIGNANT vs NORMAL WITHIN EPITHELIAL (kmeans on module scores)
###############################################################################
epi_genes <- c("EPCAM","KRT8","KRT18","KRT19","MUC1")
malig_genes <- c("MKI67","TOP2A","PCNA","MYC","CCNB1","CDK1","AXIN2","LGR5","SOX9")
normal_epi_genes <- c("KRT20","FABP1","CA1","GUCA2A","GUCA2B")

epi_genes        <- intersect(epi_genes, rownames(so))
malig_genes      <- intersect(malig_genes, rownames(so))
normal_epi_genes <- intersect(normal_epi_genes, rownames(so))

so <- AddModuleScore(so, features = list(epi_genes),        name = "EpiScore")
so <- AddModuleScore(so, features = list(malig_genes),      name = "MaligScore")
so <- AddModuleScore(so, features = list(normal_epi_genes), name = "NormalEpiScore")

epi_cells <- WhichCells(so, expression = celltype_final6 == "Epithelial")

# initialise EpiType as character to avoid factor-level NA problems
so$EpiType <- as.character(so$celltype_final6)

df <- so@meta.data[epi_cells, c("MaligScore1","NormalEpiScore1")]
df <- df[complete.cases(df), , drop = FALSE]

if (nrow(df) >= 10) {
  set.seed(1234)
  km <- kmeans(df, centers = 2, nstart = 50)
  malig_cluster <- which.max(tapply(df$MaligScore1, km$cluster, mean))
  
  # assign using rownames(df) (cell names)
  so$EpiType[rownames(df)] <- ifelse(km$cluster == malig_cluster, "Malignant", "Normal")
} else {
  warning("Too few epithelial cells with complete scores for kmeans; leaving EpiType as celltype_final6.")
}

so$EpiType <- factor(so$EpiType)

table(so$EpiType, useNA = "ifany")



# Save
saveRDS(so, file = "GSE178318_QC_Seurat_Modern_CRC_LM_3_3_samples_SingleR_only.rds")

cat("\nQC pipeline complete.\n")
