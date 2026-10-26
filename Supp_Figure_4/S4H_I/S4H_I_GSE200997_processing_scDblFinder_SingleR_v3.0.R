################################################################################
# Figure Panel : Supp S4H, S4I (step 1 of 2)
# Description  : GSE200997 scRNA-seq processing following the original publication's methods:
#                QC, cell-cycle scoring, scDblFinder doublet removal, Harmony integration and
#                SingleR annotation into six major cell types.
# Input data   : GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz
#                GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz
# Output       : GSE200997_integrated_publication_style_scDblFinder_SingleR_only.rds
#                20260226_GSE200997_tumor_MSS_only_scRNASeq_scATOMIC_UMAP_cell_types_plot.pdf
# R packages   : Seurat, dplyr, Matrix, ggplot2, clustree, DoubletFinder, scDblFinder, harmony, SingleR, SingleCellExperiment, celldex, reshape2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Run before the GSE200997 scATOMIC script (step 2).
################################################################################

###############################################################################
# scRNA-seq processing pipeline (Seurat v3.2.2)
# Follows the methods described in the original GSE200997 publication.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)      # v3.2.2
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(clustree)    # v0.4.3
  library(DoubletFinder) # v2.0.2
})

set.seed(1234)

###############################################################################
# 0) SETTINGS: metadata column names
###############################################################################
SAMPLE_COL      <- "samples"        # e.g., "Sample", "patient", "orig.ident"
COMPARTMENT_COL <- "Compartment"   # e.g., "compartment", "major_celltype"
CONDITION_COL   <- "Condition"     # e.g., "Condition" (Tumor/Normal)
MSI_COL         <- "MSI_Status"    # e.g., "MSI_Status" (MSS/MSI)

# Optional filters (only applied if columns exist)
FILTER_TUMOR_ONLY <- FALSE  # set TRUE to keep tumour samples only
FILTER_MSS_ONLY   <- FALSE  # set TRUE to keep MSS samples only

###############################################################################
# 1) LOAD DATA
###############################################################################
counts <- read.csv(
  "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz",
  row.names = 1,
  check.names = FALSE
)
counts <- as(as.matrix(counts), "dgCMatrix")

metadata <- read.csv(
  "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz",
  row.names = 1,
  check.names = FALSE
)
###############################################################################
# 2) GENE FILTER: keep genes detected in >= 3 cells (publication)
###############################################################################
genes_keep <- Matrix::rowSums(counts > 0) >= 3
counts <- counts[genes_keep, , drop = FALSE]
cat("Genes kept (>=3 cells):", nrow(counts), "\n")

###############################################################################
# 3) CREATE SEURAT OBJECT: cells with <200 detected genes filtered (publication)
###############################################################################
seu <- CreateSeuratObject(
  counts = counts,
  project = "GSE200997",
  min.cells = 0,        # already filtered genes manually
  min.features = 200    # publication
)

# Align metadata to cells
metadata <- metadata[colnames(seu), , drop = FALSE]
seu <- AddMetaData(seu, metadata = metadata)
seu
###############################################################################
# 4) OPTIONAL: Subset tumour/MSS if requested (not in provided publication text)
###############################################################################
seu <- subset(seu, subset = Condition == "Tumor" & MSI_Status == "MSS")
table(seu@meta.data$samples)

###############################################################################
# 5) QC: remove cells with >25% mitochondrial (publication)
###############################################################################
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu

# QC plots (recommended before filtering)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Apply mito cutoff from publication
seu <- subset(seu, subset = percent.mt <= 25)
cat("Cells after mito<=25%:", ncol(seu), "\n")

###############################################################################
# 6) CELL CYCLE SCORING + PHASE (publication)
###############################################################################
cc.genes <- Seurat::cc.genes.updated.2019

# Note: CellCycleScoring expects normalised data, so data are normalised here
# (and re-normalised within each sample later for integration).
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

seu <- CellCycleScoring(
  seu,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)
table(seu$Phase)

###############################################################################
# 6) DOUBLETS / MULTIPLETS REMOVAL (replacement for DoubletFinder)
#    scDblFinder labels singlet/doublet and can flag higher-order multiplets.
###############################################################################
sce <- as.SingleCellExperiment(seu)

library(scDblFinder)

# scDblFinder is run per sample (performs better with sample IDs)
#if (SAMPLE_COL %in% colnames(seu@meta.data)) {
  sce$sample_id <- seu@meta.data[[SAMPLE_COL]]
  sce <- scDblFinder(sce, samples = "sample_id")
#} else {
 # sce <- scDblFinder(sce)
#}

# Bring calls back to Seurat
seu$doublet_class <- sce$scDblFinder.class        # "singlet" / "doublet"
seu$doublet_score <- sce$scDblFinder.score
seu$doublet_type  <- if ("scDblFinder.class" %in% colnames(colData(sce))) sce$scDblFinder.class else NA

print(table(seu$doublet_class))

# Remove doublets (and any multiplets, if labelled separately)
seu <- subset(seu, subset = doublet_class == "singlet")
cat("Cells after doublet removal:", ncol(seu), "\n")

###############################################################################
# 7) DISCARD SAMPLE B-cac13 (publication)
###############################################################################
#if (SAMPLE_COL %in% colnames(seu@meta.data) && "B-cac13" %in% unique(seu@meta.data[[SAMPLE_COL]])) {
 # seu <- subset(seu, subset = get(SAMPLE_COL) != "B-cac13")
#  cat("Removed sample B-cac13. Cells now:", ncol(seu), "\n")
#} else {
 # cat("Sample B-cac13 not found (or SAMPLE_COL not set correctly). Skipping.\n")
#}

###############################################################################
# 8) INTEGRATION (Stuart et al.) WITH Harmony (publication)
###############################################################################

DefaultAssay(seu) <- "RNA"

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = VariableFeatures(seu))

###############################################################################
# 2) PCA (REQUIRED before Harmony)
###############################################################################

seu <- RunPCA(seu, npcs = 30)

###############################################################################
# 3) Harmony batch correction
###############################################################################

library(harmony)

seu <- RunHarmony(
  object = seu,
  group.by.vars = "samples"
)

Reductions(seu)

###############################################################################
# 4) Use Harmony embeddings for downstream steps
###############################################################################

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)

seu <- FindClusters(
  seu,
  resolution = c(0.2, 0.4, 0.6, 0.8, 1.0),
  algorithm = 1  # Louvain
)

#Cluster Stability Using Clustree
library(clustree)

clustree(seu@meta.data, prefix = "RNA_snn_res.")
Idents(seu) <- "RNA_snn_res.0.4"


seu <- RunTSNE(seu, reduction = "harmony", dims = 1:30)

DimPlot(seu, reduction = "tsne", group.by = "samples")

seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)
seu

Reductions(seu)
DimPlot(seu, reduction = "umap", group.by = "samples")

# Setup (choose clustering identity + define 6 marker sets)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

DefaultAssay(seu) <- "RNA"

###############################################################################
# SingleR
###############################################################################
#SingleR annotation (orthogonal evidence)
suppressPackageStartupMessages({
  library(SingleR)
  library(SingleCellExperiment)
  library(celldex)
})

# Convert to SCE for SingleR
sce <- as.SingleCellExperiment(seu)

# Reference (pick one; HumanPrimaryCellAtlas is broad, Monaco is immune-focused)
ref <- HumanPrimaryCellAtlasData()

# Use log-normalized data if present; otherwise SingleR will use assay data
pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

# Add SingleR labels to Seurat
seu$SingleR_label_main <- pred$labels
seu$SingleR_pruned     <- pred$pruned.labels

table(seu$SingleR_label_main)

DimPlot(seu, reduction = "umap", group.by = "SingleR_label_main")


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

seu$SingleR_6 <- vapply(seu$SingleR_label_main, map_singleR_to_6, character(1))
table(seu$SingleR_6, useNA = "ifany")

DimPlot(seu, reduction = "umap", group.by = "SingleR_6")

#myeloid cells
FeaturePlot(
  seu,
  features = c("LYZ", "CD68", "CSF1R"),
  reduction = "umap",
  cols = c("lightgrey", "red")
)



#Monocytes vs Macrophages vs Dendritic Cells

FeaturePlot(
  seu,
  features = c("LYZ", "FCGR3A", "CD14", "ITGAX", "CLEC9A"),
  reduction = "umap",
  cols = c("lightgrey", "blue")
  
)

#Tumour-Associated Macrophage (TAM) Markers (CRC-relevant)
FeaturePlot(
  seu,
  features = c("CD68", "MRC1", "CD163", "MSR1"),
  reduction = "umap"
)

#Pan-T Cell Markers
FeaturePlot(
  seu,
  features = c("CD3D", "CD3E", "TRAC"),
  reduction = "umap",
  cols = c("grey90", "darkred"),
  pt.size = 0.3
)

#CD4 vs CD8 T Cells
FeaturePlot(
  seu,
  features = c("CD4", "CD8A", "CD8B"),
  reduction = "umap",
  cols = c("grey90", "blue"),
  pt.size = 0.3
)

#Exhausted / Activated T Cells (Important in CRC)
FeaturePlot(
  seu,
  features = c("PDCD1", "LAG3", "TIGIT", "GZMB"),
  reduction = "umap",
  cols = c("grey90", "purple")
)

#Pan-B Cell Markers
FeaturePlot(
  seu,
  features = c("MS4A1", "CD79A", "CD79B"),
  reduction = "umap",
  cols = c("grey90", "darkgreen"),
  pt.size = 0.3
)

#Plasma Cells (Often present in CRC)
FeaturePlot(
  seu,
  features = c("MZB1", "SDC1", "XBP1"),
  reduction = "umap"
)


#Core fibroblast markers
FeaturePlot(
  seu,
  features = c("COL1A1", "COL1A2", "DCN", "LUM"),
  reduction = "umap",
  cols = c("grey90", "brown"),
  pt.size = 0.3,
  order = TRUE
)

#Activated CAF markers (CRC-relevant)
FeaturePlot(
  seu,
  features = c("ACTA2", "FAP", "PDPN"),
  reduction = "umap",
  cols = c("grey90", "darkorange")
)

#Endothelial Cells
FeaturePlot(
  seu,
  features = c("PECAM1", "VWF", "KDR", "CLDN5"),
  reduction = "umap",
  cols = c("grey90", "darkblue"),
  pt.size = 0.3,
  order = TRUE
)

#Pan-epithelial markers
FeaturePlot(
  seu,
  features = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  reduction = "umap",
  cols = c("grey90", "red"),
  pt.size = 0.3,
  order = TRUE
)

#CRC-specific epithelial markers
FeaturePlot(
  seu,
  features = c("MUC1", "MUC13", "CEACAM5"),
  reduction = "umap"
)


#Epithelial markers
epi_genes <- c("EPCAM","KRT8","KRT18","KRT19","MUC1")

#Malignant-associated genes (CRC-relevant)
malig_genes <- c(
  "MKI67","TOP2A","PCNA","MYC",
  "CCNB1","CDK1",
  "AXIN2","LGR5",
  "SOX9"
)

#Normal epithelial differentiation genes
normal_epi_genes <- c(
  "KRT20","FABP1","CA1","GUCA2A","GUCA2B"
)

#Add Module Scores
seu <- AddModuleScore(
  seu,
  features = list(epi_genes),
  name = "EpiScore"
)

seu <- AddModuleScore(
  seu,
  features = list(malig_genes),
  name = "MaligScore"
)

seu <- AddModuleScore(
  seu,
  features = list(normal_epi_genes),
  name = "NormalEpiScore"
)

#Subset Epithelial Cells First
#epi_cells <- WhichCells(seu, expression = EpiScore1 > 0)

epi_cells <- WhichCells(seu, expression = SingleR_6 == "Epithelial")

#Classify Malignant vs Normal
#We classify based on relative score dominance:
 # seu$EpiType <- "NonEpithelial"

#seu$EpiType[epi_cells] <- ifelse(
 # seu$MaligScore1[epi_cells] > seu$NormalEpiScore1[epi_cells],
#  "Malignant",
 # "Normal"
#)

#Visualise
#DimPlot(
 # seu,
#  reduction = "umap",
#  group.by = "EpiType",
 # cols = c("grey80","blue","red")
#) + theme_classic()


#More Rigorous Alternative for malignant vs normal 
#Instead of thresholding:
df <- seu@meta.data[epi_cells, c("MaligScore1","NormalEpiScore1")]
df <- df[complete.cases(df), ]  # safer than na.omit

km <- kmeans(df, centers = 2, nstart = 50)

# identify malignant cluster
malig_cluster <- which.max(
  tapply(df$MaligScore1, km$cluster, mean)
)

# reset first (important)
seu$EpiType <- as.character(seu$SingleR_6)  # initialise from SingleR labels

# assign using cell names, NOT epi_cells
seu$EpiType[rownames(df)] <- ifelse(
  km$cluster == malig_cluster,
  "Malignant",
  "Normal"
)

table(seu$EpiType)

#Visualise
DimPlot(
 seu,
  reduction = "umap",
  group.by = "EpiType"
) + theme_classic()


library(ggplot2)

df$Cluster <- factor(km$cluster)

ggplot(df, aes(x = NormalEpiScore1, y = MaligScore1, color = Cluster)) +
  geom_point(size = 0.6, alpha = 0.6) +
  theme_classic() +
  labs(
    x = "Normal Epithelial Score",
    y = "Malignant Score",
    title = "K-means Separation of Epithelial Cells"
  )

#UMAP Visualisation (Most Important)
#pdf("20260226_GSE200997_tumor_MSS_only_scRNASeq_scATOMIC_UMAP_cell_types_plot.pdf")
DimPlot(
  seu,
  group.by = "EpiType",
  reduction = "umap",
  pt.size = 0.3
) +
  theme_classic() +
  ggtitle("Malignant vs Normal Epithelial Cells")
#dev.off()

#Density Plot (Score Distribution)
library(reshape2)

df_long <- melt(
  cbind(df, EpiType = seu$EpiType[epi_cells]),
  id.vars = "EpiType"
)

#ggplot(df_long, aes(value, fill = EpiType)) +
#  geom_density(alpha = 0.4) +
 # facet_wrap(~variable, scales = "free") +
  #theme_classic()

#FeaturePlot Overlay (Sanity Check)
FeaturePlot(
  seu,
  features = c("MKI67","KRT20"),
  reduction = "umap",
  order = TRUE
)

#Important Sanity Check
#Verify cluster assignment logic:
  aggregate(df$MaligScore1, by=list(km$cluster), mean)
aggregate(df$NormalEpiScore1, by=list(km$cluster), mean)

#Cleaner Separation
FeaturePlot(seu, features = c("MaligScore1","NormalEpiScore1"))



#DotPlot by Final Cell Type
DefaultAssay(seu) <- "RNA"

DotPlot(
  seu,
  features = c( "EGFR","PTK2"),
  group.by = "EpiType"
) +
  RotatedAxis() +
  theme_classic() +
  labs(title = "PTK2 and EGFR Expression Across Major Cell Types")

###############################################################################
# 9) SAVE
###############################################################################
saveRDS(seu, file = "GSE200997_integrated_publication_style_scDblFinder_SingleR_only.rds")
cat("Saved integrated object.\n")


