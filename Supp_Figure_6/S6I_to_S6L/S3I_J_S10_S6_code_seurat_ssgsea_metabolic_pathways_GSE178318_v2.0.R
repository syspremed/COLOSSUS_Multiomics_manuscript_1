################################################################################
# Figure Panels : Supplementary Figures 5I–5L (also S5R)
# Description  : Seurat scRNA-seq analysis of public dataset GSE178318 (v2.0).
#                Geyser enrichment and cell-type expression plots for KEGG
#                alanine/aspartate/glutamate metabolism, GLS/GLUL/GOT2,
#                TCA genes, and EIF4B.
# Input data   : GSE178318 scRNA-seq Seurat object
#                c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt
# Output       : S5I_..._geyserEnrichment_GSE178318.pdf
#                S5J_..._GLS_GLUL_GOT2_GSE178318.pdf
#                S5K_..._tca_cancer_expressed_genes_GSE178318.pdf
#                S5L_..._ala_asp_GSE178318.pdf
#                S5R_..._EIF4B_GSE178318.pdf
# R packages   : escape, Seurat, SeuratObject, SingleCellExperiment, scran,
#                qusage, RColorBrewer, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026 (v2.0)
################################################################################



#ESCAPE R package to do SSGSEA

# escape: single-cell gene set scoring
# Seurat/SeuratObject: scRNA-seq/snRNA-seq processing
# scran: normalisation; qusage: read GMT files
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))



# ── Step 1: Load Seurat object ───────────────────────────────────────────────
# Set this path to the relevant dataset Seurat object.
# For S5A-D: COLOSSUS snRNA-seq (RAS-WT, 4 samples)
# For S5E-H: GSE200997 scRNA-seq
# For S5I-L: GSE178318 scRNA-seq
seurat_object1<-readRDS("GSE178318_MSS_CRC_LM_3_3_integrated_publication_style_scDblFinder_scATOMIC_epi_malignancy.Rds") #seurat object

# ── Step 2: Load KEGG gene sets ──────────────────────────────────────────────
# GMT file from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/)
#KEGG legacy gene sets
library(qusage)

kegg <- read.gmt("~/Downloads/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")

# reading scaled loadings
#mv_loadings <- read.delim("MV1_marker.txt", stringsAsFactors = FALSE)
#dim(mv_loadings)

m1 <- kegg$KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM
m2 <- kegg$KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM

#create gene list
# Define gene set list for escape scoring.
# MV1Positive = glyoxylate/dicarboxylate pathway genes
# MV1Negative = alanine/aspartate/glutamate pathway genes
gene.sets <- list(MV1Positive = m1, MV1Negative = m2)
gene.sets

#Anguraj Script
#enrichment.scores <- escape.matrix(seurat_object1, 
                                 #  gene.sets = gene.sets, 
                                  # groups = 1000, 
                                  # min.size = 5)
Idents(seurat_object1) <- "main_annot3"
ES <- runEscape(seurat_object1, 
                method = "ssGSEA",
                gene.sets = gene.sets, 
                groups = 1000, 
                min.size = 5,
                new.assay.name = "escape.ssGSEA")

Idents(ES) <- "main_annot3"



#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

FeaturePlot(ES, "MV1Positive") + 
  scale_color_gradientn(colors = colorblind_vector) + 
  theme(plot.title = element_blank(), keep.scale = "all")

FeaturePlot(ES, "MV1Negative") + 
  scale_color_gradientn(colors = colorblind_vector) + 
  theme(plot.title = element_blank(), keep.scale = "all")


pdf("S5I_20260303_geyserEnrichment_scRNAseq_Epi_split_KEGG_ALANINE_ASPARTATE_GLUTAMATE_METABOLISM_genes_GSE178318.pdf")
geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Negative")
dev.off()



# Featureplot for GLS

metab_genes <- c("GLS", "GLUL", "GOT2")


DotPlot(object = seurat_object1, features = metab_genes)


#pdf("S5J_20260303_geyserEnrichment_scRNAseq_Epi_split_GLS_GLUL_GOT2_genes_GSE178318.pdf")

DotPlot(
  seurat_object1,
  features = metab_genes,
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2.5, 2.5)
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
#dev.off()


tca_cancer_expressed_genes <- c(
  "CS", 
  "ACO2", 
  "OGDH", 
  "SUCLG1", 
  "SUCLG2", 
  "SUCLA2", 
  "SDHA", 
  "SDHB", 
  "FH"
  
)


DotPlot(object = seurat_object1, features = tca_cancer_expressed_genes)

pdf("S5K_20260303_geyserEnrichment_scRNAseq_Epi_split_tca_cancer_expressed_genes_genes_GSE178318.pdf")

DotPlot(
  seurat_object1,
  features = tca_cancer_expressed_genes,
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2.5, 2.5)
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
dev.off()


# Subset only epithelial cells
Idents(seurat_object1) <- "main_annot3"
epithelial_cells <- subset(seurat_object1, idents = "Epithelial (Malignant)")

library(patchwork)
# Step 1: Scale data for epithelial cells
#epithelial_cells <- ScaleData(epithelial_cells, features = tca_cancer_expressed_genes)

# Compute average expression
avg_expr <- AverageExpression(epithelial_cells, features = tca_cancer_expressed_genes)

# Extract as a named numeric vector
avg_vector <- as.numeric(avg_expr$RNA)
names(avg_vector) <- rownames(avg_expr$RNA)

# Create barplot
barplot(avg_vector,
        main = "Average Expression of TCA Genes in Epithelial Cells",
        las = 2,                # rotate x-axis labels
        col = "steelblue",      # color
        cex.names = 0.8,        # shrink label size
        ylab = "Average Expression")

# DotPlot with custom blue-to-red gradient
DotPlot(epithelial_cells, features = tca_cancer_expressed_genes) +
  RotatedAxis() +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  ggtitle("TCA Gene Expression in Epithelial Cells") +
  theme_minimal()

ala_asp <- c(
"GPT",
"ASL",
"ASPA",
"ASPG"
  
)


DotPlot(object = seurat_object1, features = ala_asp)


pdf("S5L_20260303_geyserEnrichment_scRNAseq_Epi_split_ala_asp_GSE178318.pdf")

DotPlot(
  seurat_object1,
  features = ala_asp,
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2.5, 2.5)
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
    title = "Alanine and Aspartate Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )

dev.off()

#Figure 2M
dp <- DotPlot(object = seurat_object1, features = c("EIF4B","CCDC107","MTRF1L","CDC26","SUMO2","RPS17","RPL17"))

# Apply blue-to-red gradient
dp + scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#pdf("S5Q_20260303_geyserEnrichment_scRNAseq_Epi_split_EIF4B_GSE178318.pdf")

DotPlot(
  seurat_object1,
  features = "EIF4B",
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2.5, 2.5)
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
    title = "EIF4B Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )
#dev.off()

crc_cells_1 <- grep("_CRC$", colnames(seurat_object1), value = TRUE)

seu_crc <- subset(seurat_object1, cells = crc_cells_1)

cat("CRC cells:", length(crc_cells_1), "\n")

DotPlot(
  seu_crc,
  features = "EIF4B",
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2.5, 2.5)
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
    title = "EIF4B Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )

lm_cells <- grep("_LM$", colnames(seurat_object1), value = TRUE)

seu_lm <- subset(seurat_object1, cells = lm_cells)

cat("LM cells:", length(lm_cells), "\n")

DotPlot(
  seu_lm,
  features = "EIF4B",
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2.5, 2.5)
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
    title = "EIF4B Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )


dp_data <- dp$data

gene <- "RPS17"
cell <- "Epithelial cells"

# Filter for specific cell types
filtered_data <- dp_data %>%
  filter(id %in% cell, features.plot %in% gene) %>%
  select(Gene = features.plot, CellType = id, PercentExpression = pct.exp)

# View the filtered table
print(filtered_data)

library(Seurat)
library(ggplot2)

p <- DotPlot(
  object = seurat_object1,
  features = "EIF4B",
  group.by = "main_annot3",
  scale = TRUE
)

p <- p +
  scale_colour_gradient2(
    low = "#deebf7",     # very light blue
    mid = "white",
    high = "#08519c",    # dark blue
    midpoint = 0,
    limits = c(-1.5, 1.0),
    name = "Average Expression"
  ) +
  scale_size(
    range = c(2, 10),
    breaks = c(32, 36, 40),
    limits = c(32, 100),
    name = "% Expressed"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(title = "EIF4B Expression Across Cell Types")

print(p)

##Figure 2M
#dp <- DotPlot(object = seurat_object1, features = c("CCDC107", "MTRF1L", "EIF4B","CDC26", "SUMO2"))
dp <- DotPlot(object = seurat_object1, features = c("EIF4B"))


# Apply blue-to-red gradient
dp + scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("20260102_COLOSSUS_EIF5B_CRIPSR_selected_snRNAseq_dotplot.pdf")

# Apply blue-to-red gradient
dp + scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

pdf("20250610_COLOSSUS_EGFR_PTK2_CRIPSR_selected_figure1_dotplot.pdf")

##Figure 1O
#dp <- DotPlot(object = seurat_object1, features = c("EGFR","PTK2"))

# Apply blue-to-red gradient
#dp + scale_color_gradientn(colors = c("blue", "white", "red")) +
 # theme(axis.text.x = element_text(angle = 45, hjust = 1))

features_of_interest <- c("EGFR", "PTK2")

dp <- DotPlot(
  object   = seurat_object1,
  features = features_of_interest,
  group.by = "main_annot3"   # or "main_annot2" if that has better labels
) +
  scale_colour_gradientn(
    colours = c("blue", "white", "red"),
    name    = "Average\nexpression"
  ) +
  scale_size(range = c(0.5, 6), name = "Percent\nexpressed") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# Show plot
pdf("20260102_COLOSSUS_scRNAseq_EGFR_PTK2_expression_cell_types_Epi_classification_plot.pdf")
print(dp)
dev.off()

pdf("20260102_COLOSSUS_scRNAseq_EGFR_PTK2_expression_cell_types_Epi_classification_feature_plot.pdf")

#feature plot
FeaturePlot(object = seurat_object1, features = 'PTK2')
dev.off()

pdf("S3I_20260102_COLOSSUS_scRNAseq_EGFR_PTK2_expression_cell_types_Epi_classification_dim_plot.pdf")

DimPlot(seurat_object1)

dev.off()

#pdf("S3J_20260303_GSE178378_tumor_MSS_CRC_LM_3_3_samples_only_scRNASeq_scATOMIC_EGFR_PTK2_dotplot.pdf")

DotPlot(
  seurat_object1,
  features = features_of_interest,
  group.by = "main_annot3",
  scale = TRUE
) +
  scale_colour_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-2.5, 2.5)
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
    title = "EIF4B Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )
#dev.off()



# Subset only epithelial cells
Idents(seurat_object1) <- "main_annot3"
epithelial_cells <- subset(seurat_object1, idents = "Epithelial (Malignant)")

# Identify variable features (if not already done)
epithelial_cells <- FindVariableFeatures(epithelial_cells)

# Scale the data
epithelial_cells <- ScaleData(epithelial_cells)

# Run PCA for dimensionality reduction
epithelial_cells <- RunPCA(epithelial_cells, features = VariableFeatures(epithelial_cells))

# Identify neighbors and clusters
epithelial_cells <- FindNeighbors(epithelial_cells, dims = 1:10)
epithelial_cells <- FindClusters(epithelial_cells, resolution = 0.5)

# Optional: Run UMAP for visualization
epithelial_cells <- RunUMAP(epithelial_cells, dims = 1:10)

#pdf("20250224_COLOSSUS_single_cell_epithelial_subclusters_Dim_Plot.pdf")
DimPlot(epithelial_cells, reduction = "umap", label = TRUE)
#dev.off()

#identify cluster markers
epi_markers <- FindAllMarkers(epithelial_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(epi_markers)

# Load dplyr for data manipulation
library(dplyr)

# Get the top 10 genes for each cluster based on average log2 fold change
top10_markers <- epi_markers %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)

# View the results
print(top10_markers)

setwd("Anguraj")

#write.table(top10_markers,"20250224_COLOSSUS_single_cell_top10markers_epithelial_subclusters.txt",sep ="\t", quote = FALSE, row.names = FALSE)

cluster0_markers <- subset(epi_markers, cluster == "0")
head(cluster0_markers, 25)  # Shows top 10 markers for cluster 0

cluster1_markers <- subset(epi_markers, cluster == "1")
head(cluster1_markers, 25)  # Shows top 10 markers for cluster 0

# Replace with your custom gene names
custom_genes <- c("CDC25A", "PRC1", "TTK","KIF23","CDCA2","BUB1B","EXO1","NCAPG2","FANCI","WDHD1","KIF11","DTL","CEP55","RACGAP1","ATAD2","ECT2","KIF18A","STIL")
FeaturePlot(epithelial_cells, features = custom_genes)
VlnPlot(epithelial_cells, features = custom_genes)

# Custom markers list should be provided as a list of gene vectors
custom_markers <- list(custom_genes)
epithelial_cells <- AddModuleScore(epithelial_cells, features = custom_markers, name = "CustomScore")

# Visualize the module score
FeaturePlot(epithelial_cells, features = "CustomScore1")


#Anguraj Script Methylation
#enrichment.scores <- escape.matrix(epithelial_cells, 
 #                                  gene.sets = gene.sets, 
  #                                 groups = 1000, 
   #                                min.size = 5)

ES <- runEscape(epithelial_cells, 
                method = "ssGSEA",
                gene.sets = gene.sets, 
                groups = 1000, 
                min.size = 5,
                new.assay.name = "escape.ssGSEA")

pdf("20250224_COLOSSUS_single_cell_epithelial_subclusters_hypo_hyper_methylation_genes_enrichment.pdf")
geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Positive")

geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Negative")
dev.off()


# reading scaled loadings
mv_loadings <- read.delim("MV1_marker.txt", stringsAsFactors = FALSE)
dim(mv_loadings)

m1 <- unique(mv_loadings$V1[mv_loadings$V3 == 1])
m2 <- unique(mv_loadings$V2[mv_loadings$V3 == 2])

#create gene list
gene.sets <- list(MV1Positive = m1,MV1Negative = m2)
gene.sets


ES_mv <- runEscape(seurat_object1, 
                method = "ssGSEA",
                gene.sets = gene.sets, 
                groups = 1000, 
                min.size = 5,
                new.assay.name = "escape.ssGSEA")
Idents(ES_mv) <- "main_annot3"

pdf("S10D_20260303_geyserEnrichment_scRNAseq_Epi_split_GSE178318_MV1_positive_genes.pdf")
geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Positive")

dev.off()

pdf("S10H_20260303_geyserEnrichment_scRNAseq_Epi_split_GSE178318_MV1_negative_genes.pdf")
geyserEnrichment(ES_mv, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Negative")
dev.off()
