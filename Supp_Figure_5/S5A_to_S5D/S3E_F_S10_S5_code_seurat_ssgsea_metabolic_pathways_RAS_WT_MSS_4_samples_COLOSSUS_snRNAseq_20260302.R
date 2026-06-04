################################################################################
# Figure Panels : Supplementary Figures 5A–5D (also generates S5P)
# Description  : Seurat snRNA-seq analysis of 4 COLOSSUS RAS-wildtype MSS
#                samples. Runs escape/ssGSEA for KEGG alanine-aspartate-
#                glutamate metabolism; generates geyser enrichment plots and
#                cell-type-stratified expression plots for GLS, GLUL, GOT2,
#                TCA genes, alanine/aspartate genes, and EIF4B.
# Input data   : COLOSSUS snRNA-seq Seurat object (RAS-WT, 4 samples)
#                c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt
# Output       : S5A_..._geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype.pdf
#                S5B_..._GLS_GLUL_GOT2_expression_cell_types.pdf
#                S5C_..._tca_cancer_expressed_genes.pdf
#                S5D_..._ala_asp_genes.pdf
#                S5P_..._EIF4B_gene.pdf
# R packages   : escape, Seurat, SeuratObject, SingleCellExperiment, scran,
#                qusage, RColorBrewer, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
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
seurat_object1<-readRDS("query_seurat_WILDTYPE_MSI_MSS_with_scATOMIC_epi_malignancy_20260302.Rds") #seurat object

write.table(table(seurat_object1$main_annot3),"20260302_RAS_wildtype_MSS_4_samples_COLOSSUS_selected_snRNAseq_table_cell_types.txt", sep ="\t", quote = FALSE, row.names = FALSE)
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



pdf("S5A_20260302_geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype_4_samples_KEGG_ALANINE_ASPARTATE_GLUTAMATE_METABOLISM_genes.pdf")
geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Negative")
dev.off()

r2 <- geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Negative")

kruskal.test(r2$data$MV1Negative~r2$data$ident)


# Featureplot for GLS

metab_genes <- c("GLS", "GLUL", "GOT2")


#pdf("S5B_20260102_COLOSSUS_scRNAseq_GLS_GLUL_expression_cell_types_flux.pdf")
DotPlot(object = seurat_object1, features = metab_genes)
#dev.off()

pdf("S5B_20260302_geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype_4_samples_expression_cell_types_Epi_classification_GLS_GLUL_GOT2_plot.pdf")

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
    title = "Glutamate Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )


dev.off()





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

#pdf("20260102_COLOSSUS_scRNAseq_TCA_selected_genes__expression_cell_types_flux.pdf")

DotPlot(object = seurat_object1, features = tca_cancer_expressed_genes)
#dev.off()

pdf("S5C_20260302_geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype_4_samples_expression_cell_types_Epi_classification_tca_cancer_expressed_genes_plot.pdf")

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
epithelial_cells <- subset(seurat_object1, idents = "Epithelial cells")

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

#pdf("20260102_COLOSSUS_scRNAseq_Alanine_Aspartate_metabolism_selected_genes__expression_snRNAseq_cell_types.pdf")

DotPlot(object = seurat_object1, features = ala_asp)
#dev.off()

pdf("S5D_20260302_geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype_4_samples_expression_cell_types_Epi_classification_ala_asp_genes_plot.pdf")

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
    title = "PTK2 and EGFR Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )


dev.off()

#Figure S5O

#pdf("S5O_20260302_geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype_4_samples_expression_cell_types_Epi_classification_EIF4B_gene_plot.pdf")

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



cripr_genes <- c("TRIM46",
                 "ZNF689",
                 "KRTCAP2",
                 "REM1",
                 "CCND2",
                 "MUC20",
                 "MRPL14",
                 "EGFR",
                 "ORC2",
                 "POLR1D",
                 "PTK2",
                 "MOB4",
                 "MARS2",
                 "NDUFB3",
                 "PPP1R2",
                 "CFLAR",
                 "PCYT1A",
                 "FNTA",
                 "TFRC",
                 "SF3B1",
                 "HSPD1",
                 "GINS1",
                 "HSPE1")

d <- DotPlot(object = seurat_object1, features = cripr_genes)

dot_data <- d$data

dot_data %>%
  dplyr::filter(features.plot %in% c("EGFR", "PTK2")) %>%
  dplyr::select(features.plot, id, pct.exp) %>%
  dplyr::arrange(features.plot, desc(pct.exp))

#pdf("20250610_COLOSSUS_EGFR_PTK2_CRIPSR_selected_figure1_dotplot.pdf")

##Figure 1O
#dp <- DotPlot(object = seurat_object1, features = c("EGFR","PTK2"))

# Apply blue-to-red gradient
#dp + scale_color_gradientn(colors = c("blue", "white", "red")) +
 # theme(axis.text.x = element_text(angle = 45, hjust = 1))

features_of_interest <- c("EGFR", "PTK2")

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

pdf("S3F_20260302_geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype_4_samples_expression_cell_types_Epi_classification_EGFR_PTK2_plot.pdf")

p

dev.off()


pdf("20260102_COLOSSUS_scRNAseq_EGFR_PTK2_expression_cell_types_Epi_classification_feature_plot.pdf")

#feature plot
FeaturePlot(object = seurat_object1, features = 'PTK2')
dev.off()

pdf("S3E_20260302_geyserEnrichment_snRNAseq_Epi_split_RAS_wildtype_4_samples_expression_cell_types_Epi_classification_dim_plot.pdf")

DimPlot(seurat_object1)

dev.off()

# reading scaled loadings
mv_loadings <- read.delim("MV1_marker.txt", stringsAsFactors = FALSE)
dim(mv_loadings)

m1 <- unique(mv_loadings$V1[mv_loadings$V3 == 1])
m2 <- unique(mv_loadings$V2[mv_loadings$V3 == 2])

#create gene list
gene.sets <- list(MV1Positive = m1,MV1Negative = m2)
gene.sets

#Anguraj Script
#enrichment.scores <- escape.matrix(seurat_object1, 
#                                  gene.sets = gene.sets, 
#                                 groups = 1000, 
#                                min.size = 5)

ES <- runEscape(seurat_object1, 
                method = "ssGSEA",
                gene.sets = gene.sets, 
                groups = 1000, 
                min.size = 5,
                new.assay.name = "escape.ssGSEA")
Idents(ES) <- "main_annot3"

pdf("S10C_20260302_geyserEnrichment_scRNAseq_RAS_wildtype_4_samples_MV1_positive_genes.pdf")
geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Positive")

dev.off()

pdf("S10D_20260302_geyserEnrichment_scRNAseq_RAS_wildtype_4_samples_MV1_negative_genes.pdf")
geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Negative")
dev.off()


# cell cycle genes
MCM_genes <- c("MCM2","MCM3","MCM4","MCM5","MCM6","MCM7")
DotPlot(object = seurat_object1, features = MCM_genes)

# Define MCM genes
mcm_genes <- c("MCM3", "MCM4", "MCM5", "MCM6", "MCM7")

# Fetch expression data (normalized data)
mcm_expr <- FetchData(seurat_object1, vars = mcm_genes)

# Add Phase info to the same data frame
mcm_expr$Phase <- seurat_object1$Phase

# Load dplyr for summarizing
library(dplyr)

# Summarize mean expression of each gene by Phase
summary_table <- mcm_expr %>%
  group_by(Phase) %>%
  summarise(across(all_of(mcm_genes), mean, .names = "mean_{.col}"))

# View result
print(summary_table)

library(Seurat)
library(ggplot2)

# Set identity to Phase
Idents(seurat_object1) <- "Phase"

# Define MCM genes
mcm_genes <- c("MCM3", "MCM4", "MCM5", "MCM6", "MCM7")

# Create DotPlot

dot_plot <- DotPlot(seurat_object1, features = mcm_genes) +
  scale_color_gradientn(colors = c("lightblue", "blue", "darkblue")) +  # custom color gradient
  theme_minimal(base_size = 16) +  # cleaner base theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    size = guide_legend(title = "% Cells Expressing"),
    color = guide_colorbar(title = "Avg Expression")
  ) +
  ggtitle("MCM Gene Expression Across Cell Cycle Phases")

# Display
pdf("20250615_COLOSSUS_cell_cycle_MCM_genes_dot_plot.pdf")
print(dot_plot)
dev.off()
