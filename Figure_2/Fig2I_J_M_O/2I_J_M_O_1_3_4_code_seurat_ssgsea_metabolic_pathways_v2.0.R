################################################################################
# Figure Panels: 2I, 2J, 2M, 2O
# Description  : Seurat scRNA-seq analysis of COLOSSUS MSS KRAS-mutant tumours.
#                Runs escape/ssGSEA on single cells for selected KEGG metabolic
#                pathways, then generates geyser enrichment plots and
#                cell-type-stratified gene expression plots (GLS, GLUL, GOT2,
#                EIF4B, TCA genes, alanine/aspartate genes).
# Input data   : query_seurat_mss_kras_with_scATOMIC_epi_malignancy.Rds
#                c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt
# Output       : 2I_..._geyserEnrichment_scRNAseq_Epi_split_KEGG_ALANINE_...pdf
#                2J_..._GLS_GLUL_expression_cell_types_flux...pdf
#                2M_..._TCGA_genes_expression_cell_types...pdf
#                2O_..._EIF4B_genes_expression_cell_types...pdf
# R packages   : escape, Seurat, SeuratObject, SingleCellExperiment, scran,
#                qusage, RColorBrewer, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026 (v2.0)
################################################################################



# ══════════════════════════════════════════════════════════════════════════════
# Figures 2I, 2J, 2M, 2O — Single-cell ssGSEA and cell-type expression plots
# ══════════════════════════════════════════════════════════════════════════════
# Strategy:
#   1. Load annotated COLOSSUS scRNA-seq Seurat object (MSS KRAS-mutant).
#   2. Run escape package to score each cell against KEGG metabolic pathways
#      (single-cell ssGSEA equivalent).
#   3. Generate geyser enrichment plots (Fig 2I) and cell-type-stratified
#      violin plots for GLS/GLUL (Fig 2J), TCGA signature genes (Fig 2M),
#      and EIF4B (Fig 2O).
#
# Prerequisites: run the scATOMIC annotation script first to generate
#   'main_annot3' cell-type labels in the Seurat object.
#
# ESCAPE R package to do SSGSEA

# escape: single-cell gene set scoring (ssGSEA equivalent for scRNA-seq)
# scran: normalisation and HVG selection
# qusage: read GMT gene set files
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))



# ── Step 1: Load annotated Seurat object ─────────────────────────────────────
# This Seurat object must already contain scATOMIC cell-type annotations
# in the 'main_annot3' metadata column (run scATOMIC script first).
seurat_object1<-readRDS("query_seurat_mss_kras_with_scATOMIC_epi_malignancy.Rds") #seurat object

# ── Step 2: Load KEGG Legacy gene sets from MSigDB GMT file ──────────────────
# The GMT file is freely available at https://www.gsea-msigdb.org/gsea/msigdb/
#KEGG legacy gene sets
library(qusage)

kegg <- read.gmt("~/Downloads/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")

# reading scaled loadings
#mv_loadings <- read.delim("MV1_marker.txt", stringsAsFactors = FALSE)
#dim(mv_loadings)

m1 <- kegg$KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM
m2 <- kegg$KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM

#create gene list
# Define a named list of gene sets to score.
# MV1Positive and MV1Negative correspond to glyoxylate/dicarboxylate and
# alanine/aspartate/glutamate KEGG pathways respectively.
gene.sets <- list(MV1Positive = m1, MV1Negative = m2)
gene.sets

#Anguraj Script

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


#pdf("20241009_geyserEnrichment_scRNAseq_KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM_genes.pdf")
#geyserEnrichment(ES, 
 #                assay = "escape.ssGSEA",
  #               gene.set = "MV1Positive")

#dev.off()

pdf("20260102_geyserEnrichment_scRNAseq_Epi_split_KEGG_ALANINE_ASPARTATE_GLUTAMATE_METABOLISM_genes.pdf")
g2<-geyserEnrichment(ES, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Negative")
dev.off()

kruskal.test(g2$data$MV1Negative~g2$data$ident)


# Featureplot for GLS

metab_genes <- c("GLS", "GLUL", "GOT2")


pdf("20260102_COLOSSUS_scRNAseq_GLS_GLUL_expression_cell_types_flux.pdf")
DotPlot(object = seurat_object1, features = metab_genes)
dev.off()

p <- DotPlot(
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

pdf("2J_20260323_COLOSSUS_scRNAseq_GLS_GLUL_expression_cell_types_flux_increased_resolution.pdf")

p

dev.off()

#tca_cycle_genes <- c("CS", "ACO2", "IDH3A", "IDH3B", "IDH3G", "IDH2", 
   #                  "OGDH", "SUCLG1", "SUCLG2", "SUCLA2", 
   #                  "SDHA", "SDHB", "SDHC", "SDHD", 
    #                 "FH", "MDH1", "MDH2")
#DotPlot(object = seurat_object1, features = tca_cycle_genes)

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

p1 <- DotPlot(
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

pdf("2J_20260323_COLOSSUS_scRNAseq_TCGA_genes_expression_cell_types_flux_increased_resolution.pdf")

p1

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

p2 <- DotPlot(
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
    title = "Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )

pdf("S4L_20260323_COLOSSUS_scRNAseq_ala_asp_genes_expression_cell_types_flux_increased_resolution.pdf")

p2

dev.off()

#Figure 2O

# Filter for specific cell types
#filtered_data <- dp_data %>%
 # filter(id %in% cell, features.plot %in% gene) %>%
  #select(Gene = features.plot, CellType = id, PercentExpression = pct.exp)

# View the filtered table
#print(filtered_data)

##Figure 2M
dp <- DotPlot(object = seurat_object1, features = c("EIF4B"))


# Apply blue-to-red gradient
dp + scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#pdf("20260102_COLOSSUS_EIF5B_CRIPSR_selected_snRNAseq_dotplot.pdf")

# Apply blue-to-red gradient
dp + scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#dev.off()

p3 <- DotPlot(
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
    title = "Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )

pdf("2O_20260323_COLOSSUS_scRNAseq_EIF4B_genes_expression_cell_types_flux_increased_resolution.pdf")

p3

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

p4 <- DotPlot(
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
    title = "Expression Across Major Cell Types",
    colour = "Average Expression",
    size = "Percent Expressed"
  )

pdf("1M_20260323_COLOSSUS_scRNAseq_EIF4B_genes_expression_cell_types_flux_increased_resolution.pdf")

p4

dev.off()

pdf("1N_20260102_COLOSSUS_scRNAseq_EGFR_PTK2_expression_cell_types_Epi_classification_feature_plot.pdf")

#feature plot
FeaturePlot(object = seurat_object1, features = 'PTK2')
dev.off()

pdf("1N_20260102_COLOSSUS_scRNAseq_EGFR_PTK2_expression_cell_types_Epi_classification_dim_plot.pdf")

DimPlot(seurat_object1)

dev.off()

library(tidyr)

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
  dplyr::filter(features.plot %in% cripr_genes) %>%
  dplyr::select(features.plot, id, pct.exp) %>%
  dplyr::arrange(features.plot, desc(pct.exp))



#figure 4
# reading scaled loadings
mv_loadings <- read.delim("MV1_marker.txt", stringsAsFactors = FALSE)
dim(mv_loadings)

m1 <- unique(mv_loadings$V1[mv_loadings$V3 == 1])
m2 <- unique(mv_loadings$V1[mv_loadings$V3 == 2])

#create gene list
gene.sets <- list(MV1Positive = m1,MV1Negative = m2)
gene.sets

#Anguraj Script
#enrichment.scores <- escape.matrix(seurat_object1, 
#                                  gene.sets = gene.sets, 
#                                 groups = 1000, 
#                                min.size = 5)

ES_mv <- runEscape(seurat_object1, 
                method = "ssGSEA",
                gene.sets = gene.sets, 
                groups = 1000, 
                min.size = 5,
                new.assay.name = "escape.ssGSEA")
Idents(ES_mv) <- "main_annot3"
#GS.hallmark <- getGeneSets(library = "H")


#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

FeaturePlot(ES_mv, "MV1Positive") + 
  scale_color_gradientn(colors = colorblind_vector) + 
  theme(plot.title = element_blank(), keep.scale = "all")

FeaturePlot(ES_mv, "MV1Negative") + 
  scale_color_gradientn(colors = colorblind_vector) + 
  theme(plot.title = element_blank(), keep.scale = "all")


pdf("4F_20260105_geyserEnrichment_scRNAseq_MV1_positive_genes.pdf")
geyserEnrichment(ES_mv, 
                 assay = "escape.ssGSEA",
                 gene.set = "MV1Positive")

dev.off()

pdf("S10B_20260105_geyserEnrichment_scRNAseq_MV1_negative_genes.pdf")
geyserEnrichment(ES_mv, 
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
  theme_minimal(base_size = 24) +  # cleaner base theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 18),
    axis.title = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    size = guide_legend(title = "% Cells Expressing"),
    color = guide_colorbar(title = "Avg Expression")
  ) +
  ggtitle("MCM Gene Expression Across Cell Cycle Phases")

# Display
print(dot_plot)



dot_plot <- DotPlot(
  seurat_object1,
  features = mcm_genes,
  dot.scale = 6   # increases overall dot size (default ~6)
) +
  scale_color_gradientn(colors = c("darkblue", "white", "red")) +
  
  #refine size scaling range
  scale_size(range = c(4, 8)) +
  
  theme_minimal(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 18),
    axis.title = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    
    # clean panel
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    #  add border around full plot
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  ) +
  guides(
    size = guide_legend(title = "% Cells Expressing"),
    color = guide_colorbar(title = "Avg Expression")
  ) +
  ggtitle("MCM Gene Expression Across Cell Cycle Phases")

pdf("4G_20250615_COLOSSUS_cell_cycle_MCM_genes_dot_plot.pdf")

print(dot_plot)
dev.off()
