################################################################################
# Figure Panels : 3A, 3B, 3C, 3E
# Description   : ELMER (Enhancer Linking by Methylation/Expression Relationships)
#                 analysis of TCGA COAD RAS-mutant tumours — HYPOMETHYLATION arm.
#                 Identifies distal CpG probes hypomethylated in tumours vs normal,
#                 links them to nearby gene expression changes, and identifies
#                 enriched transcription factor motifs at hypomethylated enhancers.
#                 Generates: (3A) volcano/scatter of hypo probes, (3B) gene
#                 expression heatmap, (3C) motif enrichment plot, (3E) probe–gene
#                 pair visualisation.
# Input data    : 20251105_TCGA_normal_solid_tissue_Methylation_full_probes.txt
#                 20251105_TCGA_tumour_solid_tissue_Methylation_full_probes.txt
#                 TCGA RNA-seq expression matrix
#                 tcga_coadread_pancanatlas_RAS_mutant_MSS_TCGA_barcodes.txt
# Output        : ELMER results in results/result/ subfolder (RDS + tables)
# R packages    : ELMER, ELMER.data, MultiAssayExperiment, SummarizedExperiment,
#                 GenomicInteractions, plyr, doParallel
# Authors       : Sadanandam Lab, Institute of Cancer Research, London
# Contact       : anguraj.sadanandam@icr.ac.uk
# Date          : March 2026 (v3.0)
# Note          : ELMER is computationally intensive. Run on a multi-core server
#                 with ≥32 GB RAM. Set cores with doParallel::registerDoParallel().
#                 Update input file names to match your local copies.
################################################################################

# ELMER: Enhancer Linking by Methylation/Expression Relationships
# doParallel enables multi-core execution — set cores below.
library(ELMER)
library(ELMER)

library(MultiAssayExperiment)
library(ELMER.data)
library(SummarizedExperiment)
library(plyr)
library(doParallel)
library(GenomicInteractions)
# get distal probes that are 2kb away from TSS on chromosome 1
#distal.probes <- get.feature.probe(
 # genome = "hg19", 
  #met.platform = "450K", 
  #rm.chr = paste0("chr",c(2:22,"X","Y"))
#)


# ── Step 1: Load methylation matrices ─────────────────────────────────────────
# Two separate files: normal colon tissue and primary tumour methylation.
# Step 1: Load existing methylation data
#normal.db <- read.delim("20250128_TCGA_normal_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
normal.db <- read.delim("20251105_TCGA_normal_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
message("Normal data dimensions: ", paste(dim(normal.db), collapse = " x "))

#tumor_data <- read.delim("20250128_TCGA_primary_tumour_Methylation_full_probes.txt", stringsAsFactors = FALSE)
tumor_data <- read.delim("20251105_TCGA_tumour_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
message("Tumor data dimensions: ", paste(dim(tumor_data), collapse = " x "))

# ── Step 2: Filter probes with complete data ───────────────────────────────────
# Remove any probes with NA values in either normal or tumour data.
# Step 2: Extract and preprocess probe IDs
message("Processing methylation probes...")

#filter for probes
na_filter <- complete.cases(normal.db) & complete.cases(tumor_data)
normal.db <- normal.db[na_filter, , drop = FALSE]
tumor_data <- tumor_data[na_filter, , drop = FALSE]

# Combine datasets
# Combine normal + tumour into a single matrix for ELMER input.
combined_data <- cbind(normal.db, tumor_data)
group_labels <- c(rep("normal", ncol(normal.db)), rep("tumor", ncol(tumor_data)))
dim(combined_data)

#Expression
expr <- read.delim("HiSeqV2",stringsAsFactors = FALSE)
dim(expr)

me_t <- match(colnames(combined_data), colnames(expr))
w_me_t <- which(!is.na(me_t))
length(w_me_t)

#expr_1 <- cbind(expr$sample, expr[,me_t[w_me_t]])
expr_1 <- expr[,me_t[w_me_t]]
rownames(expr_1) <- expr$sample
dim(expr_1)

combined_data_1 <- combined_data[,w_me_t]
rownames(combined_data_1) <- rownames(combined_data)
dim(combined_data_1)

group_labels_1 <- group_labels[w_me_t]

df <- data.frame(cbind(colnames(combined_data_1),group_labels_1))
colnames(df) <- c("sample ID","groups")
dim(df)

sample.info <- S4Vectors::DataFrame(
  primary =  df$`sample ID`, 
  sample.type = df$groups
  )



#ENSEMBL ID conversion
library(biomaRt)

# Function to convert gene symbols to Ensembl IDs
gene_symbol_to_ensembl <- function(gene_symbols) {
  # Connect to the Ensembl database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Query the Ensembl database to retrieve Ensembl IDs for the given gene symbols
  result <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = ensembl
  )
  
  # Return the result
  return(result)
}

ensembl <- gene_symbol_to_ensembl(rownames(expr_1))

m <- match(rownames(expr_1),ensembl$hgnc_symbol)
w <- which(!is.na(m))
length(w)

expr_2 <- expr_1[w,]
rownames(expr_2) <- ensembl$ensembl_gene_id[m[w]]
dim(expr_2)

#ELMER results
#102,518 enhancer probes and 91,252 mapped (because NA probes removed)


el <- read.delim("20250128_ELMER_102518_enhancer_probes.txt", stringsAsFactors = FALSE)
dim(el)

mm <- match(rownames(combined_data_1),el$Probe)
ww <- which(!is.na(mm))
length(ww)

combined_data_2 <- combined_data_1[ww,]
dim(combined_data_2)

#Load reference file with sample barcodes - KRAS mutation
reference_file <- read.delim("tcga_coadread_pancanatlas_RAS_mutant_MSS_TCGA_barcodes.txt", stringsAsFactors = FALSE)
valid_samples <- reference_file$TCGA_RAS_mutant_samples  # Extract valid sample barcodes

#Match and filter the datasets based on valid sample barcodes
valid_normal_tumour_samples <- substr(colnames(combined_data_2),1,12) %in% substr(gsub("-","\\.",valid_samples),1,12)
valid_expression_samples <- substr(colnames(expr_2),1,12) %in% substr(gsub("-","\\.",valid_samples),1,12)

# Filter combined_data_2 to keep only matched samples
combined_data_2_filtered <- combined_data_2[, valid_normal_tumour_samples, drop = FALSE]
dim(combined_data_2_filtered)

#write.table(combined_data_2_filtered,"20260331_TCGA_tumour_solid_tissue_NORMAL_Methylation_full_probes.txt", sep = "\t", quote = FALSE)

# Filter expr_2 to keep only matched samples
expr_2_filtered <- expr_2[, valid_expression_samples, drop = FALSE]
dim(expr_2_filtered)

#expr_1_filtered <- expr_1[, valid_expression_samples, drop = FALSE]
#dim(expr_1_filtered)
    
group_labels_2 <-  group_labels_1[valid_normal_tumour_samples, drop = FALSE]
    
    df <- data.frame(cbind(colnames(combined_data_2_filtered),group_labels_2))
    colnames(df) <- c("sample ID","groups")
    dim(df)
    
    sample.info <- S4Vectors::DataFrame(
      primary =  df$`sample ID`, 
      sample.type = df$groups
    )
    
mae <- createMAE(
  exp = expr_2_filtered, 
  met = combined_data_2_filtered,
  colData = sample.info,
  save = TRUE,
  linearize.exp = FALSE,
  save.filename = "mae_TCGA.rda",
  #filter.probes = distal.probes,
  met.platform = "450K",
  genome = "hg19",
  TCGA = TRUE
)

group.col <- "sample.type"
group1 <- "tumor"
group2 <-  "normal"
dir.out <- "result"
diff.dir <-  "hypo" # Search for hypo/hypermethylated probes in group 1
#diff.dir <-  "hyper" # uncomment if need to run hypermethylation

  sig.diff <- get.diff.meth(
    data = mae, 
    group.col = group.col,
    group1 = group1,
    group2 = group2,
    minSubgroupFrac = 1,
    sig.dif = 0.3,
    diff.dir = diff.dir,
    cores = 1, 
    dir.out = dir.out,
    pvalue = 0.05,
    mode = "supervised"
  )

dim(sig.diff)

## volcano_plot_probes.R

# Load packages ------------------------------------------------------------
library(ggplot2)

# 1. Read the results table -----------------------------------------------
# If your file is tab-separated (like the snippet you pasted), use read.delim:
#res <- read.delim("getMethdiff.probes.txt", stringsAsFactors = FALSE)

# If it is comma-separated instead, use:
 res <- read.csv("getMethdiff.hypo.probes.csv", stringsAsFactors = FALSE)

# 2. Define cutoffs --------------------------------------------------------
x.col  <- "tumor_Minus_normal"  # effect size column
p.cut  <- 0.01                  # FDR cutoff (adjusted p-value)
fc.cut <- 0.30                  # methylation difference cutoff

# 3. Prepare data for plotting --------------------------------------------
res$negLog10FDR <- -log10(res$adjust.p)

# Label points by significance & direction
res$group <- "Not Significant"
res$group[res$adjust.p < p.cut & res[[x.col]] >  fc.cut] <- "Hypermethylated in tumor"
res$group[res$adjust.p < p.cut & res[[x.col]] < -fc.cut] <- "Hypomethylated in tumor"

# 4. Create PDF volcano plot ----------------------------------------------
#pdf("volcanoPlot.probes.hypo.pdf", width = 6, height = 5)

ggplot(res, aes_string(x = x.col, y = "negLog10FDR", colour = "group")) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_vline(xintercept = c(-fc.cut, fc.cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(p.cut), linetype = "dashed") +
  scale_color_manual(values = c("Not Significant"          = "grey50",
                                "Hypermethylated in tumor" = "red",
                                "Hypomethylated in tumor"  = "blue")) +
  labs(
    x = expression(paste("DNA Methylation difference (", beta, "-values)")),
    y = expression(-Log[10]~"(FDR adjusted P-values)"),
    colour = ""
  ) +
  ggtitle("Volcano plot - Probes methylated in tumor vs normal") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

#dev.off()


nearGenes <- GetNearGenes(
  data = mae, 
  probes = sig.diff$probe, 
  numFlankingGenes = 20
) # 10 upstream and 10 dowstream genes


pair <- get.pair(
  data = mae,
  group.col = group.col,
  group1 = group1,
  mode = "supervised",
  group2 = group2,
  nearGenes = nearGenes,
  diff.dir = diff.dir,
  minSubgroupFrac = 1, # % of samples to use in to create groups U/M
  permu.dir = file.path(dir.out,"permu"),
  permu.size = 100, # Please set to 100000 to get significant results
  raw.pvalue = 0.05,   
  Pe = 0.001, # Please set to 0.001 to get significant results
  filter.probes = TRUE, # See preAssociationProbeFiltering function
  filter.percentage = 0.05,
  filter.portion = 0.3,
  dir.out = dir.out,
  cores = 1,
  label = diff.dir,
  correlation = "negative"
)

pair <- read.delim("getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation_2000_removed.txt", stringsAsFactors = FALSE)
dim(pair)
# Extract DNA methylation and gene expression data
methylation_data <- assay(getMet(mae))  # Methylation levels
expression_data <- assay(getExp(mae))   # Gene expression levels

# Create a function to compute correlation
compute_correlation <- function(probe, gene) {
  meth_values <- methylation_data[probe, ]  # Get methylation values
  exp_values <- expression_data[gene, ]     # Get expression values
  cor.test(meth_values, exp_values, method = "pearson")$estimate
}

# Apply function to all significant pairs
pair$Correlation <- mapply(compute_correlation, pair$Probe, pair$GeneID)

w_cor<- which(pair$Correlation  <= -0.2)
length(w_cor)

# check for missing values - this should be 0
sum(is.na(pair$Correlation))

# View results with correlation values
head(pair[, c("Probe", "GeneID", "Distance", "Correlation")])

#write.table(pair,"result/getPair.hypo.pairs.significant_with_correlation_Anguraj.txt", sep="\t", quote = FALSE)
#write.table(pair[w_cor,],"result/getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation.txt", sep="\t", quote = FALSE)


#average z-score methylation and expression
library(ggplot2)
publication_theme = theme(
  title = element_text(size = 20),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14)
)

ave_zscore_methy <- apply(methylation_data[unique(pair$Probe[w_cor]),],2,mean)
ave_zscore_expr <- apply(expression_data[unique(pair$GeneID[w_cor]),],2,mean)

fig3cd_data <- data.frame(cbind(scale(ave_zscore_expr),scale(ave_zscore_methy)))
names(fig3cd_data) = c("expression", "methylation")

fig3c_lm = coef(lm(formula = expression~methylation, data=fig3cd_data))

cor.test(fig3cd_data$expression,fig3cd_data$methylation)

plot(fig3cd_data$methylation,fig3cd_data$expression)
#-0.7130365 


#pdf("3C_20260331_TCGA_hyomethylation_cor_-0.2_4183_average_zscore_RNA_methylation.pdf")
ggplot(fig3cd_data, aes(y=expression, x=methylation, colour = "#F98400")) +
  geom_abline(intercept = fig3c_lm[1], slope=fig3c_lm[2], lty=2, lwd=1, col = "#F98400")+
  geom_point(size=3, alpha=1)+
  scale_color_manual(values = c("#F98400"))+
  theme_bw()+
  labs(y="Average RNA expression Z-score", x="Average Methylation Z-score")+
  publication_theme
#dev.off()


# Enrichment Hallmarks 
## Packages -------------------------------------------------------------
library(dplyr)
library(stringr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

fig3kl_hyper_elsevier = read.table("3K_20260408_TCGA_hypomethylation_corr_0.2_4183_genes_Enrichr_results_bar_2000_removed_MSigDB_Hallmark_2020_table.txt", sep = "\t", header = T, quote = "\"")

fig3l_df = head(fig3kl_hyper_elsevier,8)
fig3l_df$Term = factor(fig3l_df$Term, levels = rev(fig3l_df$Term))

#pdf("3K_20260408_TCGA_hyomethylation_cor_-0.2_4183_average_MsigDB_Hallmarks.pdf")

ggplot(data = fig3l_df, aes(y = Term, x = -log10(Adjusted.P.value), fill = Odds.Ratio)) +
  geom_bar(stat = "identity") +  # Removed show.legend = F so the legend is displayed
  geom_text(aes(label = Term), x = 0.05, hjust = 0, vjust = 0.5, size = 4) +
  scale_fill_gradient(low = "#EEEEEE", high = "#8866FF", na.value = NA) +
  labs(title = "MSigDB Hallmarks", fill = "Odds Ratio") +  # Adds legend title
  theme_bw() +
  publication_theme +
  theme(axis.text.y = element_blank())

#dev.off()

#cnet plot for hallmarks

## 2) Take TOP 5 terms --------------------------------------------------
## (Assumes the table is already ordered by significance/score)
fig3l_df <- head(fig3kl_hyper_elsevier, 5)

## Optional: keep the current order but reverse it for plotting
fig3l_df$Term <- factor(fig3l_df$Term, levels = rev(fig3l_df$Term))

## 3) Build gene sets and gene universe from the FULL table -------------
##    (geneSets = list of genes per term; universe = all unique genes)

gs_list <- fig3kl_hyper_elsevier %>%
  select(Term, Genes) %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = str_trim(Genes)) %>%
  filter(Genes != "") %>%
  group_by(Term) %>%
  summarise(genes = list(unique(Genes)), .groups = "drop")

geneSets <- setNames(gs_list$genes, gs_list$Term)

# gene universe = every unique gene that appeared in any term
gene_universe <- unique(unlist(geneSets))

## 4) Build the 'result' slot for the TOP 5 only ------------------------
##    (clusterProfiler-like result data.frame)

res_top <- fig3l_df %>%
  mutate(
    ID          = as.character(Term),             # internal ID
    Description = as.character(Term),
    GeneRatio   = Overlap,                        # e.g. "63/200"
    BgRatio     = Overlap,                        # not used by cnetplot; OK to recycle
    pvalue      = P.value,
    p.adjust    = Adjusted.P.value,
    qvalue      = Adjusted.P.value,
    geneID      = gsub(";", "/", Genes),          # clusterProfiler uses "/" separator
    Count       = as.numeric(sub("/.*", "", Overlap))
  ) %>%
  select(ID, Description, GeneRatio, BgRatio,
         pvalue, p.adjust, qvalue, geneID, Count)

# rownames as IDs (typical for enrichResult)
row.names(res_top) <- res_top$ID

## 5) Create enrichResult object ----------------------------------------
er_hallmark <- new(
  "enrichResult",
  result        = res_top,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  qvalueCutoff  = 1,
  gene          = gene_universe,              # genes considered in analysis
  universe      = gene_universe,              # MUST be character, not NULL
  geneSets      = geneSets,
  organism      = "UNKNOWN",
  ontology      = "MSigDB_Hallmark_2020",
  keytype       = "SYMBOL",
  readable      = FALSE
)

## 6) ClusterProfiler CNET plot (top 5 terms) ---------------------------
p_cnet <- cnetplot(
  er_hallmark,
  showCategory = 5,      # <--- only top 5 terms
  foldChange   = NULL,   # add named vector if you have gene logFC
  circular     = FALSE,
  colorEdge    = TRUE
) +
  ggtitle("CNET plot – MSigDB Hallmark (top 5 terms)")

#pdf("20260331_TCGA_hyomethylation_cor_-0.2_4183_average_MsigDB_Hallmarks_cnet+plot.pdf")
print(p_cnet)
 # dev.off()

# metabolic work bench

fig3_hyper_elsevier = read.table("S7F_Metabolomics_Workbench_Metabolites_2022_table.txt", sep = "\t", header = T, quote = "\"")

fig3_df = head(fig3_hyper_elsevier,8)
fig3_df$Term = factor(fig3_df$Term, levels = rev(fig3_df$Term))

#pdf("20260331_TCGA_hyomethylation_cor_-0.2_4183_average_Metabolomics_Workbench_Metabolites_2022.pdf")

ggplot(data=fig3_df, aes(y=Term, x=-log10(Adjusted.P.value), fill = Odds.Ratio))+
  geom_bar(stat = "identity", show.legend = T)+
  geom_text(data = fig3_df, aes(label=Term), x = 0.05, hjust=0, vjust=0.5, size=4)+
  scale_fill_gradient(low = "#EEEEEE", high = "#8866FF", na.value = NA)+
  labs(title="Metabolomics_Workbench_Metabolites_2022")+
  theme_bw()+
  publication_theme+
  theme(axis.text.y = element_blank())

#dev.off()



#MYC target - correlation

library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

# MYC target genes
myc_targets <- c(
  "CDK4","MYC","PLK4","NOP56","NOLC1","SYNCRIP",
  "MCM2","PA2G4","MCM4","HSPE1","HSPD1","NME1","PLK1"
)

# Input file
pair_file <- "getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation_2000_removed.txt"

# Output file
out_pdf <- "Top10_MYC_targets_probe_gene_barplot_reversed.pdf"

# Read file
pair_df <- read_tsv(pair_file, show_col_types = FALSE)

# Keep only MYC target genes
pair_myc <- pair_df %>%
  filter(Symbol %in% myc_targets)

# Stop if nothing found
if (nrow(pair_myc) == 0) {
  stop("No MYC target genes found in the input file.")
}

# Select top 10 strongest negative correlations
top10_myc <- pair_myc %>%
  arrange(Correlation) %>%
  slice_head(n = 10) %>%
  mutate(
    Pair = paste0(Symbol, " (", Probe, ")"),
    Pair = fct_reorder(Pair, Correlation, .desc = TRUE)
  )

# Make publication-quality bar plot
p <- ggplot(top10_myc, aes(x = Pair, y = Correlation)) +
  geom_col(width = 0.7, fill = "#2C7FB8") +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "black") +
  geom_text(
    aes(label = round(Correlation, 3)),
    hjust = 1.1,
    colour = "white",
    size = 3.8
  ) +
  labs(
    title = "Top 10 MYC-target probe–gene pairs",
    subtitle = "Ranked by strongest negative methylation-expression correlation",
    x = NULL,
    y = "Correlation"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.5)
  )

print(p)

# Save plot
ggsave(
  filename = "3J_20260331_TCGA_hyomethylation_cor_-0.2_4183_average_Myc_targets_Pearson_correlation_top_10.pdf",
  plot = p,
  width = 8,
  height = 5.5
)

# Optional: print selected rows
print(top10_myc)





#PTK2 analysis
# Filter the 'pair' data frame for rows with Symbol equal to "PTK2"
ptk2_pair <- subset(pair[w_cor,], Symbol == "PTK2")
dim(ptk2_pair)

# Check the filtered results
print(ptk2_pair)
dim(ptk2_pair)

# Load ggplot2 for plotting
library(ggplot2)

# Create a bar plot of the correlation coefficient for each probe associated with PTK2
ggplot(ptk2_pair, aes(x = reorder(Probe, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Correlation Coefficient for Probes Associated with PTK2",
       x = "Probe",
       y = "Correlation Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Dot plot
library(ggplot2)
#pdf("S7A_20260331_TCGA_hyomethylation_cor_-0.2_4183_PTK2_correlation_probes.pdf")
ggplot(ptk2_pair, aes(x = reorder(Probe, Correlation), y = Correlation, color = Correlation)) +
  geom_point(size = 4) +
  labs(title = "Dot Plot: Correlation Coefficient for Probes Associated with PTK2",
       x = "Probe",
       y = "Correlation Coefficient",
       color = "Correlation") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  coord_flip()
#dev.off()

# Plot with color gradient based on correlation
#pdf("3D_20260331_TCGA_PTK2_methylation_7_highly_correlated_plot.pdf")
ggplot(ptk2_pair, aes(x = Distance, y = Correlation, label = Probe, color = "blue")) +
  geom_point(size = 4) +
  geom_segment(data = ptk2_pair,
               aes(x = Distance, xend = Distance, y = 0, yend = Correlation),
               color = "gray30")+
  geom_text(vjust = -1, size = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  theme_minimal() +
  labs(
    title = "Correlation vs Distance Colored by Corr Value",
    x = "Distance",
    y = "Correlation",
    color = "Corr"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()


#MYC plot
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

# -----------------------------
# Input file
# -----------------------------
pair_file <- "getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation_2000_removed.txt"

# Output file
out_pdf <- "3H_20260331_TCGA_hypomethylation_cor_-0.2_MYC_distance_vs_correlation.pdf"

# -----------------------------
# Read data
# -----------------------------
pair_cor <- read_tsv(pair_file, show_col_types = FALSE)

# -----------------------------
# Filter for MYC
# -----------------------------
myc_pair <- pair_cor %>%
  filter(Symbol == "MYC") %>%
  arrange(Distance)

print(myc_pair)
dim(myc_pair)

if (nrow(myc_pair) == 0) {
  stop("No MYC rows found in the correlation file.")
}

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(myc_pair, aes(x = Distance, y = Correlation)) +
  geom_segment(
    aes(x = Distance, xend = Distance, y = 0, yend = Correlation),
    linewidth = 0.5,
    colour = "grey50"
  ) +
  geom_point(
    aes(fill = Correlation),
    shape = 21,
    size = 4.2,
    colour = "black",
    stroke = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "#D55E00"
  ) +
  ggrepel::geom_text_repel(
    aes(label = Probe),
    size = 3.4,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey70",
    max.overlaps = Inf
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Correlation"
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8),
    labels = comma
  ) +
  labs(
    title = "MYC-associated probes",
    subtitle = "Distance from MYC gene locus versus methylation-expression correlation",
    x = "Distance from MYC gene locus",
    y = "Correlation"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

print(p)

# -----------------------------
# Save
# -----------------------------
ggsave(
  filename = out_pdf,
  plot = p,
  width = 8,
  height = 5.5
)


# Identify enriched motif for significantly hypomethylated probes which 
# have putative target genes.

fix_getHocomocoTable <- function() {
  data("hocomoco.table", package = "ELMER.data", envir = environment())
  hocomoco.table
}

assignInNamespace(
  x = "getHocomocoTable",
  value = fix_getHocomocoTable,
  ns = "ELMER"
)

test_hocomoco <- ELMER:::getHocomocoTable()
print(dim(test_hocomoco))
print(head(test_hocomoco[, c("Model", "TF family", "TF subfamily")]))

enriched.motif <- get.enriched.motif(
  data = mae,
#  probes = unique(pair$Probe[w_cor]), 
  probes = unique(pair$Probe), 
  dir.out = "result", 
  label = diff.dir,
  min.incidence = 3,
  lower.OR = 1.1,
  pvalue = 0.2
)

TF <- get.TFs(
  data = mae, 
  mode = "supervised",
  group.col = group.col,
  group1 = group1,
  group2 = group2,
  enriched.motif = enriched.motif,
  dir.out = dir.out,
  cores = 1, 
  diff.dir = diff.dir,
  correlation = "positive",
  topTFper = 0.05,
  save.plots = FALSE,
  label = diff.dir
)


#Plots

gee <- which(pair$Symbol[w_cor] == "MYC")
gee

#pdf("20260408_TCGA_methylation_MYC_schematic_plot_Elmer.pdf")
schematic.plot(
  pair = pair, 
  data = mae,
  group.col = group.col,
  byProbe = pair$Probe[w_cor[gee[8]]],
  numFlankingGenes = 20,
  save = TRUE
)
#dev.off()




motif.enrichment.plot(motif.enrichment = "getMotif.hypo.motif.enrichment.csv", 
                      significant = list(OR = 2.8,lowerOR = 1.3), 
                      label = "hypo", 
                      save = FALSE)  

load("getTF.hypo.TFs.with.motif.pvalue.rda")
motif <- colnames(TF.meth.cor)
TF.rank.plot(
  motif.pvalue = TF.meth.cor, 
  motif = motif,
  TF.label = createMotifRelevantTfs("family")[motif],
  save = FALSE
) 


pair_results <- read.table("getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation_2000_removed.txt")

heatmapPairs(data = mae, 
             group.col = group.col,
             group1 = group1, 
             group2 = group2,
             pairs = pair,
           #  filename = NULL)

        filename =  "S6B_20260408_TCGA_getPair.hypo.pairs.significant_with_correlation_Anguraj_0.2_correlated_heatmap_methylation_genes.pdf")


# Load required libraries
library(ggplot2)
library(dplyr)

# Read the CSV file (modify the file path)
df <- read.csv("getMotif.hypo.motif.enrichment.csv")

# Filter motifs where OR > 2.5
#df_filtered <- df %>%
 # filter(OR > 3) %>%
#  arrange(desc(OR))  # Sort by OR in descending order

# Filter motifs where OR > 2.5
df_filtered <- df %>%
  filter(OR > 2, FDR < 0.05) %>%
  arrange(desc(OR))  # Sort by OR in descending order

# Convert motif names into a factor for ordered plotting
df_filtered$motif <- factor(df_filtered$motif, levels = df_filtered$motif[order(df_filtered$OR, decreasing = TRUE)])

# Compute -log10(FDR) for visualization
df_filtered$logFDR <- -log10(df_filtered$FDR)

# Define a pastel color palette
pastel_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#E6BAFF", "#FFC3E1")

pdf(file ="3G_20260403_TCGA_Enhancer_Motif_Enrichment_OR_G_3.pdf")
# Create the bar plot with pastel rainbow colors
ggplot(df_filtered[1:7,], aes(x = motif, y = OR, fill = logFDR)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +  # Flip axes for readability
  scale_fill_gradientn(colors = pastel_colors) +  # 🌈 Light pastel colors
  labs(title = "Enhancer Motif Enrichment (OR > 2.5, Light Colors)",
       x = "Motif",
       y = "Odds Ratio (OR)",
       fill = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
dev.off()

write.table(df_filtered,"3G_20260403_TCGA_Enhancer_Motif_Enrichment_OR_G_3_output.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Origial paper CRC probes
#ELMER_methy <- read.delim("20250128_Elmer_CRC_methylation_original_publication_results.txt", stringsAsFactors = FALSE)
#dim(ELMER_methy)

#mmm <- match(sig.diff_hyper$probe,ELMER_methy$Methy_ID)
#www <- which(!is.na(mmm))
#length(www)

# output table similar to their paper
# The function also writes a CSV file with the summary (OR, p-values, etc.)
dir_out <- "result"
label   <- "hypo"
# --------------------------------------------------------------
# 4) Construct motif summary from enriched_motif
# --------------------------------------------------------------
# Here we derive the number of probes per motif and format it as: "count(percentage)"
all_probes <- pair$Probe
motif_summary <- data.frame(
  motif       = names(enriched.motif),
  NumOfProbes = sapply(enriched.motif, length),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Probes = paste0(NumOfProbes, "(",
                    sprintf("%4.1f", NumOfProbes / length(all_probes) * 100),
                    "%)")
  )

# --------------------------------------------------------------
# 5) Retrieve enrichment statistics (e.g., enrichFold) from the summary CSV
# --------------------------------------------------------------
# The CSV file (written by get.enriched.motif) contains columns such as 'motif' and 'OR'
stats_file <- file.path(dir_out, paste0("getMotif.", label, ".motif.enrichment.csv"))
motif_stats <- read.csv(stats_file)
# Assume the odds ratio column is named 'OR'. We rename it to 'enrichFold' later.

# Merge the enrichment stats with our motif_summary.
# (We keep only the 'motif' and 'OR' column from motif_stats.)
motif_summary_final <- motif_summary %>%
  left_join(motif_stats %>% select(motif, OR), by = "motif") %>%
  rename(enrichFold = OR)

# --------------------------------------------------------------
# 6) Load the probe–gene pairs file
# --------------------------------------------------------------
# This file contains columns such as: Probe, GeneID, Symbol, Distance, Sides, Raw.p, FDR, Correlation
pairs_df   <- read.csv("getPair.hypo.pairs.significant_with_correlation_Anguraj_0.3_correlated.csv")

# --------------------------------------------------------------
# 7) Compute "# of pair" and "# of gene" for each motif
# --------------------------------------------------------------
# For each motif (from enriched_motif), filter pairs_df to count the total number of probe–gene pairs 
# and the number of unique genes.
df_counts <- lapply(names(enriched.motif), function(mtf) {
  motif_probes <- enriched.motif[[mtf]]
  sub_df <- pairs_df %>% filter(Probe %in% motif_probes)
  
  data.frame(
    motif       = mtf,
    `# of Probes` = length(unique(sub_df$Probe)),
    `# of pair` = nrow(sub_df),
    `# of gene` = length(unique(sub_df$GeneID))
  )
}) %>% bind_rows()

# --------------------------------------------------------------
# 8) Merge all the information into the final table
# --------------------------------------------------------------
final_output <- motif_summary_final %>%
  left_join(df_counts, by = "motif") %>%
  select(motif, `X..of.Probes`, `X..of.gene`, `X..of.pair`, enrichFold)

# --------------------------------------------------------------
# 9) Print and optionally save the final table
# --------------------------------------------------------------
print(final_output)
write_csv(final_output, file.path(dir_out, "final_motif_table.csv"))



# --------------------------------------------------------------
# 4) Compute the average methylation per motif across samples
# --------------------------------------------------------------
# Retrieve the methylation matrix from elmer.results.
meth_matrix <- methylation_data
  
# For each motif, compute the average methylation over its associated probes.
motif_meth <- lapply(names(enriched.motif), function(motif) {
  probes <- enriched.motif[[motif]]
  if (length(probes) == 1) {
    # If only one probe is present, use its methylation values directly.
    return(meth_matrix[probes,])
  } else {
    # Otherwise, calculate the average methylation (row-wise mean) across the probes.
    return(colMeans(meth_matrix[probes, ,drop = FALSE], na.rm = TRUE))
  }
})

# Combine the list into a data frame (samples as rows, motifs as columns)
motif_meth_df <- as.data.frame(motif_meth)
colnames(motif_meth_df) <- names(enriched.motif)

# --------------------------------------------------------------
# 5) Retrieve and subset TF expression data using tf_summary
# --------------------------------------------------------------
# expression_data$exp_matrix contains the full gene expression data.
# The tf_summary data frame has a column "top_5percent_TFs" that contains semicolon-separated TF names.
# We split these strings and take the union of all TFs to form our TF list.
tf_list <- unique(unlist(strsplit(TF$top_5percent_TFs[1], ";")))
exp_matrix <- expr_1_filtered
exp_matrix_tf <- exp_matrix[rownames(exp_matrix) %in% tf_list, ]
dim(exp_matrix_tf)
# --------------------------------------------------------------
# 6) For each motif, calculate the Pearson correlation between its average methylation and expression of the selected TFs
# --------------------------------------------------------------
cor_results <- lapply(colnames(motif_meth_df), function(motif) {
  meth_values <- motif_meth_df[[motif]]
  
  # For each TF (each row in exp_matrix_tf), compute Spearman correlation with methylation values.
  tf_cor <- apply(exp_matrix_tf, 1, function(tf_expr) {
    test <- cor.test(meth_values, tf_expr, method = "pearson")
    c(correlation = test$estimate, p.value = test$p.value)
  })
  
  # Transpose the result and convert to a data frame.
  tf_cor_df <- as.data.frame(t(tf_cor))
  tf_cor_df$TF <- rownames(tf_cor_df)
  tf_cor_df$motif <- motif
  tf_cor_df
})

cor_results_df <- do.call(rbind, cor_results)

# --------------------------------------------------------------
# 7) Rank all TFs by the strength of inverse (negative) correlation for each motif
# --------------------------------------------------------------
# TFs with more negative correlation (i.e. stronger inverse relationship) are ranked higher.
ranked_results <- cor_results_df %>%
  group_by(motif) %>%
  arrange(cor_results_df$correlation.rho) %>%   # ascending order: more negative values come first
  mutate(rank = row_number()) %>%
  ungroup()

# --------------------------------------------------------------
# 8) Output the ranked correlations
# --------------------------------------------------------------
print(ranked_results)
write.csv(ranked_results, file = file.path(dir_out, "20250224_TCGA_TF_Correlation_by_Motif.csv"), row.names = FALSE)


# --------------------------------------------------------------
# 4) Prepare methylation and expression data
# --------------------------------------------------------------
# Retrieve the methylation matrix (probes as rows, samples as columns)
meth_matrix <- methylation_data
# Retrieve the full expression matrix
exp_matrix <- expr_1_filtered

# Ensure that both matrices share the same samples and are ordered identically.
common_samples <- intersect(colnames(meth_matrix), colnames(exp_matrix))
meth_matrix <- meth_matrix[, common_samples, drop = FALSE]
dim(meth_matrix)
exp_matrix  <- exp_matrix[, common_samples, drop = FALSE]
dim(exp_matrix)

# --------------------------------------------------------------
# 5) Subset the expression matrix to only include the selected TFs
# --------------------------------------------------------------
# tf_summary$top_5percent_TFs contains semicolon-separated TF names.
library(dorothea)

# Load the human dorothea regulon data
data(dorothea_hs)

# Extract a unique list of transcription factors from the regulon data
human_tfs <- unique(dorothea_hs$tf)

# View the first few TFs
head(human_tfs)

#tf_list <- unique(unlist(strsplit(TF$top_5percent_TFs[1], ";")))
tf_list <- human_tfs
exp_matrix_tf <- exp_matrix[rownames(exp_matrix) %in% tf_list, ]

# --------------------------------------------------------------
# 6) For each motif, for each associated probe, calculate the Spearman correlation 
#    between that probe's methylation (across samples) and expression of each selected TF.
# --------------------------------------------------------------
cor_results <- lapply(names(enriched.motif), function(motif) {
  probes <- enriched.motif[[motif]]
  # For each probe, compute correlation with each TF.
  probe_cor_list <- lapply(probes, function(probe) {
    # Extract methylation values for the specific probe (vector across samples)
    meth_values <- meth_matrix[probe, ]
    
    # For each TF (each row in exp_matrix_tf), compute Spearman correlation.
    tf_cor <- apply(exp_matrix_tf, 1, function(tf_expr) {
      test <- cor.test(meth_values, tf_expr, method = "spearman")
      c(correlation = test$estimate, p.value = test$p.value)
    })
    
    # Transpose and convert the result to a data frame.
    tf_cor_df <- as.data.frame(t(tf_cor))
    tf_cor_df$TF <- rownames(tf_cor_df)
    tf_cor_df$probe <- probe
    tf_cor_df$motif <- motif
    tf_cor_df
  })
  do.call(rbind, probe_cor_list)
})
cor_results_df <- do.call(rbind, cor_results)

# --------------------------------------------------------------
# 7) Optionally, rank TFs by the strength of inverse (negative) correlation within each motif-probe combination.
# --------------------------------------------------------------
# For example, you might rank TFs for each motif and probe:
ranked_results <- cor_results_df %>%
  group_by(motif, probe) %>%
  arrange(cor_results_df$correlation.rho) %>%  # More negative correlations (stronger inverse relationships) come first.
  mutate(rank = row_number()) %>%
  ungroup()

# --------------------------------------------------------------
# 8) Output the correlation results
# --------------------------------------------------------------
print(ranked_results)
write.csv(ranked_results, file = file.path(dir_out, "TF_Correlation_by_Motif_Probe.csv"), row.names = FALSE)




# COLOSSUS CNA

CNA_amp <- read.delim("20240719_CNA_all_COLOSSUS_amplified_genes.txt", stringsAsFactors = FALSE)
dim(CNA_amp)

#hypermethylation file
pair_1 <- read.delim("getPair.hyper.pairs.significant_with_correlation_Anguraj_-0.2_correlation__correlation_2000_removed.txt", stringsAsFactors = FALSE)
dim(pair_1)

#hypomethylation file
#pair_1 <- read.delim("getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation_2000_removed.txt", stringsAsFactors = FALSE)
#dim(pair_1)

#m_amp <- match(unique(pair_1$Symbol[w_cor]),unique(CNA_amp$Genes))
m_amp <- match(unique(pair_1$Symbol),unique(CNA_amp$Genes))
w_amp <- which(!is.na(m_amp))
length(w_amp)

CNA_amp$Genes[m_amp[w_amp]]

CNA_amp[m_amp[w_amp],]

CNA_del <- read.delim("20240722_CNA_present_in_COLOSSUS_deleted_genes.txt", stringsAsFactors = FALSE)
dim(CNA_del)

#m_del <- match(unique(pair$Symbol[w_cor]),CNA_del$Genes)
m_del <- match(unique(pair_1$Symbol),unique(CNA_del$Genes))
w_del <- which(!is.na(m_del))
length(w_del)

CNA_del$Genes[m_del[w_del]]

CNA_del[m_del[w_del],]

#  gene lists 
cnv_amplification_genes <- CNA_amp$Genes
cnv_deletion_genes <- CNA_del$Genes
hypermethylation_genes <- pair_1$Symbol

# Create a named list of gene sets
gene_sets <- list(
  "CNV Amplification" = cnv_amplification_genes,
  "CNV Deletion" = cnv_deletion_genes,
  "Hypomethylation" = hypermethylation_genes
)

library(ggVennDiagram)
library(ggplot2)

# Generate the 3-set Venn diagram
#pdf("20251108_TCGA_hypermethylation_cor_-0.2_4183_CNV_methylation_venn_diagram.pdf")

#pdf("3E_20260401_TCGA_hypomethylation_cor_-0.2_4183_CNV_methylation_venn_diagram.pdf")

ggVennDiagram(gene_sets, label = "count", edge_size = 1.2) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Gradient fill
  theme_minimal() +  # Clean theme
  labs(title = "Venn Diagram of CNV Amplification, CNV Deletion & Hypomethylation") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
#dev.off()

#GLS hyper methy
g_g <- which(rownames(expr_2_filtered) == "ENSG00000115419")
#[1] 6278
m_g <- which(rownames(combined_data_2_filtered) == "cg22742965")
#[1] 9103

# Install packages if needed
# install.packages(c("ggplot2", "viridis"))

library(ggplot2)
library(dplyr)

# ----------------------------------------
# 1. Extract numeric vectors
# ----------------------------------------
x <- as.numeric(combined_data_2_filtered[m_g, ])
y <- as.numeric(expr_2_filtered[g_g, ])

df <- data.frame(
  sample = colnames(combined_data_2_filtered),
  x = x,
  y = y
)

# ----------------------------------------
# 2. Merge sample info
# ----------------------------------------
df2 <- df %>%
  inner_join(as.data.frame(sample.info),
             by = c("sample" = "primary"))

df2 <- df2 %>% filter(!is.na(x), !is.na(y))

# ----------------------------------------
# 3. Darker pastel color palette
# ----------------------------------------
dark_pastel_colors <- c(
  "#FF8A5C",  # coral (replaces pink)
  "#FFBF8A",
  "#FFF08A",
  "#8AFFB0",
  "#8AC6FF",
  "#D28AFF",
  "#A6A0FF"
)
# match number of groups
sample_levels <- unique(df2$sample.type)
color_map <- dark_pastel_colors[seq_along(sample_levels)]
names(color_map) <- sample_levels

# correlation
cor_val <- cor(df2$x, df2$y, method = "pearson")

# ----------------------------------------
# 4. Beautiful darker-pastel publication-quality plot
# ----------------------------------------
pdf("20251118_TCGA_hypomethylation_cg22742965_GLS_gene_expression_cor_plot.pdf")
ggplot(df2, aes(x = x, y = y, color = sample.type)) +
  
  # white halo for clarity
  geom_point(size = 4.5, alpha = 0.95, stroke = 0.7, color = "white") +
  geom_point(size = 3.3, alpha = 0.95, aes(color = sample.type)) +
  
  # regression line
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.5, color = "#333333") +
  
  scale_color_manual(values = color_map) +
  
  labs(
    title = "Expression Comparison by Sample Type",
    subtitle = paste0("Pearson correlation: r = ", round(cor_val, 3)),
    x = "cg22742965 methylation profile",
    y = "GLS gene expression",
    color = "Sample Type"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 20),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    axis.title    = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
dev.off()
