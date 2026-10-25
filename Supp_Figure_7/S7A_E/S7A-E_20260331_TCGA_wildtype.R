################################################################################
# Figure Panels : Supplementary Figures 7A, 7B, 7C, 7D, 7E
# Description   : ELMER-based analysis of TCGA COAD RAS-WILDTYPE tumours.
#                 Mirror of the RAS-mutant analysis (Fig 3A–K) applied to the
#                 wildtype arm for comparison.
#
#                 7A — Correlation: PTK2 gene and CpG enhancer methylation
#                      (hypomethylation, RAS-wildtype, 1333 probes / 726 genes)
#                 7B — Correlation: average hypomethylation enhancer regions
#                      and average gene expression (scatter, RAS-wildtype)
#                 7C — Correlation: MYC gene and CpG enhancer methylation
#                      (RAS-wildtype cohort)
#                 7D — GLS (cg02514318): methylation vs gene expression scatter
#                 7E — Association of top 10 hypomethylated probes with MYC
#                      target genes (RAS-wildtype cohort)
#
# Input data    : 20251105_TCGA_normal_solid_tissue_Methylation_full_probes.txt
#                 20251105_TCGA_tumour_solid_tissue_Methylation_full_probes.txt
#                 TCGA RNA-seq expression matrix (RAS-wildtype samples)
#                 tcga_coadread_pancanatlas_RAS_WILDTYPE_MSS_TCGA_barcodes.txt
# Output        : ELMER results in results/next/wildtype/  subfolder
#                 S7A–E scatter/correlation PDFs
# R packages    : ELMER, ELMER.data, MultiAssayExperiment, SummarizedExperiment,
#                 GenomicInteractions, plyr, doParallel, ggplot2, dplyr
# Authors       : Sadanandam Lab, Institute of Cancer Research, London
# Contact       : anguraj.sadanandam@icr.ac.uk
# Date          : March 2026
# Note          : Update all file paths before running. Results directory must
#                 exist at results/next/wildtype/ before ELMER writes output.
################################################################################

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


# Step 1: Load existing methylation data
#normal.db <- read.delim("20250128_TCGA_normal_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
normal.db <- read.delim("20251105_TCGA_normal_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
message("Normal data dimensions: ", paste(dim(normal.db), collapse = " x "))

#tumor_data <- read.delim("20250128_TCGA_primary_tumour_Methylation_full_probes.txt", stringsAsFactors = FALSE)
tumor_data <- read.delim("20251105_TCGA_tumour_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
message("Tumor data dimensions: ", paste(dim(tumor_data), collapse = " x "))

# Step 2: Extract and preprocess probe IDs
message("Processing methylation probes...")

#filter for probes
na_filter <- complete.cases(normal.db) & complete.cases(tumor_data)
normal.db <- normal.db[na_filter, , drop = FALSE]
tumor_data <- tumor_data[na_filter, , drop = FALSE]

# Combine datasets
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
reference_file <- read.delim("tcga_coadread_pancanatlas_RAS_wildtype_MSS_TCGA_barcodes.txt", stringsAsFactors = FALSE)
valid_samples <- reference_file$TCGA_RAS_wt_samples  # Extract valid sample barcodes

#Match and filter the datasets based on valid sample barcodes
valid_normal_tumour_samples <- substr(colnames(combined_data_2),1,12) %in% substr(gsub("-","\\.",valid_samples),1,12)
valid_expression_samples <- substr(colnames(expr_2),1,12) %in% substr(gsub("-","\\.",valid_samples),1,12)

# Filter combined_data_2 to keep only matched samples
combined_data_2_filtered <- combined_data_2[, valid_normal_tumour_samples, drop = FALSE]
dim(combined_data_2_filtered)
# Filter expr_2 to keep only matched samples
expr_2_filtered <- expr_2[, valid_expression_samples, drop = FALSE]
dim(expr_2_filtered)

expr_1_filtered <- expr_1[, valid_expression_samples, drop = FALSE]
dim(expr_1_filtered)

group_labels_2 <-  group_labels_1[valid_normal_tumour_samples, drop = FALSE]

df <- data.frame(cbind(colnames(combined_data_2_filtered),group_labels_2))
colnames(df) <- c("sample ID","groups")
dim(df)

sample.info <- S4Vectors::DataFrame(
  primary =  df$`sample ID`, 
  sample.type = df$groups
)

# Extract DNA methylation and gene expression data
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
methylation_data <- assay(getMet(mae))  # Methylation levels
expression_data <- assay(getExp(mae))   # Gene expression levels

group.col <- "sample.type"
group1 <- "tumor"
group2 <-  "normal"
dir.out <- "result"
#diff.dir <-  "hypo" # Search for hypo/hypermethylated probes in group 1
diff.dir <-  "hypo" # uncomment if need to run hypermethylation

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

dim(pair)

## --------------------------------------------------------
## INPUTS YOU ALREADY HAVE IN THE SESSION
## --------------------------------------------------------
## expr_2_filtered:   gene expression matrix (rows = Ensembl IDs, cols = samples)
## combined_data_2_filtered: methylation matrix (rows = CpG probes, cols = samples)
## sample.info:       data.frame with sample information
##   - must contain a column "primary" with sample IDs
##   - any other columns (e.g. sample.type) are fine but not required

## --------------------------------------------------------
## 0. Load packages
## --------------------------------------------------------
library(dplyr)

## --------------------------------------------------------
## 1. Read the CpG–gene list
##    File should have columns: Probe, GeneID, Symbol, Distance, Sides, Raw.p, FDR, Correlation
##    (extra columns are okay; they’ll just be carried along)
## --------------------------------------------------------

pairs_file <- "getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation_2000_removed.txt"  # <-- change to your filename

pairs_df <- read.delim(
  pairs_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

## Optional: check that the key columns exist
required_cols <- c("Probe", "GeneID")
stopifnot(all(required_cols %in% colnames(pairs_df)))

## --------------------------------------------------------
## 2. Helper: compute correlation for ONE CpG–gene pair
## --------------------------------------------------------

compute_pair_cor <- function(probe_id, gene_id, pairs_row, 
                             meth_mat, expr_mat, sample_info) {
  # Check presence
  if (!probe_id %in% rownames(meth_mat)) {
    return(
      data.frame(
        Probe      = probe_id,
        GeneID     = gene_id,
        Symbol     = pairs_row[["Symbol"]],
        Distance   = pairs_row[["Distance"]],
        Sides      = pairs_row[["Sides"]],
        cor        = NA_real_,
        p.value    = NA_real_,
        n.samples  = 0,
        note       = "Probe_not_found"
      )
    )
  }
  
  if (!gene_id %in% rownames(expr_mat)) {
    return(
      data.frame(
        Probe      = probe_id,
        GeneID     = gene_id,
        Symbol     = pairs_row[["Symbol"]],
        Distance   = pairs_row[["Distance"]],
        Sides      = pairs_row[["Sides"]],
        cor        = NA_real_,
        p.value    = NA_real_,
        n.samples  = 0,
        note       = "Gene_not_found"
      )
    )
  }
  
  # Extract numeric vectors (like in your example)
  x <- as.numeric(meth_mat[probe_id, ])
  y <- as.numeric(expr_mat[gene_id, ])
  
  df <- data.frame(
    sample = colnames(meth_mat),
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
  
  # Merge with sample.info (same as your code)
  df2 <- df %>%
    inner_join(as.data.frame(sample_info), 
               by = c("sample" = "primary")) %>%
    filter(!is.na(x), !is.na(y))
  
  n_use <- nrow(df2)
  
  if (n_use < 3) {
    # Not enough samples for correlation
    return(
      data.frame(
        Probe      = probe_id,
        GeneID     = gene_id,
        Symbol     = pairs_row[["Symbol"]],
        Distance   = pairs_row[["Distance"]],
        Sides      = pairs_row[["Sides"]],
        cor        = NA_real_,
        p.value    = NA_real_,
        n.samples  = n_use,
        note       = "Too_few_samples"
      )
    )
  }
  
  # Spearman correlation
  ct <- suppressWarnings(
    cor.test(df2$x, df2$y, method = "spearman")
  )
  
  # Spearman correlation
  cp <- suppressWarnings(
    cor.test(df2$x, df2$y, method = "pearson")
  )
  
  data.frame(
    Probe      = probe_id,
    GeneID     = gene_id,
    Symbol     = pairs_row[["Symbol"]],
    Distance   = pairs_row[["Distance"]],
    Sides      = pairs_row[["Sides"]],
    cor        = unname(ct$estimate),
    p.value    = ct$p.value,
    Correlation = unname(cp$estimate),
    n.samples  = n_use,
    note       = "OK",
    stringsAsFactors = FALSE
  )
}

## --------------------------------------------------------
## 3. Loop over all rows in the pairs table
## --------------------------------------------------------

results_list <- lapply(seq_len(nrow(pairs_df)), function(i) {
  row_i <- pairs_df[i, ]
  compute_pair_cor(
    probe_id   = row_i[["Probe"]],
    gene_id    = row_i[["GeneID"]],
    pairs_row  = row_i,
    meth_mat   = combined_data_2_filtered,
    expr_mat   = expr_2_filtered,
    sample_info = sample.info
  )
})

results_df <- bind_rows(results_list)

## --------------------------------------------------------
## 4. (Optional) adjust p-values and save
## --------------------------------------------------------

results_df$FDR.new <- p.adjust(results_df$p.value, method = "BH")

length(which(results_df$Correlation < -0.2))

wc <- which(results_df$Correlation < -0.2)

results_df_c <- results_df[wc,]
dim(results_df_c)

output_file <- "getPair.hypo.pairs.significant_with_correlation_Anguraj_corr_-0.2_hypo_methy_cor_results.txt"
write.table(
  results_df_c,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

results_df[1:5, ]


## change these to your actual filenames / paths
file1 <- "getMethdiff.hypo.probes.significant.csv"  # first table (probe, pvalue, tumor_Minus_normal, adjust.p)
file2 <- "getPair.hypo.pairs.significant_with_correlation_Anguraj_corr_-0.2_hypo_methy_cor_results.txt"  # second table (Probe, GeneID, Symbol, ...)

## read the tables (tab-separated)
df1 <- read.table(file1,
                  header = TRUE,
                  sep = ",",
                  stringsAsFactors = FALSE,
                  check.names = FALSE)
dim(df1)

df2 <- read.table(file2,
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  check.names = FALSE)

#df2 <- df2[df2$FDR < 0.05,]
dim(df2)

mu <- match(df2$Probe,df1$probe)
wu <- which(!is.na(mu))
length(wu)

## merge by probe (note different column name capitalization)
merged <- cbind(df2[wu,],df1[mu[wu],])

## look at the result
head(merged)
dim(merged)

## (optional) write merged table to a file
#write.table(merged,
#           file = "getPair.hypo.pairs.significant_with_correlation_Anguraj_corr_-0.2_hypo_methy_cor_results_merged_by_probe.txt",
#          sep = "\t",
#         quote = FALSE,
#        row.names = FALSE)


mer_ptk2 <- which(merged$Symbol == "PTK2")
mer_myc <- which(merged$Symbol == "MYC")
mer_cdk4 <- which(merged$Symbol == "CDK4")
mer_gls <- which(merged$Symbol == "GLS")
mer_nfkbia <- which(merged$Symbol == "NFKBIA")

#average expr vs methyl analysis
ave_zscore_methy <- apply(combined_data_2_filtered[unique(merged$Probe),],2,mean)
ave_zscore_expr <- apply(expr_2_filtered[unique(merged$GeneID),],2,mean)

fig3cd_data <- data.frame(cbind(scale(ave_zscore_expr),scale(ave_zscore_methy)))
names(fig3cd_data) = c("expression", "methylation")

fig3c_lm = coef(lm(formula = fig3cd_data$expression~fig3cd_data$methylation, data=fig3cd_data))

cor.test(fig3cd_data$expression,fig3cd_data$methylation)

plot(fig3cd_data$methylation,fig3cd_data$expression)
#-0.7130365 

#average z-score methylation and expression
library(ggplot2)
publication_theme = theme(
  title = element_text(size = 20),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14)
)

#pdf("20260403_TCGA_RAS_wildtype_hyomethylation_cor_-0.2_average_zscore_RNA_methylation.pdf")
ggplot(fig3cd_data, aes(y=expression, x=methylation, colour = "#F98400")) +
  geom_abline(intercept = fig3c_lm[1], slope=fig3c_lm[2], lty=2, lwd=1, col = "#F98400")+
  geom_point(size=3, alpha=1)+
  scale_color_manual(values = c("#F98400"))+
  theme_bw()+
  labs(y="Average RNA expression Z-score", x="Average Methylation Z-score")+
  publication_theme
#dev.off()

#PTK2 analysis
# Filter the 'pair' data frame for rows with Symbol equal to "PTK2"
ptk2_pair <- subset(merged, Symbol == "PTK2")
dim(ptk2_pair)

# Check the filtered results
print(ptk2_pair)
dim(ptk2_pair)


# Load libraries
library(ggplot2)
#  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = mean(data$Corr)) +

# Plot with color gradient based on correlation
pdf("20260403_TCGA_RAS_wildtype_PTK2_methylation_11_highly_correlated_plot.pdf")
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
dev.off()

#MYC analysis
# Filter the 'pair' data frame for rows with Symbol equal to "MYC"
myc_pair <- subset(merged, Symbol == "MYC")
dim(myc_pair)

# Check the filtered results
print(myc_pair)
dim(myc_pair)


# Load libraries
library(ggplot2)
#  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = mean(data$Corr)) +

# Plot with color gradient based on correlation
#pdf("20260403_TCGA_RAS_wildtype_MYC_methylation_5_highly_correlated_plot.pdf")
ggplot(myc_pair, aes(x = Distance, y = Correlation, label = Probe, color = "blue")) +
  geom_point(size = 4) +
  geom_segment(data = myc_pair,
               aes(x = Distance, xend = Distance, y = 0, yend = Correlation),
               color = "gray50")+
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

#validation
cor.test(methylation_data["cg08526705",],unlist(expr_1_filtered["MYC",]), method = "pearson")
cor.test(methylation_data["cg23399286",],unlist(expr_1_filtered["MYC",]), method = "pearson")
cor.test(methylation_data["cg13394191",],unlist(expr_1_filtered["MYC",]), method = "pearson")


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
pair_file <- "getPair.hypo.pairs.significant_with_correlation_Anguraj_corr_-0.2_hypo_methy_cor_results.txt"

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
  filename = "S7J_20260403_TCGA_wt_hyomethylation_cor_-0.2_4183_average_Myc_targets_Pearson_correlation_top_10.pdf",
  plot = p,
  width = 8,
  height = 5.5
)

# Optional: print selected rows
print(top10_myc)


# ============================================================
# Publication-quality scatter plot for probe–gene association
# Example: pair_df[69, ]
# Assumes:
#   1. pair_df exists
#   2. methylation_data = methylation matrix (rows = probes, cols = samples)
#   3. expression_data = expression matrix (rows = genes, cols = samples)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(ggpubr)
  library(scales)
})

# ------------------------------------------------------------
# 1. Select the probe-gene pair
# ------------------------------------------------------------
pair_row <- 69
pair_sel <- pair_df[pair_row, ]

probe_id <- pair_sel$Probe
gene_id  <- pair_sel$GeneID
gene_sym <- pair_sel$Symbol

# ------------------------------------------------------------
# 2. Check rownames and handle Ensembl version suffix if needed
# ------------------------------------------------------------
expr_gene_ids <- rownames(expression_data)

# direct match first
gene_match <- which(expr_gene_ids == gene_id)

# if not found, try removing Ensembl version suffix
if (length(gene_match) == 0) {
  expr_gene_ids_clean <- sub("\\..*$", "", expr_gene_ids)
  gene_match <- which(expr_gene_ids_clean == gene_id)
}

if (!(probe_id %in% rownames(methylation_data))) {
  stop(paste("Probe not found inmethylation_data:", probe_id))
}

if (length(gene_match) == 0) {
  stop(paste("Gene not found in expr_mat2:", gene_id))
}

if (length(gene_match) > 1) {
  warning("Multiple matching genes found; using the first match.")
  gene_match <- gene_match[1]
}

# ------------------------------------------------------------
# 3. Extract methylation and expression values
# ------------------------------------------------------------
meth_vec <- as.numeric(methylation_data[probe_id, ])
names(meth_vec) <- colnames(methylation_data)

expr_vec <- as.numeric(expression_data[gene_match, ])
names(expr_vec) <- colnames(expression_data)

# ------------------------------------------------------------
# 4. Match common samples
# ------------------------------------------------------------
common_samples <- intersect(names(meth_vec), names(expr_vec))

if (length(common_samples) < 3) {
  stop("Too few overlapping samples between methylation and expression matrices.")
}

plot_df <- data.frame(
  Sample       = common_samples,
  Methylation  = meth_vec[common_samples],
  Expression   = expr_vec[common_samples],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Methylation), !is.na(Expression))

# ------------------------------------------------------------
# 5. Compute Pearson correlation
# ------------------------------------------------------------
ct <- cor.test(plot_df$Methylation, plot_df$Expression, method = "pearson")

r_val <- unname(ct$estimate)
p_val <- ct$p.value
n_val <- nrow(plot_df)

label_text <- paste0(
  "Pearson r = ", sprintf("%.3f", r_val),
  "\nP = ", format.pval(p_val, digits = 2, eps = 1e-4),
  "\nn = ", n_val
)

# ------------------------------------------------------------
# 6. Build plot labels
# ------------------------------------------------------------
main_title <- paste0(gene_sym, " expression vs ", probe_id, " methylation")
sub_title  <- paste0("Gene: ", gene_id,
                     " | Distance: ", pair_sel$Distance,
                     " | Side: ", pair_sel$Sides)

x_lab <- paste0(probe_id, " methylation")
y_lab <- paste0(gene_sym, " expression")

# ------------------------------------------------------------
# 7. Publication-style theme
# ------------------------------------------------------------
pub_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.line = element_line(linewidth = 0.6),
    plot.margin = margin(10, 12, 10, 10)
  )

# ------------------------------------------------------------
# 8. Draw plot
# ------------------------------------------------------------
p <- ggplot(plot_df, aes(x = Methylation, y = Expression)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1.05, vjust = 1.2,
    size = 4.2
  ) +
  labs(
    title = main_title,
    subtitle = sub_title,
    x = x_lab,
    y = y_lab
  ) +
  pub_theme

print(p)

# ------------------------------------------------------------
# 9. Save high-resolution outputs
# ------------------------------------------------------------
ggsave(
  filename = paste0("scatter_", gene_sym, "_", probe_id, ".pdf"),
  plot = p,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = paste0("scatter_", gene_sym, "_", probe_id, ".png"),
  plot = p,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 600
)

#MYC plot
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

# -----------------------------
# Input file
# -----------------------------
# Input file
pair_file <- "getPair.hypo.pairs.significant_with_correlation_Anguraj_corr_-0.2_hypo_methy_cor_results.txt"

# Output file
out_pdf <- "Top10_MYC_targets_probe_gene_barplot_reversed.pdf"

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

