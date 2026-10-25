################################################################################
# Figure Panels : 3B, 3N, 3O
# Description   : ELMER analysis — HYPERMETHYLATION arm.
#                 Identifies CpG probes hypermethylated in RAS-mutant tumours
#                 vs normal colon, links them to downregulated nearby genes,
#                 and characterises enriched TF motifs at hypermethylated loci.
#                 Mirror analysis to the hypomethylation script (3A_B_C_E).
# Input data    : 20251105_TCGA_normal_solid_tissue_Methylation_full_probes.txt
#                 20251105_TCGA_tumour_solid_tissue_Methylation_full_probes.txt
#                 TCGA RNA-seq expression matrix
#                 tcga_coadread_pancanatlas_RAS_mutant_MSS_TCGA_barcodes.txt
# Output        : ELMER results in results/next/hyper/ subfolder
# R packages    : ELMER, ELMER.data, MultiAssayExperiment, SummarizedExperiment,
#                 GenomicInteractions, plyr, doParallel
# Authors       : Sadanandam Lab, Institute of Cancer Research, London
# Contact       : anguraj.sadanandam@icr.ac.uk
# Date          : November 2025 (v3.0)
# Note          : Run after the hypomethylation script. Results from both arms
#                 are combined for downstream Fig 3B comparison.
################################################################################

# ELMER hypermethylation arm — same pipeline as the hypo script
# but targets CpG probes with HIGHER methylation in tumour vs normal.
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
# Step 1: Load existing methylation data
#normal.db <- read.delim("20250128_TCGA_normal_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
normal.db <- read.delim("20251105_TCGA_normal_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
message("Normal data dimensions: ", paste(dim(normal.db), collapse = " x "))

#tumor_data <- read.delim("20250128_TCGA_primary_tumour_Methylation_full_probes.txt", stringsAsFactors = FALSE)
tumor_data <- read.delim("20251105_TCGA_tumour_solid_tissue_Methylation_full_probes.txt", stringsAsFactors = FALSE)
message("Tumor data dimensions: ", paste(dim(tumor_data), collapse = " x "))

# ── Step 2: Filter and combine matrices ────────────────────────────────────────
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
reference_file <- read.delim("tcga_coadread_pancanatlas_RAS_mutant_MSS_TCGA_barcodes.txt", stringsAsFactors = FALSE)
valid_samples <- reference_file$TCGA_RAS_mutant_samples  # Extract valid sample barcodes

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
#diff.dir <-  "hypo" # Search for hypo/hypermethylated probes in group 1
diff.dir <-  "hyper" # uncomment if need to run hypermethylation

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
res <- read.csv("getMethdiff.hyper.probes.csv", stringsAsFactors = FALSE)

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
pdf("volcanoPlot.probes.hyper.pdf", width = 6, height = 5)

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

dev.off()



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

# View results with correlation values
head(pair[, c("Probe", "GeneID", "Distance", "Correlation")])

write.table(pair,"result/getPair.hyper.pairs.significant_with_correlation_Anguraj.txt", sep="\t", quote = FALSE)
write.table(pair[w_cor,],"result/getPair.hyper.pairs.significant_with_correlation_Anguraj_-0.2_correlation.txt", sep="\t", quote = FALSE)


#average z-score methylation and expression
library(ggplot2)
publication_theme = theme(
  title = element_text(size = 20),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14)
)

ave_zscore_methy <- apply(methylation_data[pair$Probe[w_cor],],2,mean)
ave_zscore_expr <- apply(expression_data[pair$GeneID[w_cor],],2,mean)

fig3cd_data <- data.frame(cbind(scale(ave_zscore_expr),scale(ave_zscore_methy)))
names(fig3cd_data) = c("expression", "methylation")

fig3c_lm = coef(lm(formula = expression~methylation, data=fig3cd_data))

cor.test(fig3cd_data$expression,fig3cd_data$methylation)

plot(fig3cd_data$methylation,fig3cd_data$expression)
#-0.7130365 


pdf("20260303_TCGA_hypermethylation_cor_-0.2_915_average_zscore_RNA_methylation.pdf")
ggplot(fig3cd_data, aes(y=expression, x=methylation, colour = "#F98400")) +
  geom_abline(intercept = fig3c_lm[1], slope=fig3c_lm[2], lty=2, lwd=1, col = "#F98400")+
  geom_point(size=3, alpha=1)+
  scale_color_manual(values = c("#F98400"))+
  theme_bw()+
  labs(y="Average RNA expression Z-score", x="Average Methylation Z-score")+
  publication_theme
dev.off()


# Enrichment Hallmarks 

fig3kl_hyper_elsevier = read.table("20250224_TCGA_hypomethylation_corr_0.2_2881_genes_Enrichr_results_bar_MsigDB_Hallmark_2020_table.txt", sep = "\t", header = T, quote = "\"")

fig3l_df = head(fig3kl_hyper_elsevier,6)
fig3l_df$Term = factor(fig3l_df$Term, levels = rev(fig3l_df$Term))

pdf("20250224_TCGA_hyomethylation_cor_-0.2_2881_average_MsigDB_Hallmarks.pdf")

ggplot(data = fig3l_df, aes(y = Term, x = -log10(Adjusted.P.value), fill = Odds.Ratio)) +
  geom_bar(stat = "identity") +  # Removed show.legend = F so the legend is displayed
  geom_text(aes(label = Term), x = 0.05, hjust = 0, vjust = 0.5, size = 4) +
  scale_fill_gradient(low = "#EEEEEE", high = "#8866FF", na.value = NA, limits = c(0, 50)) +
  labs(title = "MSigDB Hallmarks", fill = "Odds Ratio") +  # Adds legend title
  theme_bw() +
  publication_theme +
  theme(axis.text.y = element_blank())

dev.off()

# metabolic work bench

fig3_hyper_elsevier = read.table("S7E_Metabolomics_Workbench_Metabolites_2022_table.txt", sep = "\t", header = T, quote = "\"")

fig3_df = head(fig3_hyper_elsevier,8)
fig3_df$Term = factor(fig3_df$Term, levels = rev(fig3_df$Term))

pdf("S7E_20260408_TCGA_hyomethylation_cor_-0.2_4140_average_Metabolomics_Workbench_Metabolites_2022.pdf")

ggplot(data=fig3_df, aes(y=Term, x=-log10(Adjusted.P.value), fill = Odds.Ratio))+
  geom_bar(stat = "identity", show.legend = T)+
  geom_text(data = fig3_df, aes(label=Term), x = 0.05, hjust=0, vjust=0.5, size=4)+
  scale_fill_gradient(low = "#EEEEEE", high = "#8866FF", na.value = NA)+
  labs(title="Metabolomics_Workbench_Metabolites_2022")+
  theme_bw()+
  publication_theme+
  theme(axis.text.y = element_blank())

dev.off()



#MYC target - lm

myc_targets <- c("CDK4","MYC","PLK4","NOP56","NOLC1","SYNCRIP","MCM2","PA2G4","MCM4","HSPE1","HSPD1","NME1", "PLK1")
expr_1_filtered_myc <- expr_1_filtered[myc_targets,]
dim(expr_1_filtered_myc)
# Subset the pair data frame to include only rows with Symbol in myc_targets
pair_myc <- pair[pair$Symbol %in% myc_targets, ]
dim(pair_myc)

methylation_data_myc <- methylation_data[pair_myc$Probe,]
dim(methylation_data_myc )

# - methylation_data: matrix/data frame with rownames as probe IDs and columns as samples.
# - expr_1_filtered_myc: matrix/data frame of gene expression with rownames as gene symbols.
# - pair_myc: data frame containing the mapping with columns "Probe" and "Symbol".
# - myc_targets: a character vector of target gene symbols.

# 1. Subset methylation_data to only include the probes in pair_myc
selected_probes <- unique(pair_myc$Probe)
meth_myc <- methylation_data[rownames(methylation_data) %in% selected_probes, ]

# 2. Prepare matrices for storing estimates and p-values.
#    Rows: probes, Columns: genes (MYC targets)
unique_probes <- unique(pair_myc$Probe)
unique_genes  <- myc_targets  # ensure these match the rownames of expr_1_filtered_myc

estimate_matrix <- matrix(NA, nrow = length(unique_probes), ncol = length(unique_genes))
pvalue_matrix   <- matrix(NA, nrow = length(unique_probes), ncol = length(unique_genes))
rownames(estimate_matrix) <- unique_probes
colnames(estimate_matrix) <- unique_genes
rownames(pvalue_matrix)   <- unique_probes
colnames(pvalue_matrix)   <- unique_genes

# 3. Loop through each probe and its associated gene(s) to perform linear regression.
for (probe in unique_probes) {
  # Identify gene(s) associated with this probe using pair_myc.
  associated_genes <- unique(pair_myc$Symbol[pair_myc$Probe == probe])
  
  for (gene in associated_genes) {
    # Check that the gene is present in the expression data.
    if (!(gene %in% rownames(expr_1_filtered_myc))) next
    
    # Extract the predictor values (gene expression) for the gene.
    expr_values <- as.numeric(expr_1_filtered_myc[gene, ])
    
    # Ensure the probe exists in the methylation dataset.
    if (probe %in% rownames(meth_myc)) {
      # Extract the response values (methylation) for the probe.
      meth_values <- as.numeric(meth_myc[probe, ])
      
      # Perform the linear regression: methylation ~ gene expression
      lm_model <- lm(meth_values ~ expr_values)
      model_summary <- summary(lm_model)
      
      # Extract the estimate and p-value for the gene expression coefficient.
      coef_estimate <- model_summary$coefficients["expr_values", "Estimate"]
      p_val         <- model_summary$coefficients["expr_values", "Pr(>|t|)"]
      
      # Store the results in the corresponding cell of the matrices.
      estimate_matrix[probe, gene] <- coef_estimate
      pvalue_matrix[probe, gene]   <- p_val
    }
  }
}

# Optional: Combine the two matrices into a list for easier access.
results <- list(Estimates = estimate_matrix, PValues = pvalue_matrix)

# Load necessary package
library(ggplot2)

# Step 1: Flatten the matrices into a data frame.
# Assume estimate_matrix and pvalue_matrix have rownames (probes) and colnames (genes).
results_df <- data.frame(
  Probe = rep(rownames(estimate_matrix), times = ncol(estimate_matrix)),
  Gene = rep(colnames(estimate_matrix), each = nrow(estimate_matrix)),
  Estimate = as.vector(estimate_matrix),
  PValue = as.vector(pvalue_matrix)
)

# Remove any rows with NA values (i.e. where no regression was performed)
results_df <- na.omit(results_df)

# Step 2: Perform multiple testing correction using the Benjamini-Hochberg method.
results_df$AdjPValue <- p.adjust(results_df$PValue, method = "BH")

# Step 3: Filter for significant associations. Here we use a threshold of 0.05.
sig_results <- subset(results_df, AdjPValue < 0.05)

# Step 4: Order by adjusted p-value and select the top 10 associations.
sig_results <- sig_results[order(sig_results$AdjPValue), ]
top10 <- head(sig_results, 10)

# Optionally, create a label combining probe and gene for plotting.
top10$Pair <- paste(top10$Probe, top10$Gene, sep = " - ")

# Step 5: Create a bar plot for the top 10 associations.

pdf("20250224_TCGA_hyomethylation_cor_-0.2_1524_average_Myc_targets_lm_top_10.pdf")
ggplot(top10, aes(x = reorder(Pair, -Estimate), y = Estimate)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for better readability
  labs(title = "Top 10 Significant Probe-Gene Associations",
       x = "Probe - Gene Pair",
       y = "Regression Estimate (Slope)") +
  theme_minimal()
dev.off()

# Optionally, you can print the top10 data frame to inspect the numbers.
print(top10)

# View dimensions and a snippet of the results
dim(estimate_matrix)
dim(pvalue_matrix)
head(estimate_matrix)
head(pvalue_matrix)

# Optionally, view the regression results for each gene-probe pair:
regression_results



# Select a significant probe-gene pair
selected_probe <- pair$Probe[1]
selected_gene <- pair$GeneID[1]

# Extract values
meth_values <- methylation_data[selected_probe, ]
exp_values <- expression_data[selected_gene, ]

# Compute correlation
cor_value <- cor(meth_values, exp_values, method = "pearson")



# Plot the scatter plot with correlation value
library(ggplot2)
ggplot(data.frame(Methylation = meth_values, Expression = exp_values), 
       aes(x = Methylation, y = Expression)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  ggtitle(paste("Correlation =", round(cor_value, 3))) +
  xlab(paste("Methylation at", selected_probe)) +
  ylab(paste("Expression of", selected_gene))

# scatter plots

scatter.plot(data = mae,
             byProbe = list(probe = c("cg13335974"), numFlankingGenes = 20), 
             category = group.col, correlation = TRUE,
             lm = TRUE, # Draw linear regression curve
             save = TRUE) 

scatter.plot(data = mae,
             byPair = list(probe = c("cg16182447"), gene = c("ENSG00000184828")), 
             category = group.col, , correlation = TRUE, save = TRUE, lm_line = TRUE) 


# Identify enriched motif for significantly hypomethylated probes which 
# have putative target genes.
enriched.motif <- get.enriched.motif(
  data = mae,
  probes = pair$Probe[w_cor], 
  dir.out = dir.out, 
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

gee <- which(pair$Symbol == "MYC")
gee

schematic.plot(
  pair = pair, 
  data = mae,
  group.col = group.col,
  byProbe = pair$Probe[gee[1]],
  numFlankingGenes = 20,
  save = FALSE
)

schematic.plot(
  pair = pair, 
  data = mae,   
  group.col = group.col, 
  byGene = pair$GeneID[510],
  save = FALSE
)

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


pair_results <- read.table("getPair.hyper.pairs.significant_with_correlation_Anguraj_0.2_correlated.txt")

heatmapPairs(data = mae, 
             group.col = group.col,
             group1 = group1, 
             group2 = group2,
             pairs = pair[w_cor,],
           #  filename = NULL)

        filename =  "20251119_TCGA_getPair.hyper.pairs.significant_with_correlation_Anguraj_0.2_correlated_heatmap_methylation_genes.pdf")


# Load required libraries
library(ggplot2)
library(dplyr)

# Read the CSV file (modify the file path)
df <- read.csv("result/getMotif.hypo.motif.enrichment.csv")

# Filter motifs where OR > 2.5
df_filtered <- df %>%
  filter(OR > 3) %>%
  arrange(desc(OR))  # Sort by OR in descending order

# Convert motif names into a factor for ordered plotting
df_filtered$motif <- factor(df_filtered$motif, levels = df_filtered$motif[order(df_filtered$OR, decreasing = TRUE)])

# Compute -log10(FDR) for visualization
df_filtered$logFDR <- -log10(df_filtered$FDR)

# Define a pastel color palette
pastel_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#E6BAFF", "#FFC3E1")

pdf(file ="20250224_TCGA_Enhancer_Motif_Enrichment_OR_G_3.pdf")
# Create the bar plot with pastel rainbow colors
ggplot(df_filtered, aes(x = motif, y = OR, fill = logFDR)) +
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
    return(meth_matrix[, probes])
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




# COLOSSUS probes

pair_results <- read.delim("getPair.hypo.pairs.significant_with_correlation_Anguraj_-0.2_correlation_2000_removed.txt")

hari_rna <- read.csv("2025_12_02_significant_correlated_RNA_Methylation_genes.csv")
dim(hari_rna)

mr <- match(hari_rna$RNA,pair_results$Symbol)
wr <- which(!is.na(mr))
length(wr)

# Identify enriched motif for significantly hypomethylated probes which 
# have putative target genes.
enriched.motif_COLOSSUS <- get.enriched.motif(
  data = mae,
  probes = pair_results$Probe[mr[wr]], 
  dir.out = dir.out, 
  label = diff.dir,
  min.incidence = 3,
  lower.OR = 1.1,
  pvalue = 0.2
)

# Load required libraries
library(ggplot2)
library(dplyr)

# Read the CSV file (modify the file path)
df <- read.csv("getMotif.hypo.motif.enrichment.csv")

# Filter motifs where OR > 2.5
df_filtered <- df %>%
  filter(OR > 2.8, FDR < 0.2) %>%
  arrange(desc(OR))  # Sort by OR in descending order

# Convert motif names into a factor for ordered plotting
df_filtered$motif <- factor(df_filtered$motif, levels = df_filtered$motif[order(df_filtered$OR, decreasing = TRUE)])

# Compute -log10(FDR) for visualization
df_filtered$logFDR <- -log10(df_filtered$FDR)

# Define a pastel color palette
pastel_colors <- c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#E6BAFF", "#FFC3E1")

pdf(file ="20250110_TCGA_Enhancer_Motif_Enrichment_OR_G_2.8.pdf")
# Create the bar plot with pastel rainbow colors
ggplot(df_filtered, aes(x = motif, y = OR, fill = logFDR)) +
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

TF <- get.TFs(
  data = mae, 
  mode = "supervised",
  group.col = group.col,
  group1 = group1,
  group2 = group2,
  enriched.motif = enriched.motif_COLOSSUS,
  dir.out = dir.out,
  cores = 1, 
  diff.dir = diff.dir,
  correlation = "positive",
  topTFper = 0.01,
  save.plots = TRUE,
  label = diff.dir
)

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
  "Hypermethylation" = hypermethylation_genes
)

library(ggVennDiagram)
library(ggplot2)

# Generate the 3-set Venn diagram
#pdf("3E_20260401_TCGA_hypermethylation_cor_-0.2_915_CNV_methylation_venn_diagram.pdf")
ggVennDiagram(gene_sets, label = "count", edge_size = 1.2) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Gradient fill
  theme_minimal() +  # Clean theme
  labs(title = "Venn Diagram of CNV Amplification, CNV Deletion & Hypermethylation") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
#dev.off()

#NKBIA hyper methy
nfk_g <- which(rownames(expr_2_filtered) == "ENSG00000100906")

nfk_m <- which(rownames(combined_data_2_filtered) == "cg05134015")


# Install packages if needed
# install.packages(c("ggplot2", "viridis"))

library(ggplot2)
library(dplyr)

# ----------------------------------------
# 1. Extract numeric vectors
# ----------------------------------------
x <- as.numeric(combined_data_2_filtered[nfk_m, ])
y <- as.numeric(expr_2_filtered[nfk_g, ])

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
pdf("2N_20260401_TCGA_hypermethylation_cg05134015_NFKBIA_gene_expression_cor_plot.pdf")
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
    x = "cg05134015 methylation profile",
    y = "NFKBIA gene expression",
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


# violin plot of 2 NFKBIA vs methylation probes

library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------
# 1. Probes of interest
# ----------------------------------------
probes_use <- c("cg13012916", "cg05134015")

probe_idx <- match(probes_use, rownames(combined_data_2_filtered))
print(data.frame(probe = probes_use, idx = probe_idx))

valid_probes <- probes_use[!is.na(probe_idx)]
valid_idx <- probe_idx[!is.na(probe_idx)]

if (length(valid_probes) == 0) {
  stop("None of the requested probes were found in rownames(combined_data_2_filtered).")
}

# ----------------------------------------
# 2. Extract methylation values
# ----------------------------------------
meth_mat <- combined_data_2_filtered[valid_idx, , drop = FALSE]

df <- as.data.frame(t(meth_mat))
df$sample <- rownames(df)

# ----------------------------------------
# 3. Merge with sample info
# ----------------------------------------
sample_info_df <- as.data.frame(sample.info)

df_long <- df %>%
  inner_join(sample_info_df, by = c("sample" = "primary")) %>%
  pivot_longer(
    cols = all_of(valid_probes),
    names_to = "Probe",
    values_to = "methylation"
  ) %>%
  mutate(
    sample.type = tolower(trimws(as.character(sample.type))),
    sample.type = case_when(
      sample.type == "normal" ~ "Normal",
      sample.type == "tumor"  ~ "Tumor",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(methylation), !is.na(sample.type))

print(table(df_long$Probe, df_long$sample.type, useNA = "ifany"))

df_long$Probe <- factor(df_long$Probe, levels = valid_probes)
df_long$sample.type <- factor(df_long$sample.type, levels = c("Normal", "Tumor"))

# ----------------------------------------
# 4. Create x-axis groups
# ----------------------------------------
df_long$group <- factor(
  paste(df_long$Probe, df_long$sample.type, sep = "_"),
  levels = c(
    "cg13012916_Normal", "cg13012916_Tumor",
    "cg05134015_Normal", "cg05134015_Tumor"
  )
)

# ----------------------------------------
# 5. Wilcoxon p-values per probe
# ----------------------------------------
pvals <- df_long %>%
  group_by(Probe) %>%
  group_modify(~{
    dat <- droplevels(.x)
    p <- wilcox.test(methylation ~ sample.type, data = dat)$p.value
    ymax <- max(dat$methylation, na.rm = TRUE)
    lab <- case_when(
      p <= 0.0001 ~ "****",
      p <= 0.001  ~ "***",
      p <= 0.01   ~ "**",
      p <= 0.05   ~ "*",
      TRUE ~ "ns"
    )
    data.frame(
      p = p,
      label = lab,
      y = ymax + 0.07
    )
  }) %>%
  ungroup()

pvals$x1 <- c(1, 3)[seq_len(nrow(pvals))]
pvals$x2 <- c(2, 4)[seq_len(nrow(pvals))]
pvals$xmid <- (pvals$x1 + pvals$x2) / 2

print(pvals)

# ----------------------------------------
# 6. Colors
# ----------------------------------------
cols <- c(
  "Normal" = "#2B6CB0",
  "Tumor"  = "#E07A4F"
)

# ----------------------------------------
# 7. Plot
# ----------------------------------------
pdf("NFKBIA_both_probes_side_by_side_single_plot.pdf", width = 8, height = 5.8)

ggplot(df_long, aes(x = group, y = methylation, fill = sample.type, color = sample.type)) +
  geom_violin(alpha = 0.18, linewidth = 1.0, trim = FALSE, width = 0.9) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.30,
    linewidth = 0.5,
    color = "grey30"
  ) +
  geom_jitter(width = 0.08, size = 2.0, alpha = 0.65) +
  geom_segment(
    data = pvals,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.6, color = "black"
  ) +
  geom_segment(
    data = pvals,
    aes(x = x1, xend = x1, y = y - 0.02, yend = y),
    inherit.aes = FALSE, linewidth = 0.6, color = "black"
  ) +
  geom_segment(
    data = pvals,
    aes(x = x2, xend = x2, y = y - 0.02, yend = y),
    inherit.aes = FALSE, linewidth = 0.6, color = "black"
  ) +
  geom_text(
    data = pvals,
    aes(x = xmid, y = y + 0.02, label = label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  scale_x_discrete(labels = c(
    "cg13012916_Normal" = "Normal",
    "cg13012916_Tumor"  = "Tumor",
    "cg05134015_Normal" = "Normal",
    "cg05134015_Tumor"  = "Tumor"
  )) +
  labs(
    title = "NFKBIA enhancer hypermethylation across two CpG probes",
    x = NULL,
    y = expression("DNA methylation ("*beta*"-value)")
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7)
  )

dev.off()