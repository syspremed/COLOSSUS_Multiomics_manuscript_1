################################################################################
# Figure Panel : 1F
# Description  : Association between significant CNA peaks and expression of genes within each
#                peak (COLOSSUS retrospective cohort, VST/ComBat-corrected RNA-seq), per cytoband.
# Input data   : 20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt
#                COLOSUSS_IDs_Matched.xlsx
#                all_lesions.conf_5.txt
#                amp_genes.conf_5.txt
#                del_genes.conf_5.txt
# Output       : 20260506_COLOSSUS_Amp_sig_genes_list.txt
#                1G_new_20260506_COLOSSUS_Figure_2F_amplification_alone.txt
#                1G_new_20260506_COLOSSUS_Figure_2F_deletion_alone.txt
#                20260507_COLOSSUS_Figure_2F_amplification_deletion_both.txt
#                1G_new_20260506_CNV_RNA_logistic_regression_plot.pdf
# R packages   : tidyverse, readxl, broom, dplyr, ggplot2, ggbeeswarm, ggrepel
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load libraries
library(tidyverse)
library(readxl)
library(broom)

# Load RNA data
rna <- read.delim(
  "20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt",
  check.names = FALSE
)


# Fix sample names
colnames(rna) <- gsub("\\.", "-", colnames(rna))

# Load clinical map
clinical_map <- read_xlsx(
  "COLOSUSS_IDs_Matched.xlsx"
)

colnames(clinical_map) <- c("COLOSSUS_ID",   "Sample")

# Match RNA samples to clinical map
matching_samples <- intersect(colnames(rna), clinical_map$COLOSSUS_ID)
#rna_matched <- round(rna[, matching_samples],3)
rna_matched <- rna[, matching_samples]
dim(rna_matched)

# ---------------------
# Load and prepare CNV data
# ---------------------

# Load CNV data
cnv_colossus <- read.delim(
  "all_lesions.conf_5.txt",
  check.names = FALSE
)

cnv_colossus <- cnv_colossus[cnv_colossus$`q values` < 0.2,]
dim(cnv_colossus)

# Clean column names
cnv_colossus_clean <- cnv_colossus[1:114, !is.na(colnames(cnv_colossus)) & colnames(cnv_colossus) != ""]

# Extract cytoband mapping
cnv_region_map <- cnv_colossus_clean %>%
  dplyr::select(`Unique Name`, Descriptor) %>%
  dplyr::filter(!is.na(`Unique Name`) & !is.na(Descriptor)) %>%
  dplyr::distinct() %>%
  dplyr::rename(cytoband = Descriptor, peak_id = `Unique Name`)

# Create long-format CNV data
cnv_long <- cnv_colossus_clean %>%
  dplyr::filter(!is.na(`Unique Name`) & `Unique Name` %in% cnv_region_map$peak_id) %>%
  dplyr::select(-Descriptor) %>%
  pivot_longer(
    cols = starts_with("sRCSI"),
    names_to = "CNV_ID",
    values_to = "CNV"
  ) %>%
  dplyr::rename(peak_id = `Unique Name`) %>%
  left_join(cnv_region_map, by = "peak_id") %>%
  dplyr::select(peak_id, cytoband, CNV_ID, CNV)

# ---------------------
# Load and reshape GISTIC gene peaks
# ---------------------

# Load amplification peak genes
amp <- read.delim(
  "amp_genes.conf_5.txt",
  header = TRUE, check.names = FALSE
)

# Clean gene table
#remove "!duplicated(colnames(amp)) &"
amp_dedup <- amp[, colnames(amp) != ""]
gene_rows <- amp_dedup[4:nrow(amp_dedup), -1]
gene_rows <- gene_rows[, !is.na(colnames(gene_rows)) & colnames(gene_rows) != ""]
regions <- colnames(gene_rows)

# Reshape to long format
amp_long <- gene_rows %>%
  dplyr::mutate(row_id = row_number()) %>%
  pivot_longer(
    cols = all_of(regions),
    names_to = "cytoband",
    values_to = "gene"
  ) %>%
  dplyr::filter(gene != "") %>%
  dplyr::mutate(gene = as.character(gene))

# Filter to genes present in RNA
amp_genes_in_rna <- amp_long %>%
  dplyr::filter(gene %in% rownames(rna_matched))

write.table(amp_genes_in_rna,"20260506_COLOSSUS_Amp_sig_genes_list.txt", sep = "\t",quote = FALSE, row.names = FALSE)

# ---------------------
# Compute median expression per cytoband per sample
# ---------------------

median_expr_by_cytoband <- amp_genes_in_rna %>%
  dplyr::group_by(cytoband) %>%
  dplyr::summarise(
    expr_matrix = list({
      genes <- unique(gene)
      sub_rna <- rna_matched[genes, , drop = FALSE]
      apply(sub_rna, 2, median, na.rm = TRUE)
    }),
    .groups = "drop"
  )

# Reshape to long format
median_expr_long <- median_expr_by_cytoband %>%
  unnest_wider(expr_matrix) %>%
  pivot_longer(
    cols = -cytoband,
    names_to = "Sample_ID",
    values_to = "Median_Expression"
  )

#not rounded
#median_expr_long$Median_Expression <- round(median_expr_long$Median_Expression,3)
median_expr_long$Median_Expression <- median_expr_long$Median_Expression

# ---------------------
# Map sample IDs and join with CNV
# ---------------------

# Get mapping from RNA IDs to CNV sample IDs
sample_map <- clinical_map %>%
  dplyr::select(COLOSSUS_ID, CNV_ID = Sample) %>%
  dplyr::filter(!is.na(COLOSSUS_ID) & !is.na(CNV_ID))

# Join CNV ID into expression data
median_expr_annot <- median_expr_long %>%
  left_join(sample_map, by = c("Sample_ID" = "COLOSSUS_ID"))

#remove space from cytoband and select only amplifications
cnv_long_filtered <- cnv_long %>%
  mutate(
    peak_id = str_trim(peak_id),
    cytoband = str_trim(cytoband)
  ) %>%
filter(str_detect(peak_id, "Amplification"))
# filter(str_detect(peak_id, "Deletion"))

head(cnv_long_filtered)

# Final join: expression + CNV by cytoband and CNV_ID
median_expr_with_cnv <- median_expr_annot %>%
  left_join(cnv_long_filtered, by = c("cytoband", "CNV_ID"))

#check
glimpse(median_expr_with_cnv)
table(is.na(median_expr_with_cnv$CNV))

# ---------------------
# Logistic regression: CNV 0 vs 1/2 ~ expression
# ---------------------

# Prepare data
logit_data <- median_expr_with_cnv %>%
  filter(!is.na(CNV) & !is.na(Median_Expression)) %>%
  mutate(CNV_binary = ifelse(CNV == 0, 0, 1))


#test
# Subset to cytoband of interest
subset_rows <- logit_data %>% 
  filter(peak_id == unique(logit_data$peak_id)[20])

# Keep every alternate row (i.e., odd-numbered rows)
subset_alt <- subset_rows[seq(1, nrow(subset_rows), by = 1), ]

# Fit logistic regression
fit <- glm(as.factor(CNV_binary) ~ Median_Expression, data = subset_alt, family = binomial)

# Show summary
summary(fit)

peak <- NULL
cyto <- NULL
estimate <- NULL
pvalue <- NULL

# Subset to cytoband of interest
for(i in 1:length(unique(logit_data$peak_id))){
  subset_rows <- logit_data %>% 
    filter(peak_id == unique(logit_data$peak_id)[i])
  
  # Keep every alternate row (i.e., odd-numbered rows)
  subset_alt <- subset_rows[seq(1, nrow(subset_rows), by = 1), ]
  
  # Fit logistic regression
  fit <- glm(CNV_binary ~ Median_Expression, data = subset_alt, family = binomial)
  
  # Show summary
  s <- summary(fit)
  peak[i] <- unique(logit_data$peak_id)[i]
  cyto[i] <- unique(logit_data$cytoband)[i]
  estimate[i] <- s$coefficients[2,1]
  pvalue[i] <- s$coefficients[2,4]
}

p.adj <- p.adjust(pvalue, method = "fdr")

data <- data.frame(cbind(peak,cyto,estimate,pvalue,p.adj))
colnames(data) <- c("Peak", "Cytoband", "Estimate", "Pvalue", "P.adj")


#write.table(data,"1G_new_20260506_COLOSSUS_Figure_2F_amplification_alone.txt", row.names = FALSE, sep = "\t", quote = FALSE)


# Load deletion gene peaks
del <- read.delim("del_genes.conf_5.txt",
                  header = TRUE, check.names = FALSE)

# Step 1: Clean up column names and remove bad columns
del_clean <- del[, !duplicated(colnames(del)) & colnames(del) != ""]

# Step 2: Remove first three rows (q-values, etc.)
del_gene_rows <- del_clean[4:nrow(del_clean), -1]

# Step 3: Remove unnamed or NA columns again
del_gene_rows <- del_gene_rows[, !is.na(colnames(del_gene_rows)) & colnames(del_gene_rows) != ""]

# Step 4: Get cytobands
del_regions <- colnames(del_gene_rows)

# Step 5: Long format of deletion gene data
del_long <- del_gene_rows %>%
  dplyr::mutate(row_id = row_number()) %>%
  pivot_longer(
    cols = all_of(del_regions),
    names_to = "cytoband",
    values_to = "gene"
  ) %>%
  dplyr::filter(gene != "") %>%
  dplyr::mutate(gene = as.character(gene))

# Step 6: Filter to genes in RNA
del_genes_in_rna <- del_long %>%
  dplyr::filter(gene %in% rownames(rna_matched))

# Step 7: Median expression by cytoband per sample
del_median_expr_by_cytoband <- del_genes_in_rna %>%
  dplyr::group_by(cytoband) %>%
  dplyr::summarise(expr_matrix = list({
    genes <- unique(gene)
    sub_rna <- rna_matched[genes, , drop = FALSE]
    apply(sub_rna, 2, median, na.rm = TRUE)
  }))

# Step 8: Unnest and tidy
del_median_expr_long <- del_median_expr_by_cytoband %>%
  unnest_wider(expr_matrix) %>%
  pivot_longer(
    cols = -cytoband,
    names_to = "Sample_ID",
    values_to = "Median_Expression"
  )

# Step 9: Map Sample_ID to CNV_ID
del_expr_annot <- del_median_expr_long %>%
  left_join(clinical_map %>% dplyr::select(COLOSSUS_ID, CNV_ID = Sample),
            by = c("Sample_ID" = "COLOSSUS_ID"))

#remove space from cytoband and select only amplifications
del_cnv_long_filtered <- cnv_long %>%
  dplyr::mutate(
    peak_id = str_trim(peak_id),
    cytoband = str_trim(cytoband)
  ) %>%
  dplyr::filter(str_detect(peak_id, "Deletion"))

head(del_cnv_long_filtered)

# Final join: expression + CNV by cytoband and CNV_ID
median_expr_with_del <- del_expr_annot  %>%
  left_join(del_cnv_long_filtered, by = c("cytoband", "CNV_ID"))

#check
glimpse(median_expr_with_del)
table(is.na(median_expr_with_del$CNV))

# ---------------------
# Logistic regression: CNV 0 vs 1/2 ~ expression
# ---------------------


# Step 10: Join with CNV values
#del_expr_with_cnv <- del_expr_annot %>%
  #left_join(cnv_long, by = c("cytoband", "CNV_ID"))

# Step 11: Binary CNV outcome
del_logit_data <- median_expr_with_del %>%
  dplyr::filter(!is.na(CNV)) %>%
  dplyr::mutate(CNV_binary = ifelse(CNV >= 1, 1, 0))  # 1 if CNV = 1/2, else 0

#test
# Subset to cytoband of interest
subset_rows_del <- del_logit_data %>% 
  dplyr::filter(peak_id == "Deletion Peak 27")

# Keep every alternate row (i.e., odd-numbered rows)
subset_alt_del <- subset_rows_del[seq(1, nrow(subset_rows_del), by = 1), ]

# Fit logistic regression
fit_del <- glm(CNV_binary ~ Median_Expression, data = subset_alt_del, family = binomial)

# Show summary
summary(fit_del)

peak_del <- NULL
cyto_del <- NULL
estimate_del <- NULL
pvalue_del <- NULL

# Subset to cytoband of interest
for(i in 1:length(unique(del_logit_data$peak_id))){
  subset_rows_del <- del_logit_data %>% 
    filter(peak_id == unique(del_logit_data$peak_id)[i])
  
  # Keep every alternate row (i.e., odd-numbered rows)
  subset_alt_del <- subset_rows_del[seq(1, nrow(subset_rows_del), by = 1), ]
  
  # Fit logistic regression
  fit <- glm(CNV_binary ~ Median_Expression, data = subset_alt_del, family = binomial)
  
  # Show summary
  s <- summary(fit)
  peak_del[i] <- unique(del_logit_data$peak_id)[i]
  cyto_del[i] <- unique(del_logit_data$cytoband)[i]
  estimate_del[i] <- s$coefficients[2,1]
  pvalue_del[i] <- s$coefficients[2,4]
}

p.adj_del <- p.adjust(pvalue_del, method = "fdr")

data_del <- data.frame(cbind(peak_del,cyto_del,estimate_del,pvalue_del,p.adj_del))
colnames(data_del) <- c("Peak", "Cytoband", "Estimate", "Pvalue", "P.adj")

#write.table(data_del,"1G_new_20260506_COLOSSUS_Figure_2F_deletion_alone.txt", row.names = FALSE, sep = "\t", quote = FALSE)

p.adj_all <- p.adjust(c(pvalue,pvalue_del), method = "fdr")
#,rep(NA,13)

data_all <- data.frame(rbind(data,data_del))
data_all <- data.frame(cbind(data_all,p.adj_all))
data_all$Estimate <- round(as.numeric(data_all$Estimate), digits = 3)

tmp_ggplot_data_rna <- data_all

#write.table(data_all,"20260507_COLOSSUS_Figure_2F_amplification_deletion_both.txt", row.names = FALSE, sep = "\t", quote = FALSE)

library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# --- Prepare plotting data ---
plot_df <- tmp_ggplot_data_rna %>%
  mutate(
    # ensure numeric columns (cbind/rbind often coerce to character)
    Estimate   = suppressWarnings(as.numeric(Estimate)),
    Pvalue     = suppressWarnings(as.numeric(Pvalue)),
    P.adj      = suppressWarnings(as.numeric(P.adj)),      # if present
    p.adj_all  = suppressWarnings(as.numeric(p.adj_all)),  # used in the plot
    # derive the columns used in the plot
    data_type  = ifelse(grepl("Amp|Ampl", Peak, ignore.case = TRUE),
                        "Amplification", "Deletion"),
    significant = p.adj_all < 0.2,
    is_amp      = grepl("Amp|Ampl", Peak, ignore.case = TRUE),
    neglog10    = -log10(p.adj_all)
  ) %>%
  filter(!is.na(Estimate), !is.na(p.adj_all))

# --- the plot (same look, but now backed by consistent columns) ---
ggplot(
  data = plot_df,
  aes(
    y   = is_amp,                       # two swarms: FALSE (Del) / TRUE (Amp)
    x   = Estimate,
    size= neglog10,                     # -(log10 FDR)
    label = Peak,
    col = significant
  )
) +
  geom_vline(xintercept = 0, lty = 2, col = "#FF000044") +
  geom_beeswarm(priority = "density", cex = 3, aes(shape = data_type)) +
  scale_color_manual(values = c("#44444488", "#FDD262AA")) +
  scale_shape_manual(values = c("circle", "circle")) +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.5, 6)) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank()
  )

library(ggrepel)

cytobands_to_label <- c("8q24.3",
                        "18q11.2",
                        "22q11.21",
                        "15q22.33",
                        "12q21.2")

#Add a NEW plot with labels (additional figure)
#pdf("1G_new_20260506_CNV_RNA_logistic_regression_plot.pdf")
ggplot(
  data = plot_df,
  aes(
    y   = is_amp,
    x   = Estimate,
    size= neglog10,
    label = Peak,
    col = significant
  )
) +
  geom_vline(xintercept = 0, lty = 2, col = "#FF000044") +
  geom_beeswarm(priority = "density", cex = 3, aes(shape = data_type)) +
  
  # ---- ADDITIONAL LABELS ----
geom_text_repel(
  data = subset(plot_df, Cytoband %in% cytobands_to_label),
  aes(label = Cytoband),
  size = 4,
  fontface = "bold",
  box.padding = 0.4,
  point.padding = 0.3,
  segment.color = "grey50",
  max.overlaps = Inf
) +
  
  scale_color_manual(values = c("#44444488", "#FDD262AA")) +
  scale_shape_manual(values = c("circle", "circle")) +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.5, 6)) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank()
  )

#dev.off()
