################################################################################
# Figure Panel : Supp S3D
# Description  : Fisher's exact test of PIK3CA mutation status vs CNA peaks (TCGA MSS RAS-mutant).
# Input data   : Validation_TCGA_sample_name_conversion.tsv
#                onco_matrix_TCGA_RAS_mt.txt
#                all_lesions.conf_5.txt
# Output       : 20260208_TCGA_RAS_mutant_CNV_binary_data.txt
#                S2D_PIK3CA_CNA_significant_only_FISHER_combined_results_data.csv
#                S2D_20260321_TCGA_RAS_mutant_PIK3CA_CNA_significant_only_FISHER_combined_results.pdf
#                S2D_20260321_TCGA_RAS_mutant_PIK3CA_CNA_significant_only_FISHER_combined_lable_results.pdf
# R packages   : dplyr, stringr, tidyr, ggplot2, ggrepel
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Fischer test for PIK3CA - TCGA RAS mutant samples 

library(dplyr)
library(stringr)

#------------------------------------------------------------
# 1. Load TCGA sample name conversion table
#------------------------------------------------------------
conv <- read.delim(
  "Validation_TCGA_sample_name_conversion.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

#------------------------------------------------------------
# 2. Keep only COAD / READ
#------------------------------------------------------------
conv <- conv %>%
  filter(Cancer_Type_.TCGA_Label. %in% c("COAD", "READ", "COAD/READ"))

#------------------------------------------------------------
# 3a. Build clean patient-level TCGA IDs
#------------------------------------------------------------
conv <- conv %>%
  mutate(
    TCGA_patient = str_replace(Patient_id, "-01.*$", ""),
    WES_sample   = str_replace(WES_data_Sample_Name, "-01.*$", ""),
    CNA_sample   = CNA_Data_Sample_Names
  )



#------------------------------------------------------------
# 3b. Expand rows where CNA_Data_Sample_Names has multiple entries
#------------------------------------------------------------
library(tidyr)
conv <- conv %>%
  separate_rows(CNA_sample, sep = ",") %>%
  mutate(
    CNA_sample = trimws(CNA_sample)
  )

#------------------------------------------------------------
# 4. Load & binarize mutation matrix (WES-based)
#------------------------------------------------------------
oncomatrix <- read.delim("onco_matrix_TCGA_RAS_mt.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)


mutation_binary <- ifelse(oncomatrix != "", 1, 0)
colnames(mutation_binary) <- sub("^e", "s", colnames(mutation_binary))  # Rename to match CNA sample IDs

#------------------------------------------------------------
# FIX mutation matrix column names (TCGA format harmonisation)
#------------------------------------------------------------

colnames(mutation_binary) <- colnames(mutation_binary) |>
  gsub("\\.", "-", x = _) |>      # dots → hyphens
  gsub("-01$", "", x = _)         # drop sample-level suffix

#------------------------------------------------------------
# 5. Load & binarize CNA data
#------------------------------------------------------------
cna <- read.delim(
  "all_lesions.conf_5.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cna <- cna[1:132, ]

# Remove problematic peaks
peaks_to_remove <- c(
  "Amplification Peak 10", "Amplification Peak 12", "Amplification Peak 1",
  "Deletion Peak 28", "Deletion Peak 60", "Amplification Peak 60",
  "Amplification Peak 45", "Deletion Peak 70", "Amplification Peak 14",
  "Amplification Peak 6", "Deletion Peak 50", "Deletion Peak 71",
  "Deletion Peak 7", "Deletion Peak 10", "Deletion Peak 49",
  "Amplification Peak 3", "Amplification Peak 34", "Deletion Peak 38"
)

clean_name <- function(x) gsub("\\s+", " ", trimws(x))
cna <- cna[!clean_name(cna$Unique.Name) %in% clean_name(peaks_to_remove), ]

# Extract CNA calls
cna_calls <- cna[, -c(1:8)]
rownames(cna_calls) <- cna$Unique.Name
#cna_calls <- cna[, grep("TCGA", colnames(cna))]
#rownames(cna_calls) <- cna$Unique.Name

#grep("DLP",colnames(cna_calls))

# Binarize CNA (amplifications only)
cna_binary <- apply(cna_calls, 2, function(x) ifelse(x > 0, 1, 0))
cna_binary <- t(cna_binary)

#------------------------------------------------------------
# 6. Map CNA array IDs → TCGA patients
#------------------------------------------------------------
# Clean CNA rownames
rownames(cna_binary) <- sub("\\.nocnv_hg19$", "", rownames(cna_binary))

# Subset CNA to samples present in conversion table
cna_binary <- cna_binary[rownames(cna_binary) %in% conv$CNA_sample, ]

tcga_map <- setNames(conv$TCGA_patient, conv$CNA_sample)
rownames(cna_binary) <- tcga_map[rownames(cna_binary)]

# Adjust for TCGA-A6-2684 sample to choose the highly variable one from Adrian
cna_binary <- cna_binary[-c(48,277,274,211,114),]

#------------------------------------------------------------
# 7. Final intersection
#------------------------------------------------------------
common_samples <- intersect(
  colnames(mutation_binary),
  rownames(cna_binary)
)

mutation_binary <- mutation_binary[, common_samples]
cna_binary      <- cna_binary[common_samples, ]

#write.table(cna_binary,"20260208_TCGA_RAS_mutant_CNV_binary_data.txt", sep = "\t", quote = FALSE)


#------------------------------------------------------------
# 8. Example: PIK3CA mutation vector
#------------------------------------------------------------
pik3ca_mut <- mutation_binary["PIK3CA", ]

table(pik3ca_mut)

pik3ca_mut <- pik3ca_mut[common_samples]
cna_binary_subset <- cna_binary[common_samples, ]

# --- Initialize storage vectors ---
chisq_stats <- numeric(ncol(cna_binary_subset))
p_vals <- numeric(ncol(cna_binary_subset))
combined_labels <- character(ncol(cna_binary_subset))
cytobands <- character(ncol(cna_binary_subset))
event_labels <- character(ncol(cna_binary_subset))

# --- Run chi-squared tests ---
for (i in seq_len(ncol(cna_binary_subset))) {
  cna_status <- cna_binary_subset[, i]
  tbl <- table(PIK3CA_mut = pik3ca_mut, CNA = cna_status)
  
  if (length(tbl) >= 2) {
    test_result <- tryCatch({
      fisher.test(tbl)
    }, error = function(e) NULL)
    
    if (!is.null(test_result)) {
      chisq_stats[i] <- test_result$estimate
      p_vals[i] <- test_result$p.value
    } else {
      chisq_stats[i] <- NA
      p_vals[i] <- NA
    }
  } else {
    chisq_stats[i] <- NA
    p_vals[i] <- NA
  }
  
  combined_labels[i] <- paste0(cna$Unique.Name[i], " (", cna$Descriptor[i], ")")
  cytobands[i] <- cna$Wide.Peak.Limits[i]
  event_labels[i] <- ifelse(grepl("Amp", cna$Unique.Name[i], ignore.case = TRUE), "Amplified",
                            ifelse(grepl("Del", cna$Unique.Name[i], ignore.case = TRUE), "Deleted", NA))
}

# --- FDR correction ---
fdr_vals <- p.adjust(p_vals, method = "fdr")

# --- Assemble results ---
results_df <- data.frame(
  Peak_Region = combined_labels,
  CNA_Status = event_labels,
  ChiSq = chisq_stats,
  PValue = p_vals,
  FDR = fdr_vals,
  stringsAsFactors = FALSE
)

# --- Final filtering and sorting ---
results_df <- results_df[!is.na(results_df$ChiSq) & !is.na(results_df$CNA_Status), ]
results_df <- results_df[order(results_df$PValue), ]

# --- Output ---
print("Top 10 CNA peaks associated with TP53 mutation:")
print(head(results_df, 30))

# --- Optional: write to file ---
# write.csv(results_df, "S2D_PIK3CA_CNA_significant_only_FISHER_combined_results_data.csv", row.names = FALSE)

library(ggplot2)
library(ggrepel)

fdr_threshold <- 0.05
yintercept_fdr <- -log10(fdr_threshold)

#pdf("S2D_20260321_TCGA_RAS_mutant_PIK3CA_CNA_significant_only_FISHER_combined_results.pdf")
ggplot(data = results_df, aes(x = ChiSq, y = -log10(FDR), col = grepl("Del", Peak_Region))) +
  geom_hline(yintercept = yintercept_fdr, lwd = 0.5, lty = 2) +
  geom_jitter(size = 5, width = 0, height = 0) +
  ggrepel::geom_text_repel(
    data = subset(results_df, FDR < fdr_threshold),
    aes(label = Peak_Region),
    col = "#000000",
    min.segment.length = 0
  ) +
  scale_color_manual(values = c("#F2AD0088", "#00A08A88")) +
  labs(
    x = "Fisher statistic",
    y = "-log10(FDR-adjusted p-value)",
    fill = "Deletion T/F"
  ) +
  theme_bw()
#  dev.off()

# plot with label the top 5 and in addition only top 5 amplification with FDR < 0.05

library(dplyr)
library(ggplot2)
library(ggrepel)

fdr_threshold <- 0.2

# Top 5 overall (by smallest FDR)
top5_overall <- results_df %>%
  filter(!is.na(FDR)) %>%
  arrange(FDR, PValue) %>%
  slice_head(n = 5)

# Top 5 significant amplifications (Amplified & FDR < 0.05)
top5_amp_sig <- results_df %>%
  filter(CNA_Status == "Amplified", !is.na(FDR), FDR < fdr_threshold) %>%
  arrange(FDR, PValue) %>%
  slice_head(n = 5)

# Union of labels (deduplicate)
label_df <- bind_rows(top5_overall, top5_amp_sig) %>%
  distinct(Peak_Region, .keep_all = TRUE)

#pdf("S2D_20260321_TCGA_RAS_mutant_PIK3CA_CNA_significant_only_FISHER_combined_lable_results.pdf")
ggplot(results_df, aes(x = ChiSq, y = -log10(FDR), col = grepl("Del", Peak_Region))) +
  geom_hline(yintercept = -log10(fdr_threshold), lwd = 0.5, lty = 2) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Peak_Region),
    color = "black",
    size = 3.5,
    min.segment.length = 0,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("#F2AD0088", "#00A08A88")) +
  labs(
    x = "Fisher statistic",
    y = "-log10(FDR-adjusted p-value)",
    color = "Deletion"
  ) +
  theme_bw()
#dev.off()


