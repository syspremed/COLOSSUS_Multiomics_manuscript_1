################################################################################
# Figure Panel : Supp S3C
# Description  : Fisher's exact test of PIK3CA mutation status vs significant CNA peaks (COLOSSUS).
# Input data   : onco_matrix.txt
#                all_lesions.conf_5.txt
# Output       : 1K_20260320_20260201_PIK3CA_CNA_FISHER_combined_results.csv
#                1K_20260320_20260201_PIK3CA_CNA_FISHER_combined_results.pdf
# R packages   : ggplot2, ggrepel
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# --- Load and binarize mutation matrix ---
oncomatrix <- read.table(
  "onco_matrix.txt",
  sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE
)

mutation_binary <- ifelse(oncomatrix != "", 1, 0)
colnames(mutation_binary) <- sub("^e", "s", colnames(mutation_binary))  # Rename to match CNA sample IDs

# --- Load and binarize CNA data ---
cna <- read.delim(
  "all_lesions.conf_5.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

cna <- cna[1:132,]
dim(cna)

# --- Remove problematic peaks before analysis ---
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

dim(cna)

# --- Continue with binarization ---
cna_calls <- cna[, grep("^sRCSI", colnames(cna))]
rownames(cna_calls) <- cna$Unique.Name

cna_binary <- apply(cna_calls, 2, function(x) ifelse(x > 0, 1, 0))
cna_binary <- t(cna_binary)

# --- Extract PIK3CA mutation status ---
pik3ca_mut <- mutation_binary["PIK3CA", ]

# --- Match sample names ---
common_samples <- intersect(names(pik3ca_mut), rownames(cna_binary))
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
print("Top 10 CNA peaks associated with PIK3CA mutation:")
print(head(results_df, 10))

# --- Optional: write to file ---
# write.csv(results_df, "1K_20260320_20260201_PIK3CA_CNA_FISHER_combined_results.csv", row.names = FALSE)

library(ggplot2)
library(ggrepel)

fdr_threshold <- 0.1
yintercept_fdr <- -log10(fdr_threshold)

#pdf("1K_20260320_20260201_PIK3CA_CNA_FISHER_combined_results.pdf")
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
 #dev.off()