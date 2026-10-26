################################################################################
# Figure Panel : 1E
# Description  : Frequency of significant GISTIC CNA peaks (q < 0.2) across the COLOSSUS
#                cohort, plotted as the percentage of samples with loss or gain/amplification.
# Input data   : all_lesions.conf_5.txt
# Output       : 20250609_COLOSSUS_CNV_frequency_Anguraj_results_plot.pdf
#                20250930_COLOSSUS_CNV_frequency_Anguraj_results_plot.pdf
#                20250930_COLOSSUS_CNV_frequency_Anguraj_results.txt
# R packages   : tidyverse, ggrepel, tidyr, dplyr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : August 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load necessary libraries
library(tidyverse)
library(ggrepel)
library(tidyr)
library(dplyr)

# --------------------------
# Input / I/O
# --------------------------
in_file  <- "all_lesions.conf_5.txt"
out_pdf  <- "20250609_COLOSSUS_CNV_frequency_Anguraj_results_plot.pdf"

# --------------------------
# Read & prep data
# --------------------------
data <- read.delim(in_file, check.names = FALSE)

# Extract sample columns (those that start with 'sRCSI')
sample_cols <- grep("^sRCSI", names(data), value = TRUE)

# Keep the first 132 GISTIC peaks with q < 0.2
keep_idx <- which(data$`q values`[1:132] < 0.2)
data_1 <- data[keep_idx, c(sample_cols, "Unique Name", "Descriptor", "q values")]

# Compute frequency of alteration across samples
data_1$Frequency <- rowSums(data_1[, sample_cols] > 0)

# percentage of samples altered
n_samples <- length(sample_cols)

data_1$Percentage <- rowSums(data_1[, sample_cols] > 0) / n_samples * 100

# CNA type from 'Unique Name' (case-insensitive), AMP vs DEL
data_1$CNA_Type <- ifelse(grepl("^amp", data_1$`Unique Name`, ignore.case = TRUE), "AMP", "DEL")

# Short label text
data_1$Label <- paste(data_1$CNA_Type, data_1$Descriptor)

# Sort by frequency for plotting and add index
data_1 <- data_1 %>% arrange(desc(Frequency)) %>% mutate(Index = row_number())

# Sort by percentage for plotting and add inde
data_1 <- data_1 %>%
  arrange(desc(Percentage)) %>%
  mutate(Index = row_number())

# --------------------------
# Choose which labels to draw (as many top ones as possible)
# Strategy: "knee" in the frequency curve + a minimum count fallback
# --------------------------
# Unique frequencies sorted high->low
#freq_sorted <- sort(unique(data_1$Frequency), decreasing = TRUE)
freq_sorted <- sort(unique(data_1$Percentage), decreasing = TRUE)

# Find "knee" = biggest drop; label all >= threshold
if (length(freq_sorted) > 1) {
  diffs <- diff(freq_sorted)
  knee_idx <- which.max(abs(diffs))
  threshold <- freq_sorted[knee_idx]  # label everything with Frequency >= this
} else {
  threshold <- ifelse(length(freq_sorted) == 1, freq_sorted[1], 0)
}

#labels_df <- data_1 %>% filter(Frequency >= threshold)
labels_df <- data_1 %>% filter(Percentage >= threshold)

# Ensure at least N labels (helps when frequencies are flat/noisy)
min_labels <- 20
if (nrow(labels_df) < min_labels) {
 # labels_df <- data_1 %>% slice_max(Frequency, n = min_labels, with_ties = TRUE)
  labels_df <- data_1 %>% slice_max(Percentage, n = min_labels, with_ties = TRUE)
}

# Optional: safety cap to avoid extreme clutter; set to Inf to truly allow all
max_labels <- 60
if (nrow(labels_df) > max_labels) {
  labels_df <- labels_df %>% slice_max(Frequency, n = max_labels, with_ties = TRUE)
}

# --------------------------
# Plot
# --------------------------


#pdf("20250930_COLOSSUS_CNV_frequency_Anguraj_results_plot.pdf")
#ggplot(data_1, aes(x = Index, y = Frequency)) +
  ggplot(data_1, aes(x = Index, y = Percentage)) +
  geom_point(aes(fill = CNA_Type), shape = 21, size = 8, color = "black", stroke = 1) +
  scale_fill_manual(values = c("AMP" = "#F4C542", "DEL" = "#76C7C0")) +
 # labs(y = "number of samples with loss or gain/amp", fill = "CNA type") +
  labs(
    y = "Percentage of samples with loss or gain/amplification",
    fill = "CNA type"
  )
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = labels_df,
    aes(label = paste(CNA_Type, Descriptor)),
    box.padding      = 0.5,     # space around text
    point.padding    = 0.3,     # space between label and point
    segment.color    = "grey50",
    segment.size     = 0.4,
    min.segment.length = 0,     # always draw line
    max.overlaps     = Inf,     # allow many labels
    force            = 2,       # stronger repulsion
    max.time         = 2,       # time budget for solving overlaps
    nudge_y          = 0.5      # small upward nudge
  )

#dev.off()

#write.table(data_1,"20250930_COLOSSUS_CNV_frequency_Anguraj_results.txt", sep ="\t", quote = FALSE,row.names = FALSE)

