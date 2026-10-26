################################################################################
# Figure Panel : Supp S2H
# Description  : CNA frequency plot for TCGA MSS RAS-mutant CRC.
# Input data   : Validation_TCGA_sample_name_conversion.tsv
#                onco_matrix_TCGA_RAS_mt.txt
#                all_lesions.conf_5.txt
# Output       : S1G_20260320_TCGA_RAS_mutant_CNV_frequency_data.txt
#                conv_cna_ras_mt.txt
#                20260208_TCGA_RAS_mutant_CNV_binary_data.txt
#                S1G_20260320_TCGA_RAS_mutant_CNV_frequency_plot.pdf
# R packages   : dplyr, stringr, tidyr, ggplot2, ggrepel
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

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
dim(cna_binary)

#conv_cna <- conv[conv$Patient_id %in% rownames(cna_binary),]
#dim(conv_cna)

#write.table(conv_cna,"conv_cna_ras_mt.txt", sep = "\t", quote = FALSE)

#m<-match(common_samples,conv$Patient_id)
#w<- which(!is.na(m))
#length(w)

#mutation_binary <- mutation_binary[, common_samples[w]]
#cna_binary      <- cna_binary[common_samples[w], ]
#dim(cna_binary)


#write.table(cna_binary,"20260208_TCGA_RAS_mutant_CNV_binary_data.txt", sep = "\t", quote = FALSE)

#------------------------------------------------------------
# 8. Separate AMP and DEL peaks
#------------------------------------------------------------

library(tidyr)
library(ggplot2)

peak_type <- data.frame(
  Peak = colnames(cna_binary),
  Type = ifelse(grepl("^Amplification", colnames(cna_binary)), "AMP", "DEL"),
  stringsAsFactors = FALSE
)

#------------------------------------------------------------
# 9. Count patients with CNA per peak
#------------------------------------------------------------
cna_binary_num <- apply(cna_binary, 2, function(x) as.numeric(as.character(x)))
cna_freq <- apply(cna_binary_num, 2, sum)

cna_freq_df <- data.frame(
  Peak = names(cna_freq),
  Patients = as.numeric(cna_freq),
  stringsAsFactors = FALSE
)

cna_freq_df <- left_join(cna_freq_df, peak_type, by = "Peak")


#------------------------------------------------------------
#10. Rank peaks by frequency (global ranking)
#------------------------------------------------------------

cna_freq_df <- cna_freq_df %>%
  arrange(desc(Patients)) %>%
  mutate(Rank = row_number())

write.table(cna_freq_df,"S1G_20260320_TCGA_RAS_mutant_CNV_frequency_data.txt", sep = "\t", quote = FALSE)


#------------------------------------------------------------
#11. Plot
#------------------------------------------------------------

ggplot(cna_freq_df, aes(x = Rank, y = Patients, colour = Type)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_colour_manual(
    values = c(
      "DEL" = "#00A087FF",  # teal
      "AMP" = "#E69F00FF"   # orange
    )
  ) +
  labs(
    title = "Frequency of CNA in\nTCGA MSS RAS mt CRC",
    x = NULL,
    y = "Patients with any loss or gain/amp",
    colour = "Type"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18)
  )


#Optional: label only the top peaks
top_peaks <- cna_freq_df %>% slice(1:6)

ggplot(cna_freq_df, aes(x = Rank, y = Patients, colour = Type)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(
    data = top_peaks,
    aes(label = Peak),
    hjust = -0.05,
    size = 5,
    colour = "black"
  ) +
  scale_colour_manual(
    values = c("DEL" = "#00A087FF", "AMP" = "#E69F00FF")
  ) +
  labs(
    title = "Frequency of CNA in\nTCGA MSS RAS mt CRC",
    x = NULL,
    y = "Patients with any loss or gain/amp",
    colour = "Type"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#with circles as outline
#pdf("S1G_20260320_TCGA_RAS_mutant_CNV_frequency_plot.pdf")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

#------------------------------------------------------------
# Prepare labels: top 3 AMP and top 3 DEL by Patients
#------------------------------------------------------------
label_df <- cna_freq_df %>%
  group_by(Type) %>%
  arrange(desc(Patients), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(
    Peak_label = Peak,
    Peak_label = str_replace(Peak_label, "^Amplification Peak\\s*", "Amp "),
    Peak_label = str_replace(Peak_label, "^Deletion Peak\\s*", "Del ")
  )

#------------------------------------------------------------
# Plot
#------------------------------------------------------------
p <- ggplot(cna_freq_df, aes(x = Rank, y = Patients, fill = Type)) +
  geom_point(
    shape = 21,
    size = 4,
    colour = "black",
    stroke = 0.3,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = Peak_label),
    colour = "black",
    size = 5,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.colour = "grey40",
    segment.size = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_fill_manual(
    values = c(
      "DEL" = "#00A087FF",
      "AMP" = "#E69F00FF"
    )
  ) +
  labs(
    title = "Frequency of CNA in\nTCGA MSS RAS mt CRC",
    x = NULL,
    y = "Patients with any loss or gain/amp",
    fill = "Type"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p
#dev.off()

