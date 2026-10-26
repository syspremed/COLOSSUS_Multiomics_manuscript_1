################################################################################
# Figure Panel : Supp S2K
# Description  : Bar plots of CNA amplification/deletion frequencies in the COLOSSUS mCRC RAS-mt
#                cohorts (one script per proxy cohort: 1, 2-B1 and 2-B2).
# Input data   : 20241012_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B1.txt
# Output       : S1K_20260320_COLOSSUS_Proxy_Cohort_2_1_CNV_analysis_percentage_table.txt
# R packages   : ggplot2, reshape2, tidyr, dplyr
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : October 2024 (updated 2026)
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)

# Read CNA calls at the selected GISTIC peaks
data <- data.frame(read.delim("20241012_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B1.txt", stringsAsFactors = FALSE))
dim(data)

# Melt the data for ggplot
melted_data <- melt(data, id.vars = "Descriptor")

# Create a count of CNA values by descriptor
cna_counts <- as.data.frame(table(melted_data$Descriptor, melted_data$value))

# Rename columns
colnames(cna_counts) <- c("Descriptor", "CNA_value", "Frequency")

# Convert CNA_value to a factor for better plotting
cna_counts$CNA_value <- as.factor(cna_counts$CNA_value)

# Normalize the frequencies to get percentages for each descriptor
cna_counts <- cna_counts %>%
  group_by(Descriptor) %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)

write.table(cna_counts,"S1K_20260320_COLOSSUS_Proxy_Cohort_2_1_CNV_analysis_percentage_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Define colours for CNA values
custom_colors <- c("0" = "#1f77b4",   # Blue for CNA value 0
                   "1" = "#ff7f0e",   # Orange for CNA value 1
                   "2" = "#2ca02c")   # Green for CNA value 2

# Plot the stacked bar plot with percentages
ggplot(cna_counts, aes(x = Descriptor, y = Percentage, fill = CNA_value)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Stacked Bar Plot of CNA Values (100% Total per Descriptor)", x = "Descriptor", y = "Percentage (%)") +
  scale_fill_manual(values = custom_colors) +  # Use custom colors
  theme_minimal()
