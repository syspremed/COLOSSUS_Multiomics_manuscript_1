################################################################################
# Figure Panel : Supp S2K
# Description  : Bar plots of CNA amplification/deletion frequencies in the COLOSSUS mCRC RAS-mt
#                cohorts (one script per proxy cohort: 1, 2-B1 and 2-B2).
# Input data   : 20241012_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B2_Barplot.txt
#                20260404_20241012_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B2_Barplot_1024_1042_removed.txt
# Output       : S1K_20260404_COLOSSUS_Proxy_Cohort_2_2_CNV_analysis_percentage_table_2_samples_removed.txt
#                S1K_20260404_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B2_Barplot_2_samples_removed.pdf
# R packages   : ggplot2, reshape2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : October 2024 (updated 2026)
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load necessary libraries
library(ggplot2)
library(reshape2)

# Read CNA calls at the selected GISTIC peaks
data <- data.frame(read.delim("20241012_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B2_Barplot.txt", stringsAsFactors = FALSE))
dim(data)

##sRCSI1041 sRCSI1042  had different profiles but same number of "1s" - so #sRCSI1042 was removed

data <- data[,-c(2,20)]
dim(data)
colnames(data)

data <- read.delim("20260404_20241012_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B2_Barplot_1024_1042_removed.txt", stringsAsFactors = FALSE)

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

# Define colours for CNA values
custom_colors <- c("0" = "#1f77b4",   # Blue for CNA value 0
                   "1" = "#ff7f0e",   # Orange for CNA value 1
                   "2" = "#2ca02c")   # Green for CNA value 2

write.table(cna_counts,"S1K_20260404_COLOSSUS_Proxy_Cohort_2_2_CNV_analysis_percentage_table_2_samples_removed.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Plot the stacked bar plot with percentages
pdf("S1K_20260404_COLOSSUS_Proxy_cohort_GISTIC_data_selected_data_Proxy cohort 2 B2_Barplot_2_samples_removed.pdf")
ggplot(cna_counts, aes(x = Descriptor, y = Percentage, fill = CNA_value)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Stacked Bar Plot of CNA Values (100% Total per Descriptor)", x = "Descriptor", y = "Percentage (%)") +
  scale_fill_manual(values = custom_colors) +  # Use custom colors
  theme_minimal()
dev.off()
