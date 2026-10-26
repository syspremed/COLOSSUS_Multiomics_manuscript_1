################################################################################
# Figure Panel : 1D
# Description  : Heatmap of key CNA peaks annotated with fraction of genome altered (FGA,
#                CINmetrics) and aneuploidy scores (ASCETS) across COLOSSUS samples.
# Input data   : 20240526_CINmetrics_COLOSSUS_output.txt
#                Fig1F_Colossus_retrospectivecohort_logR_editedforGISTIC_ASCETSoutput.txt_aneuploidy_scores.txt
#                all_lesions.conf_5.txt
# Output       : fig1g_20250930_key CNAs with FGA and aneuploidy annotations heatmap.pdf
# R packages   : readr, dplyr, pheatmap, textshape, readxl
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : September 2025
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load required libraries
library(readr)
library(dplyr)
library(pheatmap)
library(textshape)

# Step 1: Read FGA data (CINmetrics output)
fga_data <- read.delim("20240526_CINmetrics_COLOSSUS_output.txt", stringsAsFactors = FALSE)

# Check column names
# Expected columns: sample_id, tai, cna, base_segments, break_points, fga
colnames(fga_data)

# Step 1b: Read aneuploidy data from file
aneuploidy_data <- read.delim("Fig1F_Colossus_retrospectivecohort_logR_editedforGISTIC_ASCETSoutput.txt_aneuploidy_scores.txt", 
                              stringsAsFactors = FALSE)

# Check the structure of aneuploidy data
print("Aneuploidy data structure:")
print(head(aneuploidy_data))
print(colnames(aneuploidy_data))

# Step 2: Read GISTIC categorical data
gistic_data <- read.delim("all_lesions.conf_5.txt")
gistic_data <- gistic_data[1:132,]
dim(gistic_data)

# Step 3: Cytoband regions of interest
regions <- c(
  "7p11.2",    # EGFR
  "8q24.3",    # PTK2 + others
  "17q23.2",   # APPBP2
  "20q13.12",  # JPH2, TOX2
  "10q23.31",  # PTEN
  "15q21.1",   # B2M
  "16q23.1",   # MAF, WWOX
  "18q12.1",   # CDH2
  "18q21.2",   # SMAD4
  "18q21.33",  # BCL2
  "18q23"      # ZNF516
)

# Trim spaces from Descriptor values
gistic_data$Descriptor <- trimws(gistic_data$Descriptor)

# Step 4: Filter to selected regions
gistic_filtered <- gistic_data %>% 
  filter(`Descriptor` %in% regions)

# Step 5: Extract only sRCSI sample columns
sRCSI_cols <- grep("^sRCSI", colnames(gistic_filtered), value = TRUE)

# Step 6: Build matrix for heatmap
gistic_matrix <- gistic_filtered[, sRCSI_cols]
rownames(gistic_matrix) <- gistic_filtered$`Unique Name`

# Convert values to numeric
gistic_matrix <- gistic_matrix %>% 
  mutate_all(as.numeric)

# Update row names with descriptors
gistic_filtered$Descriptor[5] <- "8q24.3_D"
rownames(gistic_matrix) <- gistic_filtered$Descriptor
dim(gistic_matrix)

gistic_matrix <- gistic_matrix[-5,]
dim(gistic_matrix)

# Step 7: Create FGA and aneuploidy annotations for columns and order by FGA (highest to lowest)
# Filter FGA data to match samples in the heatmap
fga_annotation <- fga_data %>%
  filter(sample_id %in% colnames(gistic_matrix)) %>%
  arrange(desc(fga)) %>%  # Order by FGA descending (highest to lowest)
  dplyr::select(sample_id, fga)

# Add aneuploidy scores to the annotation
annotation_df <- fga_annotation %>%
  left_join(aneuploidy_data, by = c("sample_id" = "sample")) %>%
  column_to_rownames("sample_id")

# Reorder the gistic_matrix columns to match FGA ordering
sample_order <- rownames(annotation_df)
gistic_matrix <- gistic_matrix[, sample_order]

# Step 8: Create color scales for FGA and aneuploidy annotations
annotation_colors <- list(
  fga = colorRampPalette(c("lightblue", "darkblue"))(100),
  aneuploidy_score = colorRampPalette(c("lightyellow", "orange", "red"))(100)
)

# Step 9: Plot heatmap with FGA and aneuploidy ann  otations
#pdf("fig1g_20250930_key CNAs with FGA and aneuploidy annotations heatmap.pdf")
pheatmap(gistic_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "pink", "red"))(3),
         main = "CNA Heatmap with FGA and Aneuploidy Annotations",
         fontsize_row = 10,
         fontsize_col = 6,
         border_color = NA,
         legend_breaks = c(0, 1, 2),
         legend_labels = c("Neutral", "Gain", "Amplification"),
         annotation_col = annotation_df,
         annotation_colors = annotation_colors)
#dev.off()
