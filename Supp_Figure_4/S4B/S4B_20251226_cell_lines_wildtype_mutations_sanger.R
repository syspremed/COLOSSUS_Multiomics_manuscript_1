################################################################################
# Figure Panel : Supp S4B
# Description  : Oncoplot of mutations in CRC RAS-wildtype cell lines (Sanger Cell Model Passports).
# Input data   : mutations_all_20250318.csv
# Output       : S3B_20260407_CRC_CELL_LINES_RAS_MT_WILDTYPE_Sanger_categorical_oncoplot_maftools.pdf
# R packages   : tidyverse, maftools
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : December 2025
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Oncoplot using maftools for RAS wildtype cell lines
# Load required libraries
library(tidyverse)
library(maftools)

# Define cell lines - removed "HT-115",
crc_cell_lines <- c(
  "C2BBe1","COLO-320-HSR","SW1417","NCI-H508","NCI-H716","HT-29","SNU-283",
  "SNU-C1","HT55","MDST8","COLO-205","DIFi","CaR-1"
)

# Define genes in original order (will be reversed for plotting)
genes_original <- c(
  "ERBB4", "FAT3", "NRAS", "PTPRT", "KMT2C", "FBXW7", "TCF7L2", 
  "ROBO2", "NBEA", "SMAD4", "CSMD3", "AMER1", "ATM", "FAT4", 
  "LRP1B", "PIK3CA", "MUC16", "TP53", "APC", "KRAS"
)

# Reverse for display (KRAS at top)
genes_of_interest <- rev(genes_original)

# Read mutation data
mutations <- read_csv("mutations_all_20250318.csv")

# Filter for genes and cell lines of interest
mutations_filtered <- mutations %>%
  filter(gene_symbol %in% genes_of_interest,
         model_name %in% crc_cell_lines)

# Convert to MAF format
# Maftools requires specific column names
maf_data <- mutations_filtered %>%
  mutate(
    Hugo_Symbol = gene_symbol,
    Tumor_Sample_Barcode = model_name,
    # Map effect to Variant_Classification
    Variant_Classification = case_when(
      grepl("frameshift.*del", effect, ignore.case = TRUE) ~ "Frame_Shift_Del",
      grepl("frameshift.*ins", effect, ignore.case = TRUE) ~ "Frame_Shift_Ins",
      grepl("frameshift", effect, ignore.case = TRUE) ~ "Frame_Shift_Del",
      grepl("nonsense", effect, ignore.case = TRUE) ~ "Nonsense_Mutation",
      grepl("missense", effect, ignore.case = TRUE) ~ "Missense_Mutation",
      grepl("inframe.*del", effect, ignore.case = TRUE) ~ "In_Frame_Del",
      grepl("inframe.*ins", effect, ignore.case = TRUE) ~ "In_Frame_Ins",
      grepl("inframe", effect, ignore.case = TRUE) ~ "In_Frame_Del",
      grepl("splice", effect, ignore.case = TRUE) ~ "Splice_Site",
      TRUE ~ "Other"
    ),
    Variant_Type = "SNP",  # Simplified
    Chromosome = chromosome,
    Start_Position = position,
    End_Position = position,
    Reference_Allele = reference,
    Tumor_Seq_Allele2 = alternative,
    Protein_Change = protein_mutation
  ) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification, 
         Variant_Type, Chromosome, Start_Position, End_Position,
         Reference_Allele, Tumor_Seq_Allele2, Protein_Change)

# Create MAF object
maf <- read.maf(maf = maf_data, verbose = FALSE)

# Define custom colors matching the example
vc_cols <- c(
  "Missense_Mutation" = "#0066CC",      # Blue
  "Nonsense_Mutation" = "#00AA00",      # Green
  "Frame_Shift_Del" = "#33CCCC",        # Cyan
  "Frame_Shift_Ins" = "#FF6600",        # Orange
  "In_Frame_Del" = "#FFCC00",           # Yellow
  "In_Frame_Ins" = "#00CC66",           # Light green
  "Splice_Site" = "#9966CC",            # Purple
  "Multi_Hit" = "#000000"               # Black
)

# Create oncoplot - PDF version
#pdf("S3B_20260407_CRC_CELL_LINES_RAS_MT_WILDTYPE_Sanger_categorical_oncoplot_maftools.pdf", width = 16, height = 10)

oncoplot(
  maf = maf,
  genes = genes_of_interest,
  keepGeneOrder = TRUE,
  showTumorSampleBarcodes = TRUE,
  SampleNamefontSize = 0.8,
  titleFontSize = 1.2,
  legendFontSize = 1,
  colors = vc_cols,
  removeNonMutated = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  sortByAnnotation = FALSE,
  bgCol = "#E0E0E0",
  borderCol = "white",
  gene_mar = 8,
  barcode_mar = 6,
  top = length(genes_of_interest),
  fontSize = 0.8
)

#dev.off()


# Print summary
cat("\n=== MAF Summary ===\n")
print(getSampleSummary(maf))
cat("\n")
print(getGeneSummary(maf))

cat("\n=== Files Created ===\n")
cat("oncoplot_maftools.pdf\n")
