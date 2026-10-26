################################################################################
# Figure Panel : Supp S10Q (step 1)
# Description  : MCP-counter immune deconvolution of COLOSSUS St2/3 bulk RNA-seq (n = 134).
# Input data   : 20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt
# Output       : 4L_20260422_COLOSSUS_st2_3_MCPCounter_immune_1M_134_samples.txt
#                20260422_COLOSSUS_st2_3_MCPCounter_immune_1M_134_samples.txt
# R packages   : MCPcounter
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : April 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

library(MCPcounter)

expr <- read.delim("20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt", stringsAsFactors = FALSE)

# expr: genes x samples matrix
# rownames(expr) = gene symbols

mcp_res <- MCPcounter.estimate(
  expression = expr,
  featuresType = "HUGO_symbols"
)

# Output: cell populations x samples
dim(mcp_res)
head(mcp_res[,1:5])

write.table(mcp_res,"4L_20260422_COLOSSUS_st2_3_MCPCounter_immune_1M_134_samples.txt", sep = "\t", quote = FALSE)
write.table(mcp_res,"20260422_COLOSSUS_st2_3_MCPCounter_immune_1M_134_samples.txt", sep = "\t", quote = FALSE)
