################################################################################
# Figure Panel : Supp S2G
# Description  : GISTIC G-score plot for amplifications/gains and deletions (COLOSSUS St2/3 +
#                collated APD-extended cohort, n = 1,062).
# Input data   : all_lesions.conf_5.txt
#                amp_genes.conf_5.txt
#                del_genes.conf_5.txt
# Output       : supp_gistic_score_plot_from_Adrian.pdf
# R packages   : maftools
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# GISTIC overall plot for CNAs

library(maftools)

gistic_from_tom = readGistic(
  gisticAllLesionsFile = "all_lesions.conf_5.txt", 
  gisticAmpGenesFile = "amp_genes.conf_5.txt", 
  gisticDelGenesFile = "del_genes.conf_5.txt", 
  gisticScoresFile = "scores.gistic")

#pdf(file = "supp_gistic_score_plot_from_Adrian.pdf", width = 7, height = 5)
gisticChromPlot(gistic=gistic_from_tom, ref.build='hg38', txtSize = 0.001, fdrCutOff = 0.2)
#dev.off()
