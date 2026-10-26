################################################################################
# Figure Panel : 4E
# Description  : Bar plot of scaled MV1 loading coefficients for MCM3-7 proteins.
# Input data   : scaledLoadingsCoefficients.txt
# Output       : 4E_20260421_MV1_MCM_protein_scaled_loadings_barplot.pdf
# R packages   : plotly, dplyr, tidyverse, grid, shadowtext, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : April 2026
# Note         : Reads scaledLoadingsCoefficients.txt from the iPhenMap run (Fig4B_C step 1).
################################################################################

#MCM MV1 plot


mv <- data.frame(read.delim("scaledLoadingsCoefficients.txt", stringsAsFactors = FALSE))

mcm_loc <- grep("MCM",rownames(mv))

mcm <- data.frame(cbind(rownames(mv)[mcm_loc],mv$MV1[mcm_loc]))

colnames(mcm) <- c("Proteins", "MV1")

library(plotly)
library(dplyr)
library(tidyverse)
library(grid)
library(shadowtext)
p1 <- plot_ly(mcm, x = ~Proteins, y = ~MV1) %>%
  add_histogram()

library(ggplot2)
pdf("4E_20260421_MV1_MCM_protein_scaled_loadings_barplot.pdf")
p <- ggplot(mcm, aes(x = reorder(Proteins,-as.numeric(MV1)), y = (-1)*as.numeric(MV1)))+
  geom_col(aes(color = "#A8BAC4"), width = 0.7) + coord_flip()
p
dev.off()
