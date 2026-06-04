################################################################################
# Figure Panel : 2D
# Description  : Pathway impact scatter plot (MetaboAnalyst results).
#                Plots –log10(FDR) vs. pathway impact score for KEGG pathways
#                enriched in COLOSSUS RAS-mutant patient metabolomics.
# Input data   : results/20251201_.../pathway_results.csv
#                (MetaboAnalyst output – must be run first via web or R API)
# Output       : 2D_..._metaboanalyst_replot_patient_HMDB.pdf
# R packages   : ggplot2, scales
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : December 2025
# Note         : Set the input CSV path to your local MetaboAnalyst output.
################################################################################

# Required packages: ggplot2 for plotting, scales for axis rescaling.
library(ggplot2)
library(scales)

# Load MetaboAnalyst pathway analysis results.
# This CSV is the output of running MetaboAnalyst pathway analysis
# on the COLOSSUS RAS-mutant metabolomics data.
df <- read.csv("pathway_results.csv")

# Use FDR-adjusted p-value (recommended) rather than raw p-value.
# Change to df$Raw.p if raw p-values are needed for sensitivity checks.
df$logP <- -log10(df$FDR)

# FDR significance threshold line drawn at 0.05 (dotted red horizontal line).
cut_FDR_001 <- -log10(0.05)

# ── Colour palette ──────────────────────────────────────────────────────────
# Gradient from light yellow (low enrichment) to red (high enrichment).
# CUSTOM COLOUR PALETTE
# white → light yellow → medium yellow → orange → red
cols <- c(
  "#FFFF80",  # light yellow
  "#FFEB3B",  # medium yellow
  "#FFA500",  # orange
  "#FF0000"   # red
)

# Build the scatter plot:
# x-axis = pathway impact score (fraction of matched metabolites that are hubs)
# y-axis = -log10(FDR) = statistical enrichment significance
# Point size and colour both encode enrichment strength.
p <- ggplot(df, aes(x = Impact, y = logP)) +
  geom_point(aes(size = Impact, colour = logP)) +
  
  # Colour scale: –log10(p)=0 (p=1) = white
  scale_colour_gradientn(
    colours = cols,
    values = rescale(c(0, max(df$logP)*0.1, max(df$logP)*0.3, max(df$logP))),
    name = "-log10(p)"
  ) +
  
  # Size scale for circles corresponding to Impact
  scale_size_continuous(
    range = c(0.5, 8),
    name = "Impact (circle size)"
  ) +
  
  # Add horizontal dotted red line at FDR < 0.01
 geom_hline(yintercept = cut_FDR_001,
           linetype = "dotted",
          colour = "red",
           linewidth = 1) +
  
  
  theme_bw(base_size = 14) +
  theme(
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.8, "cm")
  ) +
  labs(
    x = "Pathway Impact",
    y = "-log10(FDR)"
  )

# Save figure to PDF. Uncomment to write file.
# pdf("20251208_metabolites_Anguraj_metablite_HMDB_list_COLOSSUS_RAS_mutant_patient_results_Metaboanalyst.pdf")
p
#dev.off()
