################################################################################
# Figure Panel : Supp S12A
# Description  : Scatter plot of scaled iPhenMap loadings of known metabolites on MV2 vs MV6
#                (draft MV numbering); metabolites beyond |loading| > 3 on either axis are
#                highlighted and labelled.
# Input data   : scaledLoadingsCoefficients.txt
# Output       : MV2_vs_MV6_scatter_highlighted.pdf
#                MV2_vs_MV6_scatter_highlighted.png
# R packages   : ggplot2, dplyr, ggrepel
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : April 2026
# Note         : Panel letter kept from the draft (11A); draft MV numbering in the code.
################################################################################

## ============================================================
## MV2 vs MV6 scatter plot with neat labels
## - full axis range
## - metabolites outside cutoff highlighted
## - metabolites within cutoff shown in grey and not labelled
## ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

## ------------------------------------------------------------
## 1. Read scaled loadings
## ------------------------------------------------------------
loadings <- read.delim(
  "scaledLoadingsCoefficients.txt",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## ------------------------------------------------------------
## 2. Build plotting data frame
## ------------------------------------------------------------
plot_df <- data.frame(
  Metabolite = rownames(loadings),
  MV2 = as.numeric(loadings$MV2),
  MV6 = as.numeric(loadings$MV6),
  stringsAsFactors = FALSE
)

## keep metabolites with both MV2 and MV6 available
plot_df <- plot_df %>%
  filter(!is.na(MV2), !is.na(MV6))

## optionally remove unknown metabolites
plot_df <- plot_df %>%
  filter(!grepl("unknown", Metabolite, ignore.case = TRUE))

## ------------------------------------------------------------
## 3. Clean selected metabolite names
## ------------------------------------------------------------
rename_map <- c(
  "X2...2S.4aR.8aS..2.Hydroxy.4a.methyl.8.methylenedecahydro.2.naphthalenyl.acrylic_acid_250.16_1.658" = "HMMNA",
  "X.S..2-Methylbutanoic acid" = "(S)-2-Methylbutanoic acid",
  "X...9.10..DiHOME" = "9,10-DiHOME"
)

plot_df$Metabolite <- dplyr::recode(plot_df$Metabolite, !!!rename_map)

## ------------------------------------------------------------
## 4. Define highlight and label sets
## highlight anything outside cutoff on either axis
## ------------------------------------------------------------
cutoff <- 3

plot_df <- plot_df %>%
  mutate(
    highlight = abs(MV2) > cutoff | abs(MV6) > cutoff
  )

label_df <- plot_df %>%
  filter(highlight)

## ------------------------------------------------------------
## 5. Plot
## axes extend automatically to full data range
## ------------------------------------------------------------
p <- ggplot(plot_df, aes(x = MV2, y = MV6)) +
  
  ## background metabolites
  geom_point(
    color = "grey75",
    size = 2,
    alpha = 0.7
  ) +
  
  ## highlighted metabolites
  geom_point(
    data = label_df,
    color = "#18B7B0",
    fill = "#18B7B0",
    shape = 21,
    stroke = 0.4,
    size = 3
  ) +
  
  ## cutoff lines
  geom_vline(
    xintercept = c(-cutoff, cutoff),
    linetype = "dashed",
    color = "red",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = c(-cutoff, cutoff),
    linetype = "dashed",
    color = "red",
    linewidth = 0.6
  ) +
  
  ## labels only for highlighted metabolites
  geom_text_repel(
    data = label_df,
    aes(label = Metabolite),
    size = 3,
    box.padding = 0.35,
    point.padding = 0.20,
    force = 1.2,
    force_pull = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.3,
    seed = 123
  ) +
  
  labs(
    x = "Scaled MV2 contribution",
    y = "Scaled MV6 contribution"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

print(p)

## ------------------------------------------------------------
## 6. Save
## ------------------------------------------------------------
ggsave(
  filename = "MV2_vs_MV6_scatter_highlighted.pdf",
  plot = p,
  width = 7,
  height = 6
)

ggsave(
  filename = "MV2_vs_MV6_scatter_highlighted.png",
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)