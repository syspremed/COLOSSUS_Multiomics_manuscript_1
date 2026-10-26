################################################################################
# Figure Panel : 5G
# Description  : Lollipop plot of ENCODE and ChEA consensus TF (ChIP-X) enrichment (Enrichr) for
#                the validated positive-interaction protein set.
# Input data   : 5G_ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_table.txt
# Output       : TF_enrichment_OR_plot.png
#                TF_enrichment_OR_plot.pdf
#                TF_enrichment_CS_plot.png
#                TF_enrichment_CS_plot.pdf
# R packages   : (base R)
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : Input is the Enrichr results table exported from the Enrichr website.
################################################################################

## ============================================================
## ENCODE TF ChIP-seq Enrichment — Publication Figure
## Shows Odds Ratio (X) vs Term (Y), coloured by Adjusted P-value
## Filtered to Adjusted P-value < 0.2
## ============================================================

# ----- 0. Packages ------------------------------------------
pkgs <- c("ggplot2", "dplyr", "stringr", "scales", "ggrepel", "forcats")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))

# ----- 1. Load data -----------------------------------------
df <- read.delim(
  "5G_ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_table.txt",
  stringsAsFactors = FALSE, check.names = FALSE
)


# ----- 2. Parse & filter ------------------------------------
df <- df %>%
  filter(!is.na(Genes) & Genes != "") %>%
  mutate(
    Term       = as.character(Term),
    adjp       = `Adjusted P-value`,
    pval       = `P-value`,
    odds_ratio = as.numeric(`Odds Ratio`),
    comb_score = as.numeric(`Combined Score`),
    overlap_n  = as.integer(str_extract(Overlap, "^[0-9]+")),
    set_size   = as.integer(str_extract(Overlap, "[0-9]+$")),
    gene_ratio = overlap_n / set_size,
    neg_log_p  = -log10(adjp),
    # Clean TF label: strip cell line suffix after underscore
    # ENCODE format: "TF_CELLLINE" — keep TF name only for readability
    TF_name    = str_extract(Term, "^[^_]+"),
    cell_line  = str_extract(Term, "(?<=_).+$"),
    term_label = ifelse(
      !is.na(cell_line),
      paste0(TF_name, "  (", cell_line, ")"),
      Term
    )
  ) %>%
  filter(adjp < 0.2) %>%
  slice_min(adjp, n = 10) %>%          # top 10 most significant
  arrange(desc(odds_ratio))

cat("Top 10 terms (adj. p < 0.2):", nrow(df), "\n")
print(df[, c("Term", "adjp", "odds_ratio", "overlap_n")])

# Order factor by odds ratio (highest at top)
df <- df %>%
  mutate(term_label = fct_reorder(term_label, odds_ratio))

# Significance category for shape
df <- df %>%
  mutate(sig_cat = case_when(
    adjp < 0.05  ~ "adj.p < 0.05",
    adjp < 0.10  ~ "adj.p < 0.10",
    adjp < 0.20  ~ "adj.p < 0.20"
  ))

# ----- 3. Plot ----------------------------------------------
# Colour palette: 3-level significance
sig_colours <- c(
  "adj.p < 0.05" = "#08306B",
  "adj.p < 0.10" = "#2171B5",
  "adj.p < 0.20" = "#6BAED6"
)

sig_shapes <- c(
  "adj.p < 0.05" = 16,
  "adj.p < 0.10" = 17,
  "adj.p < 0.20" = 15
)

p <- ggplot(df, aes(x = odds_ratio, y = term_label)) +
  
  # Shaded significance band: OR > 1 region
  annotate("rect",
           xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#F0F7FF", alpha = 0.6) +
  
  # OR = 1 reference line
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  
  # Horizontal grid lines
  geom_hline(aes(yintercept = as.numeric(term_label)),
             colour = "grey92", linewidth = 0.35) +
  
  # Lollipop stems
  geom_segment(
    aes(x = 1, xend = odds_ratio,
        y = term_label, yend = term_label,
        colour = sig_cat),
    linewidth = 0.6, alpha = 0.7
  ) +
  
  # Dots
  geom_point(
    aes(size    = overlap_n,
        colour  = sig_cat,
        shape   = sig_cat),
    alpha = 0.95
  ) +
  
  # Odds ratio value labels at end of stem
  geom_text(
    aes(label = sprintf("%.1f", odds_ratio),
        x     = odds_ratio),
    hjust  = -0.3,
    size   = 2.8,
    colour = "grey25",
    fontface = "plain"
  ) +
  
  scale_colour_manual(
    name   = "Significance",
    values = sig_colours,
    guide  = guide_legend(
      title.position = "top",
      override.aes   = list(size = 4, shape = 16)
    )
  ) +
  
  scale_shape_manual(
    name   = "Significance",
    values = sig_shapes
  ) +
  
  scale_size_continuous(
    name   = "Gene\ncount",
    range  = c(3, 10),
    breaks = pretty(df$overlap_n, n = 4),
    guide  = guide_legend(
      title.position = "top",
      override.aes   = list(colour = "grey40", shape = 16)
    )
  ) +
  
  scale_x_continuous(
    name   = "Odds Ratio",
    expand = expansion(mult = c(0.02, 0.18)),
    breaks = pretty(c(1, df$odds_ratio), n = 6)
  ) +
  
  labs(
    title    = "Transcription Factor Enrichment",
    subtitle = paste0(
      "ENCODE TF ChIP-seq 2015  |  top 10 by adj. p  |  adj. p < 0.2  |  n = ", nrow(df), " terms\n",
      "Dot size = gene overlap count  ·  OR > 1 region shaded"
    ),
    y = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 10, colour = "grey15", hjust = 1),
    axis.text.x        = element_text(size = 9,  colour = "grey30"),
    axis.title.x       = element_text(size = 10, face = "bold",
                                      margin = margin(t = 8)),
    plot.title         = element_text(face = "bold", size = 14,
                                      hjust = 0, margin = margin(b = 3)),
    plot.subtitle      = element_text(size = 9, colour = "grey40",
                                      hjust = 0, margin = margin(b = 12)),
    plot.margin        = margin(16, 20, 16, 16),
    legend.position    = "right",
    legend.box         = "vertical",
    legend.spacing.y   = unit(0.5, "cm"),
    legend.title       = element_text(size = 9,  face = "bold"),
    legend.text        = element_text(size = 8.5),
    legend.key         = element_blank()
  )

# ----- 4. Save ----------------------------------------------
plot_h  <- max(5, nrow(df) * 0.42 + 2.5)

out_png <- "TF_enrichment_OR_plot.png"
out_pdf <- "TF_enrichment_OR_plot.pdf"

ggsave(out_png, plot = p, width = 11, height = plot_h, dpi = 220, bg = "white")
ggsave(out_pdf, plot = p, width = 11, height = plot_h, bg = "white")

message("PNG saved -> ", out_png)
message("PDF saved -> ", out_pdf)

## ---- Optional: Combined Score version ----------------------
## Same plot using comb_score instead of odds_ratio
## Combined Score = log(p) * z-score (Enrichr's composite metric)

p2 <- p %+%
  (df %>% mutate(term_label = fct_reorder(term_label, comb_score))) +
  aes(x = comb_score) +
  geom_segment(aes(x = 0, xend = comb_score,
                   y = term_label, yend = term_label,
                   colour = sig_cat),
               linewidth = 0.6, alpha = 0.7) +
  geom_text(aes(label = sprintf("%.1f", comb_score), x = comb_score),
            hjust = -0.3, size = 2.8, colour = "grey25") +
  annotate("rect", xmin = 0, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "#F0F7FF", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  scale_x_continuous(
    name   = "Combined Score  (log p × z-score)",
    expand = expansion(mult = c(0.02, 0.18)),
    breaks = pretty(c(0, df$comb_score), n = 6)
  ) +
  labs(
    title    = "Transcription Factor Enrichment",
    subtitle = paste0(
      "ENCODE TF ChIP-seq 2015  |  top 10 by adj. p  |  adj. p < 0.2  |  n = ", nrow(df), " terms\n",
      "Dot size = gene overlap count  ·  Combined Score = log(p) × z-score"
    )
  )

ggsave("TF_enrichment_CS_plot.png", plot = p2,
       width = 11, height = plot_h, dpi = 220, bg = "white")
ggsave("TF_enrichment_CS_plot.pdf", plot = p2,
       width = 11, height = plot_h, bg = "white")

message("Combined Score plots also saved.")

