################################################################################
# Figure Panel : 5F
# Description  : Lollipop plot of MSigDB Hallmark enrichment (Enrichr) for proteins positively
#                associated with validated protein-metabolite interactions.
# Input data   : 5F_MSigDB_Hallmark_2020_table_selected.txt
# Output       : 5F_Msigdb_hallmark_dotplot.png
#                5F_Msigdb_hallmark_dotplot.pdf
# R packages   : scales
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
# Note         : Input is the Enrichr results table exported from the Enrichr website.
################################################################################

## ============================================================
## MSigDB Hallmark Enrichment Dot Plot
## Terms filtered: Adjusted P-value < 0.2
## X-axis : Gene Ratio (overlap / gene set size)
## Dot size: Gene count (overlap numerator)
## Dot colour: Adjusted P-value (dark = more significant)
## ============================================================

# ----- 0. Packages ------------------------------------------
pkgs <- c("ggplot2", "dplyr", "stringr", "scales")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))

# ----- 1. Load & filter data --------------------------------
df <- read.delim("5F_MSigDB_Hallmark_2020_table_selected.txt",
                 stringsAsFactors = FALSE, check.names = FALSE)

df <- df %>%
  filter(!is.na(Genes) & Genes != "") %>%
  mutate(
    Term      = as.character(Term),
    adjp      = `Adjusted P-value`,
    # Parse "6/200" → numerator and denominator
    overlap_n = as.integer(str_extract(Overlap, "^[0-9]+")),
    set_size  = as.integer(str_extract(Overlap, "[0-9]+$")),
    gene_ratio = overlap_n / set_size,
    neg_log_p  = -log10(adjp),
    # Clean term label
    term_label = str_wrap(paste0("HALLMARK_", toupper(Term)), width = 35)
  ) %>%
  filter(adjp < 0.2) %>%
  arrange(gene_ratio)                      # order for factor levels

# Order terms by gene ratio for a clean vertical layout
df$term_label <- factor(df$term_label, levels = df$term_label)

cat("Terms included (adj. p < 0.2):", nrow(df), "\n")
print(df[, c("Term", "adjp", "overlap_n", "gene_ratio")])


# ----- 2. Dot plot ------------------------------------------
p <- ggplot(df, aes(x = gene_ratio, y = term_label)) +
  
  # Reference grid lines
  geom_vline(xintercept = seq(0, max(df$gene_ratio) * 1.2, by = 0.01),
             colour = "grey93", linewidth = 0.35) +
  
  # Lollipop stem
  geom_segment(aes(x = 0, xend = gene_ratio,
                   y = term_label, yend = term_label),
               colour = "grey75", linewidth = 0.5, linetype = "solid") +
  
  # Dots — size = gene count, colour = adjusted p-value
  geom_point(aes(size = overlap_n, colour = adjp), alpha = 0.92) +
  
  # Colour scale: low p = dark blue, high p = light yellow
  scale_colour_gradientn(
    name   = "Adjusted\np-value",
    colours = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#FFFFD4"),
    limits  = c(0, 0.05),
    breaks  = c(0, 0.01, 0.02, 0.03, 0.04, 0.05),
    labels  = c("0", "0.01", "0.02", "0.03", "0.04", "0.05"),
    oob     = scales::squish,
    guide   = guide_colourbar(
      barwidth = 0.8,
      barheight = 6,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  # Size scale: gene overlap count
  scale_size_continuous(
    name   = "Gene\ncount",
    range  = c(3, 11),
    breaks = pretty(df$overlap_n, n = 4),
    guide  = guide_legend(title.position = "top", title.hjust = 0.5,
                          override.aes = list(colour = "grey40"))
  ) +
  
  scale_x_continuous(
    name   = "Gene Ratio (overlap / gene set size)",
    expand = expansion(mult = c(0.01, 0.12)),
    labels = label_number(accuracy = 0.001)
  ) +
  
  # Significance threshold lines
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  
  # Annotate significance brackets on colour bar
  #annotate("text", x = max(df$gene_ratio) * 1.02, y = nrow(df) * 0.15,
   #        label = "adj.p\n< 0.05", size = 2.6, colour = "#08519C",
    #       hjust = 0, fontface = "italic") +
  #annotate("text", x = max(df$gene_ratio) * 1.02, y = nrow(df) * 0.5,
   #        label = "adj.p\n< 0.15", size = 2.6, colour = "#4292C6",
  #         hjust = 0, fontface = "italic") +
  
  labs(
    title    = "MSigDB Hallmark Pathway Enrichment",
    subtitle = paste0("Adjusted p-value < 0.2  |  n = ", nrow(df), " terms"),
    y        = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y  = element_text(size = 10, colour = "grey15", hjust = 1),
    axis.text.x  = element_text(size = 9,  colour = "grey30"),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    plot.title   = element_text(face = "bold", size = 14,
                                hjust = 0, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9.5, colour = "grey40",
                                 hjust = 0, margin = margin(b = 12)),
    plot.margin  = margin(16, 16, 16, 16),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.6, "cm"),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )

# ----- 3. Save ----------------------------------------------
out_png <- "5F_Msigdb_hallmark_dotplot.png"
out_pdf <- "5F_Msigdb_hallmark_dotplot.pdf"

# Dynamic height: ~0.45 inches per term, minimum 5 inches
plot_h <- max(5, nrow(df) * 0.45 + 2)

ggsave(out_png, plot = p, width = 10, height = plot_h, dpi = 220, bg = "white")
ggsave(out_pdf, plot = p, width = 10, height = plot_h, bg = "white")

message("PNG saved -> ", out_png)
message("PDF saved -> ", out_pdf)
