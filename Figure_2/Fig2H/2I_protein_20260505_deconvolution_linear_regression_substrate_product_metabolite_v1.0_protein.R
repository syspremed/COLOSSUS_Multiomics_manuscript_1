################################################################################
# Figure Panel : 2H (protein supplementary version)
# Description  : Same deconvolution-regression pipeline as 2H but using
#                proteomics (RPPA) data instead of RNA-seq as the bulk input.
# Input data   : (1) Seurat object (same as 2H)
#                (2) COLOSSUS RPPA proteomics matrix
#                (3) Metabolomics data
# Output       : Protein-based regression coefficient plots
# R packages   : Seurat, MuSiC, SingleCellExperiment, SummarizedExperiment,
#                dplyr, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : May 2026
################################################################################

## =========================================================
## Metabolite–gene/protein linear regression plot
## L-glutamine / Glutamate vs GOT2, GLS, GLUL
## =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(broom)

## -----------------------------
## 1. Read data
## -----------------------------

metab <- read.delim(
  "20200902_ColossusResults_Metabolomics_2_tissue_weight_Centre_run_batch_corrected_data.txt",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

protein_exp <- read.delim(
  "20200928_metabolomics_median_replicates_data_80_Non_zero_sample_median_centred_imputed_MARS_combat_lysis_batch.txt",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = 1
)

## -----------------------------
## 2. Clean metabolite matrix
## -----------------------------

rownames(metab) <- metab$Metabolites
metab_mat <- metab[, setdiff(colnames(metab), "Metabolites"), drop = FALSE]

## -----------------------------
## 3. Harmonise sample IDs
## -----------------------------

clean_protein_id <- function(x) {
  sub("_FP_.*$", "", x)
}

colnames(protein_exp) <- clean_protein_id(colnames(protein_exp))

common_samples <- intersect(colnames(metab_mat), colnames(protein_exp))

cat("Metabolomics samples:", ncol(metab_mat), "\n")
cat("Proteomics samples:", ncol(protein_exp), "\n")
cat("Matched samples:", length(common_samples), "\n")

if (length(common_samples) < 10) {
  stop("Too few matched samples. Please check sample ID formatting.")
}

metab_mat <- metab_mat[, common_samples, drop = FALSE]
protein_exp <- protein_exp[, common_samples, drop = FALSE]

## -----------------------------
## 4. Select metabolites
## -----------------------------

metab_lookup <- c(
  "L-glutamine" = "L-Glutamine_146.07_14.16",
  "Glutamate"   = "L-Glutamic acid_147.05_13.23"
)

missing_metabs <- setdiff(metab_lookup, rownames(metab_mat))
if (length(missing_metabs) > 0) {
  stop("Missing metabolite rows: ", paste(missing_metabs, collapse = ", "))
}

print(metab_lookup)

## -----------------------------
## 5. Select proteins/genes
## -----------------------------

genes_to_use <- c("GLUL", "GLS", "GOT2")

find_protein_row <- function(gene, protein_names) {
  hit <- grep(paste0("(^|_)", gene, "($|;|_)"), protein_names, value = TRUE)
  
  if (length(hit) == 0) {
    stop("No protein row found for gene: ", gene)
  }
  
  if (length(hit) > 1) {
    message("Multiple protein rows for ", gene, ": using ", hit[1])
  }
  
  hit[1]
}

protein_lookup <- setNames(
  sapply(genes_to_use, find_protein_row, protein_names = rownames(protein_exp)),
  genes_to_use
)

print(protein_lookup)

## -----------------------------
## 6. Build long analysis table
## -----------------------------

plot_data <- expand.grid(
  metabolite_label = names(metab_lookup),
  gene = names(protein_lookup),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  do({
    met_name <- metab_lookup[.$metabolite_label]
    prot_name <- protein_lookup[.$gene]
    
    data.frame(
      sample = common_samples,
      metabolite_label = .$metabolite_label,
      gene = .$gene,
      metabolite_value = as.numeric(metab_mat[met_name, common_samples]),
      protein_value = as.numeric(protein_exp[prot_name, common_samples])
    )
  }) %>%
  ungroup()

## -----------------------------
## 7. Linear regression
## protein/gene abundance ~ metabolite
## -----------------------------

reg_results <- plot_data %>%
  group_by(metabolite_label, gene) %>%
  do({
    df <- .
    df <- df[complete.cases(df[, c("metabolite_value", "protein_value")]), ]
    
    if (nrow(df) < 10) {
      return(data.frame(
        n = nrow(df),
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = NA_real_
      ))
    }
    
    fit <- lm(protein_value ~ metabolite_value, data = df)
    
    slope <- broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == "metabolite_value")
    
    data.frame(
      n = nrow(df),
      estimate = slope$estimate,
      conf.low = slope$conf.low,
      conf.high = slope$conf.high,
      p.value = slope$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_label = ifelse(
      is.na(p.value),
      "NA",
      ifelse(
        p.value < 0.05,
        paste0("* p = ", sprintf("%.3f", p.value)),
        paste0("ns p = ", sprintf("%.3f", p.value))
      )
    ),
    gene = factor(gene, levels = rev(c("GLUL", "GLS", "GOT2"))),
    metabolite_label = factor(metabolite_label, levels = c("L-glutamine", "Glutamate"))
  )

print(reg_results)

## -----------------------------
## 8. Plot
## -----------------------------

gene_cols <- c(
  "GOT2" = "#0072B2",
  "GLS"  = "#D55E00",
  "GLUL" = "#009E73"
)

p <- ggplot(
  reg_results,
  aes(
    x = estimate,
    y = gene,
    colour = gene
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0,
    linewidth = 1.2,
    na.rm = TRUE
  ) +
  geom_point(size = 5, na.rm = TRUE) +
  geom_text(
    aes(
      label = p_label,
      x = ifelse(estimate >= 0, conf.high + 0.05, estimate + 0.05)
    ),
    hjust = 0,
    size = 5,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  facet_wrap(~ metabolite_label, nrow = 1) +
  scale_colour_manual(values = gene_cols) +
  labs(
    title = "Metabolite–gene linear regression\nanalysis from glutamate metabolism",
    x = "Regression coefficient (estimate ± 95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    strip.background = element_rect(fill = "grey85", colour = "black"),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.y = element_text(face = "italic", size = 20),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p)

## -----------------------------
## 9. Save outputs
## -----------------------------

out_dir <- "Glu_new"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  reg_results,
  file = file.path(out_dir, "2I_new_20260505_glutamate_metabolism_metabolite_gene_regression_results.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  filename = file.path(out_dir, "2I_new_20260505_glutamate_metabolism_metabolite_gene_regression_plot.pdf"),
  plot = p,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(out_dir, "2I_new_20260505_glutamate_metabolism_metabolite_gene_regression_plot.png"),
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)
