################################################################################
# Figure Panel : Supp S3B
# Description  : TCGA MSS RAS-wildtype: logistic regression of mutation status vs FGA with
#                75/25 resampling (1,000 iterations) across multiple random seeds.
# Input data   : Validation_TCGA_RASWT_oncomatrix.txt
#                20250515_TCGA_immune_profiles.csv
# Output       : Supplementary_20260210_TCGA_RAS_wildtype_FGA_mutation_resampling_1000_reseed_100_times.pdf
#                S2B_Supplementary_20260321_TCGA_RAS_wildtype_FGA_mutation_resampling_1000_reseed_100_times_publication_quality_figure.pdf
# R packages   : dplyr, janitor, tidyverse, ggplot2, scales
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : February 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# ================================
# TCGA RAS Wildtype — Logistic regression with 75/25 resampling
# ================================

library(dplyr)
library(janitor)
library(tidyverse)
library("wesanderson")
library(ggplot2)

# 1. Load and process mutation matrix
#oncomatrix <- read.delim("onco_matrix_TCGA_RAS_wt.txt",
#                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

#with KRAS and NRAS
oncomatrix = read.table("Validation_TCGA_RASWT_oncomatrix.txt", sep = "\t", header = T) 

# Replace NA with empty string
oncomatrix[is.na(oncomatrix)] <- ""

# Binarize: 1 if non-empty mutation entry, else 0
mutation_binary <- ifelse(oncomatrix != "", 1, 0)

# 2. Load FGA data (from clinical file)
clinical_new <- read_csv("20250515_TCGA_immune_profiles.csv") %>%
  clean_names()

fga_data <- data.frame(cbind(clinical_new$tcga_participant_barcode, clinical_new$fraction_altered))
colnames(fga_data) <- c("Sample", "fga")

# Clean sample IDs
colnames(mutation_binary) <- substr(gsub("\\.", "-", colnames(mutation_binary)), 1, 12)
fga_data$Sample <- substr(fga_data$Sample, 1, 12)

# 3. Align by common samples
common_samples <- intersect(colnames(mutation_binary), fga_data$Sample)
mutation_mat_common <- mutation_binary[, common_samples]
fga_common <- fga_data[fga_data$Sample %in% common_samples, ]
rownames(fga_common) <- fga_common$Sample
fga_common <- fga_common[common_samples, ]

# Transpose mutation matrix: samples × genes
mut_mat <- t(mutation_mat_common)

# ======================================
# 4. Logistic regression with resampling
#    + 100 random seeds
# ======================================

set.seed(999)                      # reproducible seed list
seeds_to_run <- sample(100:99999, 100)

n_iter <- 1000
train_frac <- 0.75

seed_results <- vector("list", length(seeds_to_run))
names(seed_results) <- paste0("seed_", seeds_to_run)

for (s in seeds_to_run) {
  
  cat("\n============================\n")
  cat("Running SEED:", s, "\n")
  cat("============================\n")
  
  set.seed(s)
  resample_results <- list()
  
  for (iter in 1:n_iter) {
    cat("Seed:", s, " Iteration:", iter, "\n")
    
    train_idx <- sample(seq_len(nrow(mut_mat)),
                        size = floor(train_frac * nrow(mut_mat)))
    fga_train <- fga_common[train_idx, "fga", drop = FALSE]
    
    iter_results <- list()
    
    for (gene in colnames(mut_mat)) {
      df_train <- data.frame(
        mutation = as.factor(mut_mat[train_idx, gene]),
        fga = as.numeric(fga_train$fga)
      )
      
      if (length(unique(df_train$mutation)) > 1 &&
          length(unique(df_train$fga)) > 1) {
        
        model <- glm(mutation ~ fga,
                     data = df_train,
                     family = binomial)
        
        sm <- summary(model)
        
        if ("fga" %in% rownames(sm$coefficients)) {
          coef_fga <- sm$coefficients["fga", ]
          iter_results[[gene]] <- data.frame(
            Gene = gene,
            Estimate = coef_fga["Estimate"],
            P.value = coef_fga["Pr(>|z|)"],
            Iteration = iter,
            Seed = s
          )
        }
      }
    }
    
    resample_results[[iter]] <- bind_rows(iter_results)
  }
  
  all_iter_results <- bind_rows(resample_results)
  
  reappearance_summary <- all_iter_results %>%
    group_by(Gene) %>%
    summarise(
      n_sig = sum(P.value < 0.05, na.rm = TRUE),
      n_tested = n(),
      prop_sig = n_sig / n_tested,
      median_estimate = median(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(prop_sig))
  
  seed_results[[paste0("seed_", s)]] <- list(
    all_iter_results = all_iter_results,
    reappearance_summary = reappearance_summary
  )
  
  cat("Top recurrent genes (seed", s, "):\n")
  print(head(reappearance_summary[reappearance_summary$prop_sig > 0.1, ]))
}

all_seeds_summary <- bind_rows(
  lapply(seed_results, function(x) x$reappearance_summary),
  .id = "SeedRun"
)

final_across_seed_summary <- all_seeds_summary %>%
  group_by(Gene) %>%
  summarise(
    mean_prop_sig = mean(prop_sig, na.rm = TRUE),
    median_prop_sig = median(prop_sig, na.rm = TRUE),
    min_prop_sig = min(prop_sig, na.rm = TRUE),
    max_prop_sig = max(prop_sig, na.rm = TRUE),
    mean_median_estimate = mean(median_estimate, na.rm = TRUE),
    n_seeds = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_prop_sig))

print(head(final_across_seed_summary))

reappearance_summary <- final_across_seed_summary %>%
  mutate(color_group = ifelse(mean_prop_sig > 0.1, "recurrent", "rare"))

plot_df <- final_across_seed_summary %>%
  mutate(
    color_group = ifelse(mean_prop_sig > 0.1, "recurrent", "rare")
  )
colour_wes_asteroidcity1_5 <- wes_palette("AsteroidCity1")

#pdf("Supplementary_20260210_TCGA_RAS_wildtype_FGA_mutation_resampling_1000_reseed_100_times.pdf")
ggplot(plot_df,
       aes(x = reorder(Gene, mean_prop_sig),
           y = mean_median_estimate,
           size = mean_prop_sig,
           color = color_group)) +
  geom_point(alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "recurrent" = colour_wes_asteroidcity1_5[3],  # dark red
      "rare" = "#F4A6A6"                            # light red
    )
  ) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +
  labs(
    x = "Gene",
    y = "Mean median estimate (logit scale)",
    size = "Mean prop. sig (p < 0.05)",
    color = "Reappearance\n(across seeds)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
#dev.off()

#Publication quality figures
library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 1. Ensure numeric + order genes by mean_prop_sig (increasing)
# ------------------------------------------------------------

plot_df <- plot_df %>%
  mutate(mean_prop_sig = as.numeric(mean_prop_sig)) %>%
  arrange(mean_prop_sig) %>%
  mutate(Gene = factor(Gene, levels = unique(Gene)))

# ------------------------------------------------------------
# 2. Tight legend scaling
# ------------------------------------------------------------

range_prop <- range(plot_df$mean_prop_sig, na.rm = TRUE)

tight_breaks <- signif(
  seq(range_prop[1], range_prop[2], length.out = 4),
  2
)

# ------------------------------------------------------------
# 3. Publication-quality plot
# ------------------------------------------------------------

#pdf("S2B_Supplementary_20260321_TCGA_RAS_wildtype_FGA_mutation_resampling_1000_reseed_100_times_publication_quality_figure.pdf",
#    width = 11,
#    height = 6)

ggplot(plot_df,
       aes(x = Gene,
           y = mean_median_estimate,
           size = mean_prop_sig,
           fill = color_group)) +
  
  geom_point(shape = 21,
             colour = "black",
             stroke = 0.4,
             alpha = 0.95) +
  
  geom_hline(yintercept = 0,
             linewidth = 0.5,
             colour = "grey30") +
  
  scale_size_area(
    max_size = 14,
    limits = range_prop,
    breaks = tight_breaks,
    labels = number_format(accuracy = 0.01)
  ) +
  
  scale_fill_manual(values = c(
    "recurrent" = colour_wes_asteroidcity1_5[3],
    "rare"      = "#F4A6A6"
  )) +
  
  scale_y_continuous(
    breaks = pretty_breaks(n = 6)
  ) +
  
  theme_bw(base_size = 16) +
  
  labs(
    x = "Gene",
    y = "Mean Median Estimate (logit scale)",
    size = "Mean Proportion Significant (p < 0.05)",
    fill = "Recurrence Status\n(across seeds)"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

dev.off()

