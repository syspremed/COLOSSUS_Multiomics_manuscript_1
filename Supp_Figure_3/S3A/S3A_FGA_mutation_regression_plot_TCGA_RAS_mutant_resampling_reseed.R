################################################################################
# Figure Panel : Supp S3A
# Description  : TCGA MSS RAS-mutant: logistic regression of mutation status vs FGA with 75/25
#                resampling (1,000 iterations) across multiple random seeds.
# Input data   : onco_matrix_TCGA_RAS_mt.txt
#                20250515_TCGA_immune_profiles.csv
# Output       : Supplementary_20260202_TCGA_RAS_mutant_FGA_mutation_resampling_all_iterations.csv
#                Supplementary_20260202_TCGA_RAS_mutant__FGA_mutation_resampling_reappearance_summary.csv
#                Supplementary_20260202_TCGA_RAS_mutant_FGA_mutation_resampling_1000_reseed_100_times.pdf
#                S2A_Supplementary_20260321_TCGA_RAS_mutant_FGA_mutation_resampling_1000_reseed_100_times_quality_figure.pdf
# R packages   : dplyr, janitor, tidyverse, ggplot2, scales
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : February 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# ================================
# Logistic regression with 75/25 train-test resampling
# ================================

library(dplyr)
library(janitor)
library(tidyverse)
library("wesanderson")
library(ggplot2)

# 1. Load and process mutation matrix
oncomatrix <- read.delim("onco_matrix_TCGA_RAS_mt.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# Replace NA with empty string
oncomatrix[is.na(oncomatrix)] <- ""

# Binarize: 1 if non-empty mutation entry, else 0
mutation_binary <- ifelse(oncomatrix != "", 1, 0)

# 2. Load FGA data (from clinical file)
clinical_new <- read_csv("20250515_TCGA_immune_profiles.csv") %>%
  clean_names()

fga_data <- data.frame(cbind(clinical_new$tcga_participant_barcode, clinical_new$fraction_altered))
colnames(fga_data) <- c("Sample", "fga")

# Clean sample IDs (trim to 12 characters)
colnames(mutation_binary) <- substr(gsub("\\.", "-", colnames(mutation_binary)), 1, 12)
fga_data$Sample <- substr(fga_data$Sample, 1, 12)

# 3. Align by common samples
common_samples <- intersect(colnames(mutation_binary), fga_data$Sample)
mutation_mat_common <- mutation_binary[, common_samples]
fga_common <- fga_data[fga_data$Sample %in% common_samples, ]
rownames(fga_common) <- fga_common$Sample
fga_common <- fga_common[common_samples, ]  # ensure same order

# Transpose mutation matrix: samples × genes
mut_mat <- t(mutation_mat_common)

# ======================================
# 4. Logistic regression with resampling
# ======================================

set.seed(1234)
n_iter <- 1000   # number of resampling iterations
train_frac <- 0.75
resample_results <- list()

# -------------------------------
# Multiple random-seed loop
# -------------------------------
set.seed(999)                      # reproducible seed list
seeds_to_run <- sample(100:99999, 100)

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
    
    train_idx <- sample(seq_len(nrow(mut_mat)), size = floor(train_frac * nrow(mut_mat)))
    fga_train <- fga_common[train_idx, "fga", drop = FALSE]
    
    iter_results <- list()
    
    for (gene in colnames(mut_mat)) {
      df_train <- data.frame(
        mutation = as.factor(mut_mat[train_idx, gene]),
        fga = as.numeric(fga_train$fga)
      )
      
      # Only run if mutation has variability
      if (length(unique(df_train$mutation)) > 1 && length(unique(df_train$fga)) > 1) {
        model <- glm(mutation ~ fga, data = df_train, family = binomial)
        summary_model <- summary(model)
        
        if ("fga" %in% rownames(summary_model$coefficients)) {
          coef_fga <- summary_model$coefficients["fga", ]
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
  
  # Combine all iterations (per seed)
  all_iter_results <- bind_rows(resample_results)
  
  # Summarize by gene: how often p < 0.05 (per seed)
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
  
  # Store outputs per seed
  seed_results[[paste0("seed_", s)]] <- list(
    all_iter_results = all_iter_results,
    reappearance_summary = reappearance_summary
  )
  
  # Print top recurring genes (per seed)
  print(head(reappearance_summary[reappearance_summary$prop_sig > 0.1, ]))
}

# -------------------------------
# OPTIONAL: final across-seed summary (recommended)
# -------------------------------
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

# Use this as the final result object if desired:
# reappearance_summary <- final_across_seed_summary

# ======================================
# 5. Export results
# ======================================

#write.csv(all_iter_results,
#         "Supplementary_20260202_TCGA_RAS_mutant_FGA_mutation_resampling_all_iterations.csv",
#        row.names = FALSE)

#write.csv(reappearance_summary,
#        "Supplementary_20260202_TCGA_RAS_mutant__FGA_mutation_resampling_reappearance_summary.csv",
#         row.names = FALSE)

# Print top recurring genes
# (Across-seed top genes: print(head(final_across_seed_summary[final_across_seed_summary$mean_prop_sig > 0.1, ])))
print(head(reappearance_summary[reappearance_summary$prop_sig > 0.1, ]))

# ======================================
# 6. Visualization (dark vs light red)
# ======================================

colour_wes_asteroidcity1_5 <- wes_palette("AsteroidCity1")

reappearance_summary$color_group <- ifelse(reappearance_summary$prop_sig > 0.1,
                                           "recurrent", "rare")

#pdf("Supplementary_20260202_TCGA_RAS_mutant_FGA_mutation_resampling_1000_reseed_100_times.pdf")
ggplot(reappearance_summary,
       aes(x = reorder(Gene, prop_sig),
           y = median_estimate,
           size = prop_sig,
           color = color_group)) +
  geom_point(alpha = 0.9) +
  geom_hline(yintercept = 0, lty = 1, lwd = 0.5) +
  scale_color_manual(values = c(
    "recurrent" = colour_wes_asteroidcity1_5[3],  # dark Wes Anderson red
    "rare" = "#F4A6A6"                            # light red for low reappearance
  )) +
  theme_bw() +
  labs(
    x = "Gene",
    y = "Median Estimate (logit scale)",
    size = "Prop. sig (<0.05)",
    color = "Reappearance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
#dev.off()


#Publication quality figure
library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 1. Ensure prop_sig is numeric and order genes (increasing)
# ------------------------------------------------------------

reappearance_summary_final <- reappearance_summary %>%
  mutate(prop_sig = as.numeric(prop_sig)) %>%
  arrange(prop_sig) %>%                          # increasing order
  mutate(Gene = factor(Gene, levels = unique(Gene)))

# ------------------------------------------------------------
# 2. Compute tight legend range and breaks
# ------------------------------------------------------------

range_prop <- range(reappearance_summary_final$prop_sig, na.rm = TRUE)

tight_breaks <- signif(
  seq(range_prop[1], range_prop[2], length.out = 4),
  2
)

# ------------------------------------------------------------
# 3. Publication-quality bubble plot
# ------------------------------------------------------------

#pdf("S2A_Supplementary_20260321_TCGA_RAS_mutant_FGA_mutation_resampling_1000_reseed_100_times_quality_figure.pdf")

ggplot(reappearance_summary_final,
       aes(x = Gene,
           y = median_estimate,
           size = prop_sig,
           fill = color_group)) +
  
  geom_point(shape = 21,
             colour = "black",
             stroke = 0.4,
             alpha = 0.95) +
  
  geom_hline(yintercept = 0,
             linewidth = 0.4,
             colour = "grey30") +
  
  # ---- BETTER SIZE SCALING + CLEANER LEGEND ----
scale_size_continuous(
  range = c(4, 16),
  breaks = seq(0, 1, by = 0.2),     # evenly spaced legend
  limits = c(0, 1),                 # force full scale
  labels = scales::number_format(accuracy = 0.1)
) +
  
  scale_fill_manual(values = c(
    "recurrent" = colour_wes_asteroidcity1_5[3],
    "rare" = "#F4A6A6"
  )) +
  
  # ---- EQUALLY SPACED GRID LINES ----
scale_y_continuous(
  breaks = scales::pretty_breaks(n = 6)
) +
  
  theme_bw(base_size = 16) +
  
  labs(
    x = "Gene (ranked by median estimate)",
    y = "Median Estimate (logit scale)",
    size = "Prop. significant",
    fill = "Reappearance"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    
    # Mild, equally spaced grid lines
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )
dev.off()
