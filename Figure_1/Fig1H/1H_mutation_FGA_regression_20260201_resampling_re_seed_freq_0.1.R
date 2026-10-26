################################################################################
# Figure Panel : 1H
# Description  : Logistic regression of mutation status vs FGA with 75/25 resampling
#                (1,000 iterations) repeated across 100 random seeds; genes mutated in >= 10% of samples.
# Input data   : onco_matrix.txt
#                20240526_CINmetrics_COLOSSUS_output.txt
# Output       : Supplementary_for_fig1J_20260202_mutation_FGA_logistic__resampling_re_seed_freq_0.1_FDR_order_results.pdf
#                Supplementary_for_fig1J_20260202_mutation_FGA_logistic__resampling_re_seed_freq_0.1_results.pdf
#                1J_20260320_mutation_FGA_logistic__resampling_re_seed_freq_0.1_results_publication_quality_figure.pdf
# R packages   : dplyr, ggplot2, scales
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : February 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load necessary library
library(dplyr)

set.seed(248)

# 1. Load and process mutation matrix
oncomatrix <- read.table("onco_matrix.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# Binarize: 1 if non-empty mutation entry, else 0
mutation_binary <- ifelse(oncomatrix != "", 1, 0)

# Fix sample names: change "eRCSI" to "sRCSI"
colnames(mutation_binary) <- sub("^e", "s", colnames(mutation_binary))

# 2. Load FGA data
fga_data <- read.delim("20240526_CINmetrics_COLOSSUS_output.txt",
                       stringsAsFactors = FALSE)

colnames(fga_data)[colnames(fga_data) == "sample_id"] <- "Sample"
fga_data <- fga_data[, c("Sample", "fga")]

# 3. Align by common samples
common_samples <- intersect(colnames(mutation_binary), fga_data$Sample)
mutation_mat_common <- mutation_binary[, common_samples]
fga_common <- fga_data[fga_data$Sample %in% common_samples, ]
rownames(fga_common) <- fga_common$Sample
fga_common <- fga_common[common_samples, ]
mut_mat <- t(mutation_mat_common)

### NEW SECTION: parameters for resampling
# set.seed(932)
n_iter <- 1000  # number of resampling iterations
train_frac <- 0.75

### Run the full analysis across 100 random seeds
set.seed(999)                 # for reproducibility of the seed list
seeds_to_run <- sample(100:99999, 100)   # 100 random seeds

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
    test_idx <- setdiff(seq_len(nrow(mut_mat)), train_idx)
    
    fga_train <- fga_common[train_idx, "fga", drop = FALSE]
    
    iter_results <- list()
    
    for (gene in colnames(mut_mat)) {
      df_train <- data.frame(
        mutation = mut_mat[train_idx, gene],
        fga = fga_train$fga
      )
      
      if (length(unique(df_train$mutation)) > 1 &&
          length(unique(df_train$fga)) > 1) {
        
        model <- glm(mutation ~ fga, data = df_train, family = binomial)
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
}

### OPTIONAL PLOT (dark vs light red by prop_sig threshold)
library("wesanderson")
library(ggplot2)

colour_wes_asteroidcity1_5 <- wes_palette("AsteroidCity1")

reappearance_summary$color_group <- ifelse(
  reappearance_summary$prop_sig > 0.1,
  "recurrent",
  "rare"
)

#pdf("Supplementary_for_fig1J_20260202_mutation_FGA_logistic__resampling_re_seed_freq_0.1_FDR_order_results.pdf")

ggplot(reappearance_summary,
       aes(x = reorder(Gene, prop_sig),
           y = median_estimate,
           size = prop_sig,
           color = color_group)) +
  geom_point(alpha = 0.9) +
  geom_hline(yintercept = 0, lty = 1, lwd = 0.5) +
  scale_color_manual(values = c(
    "recurrent" = colour_wes_asteroidcity1_5[3],  # dark Wes color (e.g. coral/red tone)
    "rare" = "#F4A6A6"  # light red
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

reappearance_summary


# Combine all iteration-level results across ALL seeds
all_seeds_iter_results <- bind_rows(
  lapply(seed_results, function(x) x$all_iter_results)
)

# Final summary across ALL seeds + iterations
reappearance_summary_final <- all_seeds_iter_results %>%
  group_by(Gene) %>%
  summarise(
    n_sig = sum(P.value < 0.05, na.rm = TRUE),
    n_tested = sum(!is.na(P.value)),
    prop_sig = n_sig / n_tested,
    median_estimate = median(Estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_sig)) %>%
  mutate(color_group = ifelse(prop_sig > 0.1, "recurrent", "rare"))

reappearance_summary_final <- reappearance_summary_final %>%
  arrange(abs(median_estimate)) %>%
  mutate(Gene = factor(Gene, levels = Gene))


#pdf("Supplementary_for_fig1J_20260202_mutation_FGA_logistic__resampling_re_seed_freq_0.1_results.pdf")
ggplot(reappearance_summary_final,
       aes(x = Gene,
           y = median_estimate,
           size = prop_sig,
           color = color_group)) +
  geom_point(alpha = 0.9) +
  geom_hline(yintercept = 0, lty = 1, lwd = 0.5) +
  scale_color_manual(values = c(
    "recurrent" = colour_wes_asteroidcity1_5[3],
    "rare" = "#F4A6A6"
  )) +
  theme_bw() +
  labs(
    x = "Gene (ranked by median estimate)",
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

reappearance_summary_final <- reappearance_summary_final %>%
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

#pdf("1J_20260320_mutation_FGA_logistic__resampling_re_seed_freq_0.1_results_publication_quality_figure.pdf")

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
#dev.off()
