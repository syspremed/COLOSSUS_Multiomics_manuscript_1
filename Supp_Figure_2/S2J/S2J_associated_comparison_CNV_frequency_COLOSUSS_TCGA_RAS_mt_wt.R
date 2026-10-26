################################################################################
# Figure Panel : Supp S2J
# Description  : Comparison of significant CNA frequencies between COLOSSUS St2/3 MSS RAS-mt,
#                TCGA MSS RAS-mt and TCGA MSS RAS-wt cohorts.
# Input data   : S1J_COLOSSUS_TCGA_RAS_mt_wt_CNA_frequency_table.txt
# Output       : S1J_20260320_COLOSSUS_TCGA_RAS_mt_wt_CNV_frequency_comparison_plot.pdf
# R packages   : tidyverse
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# CNA peak comparison: TCGA RAS-mutant vs RAS-wildtype vs COLOSSUS
#  - convert to long format with a single Peak column
#  - plot grouped points per peak
#  - facet by Amplification vs Deletion


library(tidyverse)

df_raw <- read.delim(
 "S1J_COLOSSUS_TCGA_RAS_mt_wt_CNA_frequency_table.txt",
  stringsAsFactors = FALSE
)

# Rename clearly
colnames(df_raw) <- c(
  "Peak_mt", "Patients_mt",
  "Peak_wt", "Patients_wt",
  "Peak_col", "Patients_col"
)

# Build clean wide table
df_clean <- df_raw %>%
  select(
    Peak = Peak_mt,
    TCGA.RAS.mt = Patients_mt,
    TCGA.RAS.wt = Patients_wt,
    COLOSSUS    = Patients_col
  )

# Now pivot safely
df_long <- df_clean %>%
  pivot_longer(
    cols = c(TCGA.RAS.mt, TCGA.RAS.wt, COLOSSUS),
    names_to = "Cohort",
    values_to = "Patients"
  ) %>%
  mutate(
    Type = ifelse(str_detect(Peak, "Amplification"), "AMP", "DEL"),
    Peak_num = as.numeric(str_extract(Peak, "\\d+"))
  )

#Publication-Quality Plot
#pdf("S1J_20260320_COLOSSUS_TCGA_RAS_mt_wt_CNV_frequency_comparison_plot.pdf")
ggplot(df_long,
       aes(y = reorder(Peak, Patients),
           x = Patients,
           colour = Cohort)) +
  geom_point(size = 2.8) +
  geom_line(aes(group = Peak),
            colour = "grey75",
            alpha = 0.5) +
  facet_wrap(~Type, scales = "free_y") +
  scale_colour_manual(values = c(
    "TCGA.RAS.mt" = "#D55E00",
    "TCGA.RAS.wt" = "#0072B2",
    "COLOSSUS"    = "#009E73"
  )) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 7)
  ) +
  labs(
    y = "CNA Peak",
    x = "Number of Patients",
    colour = "Cohort"
  )
#dev.off()
