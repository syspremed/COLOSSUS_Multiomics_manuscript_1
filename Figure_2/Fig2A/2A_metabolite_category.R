################################################################################
# Figure Panel : 2A
# Description  : Automated classification of 207 detected metabolites into
#                biochemical categories (Fatty Acids, Amino Acids, Organic Acids,
#                Nucleotides, Sugars, Sulfonic Acids, Miscellaneous) and
#                generation of a publication-ready pie chart.
#                Steps: load reference category table → match metabolite names →
#                apply manual overrides for unmatched entries → compute
#                per-category counts/percentages → plot and export figure.
# Input data   : metabolite_207_full_category_mapping.xlsx
#                (reference table mapping metabolite names to biochemical classes)
# Output       : 1A_metabolite_207_category_pie_chart.pdf / .svg / .png
#                metabolite_207_full_category_mapping.csv
#                metabolite_207_category_summary.csv
# R packages   : readxl, readr, dplyr, tidyr, stringr, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : 2026
# Note         : Place metabolite_207_full_category_mapping.xlsx in your
#                working directory before running.
################################################################################

# ============================================================
# Automate metabolite classification + recreate pie chart
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readr)

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

# Reference category file used for the original pie chart
category_file <- "metabolite_207_full_category_mapping.xlsx"

# File containing your 207 metabolites, one per line:
# e.g. Glyoxylic acid_74_4.82
#metabolite_file <- "metabolites_207.txt"

# ------------------------------------------------------------
# 2. Read category reference table
# ------------------------------------------------------------

category_ref <- read_excel(category_file)

# The Excel file uses merged cells for category names — forward-fill populates
# each row with the correct category after reading into R.
category_ref <- category_ref %>%
  mutate(Category = na_if(Category, "")) %>%
  fill(Category, .direction = "down")

# Rename to a standard column name. Values should be clean metabolite names
# (e.g. 'Glyoxylic acid', 'Glycine') without mass/RT suffixes.
category_ref <- category_ref %>%
  rename(Metabolite_Name = `Metabolite_Name`) %>%
  mutate(Metabolite_Name = str_trim(Metabolite_Name)) %>%
  distinct(Metabolite_Name, .keep_all = TRUE)

# ------------------------------------------------------------
# 3. Read 207 metabolites
# ------------------------------------------------------------

#metab_raw <- read_lines(metabolite_file)

metab_raw <- category_ref$Input_Metabolite

metab_df <- tibble(Input_Metabolite = metab_raw) %>%
  mutate(
    # Extract RT (last field after final underscore) and MW (second-to-last field).
    Retention_Time = str_extract(Input_Metabolite, "[^_]+$"),
    Molecular_Weight = str_extract(str_remove(Input_Metabolite, "_[^_]+$"), "[^_]+$"),
    Metabolite_Name = str_remove(Input_Metabolite, paste0("_", Molecular_Weight, "_", Retention_Time, "$"))
  )

# Fix odd encodings / leading underscores
metab_df <- metab_df %>%
  mutate(
    Metabolite_Name = str_replace(Metabolite_Name, "^_2-cis-Hexadecenoic acid$", "2-cis-Hexadecenoic acid"),
    Metabolite_Name = str_replace(Metabolite_Name, "\\\\xb1", "±"),
    Metabolite_Name = str_trim(Metabolite_Name)
  )

# ------------------------------------------------------------
# 4. Manual overrides for metabolites not matched in the reference table
#    These were manually curated based on HMDB/PubChem classification. for metabolites not in reference table
# ------------------------------------------------------------

manual_map <- tribble(
  ~Metabolite_Name, ~Manual_Category,
  "(S)-2-Methylbutanoic acid", "Organic Acids",
  "(R)-3-Hydroxyisobutyric acid", "Organic Acids",
  "(2E,4E)-2,4-Hexadienoic acid", "Fatty Acids",
  "Iminodiacetic acid", "Organic Acids",
  "Salicylic acid", "Miscellaneous",
  "3-Hexenedioic acid", "Organic Acids",
  "3-Oxoglutaric acid", "Organic Acids",
  "Adipic acid", "Organic Acids",
  "2,6-Dihydroxybenzoic acid", "Miscellaneous",
  "Orotic acid", "Nucleotides",
  "2-Methylbutyrylglycine", "Amino Acids",
  "m-Coumaric acid", "Miscellaneous",
  "Phthalic acid", "Organic Acids",
  "Quinolinic acid", "Organic Acids",
  "cis-Aconitic acid", "Organic Acids",
  "Hippuric acid", "Organic Acids",
  "Capryloylglycine", "Amino Acids",
  "(R)-2-Benzylsuccinate", "Organic Acids",
  "1,11-Undecanedicarboxylic acid", "Organic Acids",
  "(±)9-HpODE", "Fatty Acids",
  "(+/-)-C75", "Miscellaneous",
  "2-Naphthalenesulfonic acid", "Sulfonic Acids"
)

# ------------------------------------------------------------
# 5. Classify metabolites
# ------------------------------------------------------------

classified_df <- metab_df %>%
  left_join(
    category_ref %>% select(Metabolite_Name, Category),
    by = "Metabolite_Name"
  ) %>%
  left_join(manual_map, by = "Metabolite_Name") %>%
  mutate(
    # Merge: prefer original table category, then manual override, then mark UNMAPPED.
    Category = coalesce(Category, Manual_Category, "UNMAPPED"),
    Mapping_Source = case_when(
      !is.na(Category) & Category != "UNMAPPED" & is.na(Manual_Category) ~ "Original category table",
      !is.na(Manual_Category) ~ "Manual extension",
      TRUE ~ "Unmapped"
    )
  ) %>%
  select(
    Input_Metabolite,
    Metabolite_Name,
    Molecular_Weight,
    Retention_Time,
    Category,
    Mapping_Source
  )

# Save full mapping
write_csv(classified_df, "metabolite_207_full_category_mapping.csv")

# ------------------------------------------------------------
# Any metabolites still UNMAPPED after manual curation are assigned to Fatty Acids
# based on structural similarity review.
# Reassign unmapped metabolites to "Fatty Acids"
# ------------------------------------------------------------

classified_df <- classified_df %>%
  mutate(
    Category = ifelse(Category == "UNMAPPED",
                      "Fatty Acids",
                      Category),
    
    Mapping_Source = ifelse(Mapping_Source == "Unmapped",
                            "Reassigned to Fatty Acids",
                            Mapping_Source)
  )

# Optional: save metabolites that were originally unmapped
classified_df %>%
  filter(Mapping_Source == "Reassigned to Fatty Acids") %>%
  write_csv("metabolites_reassigned_to_fatty_acids.csv")

# ------------------------------------------------------------
# 6. Category summary
# ------------------------------------------------------------

category_order <- c(
  "Fatty Acids",
  "Amino Acids",
  "Organic Acids",
  "Miscellaneous",
  "Nucleotides",
  "Sugars",
  "Sulfonic Acids"
)

summary_df <- classified_df %>%
  filter(Category != "UNMAPPED") %>%
  count(Category, name = "Count") %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category) %>%
  mutate(
    Percent = Count / sum(Count) * 100,
    pct_label = sprintf("%.1f%%", Percent),
    Label = case_when(
      Category == "Fatty Acids" ~ "Fatty\nAcids",
      Category == "Amino Acids" ~ "Amino\nAcids",
      Category == "Organic Acids" ~ "Organic\nAcids",
      Category == "Sulfonic Acids" ~ "Sulfonic\nAcids",
      TRUE ~ as.character(Category)
    ),
    ymax = cumsum(Count),
    ymin = lag(ymax, default = 0),
    ypos = (ymax + ymin) / 2
  )

write_csv(summary_df, "metabolite_207_category_summary.csv")

# ------------------------------------------------------------
# 7. Recreate pie chart
# ------------------------------------------------------------

fill_cols <- c(
  "Fatty Acids"    = "#9bbcd0",
  "Amino Acids"    = "#2c7fb8",
  "Organic Acids"  = "#a8c97f",
  "Miscellaneous"  = "#38a547",
  "Nucleotides"    = "#f2a7a7",
  "Sugars"         = "#f31621",
  "Sulfonic Acids" = "#efb45a"
)

p <- ggplot(summary_df, aes(x = 1, y = Count, fill = Category)) +
  geom_col(width = 1, color = "black", linewidth = 0.5) +
  coord_polar(theta = "y") +
  
  # Inner labels: percentage of total for each slice.
  geom_text(
    aes(y = ypos, label = pct_label),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  
  # Outer labels: category name positioned just outside the pie radius.
  geom_text(
    aes(x = 1.28, y = ypos, label = Label),
    size = 6,
    hjust = 0,
    lineheight = 0.9,
    color = "black"
  ) +
  
  scale_fill_manual(values = fill_cols) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "grey85", color = NA),
    panel.background = element_rect(fill = "grey85", color = NA),
    plot.margin = margin(20, 60, 20, 20)
  ) +
  xlim(0.5, 1.65)

print(p)

# ------------------------------------------------------------
# 8. Save figure
# ------------------------------------------------------------

ggsave(
  "1A_metabolite_207_category_pie_chart.svg",
  p,
  width = 10,
  height = 7.5,
  bg = "grey85"
)

ggsave(
  "1A_metabolite_207_category_pie_chart.pdf",
  p,
  width = 10,
  height = 7.5,
  bg = "grey85"
)

ggsave(
  "1A_metabolite_207_category_pie_chart.png",
  p,
  width = 10,
  height = 7.5,
  dpi = 300,
  bg = "grey85"
)

# ------------------------------------------------------------
# 9. Print summary
# ------------------------------------------------------------

print(summary_df %>% select(Category, Count, Percent))