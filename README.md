# Code Repository — MSS RAS-Mutant CRC Metabolic Reprogramming
### Sadanandam Lab | Cancer Cell Submission

This repository contains all analysis scripts used to generate **Figure 2**, **Supplementary Figure 4**, and **Supplementary Figure 5** of the manuscript:

> *Integrated metabolic and transcriptomic profiling reveals alanine–aspartate–glutamate pathway reprogramming in MSS RAS-mutant colorectal cancer*

---

## Repository Structure

```
publication_scripts/
├── Figure_2/
│   ├── Fig2A/          Metabolite category pie chart (207 metabolites)
│   ├── Fig2B/          Metabolite–metabolite correlation heatmap
│   ├── Fig2D/          MetaboAnalyst pathway impact scatter
│   ├── Fig2E/          Metabolomics vs RNA-seq ssGSEA concordance
│   ├── Fig2F/          Metabolomics vs proteomics ssGSEA concordance
│   ├── Fig2H/          MuSiC deconvolution + regression (glutamine–glutamate)
│   ├── Fig2I_J_M_O/    Seurat scRNA-seq enrichment + cell-type expression
│   ├── Fig2K_L/        Glutamine–glutamate scatter and violin plots
│   ├── Fig2P/          (results only — no dedicated script)
│   └── Fig2Q/          Metabolite–RPPA protein correlation
├── Supp_Figure_4/
│   ├── S4A/            Pathway scoring – COLOSSUS RNA-seq (singscore + GSVA)
│   ├── S4B/            Public COAD cohort (GSE89076)
│   ├── S4C/            TCGA – RAS-mutant
│   ├── S4D/            TCGA – RAS-wildtype
│   ├── S4E/            Clinical proxy cohort
│   ├── S4F/            PDX (LMX/BASALE)
│   ├── S4G/            Organoids – RAS-mutant (all)
│   ├── S4H/            Organoids – RAS-wildtype (all)
│   ├── S4I/            Organoids – 9 selected RAS-mutant lines
│   └── S4J/            Organoids – 3 selected RAS-wildtype lines
└── Supp_Figure_5/
    ├── S5A_to_S5D/     COLOSSUS snRNA-seq (RAS-WT, 4 samples)
    ├── S5E_to_S5H/     GSE200997 scRNA-seq validation
    ├── S5I_to_S5L/     GSE178318 scRNA-seq validation
    ├── S5M_N/          GSE89076 glutamine–glutamate validation
    └── S5P_Q_R/        EIF4B expression (all three cohorts)
```

---

## Figure–Script Quick Reference

| Panel | Script | Description |
|-------|--------|-------------|
| **2A** | `Fig2A/2A_metabolite_category.R` | Metabolite biochemical category classification + pie chart (n=207) |
| **2B** | `Fig2B/2B_metabolite_metabolite_correlation_heatmap.Rmd` | Pairwise metabolite Pearson correlation heatmap (k=4 clusters) |
| **2D** | `Fig2D/2D_..._metaboanalyst_replot_patient_HMDB.R` | MetaboAnalyst pathway impact scatter (–log10 FDR vs impact) |
| **2E** | `Fig2E/2E_..._ssGSEA_mean_metabolanlyst_FDR_top_pathways_figures.R` | Metabolomics vs RNA-seq ssGSEA — scatter + boxplot |
| **2E** | `Fig2E/2E_S4_S5N_..._GSVA--ssgseaParam_Anguraj_modified.R` | Multi-cohort GSVA/ssGSEA (also generates S4, S5N panels) |
| **2F** | `Fig2F/2F_..._ssGSEA_protein_mean_metabolanlyst_FDR_top_pathways_figures.R` | Metabolomics vs proteomics ssGSEA — scatter + boxplot |
| **2H** | `Fig2H/2H_..._deconvolution_..._expr_corrected.R` | MuSiC cell-type deconvolution + regression (RNA-seq) |
| **2H** | `Fig2H/2H2_protein_..._protein.R` | MuSiC deconvolution + regression (proteomics) |
| **2I** | `Fig2I_J_M_O/2I_J_M_O_..._seurat_ssgsea_metabolic_pathways_v2.0.R` | Geyser enrichment plot (KEGG Ala/Asp/Glu, scRNA-seq) |
| **2J** | same script | GLS, GLUL cell-type expression violin |
| **2M** | same script | TCGA-selected gene expression by cell type |
| **2O** | same script | EIF4B expression by cell type |
| **2I** (cell annotation) | `Fig2I_J_M_O/2I_J_M_O_..._seurat_scATOMIC_Anguraj.R` | scATOMIC cell-type annotation (run before ssGSEA script) |
| **2K** | `Fig2K_L/2K_L_glutamine_glutamate_scatter_violin_plots.Rmd` | Glutamine vs glutamate Pearson scatter |
| **2L** | same script | Glutamine vs glutamate violin |
| **2Q** | `Fig2Q/2Q_..._COLOSSUS_metabolites_RPPA_correlation_analysis.Rmd` | Metabolite vs RPPA protein scatter + regression |
| **S4A** | `Supp_Figure_4/S4A/S4A_new_script_singscore_analysis.R` | singscore pathway scoring – COLOSSUS RNA-seq |
| **S4A** (alt) | `Supp_Figure_4/S4A/S4A_backup_script_..._GSVA_ssGSEA.R` | GSVA alternative for S4A (methods comparison) |
| **S4B** | `Supp_Figure_4/S4B/S4B_..._Public_COAD_metabolism_ssGSEA.R` | Public COAD cohort validation (GSE89076) |
| **S4C** | `Supp_Figure_4/S4C/S4C_..._TCGA_RAS_mutant_GSVA_ssGSEA.R` | TCGA RAS-mutant ssGSEA |
| **S4D** | `Supp_Figure_4/S4D/S4D_..._TCGA_RAS_wildtype_GSVA_ssGSEA.R` | TCGA RAS-wildtype ssGSEA |
| **S4E** | `Supp_Figure_4/S4E/S4E_..._clinical_Proxy_RAS_mutant.R` | Clinical proxy cohort ssGSEA |
| **S4F** | `Supp_Figure_4/S4F/S4F_..._LMX_BASALE_all_samples.R` | PDX (LMX/BASALE) ssGSEA |
| **S4G** | `Supp_Figure_4/S4G/S4G_..._Organoids_RAS_mutant.R` | Organoids RAS-mutant (all lines) |
| **S4H** | `Supp_Figure_4/S4H/S4H_..._Organoids_RAS_wildtype.R` | Organoids RAS-wildtype (all lines) |
| **S4I** | `Supp_Figure_4/S4I/S4I_..._Organoids_RAS_mutant.R` | Organoids – 9 selected RAS-mutant lines |
| **S4I/J** (QC) | `Supp_Figure_4/S4I/S4I_J_..._Organoid_metabolomics_QC.R` | Organoid metabolomics QC & normalisation (run first) |
| **S4J** | `Supp_Figure_4/S4J/S4J_..._Organoids_RAS_wildtype.R` | Organoids – 3 selected RAS-wildtype lines |
| **S5A–D** | `Supp_Figure_5/S5A_to_S5D/S3E_F_S10_S5_..._COLOSSUS_snRNAseq.R` | COLOSSUS snRNA-seq (RAS-WT) |
| **S5E–H** | `Supp_Figure_5/S5E_to_S5H/` (same script, GSE200997 input) | GSE200997 validation cohort |
| **S5I–L** | `Supp_Figure_5/S5I_to_S5L/S3I_J_S10_S5_..._GSE178318_v2.0.R` | GSE178318 validation cohort |
| **S5M, S5N** | `Supp_Figure_5/S5M_N/S5M_N_..._GSE89076_glutamine_glutamate.R` | GSE89076 Glu/Gln axis validation |
| **S5P, Q, R** | `Supp_Figure_5/S5P_Q_R/` | EIF4B expression – all three sc/snRNA-seq cohorts |

---

## Requirements

### R version
Developed and tested with **R ≥ 4.2**. 

### R Packages

**CRAN packages**
```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "ggrepel", "cowplot",
  "pheatmap", "RColorBrewer", "corrplot", "scales", "ggpubr",
  "readxl", "readr", "stringr"   # required for Fig 2A
))
```

**Bioconductor packages**
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GSVA",          # ssGSEA / GSVA
  "GSEABase",      # read .gmt files
  "singscore",     # single-sample scoring (S4A)
  "escape",        # single-cell ssGSEA (Fig2I, S5 panels)
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scran",
  "MuSiC"         # cell-type deconvolution (Fig2H)
))
```

**Other**
```r
install.packages("qusage")   # read .gmt files for single-cell scripts
install.packages("Seurat")   # scRNA-seq analysis
```

---

## Input Data

| Data type | Source | Access |
|-----------|--------|--------|
| COLOSSUS metabolomics | Batch-corrected in-house dataset | Available from authors on request |
| COLOSSUS RNA-seq | Combat-corrected VST matrix | Available from authors on request |
| COLOSSUS RPPA | Normalized proteomics data | Available from authors on request |
| COLOSSUS scRNA-seq | Seurat object (MSS KRAS-mutant) | Available from authors on request |
| COLOSSUS snRNA-seq | Seurat object (RAS-WT, 4 samples) | Available from authors on request |
| KEGG gene sets | `c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt` | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) |
| TCGA COAD | HiSeqV2 expression | [GDC Data Portal](https://portal.gdc.cancer.gov/) |
| GSE89076 | Agilent 8×60K microarray + metabolomics | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89076) |
| GSE200997 | scRNA-seq | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997) |
| GSE178318 | scRNA-seq | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178318) |
| GSE204805 | PDX/organoid RNA-seq | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204805) |

> **Note on file paths:** All scripts use absolute file paths. Before running, update the paths at the top of each script to match your local directory structure. Paths follow the pattern `~/SPM Dropbox/SysPreMed/COLOSSUS/...`.

---

## Recommended Execution Order

For each figure, scripts should be run in this order:

### Figure 2
1. `Fig2B` — metabolite correlation heatmap (standalone)
2. `Fig2D` — MetaboAnalyst scatter (requires MetaboAnalyst output CSV)
3. `Fig2E` — ssGSEA RNA concordance (requires RNA-seq matrix + .gmt)
4. `Fig2F` — ssGSEA protein concordance (requires proteomics matrix + .gmt)
5. `Fig2H` — MuSiC deconvolution (requires scRNA-seq Seurat object)
6. `Fig2I_J_M_O` — run scATOMIC script first, then ssGSEA script
7. `Fig2K_L` — glutamine/glutamate plots (standalone)
8. `Fig2Q` — RPPA correlation (standalone)

### Supplementary Figure 4
1. `S4I/S4I_J_..._Organoid_metabolomics_QC.R` — run this first
2. All other S4 scripts are independent of each other

### Supplementary Figure 5
- All S5 scripts are independent; set the Seurat object input path to the relevant cohort.

---

## Session Information

To record your R session for reproducibility, add this at the end of any script:
```r
sessionInfo()
```
Or use `renv::snapshot()` to capture the full package lock file.

---

## Authors

- **Anguraj Sadanandam** (Group Leader and lead bioinformatics analysis)
- **Rachel** (PhD student and additional analysis and visualisation)

Institute for Cancer Research, London / Systems & Precision Medicine Group

---

## Citation

If you use these scripts, please cite:

> [Citation to be added upon acceptance]

---

## License

Code is released under the **ICR License**. Data availability is subject to the data sharing policies of the respective cohorts (see Input Data table above).
