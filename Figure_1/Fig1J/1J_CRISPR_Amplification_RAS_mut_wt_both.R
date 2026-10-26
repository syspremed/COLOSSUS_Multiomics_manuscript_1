################################################################################
# Figure Panel : 1J
# Description  : DepMap CRISPR gene effect of genes in amplified/deleted CNA peaks across CRC
#                MSS RAS-mutant and RAS-wildtype cell lines.
# Input data   : CRISPRGeneEffect.csv
#                amp_genes.conf_5.txt
#                all_lesions.conf_5.txt
#                del_genes.conf_5.txt
#                Fig1L_amplified_genes_for_depmap.csv
#                Fig1L_deleted_genes_for_depmap.csv
# Output       : 20260407_filtered_CRISPR_CRC_MSS_RAS_mt_wt_both_data.csv
# R packages   : dplyr, stringr, tidyr, ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : April 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

# Load necessary library
library(dplyr)
library(stringr)
library(tidyr)

# Read DepMap CRISPR gene effect data
data <- read.csv("CRISPRGeneEffect.csv", header=TRUE, row.names=1)

#RAS mutant cell lines
crc_cell_lines <- c(
  "LS-513", "CL-11","LS1034","COLO-678","T84","SKCO1","NCI-H747","SW837","HCC56","SW1463","SW1116",
  "LS-123","SNU-61","RCM-1", "SW620", "SW948","CL-40","SNU-81","HCC2998","SW626", "SNU-1033", "SW403",
  "C2BBE1",	"COLO320",	"SW1417",	"NCIH508",	"NCIH716",	"HT29",	"SNU283",	"SNUC1",	"HT55",	"MDST8",	"HT115",	"COLO205",	"DIFI",	"CAR1")


# Define the target ACH IDs
target_ids <- c("ACH-000007", "ACH-000249", "ACH-000252", "ACH-000350",
                "ACH-000381", "ACH-000400", "ACH-000403", "ACH-000421", "ACH-000467", "ACH-000470",
                "ACH-000489", "ACH-000501", "ACH-000532", "ACH-000565", "ACH-000651",
                "ACH-000680", "ACH-000798", "ACH-000991", "ACH-001081", 
                "ACH-001399", "ACH-000286", "ACH-000820", "ACH-000009",	"ACH-000202",	"ACH-000236",	"ACH-000360",	
                "ACH-000491",	"ACH-000552",	"ACH-000708",	"ACH-000722",	"ACH-000926",	"ACH-000935",	"ACH-000986",	
                "ACH-001039",	"ACH-002233",	"ACH-002345")
 
# Filter the data
filtered_data <- data[rownames(data) %in% target_ids, ]

# View the filtered data
print(filtered_data)

# Optionally, save the filtered data to a new file
write.csv(filtered_data, "20260407_filtered_CRISPR_CRC_MSS_RAS_mt_wt_both_data.csv", row.names=TRUE)

crispr_data <- filtered_data

amp_raw = read.table("amp_genes.conf_5.txt", sep = "\t", header = T)
amp_raw$cytoband=NULL
amp_raw$X = NULL
amp_raw_residual_significant = amp_raw[,as.numeric(amp_raw[1,])<0.2]
amp_raw_residual_significant[3,] = str_replace(str_replace(amp_raw_residual_significant[3,], ":", "_"), "-", "_")

cnv_colossus = read.table("all_lesions.conf_5.txt", sep = "\t", header = T)
cnv_colossus = cnv_colossus[!grepl(pattern = "CN values", x = cnv_colossus$Unique.Name),]
rownames(cnv_colossus) = cnv_colossus$Unique.Name

tmp_cnv_colossus_amp = cnv_colossus[grepl("Amp",cnv_colossus$Unique.Name), c("Wide.Peak.Limits", colnames(cnv_colossus)[grepl("sRCSI",colnames(cnv_colossus))])]
tmp_cnv_colossus_amp$Wide.Peak.Limits = str_replace(str_replace(str_remove_all(str_remove_all(tmp_cnv_colossus_amp$Wide.Peak.Limits, "\\(.+\\)"), " "), ":", "_"), "-", "_")
rownames(tmp_cnv_colossus_amp) = tmp_cnv_colossus_amp$Wide.Peak.Limits
tmp_cnv_colossus_amp$Wide.Peak.Limits = NULL
tmp_cnv_colossus_amp[tmp_cnv_colossus_amp==2]=1
tmp_cnv_colossus_amp_025 = tmp_cnv_colossus_amp[(rowSums(tmp_cnv_colossus_amp)/ncol(tmp_cnv_colossus_amp))>0.25,]

amp_genes_list = unique(unlist(amp_raw_residual_significant[4:nrow(amp_raw_residual_significant),amp_raw_residual_significant[3,] %in% rownames(tmp_cnv_colossus_amp_025)]))
amp_genes_list = amp_genes_list[amp_genes_list != ""]
amp_genes_list = sort(amp_genes_list[!grepl("\\[", amp_genes_list)])

del_raw = read.table("del_genes.conf_5.txt", sep = "\t", header = T)
del_raw$cytoband=NULL
del_raw$X = NULL
del_raw_residual_significant = del_raw[,as.numeric(del_raw[1,])<0.1]
del_raw_residual_significant[3,] = str_replace(str_replace(del_raw_residual_significant[3,], ":", "_"), "-", "_")

tmp_cnv_colossus_del = cnv_colossus[grepl("Del",cnv_colossus$Unique.Name), c("Wide.Peak.Limits", colnames(cnv_colossus)[grepl("sRCSI",colnames(cnv_colossus))])]
tmp_cnv_colossus_del$Wide.Peak.Limits = str_replace(str_replace(str_remove_all(str_remove_all(tmp_cnv_colossus_del$Wide.Peak.Limits, "\\(.+\\)"), " "), ":", "_"), "-", "_")
rownames(tmp_cnv_colossus_del) = tmp_cnv_colossus_del$Wide.Peak.Limits
tmp_cnv_colossus_del$Wide.Peak.Limits = NULL
tmp_cnv_colossus_del[tmp_cnv_colossus_del==2]=1
tmp_cnv_colossus_del_025 = tmp_cnv_colossus_del[(rowSums(tmp_cnv_colossus_del)/ncol(tmp_cnv_colossus_del))>0.25,]

del_genes_list = unique(unlist(del_raw_residual_significant[4:nrow(del_raw_residual_significant),del_raw_residual_significant[3,] %in% rownames(tmp_cnv_colossus_del_025)]))
del_genes_list = del_genes_list[del_genes_list != ""]
del_genes_list = sort(del_genes_list[!grepl("\\[", del_genes_list)])

crispr_amp_genes = read.table("Fig1L_amplified_genes_for_depmap.csv", sep = ",", header = T)
crispr_del_genes = read.table("Fig1L_deleted_genes_for_depmap.csv", sep = ",", header = T)

# Extract strings before ".."
colnames(crispr_data) <- sub("\\.\\..*", "", colnames(crispr_data))

# Check result
head(colnames(crispr_data))

crispr_data_amp = crispr_data[,colnames(crispr_data) %in% amp_genes_list]
crispr_data_del = crispr_data[,colnames(crispr_data) %in% del_genes_list]

crispr_data_amp = crispr_data_amp[,apply(crispr_data_amp, 2, mean, na.rm=T) < (-0.25)]
crispr_data_del = crispr_data_del[,apply(crispr_data_del, 2, mean, na.rm=T) < (-0.25)]

tmp_crispr_sort_order = c(names(sort(apply(crispr_data_amp, 2, mean, na.rm=T), decreasing = T)),
                          names(sort(apply(crispr_data_del, 2, mean, na.rm=T), decreasing = T)))
tmp_crispr_colour = c(ifelse(sort(apply(crispr_data_amp, 2, mean, na.rm=T), decreasing = T)<(-1), "#DD2222", "#FFAA00"),
                      ifelse(sort(apply(crispr_data_del, 2, mean, na.rm=T), decreasing = T)<(-1), "#DD2222", "#FFAA00"))

crispr_data_amp_long = pivot_longer(crispr_data_amp, cols = everything(), names_to = "gene", values_to = "CRISPR_score")
crispr_data_amp_long$CNA_type = "Amplified"
crispr_data_del_long = pivot_longer(crispr_data_del, cols = everything(), names_to = "gene", values_to = "CRISPR_score")
crispr_data_del_long$CNA_type = "Deleted"

crispr_data_both_long = bind_rows(crispr_data_amp_long, crispr_data_del_long)

crispr_data_both_long$gene = factor(crispr_data_both_long$gene, levels=tmp_crispr_sort_order)

library(ggplot2)

publication_theme = theme(
  title = element_text(size = 20),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14)
)
library(ggplot2)
ggplot(crispr_data_both_long, aes(x=gene, y=CRISPR_score))+
  geom_hline(yintercept = 0, lty=2, col="#00000077")+
  geom_hline(yintercept = (-0.25), lty=2, col="#FF880077")+
  geom_hline(yintercept = (-1), lty=2, col="#FF110077")+
  geom_boxplot(col="#000000", fill=tmp_crispr_colour, outlier.colour = "#00000077")+
  facet_grid(~CNA_type, scales = "free_x", space = "free_x")+
  labs(y="CRISPR score", x="Gene")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  publication_theme+
  theme(axis.text.x = element_text(size = 10, face = "bold"), panel.grid.minor.y = element_blank())

