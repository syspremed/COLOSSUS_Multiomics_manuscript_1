################################################################################
# Figure Panel : 1G
# Description  : KEGG pathway enrichment of genes in significant amplified and deleted CNA peaks.
# Input data   : 1H_20260320_20240719_CNA_all_COLOSSUS_amplified_genes_only_sig_0.2_KEGG_2026_table-2.txt
#                1H_20260320_20240722_CNA_present_in_COLOSSUS_deleted_genes_only_sig_0.2_KEGG_2026_table-2.txt
# Output       : 2H_20260320_CNA_AMP_DEL_signifiant_only_KEGG_enrichment.pdf
# R packages   : ggplot2
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : March 2026
# Note         : Some output file names in the code keep draft (pre-publication) panel labels.
################################################################################

#Figure 1I

library(ggplot2)

publication_theme = theme(
  title = element_text(size = 20),
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14)
)

kegg_amp_table_for_plot = read.delim("1H_20260320_20240719_CNA_all_COLOSSUS_amplified_genes_only_sig_0.2_KEGG_2026_table-2.txt", sep = "\t", header = T)
kegg_del_table_for_plot = read.delim("1H_20260320_20240722_CNA_present_in_COLOSSUS_deleted_genes_only_sig_0.2_KEGG_2026_table-2.txt", sep = "\t", header = T)

kegg_amp_table_for_plot$cnv_type = "AMP"
kegg_del_table_for_plot$cnv_type = "DEL"

kegg_both_table_for_plot = bind_rows(kegg_amp_table_for_plot, kegg_del_table_for_plot)

# ggplot(data=kegg_amp_table_for_plot, aes(x=Overlap/Geneset, y=Label, fill=FDR))+
#   geom_bar(stat="identity")+
#   scale_fill_continuous(type = "viridis")+
#   theme_bw()
# 
# ggplot(data=kegg_del_table_for_plot, aes(x=Overlap/Geneset, y=Label, fill=FDR))+
#   geom_bar(stat="identity")+
#   scale_fill_continuous(type = "viridis")+
#   theme_bw()

#kegg_both_table_for_plot$Label = str_replace_all(str_replace_all(kegg_both_table_for_plot$Label, "KEGG_", ""), "_", " ")

kegg_both_table_for_plot$Label <- kegg_both_table_for_plot$Term
kegg_both_table_for_plot$gene_ratio <- sapply(kegg_both_table_for_plot$Overlap, function(x) eval(parse(text = x)))

del <- which(kegg_both_table_for_plot$cnv_type == "DEL")

kegg_both_table_for_plot_10_amp <- kegg_both_table_for_plot[1:10,]
kegg_both_table_for_plot_10_del <- kegg_both_table_for_plot[del[1:10],]
kegg_both_table_for_plot_10_amp_del <- rbind(kegg_both_table_for_plot_10_amp,kegg_both_table_for_plot_10_del)

pdf("2H_20260320_CNA_AMP_DEL_signifiant_only_KEGG_enrichment.pdf")
ggplot(data=kegg_both_table_for_plot_10_amp_del, aes(x=gene_ratio, y=Label, fill=Adjusted.P.value))+
  geom_bar(stat="identity")+
  scale_fill_continuous(type = "viridis")+
  facet_grid(~cnv_type)+
  labs(y="", x="Proportion of geneset altered")+
  theme_bw()+publication_theme
dev.off()