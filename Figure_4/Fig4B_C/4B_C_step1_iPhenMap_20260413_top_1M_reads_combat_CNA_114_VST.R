################################################################################
# Figure Panel : 4B, 4C (iPhenMap model; also inputs for 4D, 4E)
# Description  : Runs iPhenMap semi-supervised multiomics integration (transcriptome, proteome,
#                metabolome and 114 GISTIC CNA regions) with clinical covariates on the
#                COLOSSUS St2/3 MSS RAS-mt training cohort (n = 98). Writes proportion of
#                variance, scaled loadings and scaled regression coefficients for each MV.
# Input data   : all_lesions.conf_5.txt
#                MV_clustering.txt
# Output       : checking_20200526-colossus_data_prognosis_calculated_AS.txt
#                20201014_20200604_colossus_restro_copynumber_matched_clinical_data_after_correction_centre_sequencing_lane.txt
#                20200607_colossus_retro_protein_data_log2_median_normalized_imputed_cohort_lysis_batch_corrected.txt
#                20200721_colossus_retro_metabolomics_data_log2_weight_cohort_run_batch_corrected.txt
#                20200717_Colossus_Normalisation_TissueWeight.txt
#                20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt
#                _up_regulated_genes.pdf
#                _down_regulated_genes.pdf
#                CNA_114.txt
# R packages   : cluster, survminer, fpc, ggplot2, scales, plotrix, psych, mclust, survival, digest, gridExtra, reshape2, devtools, RColorBrewer, plyr, dplyr, easyGgplot2, lattice, cowplot, ggrepel, ComplexHeatmap, circlize, colorspace, GetoptLong, Hmisc, MASS, MCMCpack, ConsensusClusterPlus, NMF, compiler, fastDummies, tidyverse, siggenes, pamr, glmnet, caret, parallel, pROC, foreach, doMC, e1071, extrafont, BiocManager, mixtools, AnnotationDbi
# Authors      : Sadanandam Lab, Institute of Cancer Research, London
# Contact      : anguraj.sadanandam@icr.ac.uk
# Date         : April 2026
# Note         : Step 1 of Fig4B_C. Some output file names in the code keep draft labels.
################################################################################

library(cluster)
library(survminer)
library(fpc)
library(ggplot2)
library(scales)
library(plotrix)
library(psych)
library(mclust)
library(survival)
library(digest)
library(gridExtra)
library(reshape2)
library(devtools)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(easyGgplot2)
library(lattice)
library(fpc)
library(cowplot)
library(ggrepel) 
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(Hmisc)
library(MASS)
library(MCMCpack)
require(ConsensusClusterPlus)
require(NMF)
library(compiler) 
library(dplyr)
library(fastDummies)
library(tidyverse)
library(compiler)
library(siggenes)
library(pamr)
library(glmnet)
library(caret)
library(parallel)
library(pROC)  
library(foreach)
library(doMC)
library(e1071)

## Theme for ggplot
theme <- theme(strip.text.y = element_text(),#rotate strip text to horizontal
               axis.text = element_text(colour = "black", family="Helvetica", size=18),
               axis.title.x = element_text(colour = "black", family="Helvetica",face="bold", size=20),
               axis.title.y = element_text(colour = "black", family="Helvetica",face="bold", size=20),
               axis.text.x = element_text(colour = "black", family="Helvetica", size=18),#remove y axis labels
               legend.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "white"),
               legend.title = element_text(colour = "black", family="Helvetica",face="bold", size=20),
               legend.text = element_text(colour="black", family="Helvetica", size = 18, face = "bold"),
               strip.text.x = element_text(colour = "black", family="Helvetica", size = 18),
               plot.title = element_text(size=20, family="Helvetica", face="bold"),
               panel.border=element_rect(colour="black",size=2)) #format title

extrafont::loadfonts()


#source('20201010_requiredFunctions_iphenmapPackage_v3_modified_with_comments.R')

filename1="checking_20200526-colossus_data_prognosis_calculated_AS.txt"
filename2="20201014_20200604_colossus_restro_copynumber_matched_clinical_data_after_correction_centre_sequencing_lane.txt"
filename3="20200607_colossus_retro_protein_data_log2_median_normalized_imputed_cohort_lysis_batch_corrected.txt"
filename4="20200721_colossus_retro_metabolomics_data_log2_weight_cohort_run_batch_corrected.txt"
filename44="20200717_Colossus_Normalisation_TissueWeight.txt"
#filename5='Results/RNAseq/preprocess/2020730_colossus_retro_rnaSeq_filtered_deseq2_normlalized_batch_corrected_tissue_weight.txt'
filename5="20260418_20200705_COLOSSUS - Task 4.1-RNA Lexogen_retrospective cohort_June2020_Exp_data_names_corrected_DESeq2_VST_filtered_1M_reads_percent_3_5_Combat_Centre_Lane_1M_reads_batches_corrected_no_reads_involved.txt"

clindata=read.delim(filename1,sep='\t',row.names = 2); dim(clindata)
cnadata=read.delim(filename2,sep='\t',row.names = 1); dim(cnadata)
protdata=read.delim(filename3,sep='\t',row.names = 1); dim(protdata)
metdata=read.delim(filename4,sep='\t',row.names = 1); dim(metdata)
mInfo=read.delim(filename44,sep='\t',row.names = 1); dim(mInfo)
rseqdata=read.delim(filename5,sep='\t',row.names = 1); dim(rseqdata)

## cnadata: matrix/data.frame with rownames like 1p34.2A, 1q22A, 1p36.31D, ...
## peaks:   GISTIC peak table (the second file)

peaks <- read.delim(
  "all_lesions.conf_5.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

peaks <- peaks[1:132, ]

## ============================================================
## Goal: keep ALL peaks with q < 0.2 INCLUDING duplicates
## from BOTH peaks_main and cnadata
## ============================================================

## 1. Clean inputs
colnames(peaks) <- trimws(colnames(peaks))
rownames(cnadata) <- trimws(rownames(cnadata))

## 2. Remove "- CN values"
peaks_main <- peaks[
  !grepl("-\\s*CN values$", peaks$`Unique Name`),
  , drop = FALSE
]

## 3. Filter q < 0.2
peaks_main$q_numeric <- suppressWarnings(as.numeric(peaks_main$`q values`))

peaks_main <- peaks_main[
  !is.na(peaks_main$q_numeric) & peaks_main$q_numeric < 0.2,
  , drop = FALSE
]

cat("Rows after q<0.2 filter:", nrow(peaks_main), "\n")

## 4. Build IDs
peaks_main$Descriptor_clean <- trimws(as.character(peaks_main$Descriptor))

peaks_main$peak_id <- ifelse(
  grepl("^Amplification Peak", peaks_main$`Unique Name`, ignore.case = TRUE),
  paste0(peaks_main$Descriptor_clean, "A"),
  ifelse(
    grepl("^Deletion Peak", peaks_main$`Unique Name`, ignore.case = TRUE),
    paste0(peaks_main$Descriptor_clean, "D"),
    NA_character_
  )
)

peaks_main <- peaks_main[!is.na(peaks_main$peak_id), , drop = FALSE]
a<-which(duplicated(peaks_main$peak_id) == TRUE)
peaks_main$peak_id[a-1] <- paste(peaks_main$peak_id[a],"-1",sep="")
peaks_main$peak_id[a] <- paste(peaks_main$peak_id[a],"-2",sep="")


## 5. Convert cnadata to data.frame with peak_id column
cnadata_df <- data.frame(
  peak_id = rownames(cnadata),
  cnadata,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dim(cnadata_df)

## 6. Add row indices so duplicates are preserved visibly
peaks_main$peaks_row_id <- seq_len(nrow(peaks_main))
cnadata_df$cnadata_row_id <- seq_len(nrow(cnadata_df))

mt <-match(cnadata_df$peak_id,peaks_main$peak_id)
wt <- which(!is.na(mt))
length(wt)

cnadata_df_m <- cnadata_df[wt,]
dim(cnadata_df_m)

## optional: inspect duplicate counts
#table(peaks_main$peak_id)
#table(cnadata_df$peak_id)


cnadata <- cnadata_df_m[,c(-1,-184)]
dim(cnadata)

#write.table(cnadata, "CNA_114.txt", sep = "\t", quote = FALSE)
 
#### process clinical data  #### 
#sdata=clindata[clindata$center=='VHIO',]; dim(sdata)
sdata=clindata; dim(sdata)
cdata = dplyr::select(sdata, -center); dim(cdata)
cdata <- cdata %>%  
  #mutate_if(funs(anyNA(.) & !is.numeric(.)), funs(replace(., is.na(.), 'Unknown'))) %>%
  dplyr::select(
    gender, 
   ECOG_at_diagnosis,
    histologic_type, 
  #  histologic_grade,
    lymphovascular_invasion,
    perineural_invasion,
    #bowel_obstruction_at_diagnosis, # removed since it has only one patient with event
    tumor_location,   # created a dummy variable here
   # positive_nodes,
    #total_lymph_nodes,
    stage_recodif,
    #Adjuvant_chemotherapy,
    status_relapse,
    status_fup,
    RFS_months,
    OS_months
  ) %>%  
  mutate(#histologic_grade = as.numeric(histologic_grade),
         tumor_location = fct_collapse(tumor_location, left = c('left',
                                                                'rectum'),
                                       right = 'right'))

cdat <- cdata %>%
  mutate_if(is.character, as.factor) %>%  
  fastDummies::dummy_cols(remove_most_frequent_dummy = TRUE) %>%
  select_if(is.numeric)  
rownames(cdat)=rownames(sdata) 


cdat <- cdat[complete.cases(cdat),]
dim(cdat)

#### match metabol data  #### 
m0=match(rownames(cdat),rownames(metdata))
w0=which(!is.na(m0)); length(w0)
metdat=metdata[m0[w0],]; dim(metdat)
cda=cdat[w0,]; dim(cda)

#### match copy number data  #### 
m1=match(rownames(cda),gsub('[.]','-',colnames(cnadata)))
w1=which(!is.na(m1))
mcna=cnadata[,m1[w1]]; dim(mcna)
mcdat=cda[w1,]; dim(mcdat)
metda=metdat[w1,]; dim(metda)

#### match protein data  #### 
m2=match(rownames(mcdat),gsub('[.]','-',colnames(protdata)))
w2=which(!is.na(m2))
prot=t(protdata[,m2[w2]]); dim(prot)
clin=mcdat[w2,]; dim(clin)
cnat=t(mcna[,w2]); dim(cnat)
metd=metda[w2,]; dim(metd)

#### match rnaseq data  #### 
m3=match(rownames(clin),gsub('[.]','-',colnames(rseqdata)))
w3=which(!is.na(m3))
rdat=t(rseqdata[,m3[w3]]); dim(rdat)
prot1=prot[w3,]; dim(prot1)
clin1=clin[w3,]; dim(clin1)
cnat1=cnat[w3,]; dim(cnat1)
metd1=metd[w3,]; dim(metd1)


#### protein selection  #### 
varp=apply(prot1,2,mad)
pro=prot1[,varp>1]; dim(pro) 
propca=prcomp(pro,scale=F); 
pairs(propca$x[,1:4],pch=16,col=1,cex=1)

#### mets selection  #### 
varm=apply(metd1,2,mad)
met=metd1[,varm>1]; dim(met) 
metpca=prcomp(met,scale=F); 
pairs(metpca$x[,1:4],pch=16,col=1,cex=1)

#### cna selection  #### 
cna=cnat1; dim(cna)
cnapca=prcomp(cna); pairs(cnapca$x[,1:4],pch=16,col=1,cex=1)

#### gene selection  ####
varg=apply(rdat,2,mad)
rseq=rdat[,varg>1]; dim(rseq) 
rseqpca=prcomp(rseq,scale=F); 
pairs(rseqpca$x[,1:4],pch=16,col=1,cex=1)

#### run iPhenMap  #### 

X = cbind(dplyr::select(clin1, c(-RFS_months,-OS_months,-status_fup_death))); dim(X)


D=list(rna=rseq,protein=pro,metabol=met,cna=cna)

#no CNA
#D=list(rna=rseq,protein=pro,metabol=met)


Vp=sum(unlist(lapply(D,ncol)))/(5*length(D))
C2p=5

set.seed(131)

iphen(D=D,
      X,
      q=10,
      g=nrow(X),
      c=1,d=1,
      nu0=.01,Sigo=.01,
      s2o=.01,nu0k=.01,
      mcmc=10000,
      noiseMod='unique',
      scaleType='unit',
      iniType='EVD',
      time=NULL,
      status=NULL
) 

# MV analysis
mv <- read.delim(
  "MV_clustering.txt",
  stringsAsFactors = FALSE
)
dim(mv)

mv_ids_hyphen <- gsub("\\.", "-", mv$samples)
common_ids <- intersect(mv_ids_hyphen, rownames(clin1))

clin2 <- clin1[common_ids, , drop = FALSE]
dim(clin2)

m_mv <- match(common_ids, mv_ids_hyphen)
mv1 <- mv[m_mv, , drop = FALSE]
dim(mv1)

summary(lm(mv1$MV3~clin2$perineural_invasion_yes))

