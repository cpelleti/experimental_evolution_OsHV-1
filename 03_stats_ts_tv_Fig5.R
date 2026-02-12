#!/usr/bin/Rscript
##title: "Plot transition and transversions experimental evolution"
##author: "Pelletier Camille"
#date: 20/09/2022

'%!in%' <- function(x,y)!('%in%'(x,y))
options( "digits"=3)

###################################### Function to install / load package if it's not here ###############################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(unique(
  c( "tidyverse", "dplyr","patchwork", "ggnewscale") # plot
) )

addline_format <- function(x,...){
  gsub('//s','/n',x)
}

#------------------------------------------------------------------import datas-------------------------------------------------------------------------------
subst_table<- read.csv2("./data_frame_allelic_frequencies_indiduals.csv")

#------------------------------------------------------------------prepare datas------------------------------------------------------------------------------

subst_table<- subst_table[subst_table$TYPE == "snp",]

#determine if snps are ancestral or de novo
# considered as ancestral if more than one individual carry the snp at G0 and considered as de novo mutations if less than one individual carry the snp at G0 
for (i in 1:nrow(subst_table)){
  
  line_i <- subst_table[i,]
  
  ref_all <- line_i$REF
  alt_all <- line_i$ALT
  
  subst_table$GT3_G0[i] <- sum(line_i[c("GT3_G0_R1_1", "GT3_G0_R1_4", "GT3_G0_R1_6","GT3_G0_R2_1", "GT3_G0_R2_3", "GT3_G0_R2_6")] > 0)
  subst_table$GT5_G0[i] <- sum(line_i[c("GT5_G0_R1_3", "GT5_G0_R1_4", "GT5_G0_R1_6", "GT5_G0_R2_2", "GT5_G0_R2_3", "GT5_G0_R2_4")] > 0)
  subst_table$GT9_G0[i] <- sum(line_i[c("GT9_G0_R1_1", "GT9_G0_R1_3", "GT9_G0_R1_6", "GT9_G0_R2_2", "GT9_G0_R2_3", "GT9_G0_R2_6")] > 0)
  
  if (subst_table$GT3_G0[i] == 0 | subst_table$GT3_G0[i] == 1 ){
    subst_table$GT3_G0_ancestral[i] <- "denovo"
  } else if (subst_table$GT3_G0[i] >1 ){
    subst_table$GT3_G0_ancestral[i] <- "ancestral"
  } 
  
  if (subst_table$GT5_G0[i] == 0 | subst_table$GT5_G0[i] == 1 ){
    subst_table$GT5_G0_ancestral[i] <- "denovo"
  } else if (subst_table$GT5_G0[i] >1 ){
    subst_table$GT5_G0_ancestral[i] <- "ancestral"
  } 
  
  if (subst_table$GT9_G0[i] == 0 | subst_table$GT9_G0[i] == 1 ){
    subst_table$GT9_G0_ancestral[i] <- "denovo"
  } else if (subst_table$GT9_G0[i] >1 ){
    subst_table$GT9_G0_ancestral[i] <- "ancestral"
  } 
  
  
  if ((ref_all == "A" && alt_all == "G") || (ref_all == "G" && alt_all == "A") || (ref_all == "C" && alt_all == "T") || (ref_all == "T" && alt_all == "C")){
    subst_table$tstv[i] <- "transition"
  } else {
    subst_table$tstv[i] <- "transversion"
  }
  
  }

write.csv(subst_table,)
# create a subtable for ancestral/de novo snps 
ancestral_GT3_snps <- subst_table[subst_table$GT3_G0_ancestral == "ancestral", ] %>% select(UniqueID:GT3_G5_R2_7,GT3_G0,GT3_G0_ancestral,tstv)
denovo_GT3_snps <- subst_table[subst_table$GT3_G0_ancestral == "denovo", ]%>% select(UniqueID:GT3_G5_R2_7,GT3_G0,GT3_G0_ancestral,tstv)
ancestral_GT5_snps <- subst_table[subst_table$GT5_G0_ancestral == "ancestral", ]%>% select(UniqueID:TYPE,GT5_G0_R1_3:GT5_G5_R2_7,GT5_G0, GT5_G0_ancestral,tstv)
denovo_GT5_snps <- subst_table[subst_table$GT5_G0_ancestral == "denovo", ]%>% select(UniqueID:TYPE,GT5_G0_R1_3:GT5_G5_R2_7,GT5_G0, GT5_G0_ancestral,tstv)
ancestral_GT9_snps <- subst_table[subst_table$GT9_G0_ancestral == "ancestral", ]%>% select(UniqueID:TYPE,GT9_G0_R1_1:GT9_G5_R2_7,GT9_G0, GT9_G0_ancestral,tstv)
denovo_GT9_snps <- subst_table[subst_table$GT9_G0_ancestral == "denovo", ]%>% select(UniqueID:TYPE,GT9_G0_R1_1:GT9_G5_R2_7,GT9_G0, GT9_G0_ancestral,tstv)


##create a dataframe of substitution number per library
mat_tstv_ancestral_GT3 <- matrix(data=NA,nrow=144,ncol=5)
mat_tstv_ancestral_GT5 <- matrix(data=NA,nrow=180,ncol=5)
mat_tstv_ancestral_GT9 <- matrix(data=NA,nrow=252,ncol=5)
mat_tstv_denovo_GT3 <- matrix(data=NA,nrow=144,ncol=5)
mat_tstv_denovo_GT5 <- matrix(data=NA,nrow=180,ncol=5)
mat_tstv_denovo_GT9 <- matrix(data=NA,nrow=252,ncol=5)

counter <- 0
cp_names<- 0
for (j in 6:(ncol(ancestral_GT3_snps)-3)){
  
  ##colnames 
  mat_names <- colnames(ancestral_GT3_snps)[6:29]
  cp_names <- cp_names+1
  #initializing substitution type

  subst_df_j <- cbind(ancestral_GT3_snps[,1:5],as.numeric(ancestral_GT3_snps[,j]))
  
  AT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "A" ),6] >0))
  AG_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "G" ) | (subst_df_j$REF == "G" & subst_df_j$ALT == "A" ),6] >0))  
  AC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "A" ),6] >0))   
  GT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "G" ),6] >0)) 
  GC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "G" ),6] >0))
  TC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "T" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "T" ),6] >0))   
  
  
  mat_tstv_ancestral_GT3[(counter+1):(counter+6),1] <- mat_names[cp_names]
  mat_tstv_ancestral_GT3[counter+1,2] <- "A|C"
  mat_tstv_ancestral_GT3[counter+1,3] <- AC_snps
  mat_tstv_ancestral_GT3[counter+1,4] <- "#273572"
  mat_tstv_ancestral_GT3[counter+1,5] <-"Transversion"
  
  mat_tstv_ancestral_GT3[counter+2,2] <- "A|T"
  mat_tstv_ancestral_GT3[counter+2,3] <- AT_snps
  mat_tstv_ancestral_GT3[counter+2,4] <- "#6f91cb"
  mat_tstv_ancestral_GT3[counter+2,5] <-"Transversion"
  
  mat_tstv_ancestral_GT3[counter+3,2] <- "C|G"
  mat_tstv_ancestral_GT3[counter+3,3] <- GC_snps
  mat_tstv_ancestral_GT3[counter+3,4] <- "#8bd0df"
  mat_tstv_ancestral_GT3[counter+3,5] <-"Transversion"
 
  
  mat_tstv_ancestral_GT3[counter+4,2] <- "G|T"
  mat_tstv_ancestral_GT3[counter+4,3] <- GT_snps
  mat_tstv_ancestral_GT3[counter+4,4] <- "#a6d4ba"
  mat_tstv_ancestral_GT3[counter+4,5] <-"Transversion"
  
  mat_tstv_ancestral_GT3[counter+5,2] <- "A|G"
  mat_tstv_ancestral_GT3[counter+5,3] <- AG_snps
  mat_tstv_ancestral_GT3[counter+5,4] <- "#cc1236"
  mat_tstv_ancestral_GT3[counter+5,5] <-"Transition"
  
  mat_tstv_ancestral_GT3[counter+6,2] <- "C|T"
  mat_tstv_ancestral_GT3[counter+6,3] <- TC_snps
  mat_tstv_ancestral_GT3[counter+6,4] <- "#f39669"
  mat_tstv_ancestral_GT3[counter+6,5] <-"Transition"
  
  counter <- counter+6
  
}

counter <- 0
cp_names<- 0
for (j in 6:(ncol(denovo_GT3_snps)-3)){
  
  ##colnames 
  mat_names <- colnames(denovo_GT3_snps)[6:29]
  cp_names <- cp_names+1
  #initializing substitution type
 
  subst_df_j <- cbind(denovo_GT3_snps[,1:5],as.numeric(denovo_GT3_snps[,j]))
  
  AT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "A" ),6] >0))
  AG_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "G" ) | (subst_df_j$REF == "G" & subst_df_j$ALT == "A" ),6] >0))  
  AC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "A" ),6] >0))   
  GT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "G" ),6] >0)) 
  GC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "G" ),6] >0))
  TC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "T" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "T" ),6] >0))   
  
  
  mat_tstv_denovo_GT3[(counter+1):(counter+6),1] <- mat_names[cp_names]
  mat_tstv_denovo_GT3[counter+1,2] <- "A|C"
  mat_tstv_denovo_GT3[counter+1,3] <- AC_snps
  mat_tstv_denovo_GT3[counter+1,4] <- "#273572"
  mat_tstv_denovo_GT3[counter+1,5] <-"Transversion"
  
  mat_tstv_denovo_GT3[counter+2,2] <- "A|T"
  mat_tstv_denovo_GT3[counter+2,3] <- AT_snps
  mat_tstv_denovo_GT3[counter+2,4] <- "#6f91cb"
  mat_tstv_denovo_GT3[counter+2,5] <-"Transversion"
  
  mat_tstv_denovo_GT3[counter+3,2] <- "C|G"
  mat_tstv_denovo_GT3[counter+3,3] <- GC_snps
  mat_tstv_denovo_GT3[counter+3,4] <- "#8bd0df"
  mat_tstv_denovo_GT3[counter+3,5] <-"Transversion"
  
  
  mat_tstv_denovo_GT3[counter+4,2] <- "G|T"
  mat_tstv_denovo_GT3[counter+4,3] <- GT_snps
  mat_tstv_denovo_GT3[counter+4,4] <- "#a6d4ba"
  mat_tstv_denovo_GT3[counter+4,5] <-"Transversion"
  
  mat_tstv_denovo_GT3[counter+5,2] <- "A|G"
  mat_tstv_denovo_GT3[counter+5,3] <- AG_snps
  mat_tstv_denovo_GT3[counter+5,4] <- "#cc1236"
  mat_tstv_denovo_GT3[counter+5,5] <-"Transition"
  
  mat_tstv_denovo_GT3[counter+6,2] <- "C|T"
  mat_tstv_denovo_GT3[counter+6,3] <- TC_snps
  mat_tstv_denovo_GT3[counter+6,4] <- "#f39669"
  mat_tstv_denovo_GT3[counter+6,5] <-"Transition"
  
  counter <- counter+6
  
}

counter <- 0
cp_names<- 0
for (j in 6:(ncol(ancestral_GT5_snps)-3)){
  
  ##colnames 
  mat_names <- colnames(ancestral_GT5_snps)[6:35]
  cp_names <- cp_names+1
  #initializing substitution type
  
  subst_df_j <- cbind(ancestral_GT5_snps[,1:5],as.numeric(ancestral_GT5_snps[,j]))
  
  AT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "A" ),6] >0))
  AG_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "G" ) | (subst_df_j$REF == "G" & subst_df_j$ALT == "A" ),6] >0))  
  AC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "A" ),6] >0))   
  GT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "G" ),6] >0)) 
  GC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "G" ),6] >0))
  TC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "T" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "T" ),6] >0))   
  
  
  mat_tstv_ancestral_GT5[(counter+1):(counter+6),1] <- mat_names[cp_names]
  mat_tstv_ancestral_GT5[counter+1,2] <- "A|C"
  mat_tstv_ancestral_GT5[counter+1,3] <- AC_snps
  mat_tstv_ancestral_GT5[counter+1,4] <- "#273572"
  mat_tstv_ancestral_GT5[counter+1,5] <-"Transversion"
  
  mat_tstv_ancestral_GT5[counter+2,2] <- "A|T"
  mat_tstv_ancestral_GT5[counter+2,3] <- AT_snps
  mat_tstv_ancestral_GT5[counter+2,4] <- "#6f91cb"
  mat_tstv_ancestral_GT5[counter+2,5] <-"Transversion"
  
  mat_tstv_ancestral_GT5[counter+3,2] <- "C|G"
  mat_tstv_ancestral_GT5[counter+3,3] <- GC_snps
  mat_tstv_ancestral_GT5[counter+3,4] <- "#8bd0df"
  mat_tstv_ancestral_GT5[counter+3,5] <-"Transversion"
  
  
  mat_tstv_ancestral_GT5[counter+4,2] <- "G|T"
  mat_tstv_ancestral_GT5[counter+4,3] <- GT_snps
  mat_tstv_ancestral_GT5[counter+4,4] <- "#a6d4ba"
  mat_tstv_ancestral_GT5[counter+4,5] <-"Transversion"
  
  mat_tstv_ancestral_GT5[counter+5,2] <- "A|G"
  mat_tstv_ancestral_GT5[counter+5,3] <- AG_snps
  mat_tstv_ancestral_GT5[counter+5,4] <- "#cc1236"
  mat_tstv_ancestral_GT5[counter+5,5] <-"Transition"
  
  mat_tstv_ancestral_GT5[counter+6,2] <- "C|T"
  mat_tstv_ancestral_GT5[counter+6,3] <- TC_snps
  mat_tstv_ancestral_GT5[counter+6,4] <- "#f39669"
  mat_tstv_ancestral_GT5[counter+6,5] <-"Transition"
  
  counter <- counter+6
  
}


counter <- 0
cp_names<- 0
for (j in 6:(ncol(denovo_GT5_snps)-3)){
  
  ##colnames 
  mat_names <- colnames(denovo_GT5_snps)[6:35]
  cp_names <- cp_names+1
  #initializing substitution type
  
  subst_df_j <- cbind(denovo_GT5_snps[,1:5],as.numeric(denovo_GT5_snps[,j]))
  
  AT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "A" ),6] >0))
  AG_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "G" ) | (subst_df_j$REF == "G" & subst_df_j$ALT == "A" ),6] >0))  
  AC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "A" ),6] >0))   
  GT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "G" ),6] >0)) 
  GC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "G" ),6] >0))
  TC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "T" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "T" ),6] >0))   
  
  
  mat_tstv_denovo_GT5[(counter+1):(counter+6),1] <- mat_names[cp_names]
  mat_tstv_denovo_GT5[counter+1,2] <- "A|C"
  mat_tstv_denovo_GT5[counter+1,3] <- AC_snps
  mat_tstv_denovo_GT5[counter+1,4] <- "#273572"
  mat_tstv_denovo_GT5[counter+1,5] <-"Transversion"
  
  mat_tstv_denovo_GT5[counter+2,2] <- "A|T"
  mat_tstv_denovo_GT5[counter+2,3] <- AT_snps
  mat_tstv_denovo_GT5[counter+2,4] <- "#6f91cb"
  mat_tstv_denovo_GT5[counter+2,5] <-"Transversion"
  
  mat_tstv_denovo_GT5[counter+3,2] <- "C|G"
  mat_tstv_denovo_GT5[counter+3,3] <- GC_snps
  mat_tstv_denovo_GT5[counter+3,4] <- "#8bd0df"
  mat_tstv_denovo_GT5[counter+3,5] <-"Transversion"
  
  
  mat_tstv_denovo_GT5[counter+4,2] <- "G|T"
  mat_tstv_denovo_GT5[counter+4,3] <- GT_snps
  mat_tstv_denovo_GT5[counter+4,4] <- "#a6d4ba"
  mat_tstv_denovo_GT5[counter+4,5] <-"Transversion"
  
  mat_tstv_denovo_GT5[counter+5,2] <- "A|G"
  mat_tstv_denovo_GT5[counter+5,3] <- AG_snps
  mat_tstv_denovo_GT5[counter+5,4] <- "#cc1236"
  mat_tstv_denovo_GT5[counter+5,5] <-"Transition"
  
  mat_tstv_denovo_GT5[counter+6,2] <- "C|T"
  mat_tstv_denovo_GT5[counter+6,3] <- TC_snps
  mat_tstv_denovo_GT5[counter+6,4] <- "#f39669"
  mat_tstv_denovo_GT5[counter+6,5] <-"Transition"
  
  counter <- counter+6
  
}

counter <- 0
cp_names<- 0
for (j in 6:(ncol(ancestral_GT9_snps)-3)){
  
  ##colnames 
  mat_names <- colnames(ancestral_GT9_snps)[6:47]
  cp_names <- cp_names+1
  #initializing substitution type
  
  subst_df_j <- cbind(ancestral_GT9_snps[,1:5],as.numeric(ancestral_GT9_snps[,j]))
  
  AT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "A" ),6] >0))
  AG_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "G" ) | (subst_df_j$REF == "G" & subst_df_j$ALT == "A" ),6] >0))  
  AC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "A" ),6] >0))   
  GT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "G" ),6] >0)) 
  GC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "G" ),6] >0))
  TC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "T" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "T" ),6] >0))   
  
  
  mat_tstv_ancestral_GT9[(counter+1):(counter+6),1] <- mat_names[cp_names]
  mat_tstv_ancestral_GT9[counter+1,2] <- "A|C"
  mat_tstv_ancestral_GT9[counter+1,3] <- AC_snps
  mat_tstv_ancestral_GT9[counter+1,4] <- "#273572"
  mat_tstv_ancestral_GT9[counter+1,5] <-"Transversion"
  
  mat_tstv_ancestral_GT9[counter+2,2] <- "A|T"
  mat_tstv_ancestral_GT9[counter+2,3] <- AT_snps
  mat_tstv_ancestral_GT9[counter+2,4] <- "#6f91cb"
  mat_tstv_ancestral_GT9[counter+2,5] <-"Transversion"
  
  mat_tstv_ancestral_GT9[counter+3,2] <- "C|G"
  mat_tstv_ancestral_GT9[counter+3,3] <- GC_snps
  mat_tstv_ancestral_GT9[counter+3,4] <- "#8bd0df"
  mat_tstv_ancestral_GT9[counter+3,5] <-"Transversion"
  
  
  mat_tstv_ancestral_GT9[counter+4,2] <- "G|T"
  mat_tstv_ancestral_GT9[counter+4,3] <- GT_snps
  mat_tstv_ancestral_GT9[counter+4,4] <- "#a6d4ba"
  mat_tstv_ancestral_GT9[counter+4,5] <-"Transversion"
  
  mat_tstv_ancestral_GT9[counter+5,2] <- "A|G"
  mat_tstv_ancestral_GT9[counter+5,3] <- AG_snps
  mat_tstv_ancestral_GT9[counter+5,4] <- "#cc1236"
  mat_tstv_ancestral_GT9[counter+5,5] <-"Transition"
  
  mat_tstv_ancestral_GT9[counter+6,2] <- "C|T"
  mat_tstv_ancestral_GT9[counter+6,3] <- TC_snps
  mat_tstv_ancestral_GT9[counter+6,4] <- "#f39669"
  mat_tstv_ancestral_GT9[counter+6,5] <-"Transition"
  
  counter <- counter+6
  
}

counter <- 0
cp_names<- 0
for (j in 6:(ncol(denovo_GT9_snps)-3)){
  
  ##colnames 
  mat_names <- colnames(denovo_GT9_snps)[6:47]
  cp_names <- cp_names+1
  #initializing substitution type
  
  subst_df_j <- cbind(denovo_GT9_snps[,1:5],as.numeric(denovo_GT9_snps[,j]))
  
  AT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "A" ),6] >0))
  AG_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "G" ) | (subst_df_j$REF == "G" & subst_df_j$ALT == "A" ),6] >0))  
  AC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "A" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "A" ),6] >0))   
  GT_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "T" ) | (subst_df_j$REF == "T" & subst_df_j$ALT == "G" ),6] >0)) 
  GC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "G" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "G" ),6] >0))
  TC_snps <- length(which(subst_df_j[ (subst_df_j$REF == "T" & subst_df_j$ALT == "C" ) | (subst_df_j$REF == "C" & subst_df_j$ALT == "T" ),6] >0))   
  
  
  mat_tstv_denovo_GT9[(counter+1):(counter+6),1] <- mat_names[cp_names]
  mat_tstv_denovo_GT9[counter+1,2] <- "A|C"
  mat_tstv_denovo_GT9[counter+1,3] <- AC_snps
  mat_tstv_denovo_GT9[counter+1,4] <- "#273572"
  mat_tstv_denovo_GT9[counter+1,5] <-"Transversion"
  
  mat_tstv_denovo_GT9[counter+2,2] <- "A|T"
  mat_tstv_denovo_GT9[counter+2,3] <- AT_snps
  mat_tstv_denovo_GT9[counter+2,4] <- "#6f91cb"
  mat_tstv_denovo_GT9[counter+2,5] <-"Transversion"
  
  mat_tstv_denovo_GT9[counter+3,2] <- "C|G"
  mat_tstv_denovo_GT9[counter+3,3] <- GC_snps
  mat_tstv_denovo_GT9[counter+3,4] <- "#8bd0df"
  mat_tstv_denovo_GT9[counter+3,5] <-"Transversion"
  
  
  mat_tstv_denovo_GT9[counter+4,2] <- "G|T"
  mat_tstv_denovo_GT9[counter+4,3] <- GT_snps
  mat_tstv_denovo_GT9[counter+4,4] <- "#a6d4ba"
  mat_tstv_denovo_GT9[counter+4,5] <-"Transversion"
  
  mat_tstv_denovo_GT9[counter+5,2] <- "A|G"
  mat_tstv_denovo_GT9[counter+5,3] <- AG_snps
  mat_tstv_denovo_GT9[counter+5,4] <- "#cc1236"
  mat_tstv_denovo_GT9[counter+5,5] <-"Transition"
  
  mat_tstv_denovo_GT9[counter+6,2] <- "C|T"
  mat_tstv_denovo_GT9[counter+6,3] <- TC_snps
  mat_tstv_denovo_GT9[counter+6,4] <- "#f39669"
  mat_tstv_denovo_GT9[counter+6,5] <-"Transition"
  
  counter <- counter+6
  
}

mat_tstv_ancestral <- as.data.frame(rbind(mat_tstv_ancestral_GT3, mat_tstv_ancestral_GT5, mat_tstv_ancestral_GT9))  
mat_tstv_denovo <- as.data.frame(rbind(mat_tstv_denovo_GT3,mat_tstv_denovo_GT5,mat_tstv_denovo_GT9))

colnames(mat_tstv_ancestral) <- c("Individual",	"Substitution",	"Number",	"Color",	"TSTV")
colnames(mat_tstv_denovo) <- c("Individual",	"Substitution",	"Number",	"Color",	"TSTV")

write.csv2(mat_tstv_ancestral, "./mat_subst_type_ancestral.csv")
write.csv2(mat_tstv_denovo, "./mat_subst_type_denovo.csv")



subst_type_denovo <- mat_tstv_denovo %>% tidyr::separate(Individual, c("Family","Generation","Replicate","Ind"), sep =c("_") ,remove = FALSE) %>%
  dplyr::mutate(Family_generation= glue::glue("{Family}_{Generation}"))  %>%
  dplyr::mutate(Family_generation_subs= glue::glue("{Family}_{Generation}_{Substitution}")) %>%  select(-Replicate, -Ind, -Family, -Generation)

subst_type_ancestral <- mat_tstv_ancestral %>% tidyr::separate(Individual, c("Family","Generation","Replicate","Ind"), sep =c("_") ,remove = FALSE) %>%
  dplyr::mutate(Family_generation= glue::glue("{Family}_{Generation}"))  %>%
  dplyr::mutate(Family_generation_subs= glue::glue("{Family}_{Generation}_{Substitution}")) %>%  select(-Replicate, -Ind, -Family, -Generation)

list_FGS <- c(subst_type_ancestral$Family_generation_subs, subst_type_denovo$Family_generation_subs) %>% sort() %>% unique()

#create a matrix of mean number for each generations and family 
mat_moy_tstv_ancestral <- matrix(data=NA,nrow=length(list_FGS),ncol=5)
mat_moy_tstv_denovo <- matrix(data=NA,nrow=length(list_FGS),ncol=5)
colnames(mat_moy_tstv_ancestral) <- c("Family_generation","Substitution", "Number","Color","TSTV")
colnames(mat_moy_tstv_denovo) <- c("Family_generation","Substitution", "Number","Color","TSTV")

i=0
for (i in 1:length(list_FGS)){
  FGS <- list_FGS[i]
   subsi<- subst_type_ancestral[subst_type_ancestral$Family_generation_subs == FGS, ]
  
   mat_moy_tstv_ancestral[i,1] <- subsi$Family_generation[1]
   mat_moy_tstv_ancestral[i,2] <- subsi$Substitution[1]
   mat_moy_tstv_ancestral[i,4] <- subsi$Color[1]
   mat_moy_tstv_ancestral[i,5] <- subsi$TSTV[1]
   
   mean_tstv <- as.numeric(mean(as.numeric(subsi$Number)))
   mat_moy_tstv_ancestral[i,3] <- as.numeric(mean_tstv)
   
}

i=0
for (i in 1:length(list_FGS)){
  FGS <- list_FGS[i]
  subsi<- subst_type_denovo[subst_type_denovo$Family_generation_subs == FGS, ]
  
  mat_moy_tstv_denovo[i,1] <- subsi$Family_generation[1]
  mat_moy_tstv_denovo[i,2] <- subsi$Substitution[1]
  mat_moy_tstv_denovo[i,4] <- subsi$Color[1]
  mat_moy_tstv_denovo[i,5] <- subsi$TSTV[1]
  
  mean_tstv <- as.numeric(mean(as.numeric(subsi$Number)))
  mat_moy_tstv_denovo[i,3] <- as.numeric(mean_tstv)
  
}

df_moy_tstv_denovo <- as.data.frame(mat_moy_tstv_denovo)
df_moy_tstv_denovo$Number <- as.numeric(df_moy_tstv_denovo$Number)
df_moy_tstv_denovo$Family_generation <- gsub("_G5$", "_G05", df_moy_tstv_denovo$Family_generation)
df_moy_tstv_denovo <- df_moy_tstv_denovo[order(df_moy_tstv_denovo$Family_generation), ]
write.csv2(df_moy_tstv_denovo, "./mean_subst_type_denovo.csv")
df_moy_tstv_ancestral <- as.data.frame(mat_moy_tstv_ancestral)
df_moy_tstv_ancestral$Number <- as.numeric(df_moy_tstv_ancestral$Number)
df_moy_tstv_ancestral$Family_generation <- gsub("_G5$", "_G05", df_moy_tstv_ancestral$Family_generation)
df_moy_tstv_ancestral <- df_moy_tstv_ancestral[order(df_moy_tstv_ancestral$Family_generation), ]
write.csv2(df_moy_tstv_ancestral, "./mean_subst_type_ancestral.csv")

x_labels<-  sub(".*_", "",df_moy_tstv_ancestral$Family_generation)
rects <- data.frame(
  xmin = c(0.5, 4.5, 9.5),  
  xmax = c(4.5, 9.5, 16.5),
  ymin = -Inf,
  ymax = Inf,
  Group = c("GT3", "GT5", "GT9")
)
labels_rects <- data.frame(
  x = c(2.5, 7, 13),                   
  label = c("FA", "C", "NFA") 
)
minmax_labels_ancestral <- df_moy_tstv_ancestral %>%
  group_by(Family_generation) %>%
  summarize(
    x = unique(Family_generation),
    y = sum(Number) + 10,                    
    label = paste0(round(sum(Number[TSTV=="Transition"])), ":", round(sum(Number[TSTV=="Transversion"])))
  )
minmax_labels_denovo <- df_moy_tstv_denovo %>%
  group_by(Family_generation) %>%
  summarize(
    x = unique(Family_generation),
    y = sum(Number) + 2,                    
    label = paste0(round(sum(Number[TSTV=="Transition"])), ":", round(sum(Number[TSTV=="Transversion"])))
  )
p_moy_tstv_ancestral<-ggplot(data = df_moy_tstv_ancestral, aes(x = Family_generation, y = Number)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Group),
            inherit.aes = FALSE, alpha = 0.1) +
  scale_fill_manual(values = c("GT3" = "#688ec8",  
                               "GT5" = "#deb180",  
                               "GT9" = "#8eba6b"))+
  scale_x_discrete(labels = x_labels)+
  new_scale_fill() +
  geom_col(aes(fill = Substitution, group=TSTV), width = 0.7)+
  labs(x = "", y = "Mean number of ancestral substitutions", title = "") +
  scale_fill_manual(values=df_moy_tstv_ancestral$Color) + 
  geom_text(data = labels_rects, aes(x = x, y = 350, label = label), inherit.aes = FALSE, 
            size = 5, fontface = "bold") +
  geom_text(
    data = minmax_labels_ancestral,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold")+
  theme_classic()+
  theme(legend.position = "none")

#------------------------------------------------------------------plot ancestral or de novo substitution datas------------------------------------------------------------------------------
p_moy_tstv_denovo<-ggplot(data = df_moy_tstv_denovo, aes(x = Family_generation, y = Number)) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Group),
            inherit.aes = FALSE, alpha = 0.1) +
  scale_fill_manual(values = c("GT3" = "#688ec8",  
                               "GT5" = "#deb180",  
                               "GT9" = "#8eba6b"))+
  scale_x_discrete(labels = x_labels)+
  labs(x = "", y = "Mean number of de novo substitutions", title = "") +
  new_scale_fill() +
  geom_col(aes(fill = Substitution, group=TSTV), width = 0.7)+
  scale_fill_manual(values=df_moy_tstv_denovo$Color) + 
  geom_text(data = labels_rects, aes(x = x, y = 55, label = label), inherit.aes = FALSE, 
            size = 5, fontface = "bold") +
  geom_text(
    data = minmax_labels_denovo,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold")+
  theme_classic()



i=0
list_FGS <- df_moy_tstv_denovo$Family_generation %>% sort() %>% unique()
agg_df_final <- tibble()
for (i in 1:length(list_FGS)){
   list_FGS_i <- list_FGS[i]
  subsi<- df_moy_tstv_denovo[df_moy_tstv_denovo$Family_generation ==  list_FGS_i,]
  split_subsi <-subsi%>% tidyr::separate(Substitution, c("REF","ALT"), sep =c(">") ,remove = FALSE) 
  
  split_subsi$Pair <- apply(split_subsi[, c("REF", "ALT")], 1, function(x) paste(sort(x), collapse = "|"))
  agg_df <- aggregate(Number ~ Pair, data = split_subsi, sum)
  agg_df$Family_generation <- list_FGS_i
  
  
  agg_df <- split_subsi %>%
    group_by(Family_generation, Pair) %>%
    summarise(
      Total_Number = sum(Number),
      TSTV = names(sort(table(TSTV), decreasing = TRUE))[1], # major type
      .groups = "drop"
    )
  
  unique_pairs <- sort(unique(agg_df$Pair))
  palette <- c("#082684","#CC0033", "#6699FF", "#66FFFF","#FF9966","#99FFCC")
  
  # assign a unique color to each pair
  pair_colors <- setNames(palette, unique_pairs)
  agg_df$Color <- pair_colors[agg_df$Pair]
  agg_df_final <- rbind(agg_df_final, agg_df)
  
  
  
  
}

agg_df_final <- as.data.frame(agg_df_final)
agg_df_final$Total_Number <- as.numeric(agg_df_final$Total_Number)

#------------------------------------------------------------------plot global datas------------------------------------------------------------------------------



p_global_subst <- ggplot(data = agg_df_final, aes(x = Family_generation, y = Total_Number)) +
  geom_col(aes(fill = Pair, group=TSTV), width = 0.7)+
  scale_fill_manual(values=agg_df_final$Color) + theme_classic()



#-----------------------------------------Compute ts:tv ratio ------------------------------------------------


family_generation_names <- mat_moy_tstv_ancestral[,1] %>% sort() %>% unique()

mat_tstv_ratio_ancestral<- matrix(data=NA,nrow = length(family_generation_names), ncol=4)
colnames(mat_tstv_ratio_ancestral) <- c("Fam_Gen","ts","tv","ratio")

for (i in 1:length(family_generation_names)){
  lines.tmp<- mat_moy_tstv_ancestral[grep(family_generation_names[i],mat_moy_tstv_ancestral),]
  ts_number <- sum(as.numeric(lines.tmp[grep("Transition",lines.tmp[,5]),3]))
  tv_number <- sum(as.numeric(lines.tmp[grep("Transversion",lines.tmp[,5]),3]))
  ratio<- ts_number/tv_number
  
  mat_tstv_ratio_ancestral[i,1] <- lines.tmp[1,1]
  mat_tstv_ratio_ancestral[i,2] <- ts_number
  mat_tstv_ratio_ancestral[i,3] <- tv_number
  mat_tstv_ratio_ancestral[i,4] <- ratio
}

mat_tstv_ratio_ancestral[,1] <- gsub("_G5$", "_G05", mat_tstv_ratio_ancestral[,1])
mat_tstv_ratio_ancestral <- mat_tstv_ratio_ancestral[order(mat_tstv_ratio_ancestral[,1]), ]
mat_tstv_ratio_ancestral <- as.data.frame(mat_tstv_ratio_ancestral) %>%
  mutate(
    Group = sub("_.*", "", Fam_Gen),   # garde GT3 / GT5 / GT9
    Fam_Gen = factor(Fam_Gen, levels = Fam_Gen)   # garde ordre d’origine
  )
mat_tstv_ratio_ancestral$ratio <- as.numeric(mat_tstv_ratio_ancestral$ratio)

mat_tstv_ratio_denovo<- matrix(data=NA,nrow = length(family_generation_names), ncol=4)
colnames(mat_tstv_ratio_denovo) <- c("Fam_Gen","ts","tv","ratio")

for (i in 1:length(family_generation_names)){
  lines.tmp<- mat_moy_tstv_denovo[grep(family_generation_names[i],mat_moy_tstv_denovo),]
  ts_number <- sum(as.numeric(lines.tmp[grep("Transition",lines.tmp[,5]),3]))
  tv_number <- sum(as.numeric(lines.tmp[grep("Transversion",lines.tmp[,5]),3]))
  ratio<- ts_number/tv_number
  
  mat_tstv_ratio_denovo[i,1] <- lines.tmp[1,1]
  mat_tstv_ratio_denovo[i,2] <- ts_number
  mat_tstv_ratio_denovo[i,3] <- tv_number
  mat_tstv_ratio_denovo[i,4] <- ratio
}

mat_tstv_ratio_denovo[,1] <- gsub("_G5$", "_G05", mat_tstv_ratio_denovo[,1])
mat_tstv_ratio_denovo <- mat_tstv_ratio_denovo[order(mat_tstv_ratio_denovo[,1]), ]
mat_tstv_ratio_denovo <- as.data.frame(mat_tstv_ratio_denovo) %>%
  mutate(
    Group = sub("_.*", "", Fam_Gen),   # garde GT3 / GT5 / GT9
    Fam_Gen = factor(Fam_Gen, levels = Fam_Gen)   # garde ordre d’origine
  )
mat_tstv_ratio_denovo$ratio <- as.numeric(mat_tstv_ratio_denovo$ratio)


x_labels<-  sub(".*_", "",mat_tstv_ratio_ancestral$Fam_Gen)

#------------------------------------------------------------------plot ts:tv ratio------------------------------------------------------------------------------


p_tstv_ancestral<- ggplot(mat_tstv_ratio_ancestral, aes(x = Fam_Gen, y = ratio, group=Group, color=Group))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "solid", linewidth = 1) +
  geom_point(size = 3) +
  geom_line(linewidth = 2) + 
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Group),
            inherit.aes = FALSE, alpha = 0.1) +
  scale_fill_manual(values = c("GT3" = "#688ec8",  
                                "GT5" = "#deb180",  
                                "GT9" = "#8eba6b"))+
  scale_color_manual(values = c("GT3" = "#8faad9",  
                                "GT5" = "#fbd09c",  
                                "GT9" = "#abd18f")) + 
  scale_y_continuous(limits = c(0, 1.4), breaks = seq(0, 1.5, by = 0.2))+
  scale_x_discrete(labels = x_labels)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "ts:tv ratio", title = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey", linetype = "solid"),  # lignes horizontales
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none")


p_tstv_denovo<- ggplot(mat_tstv_ratio_denovo, aes(x = Fam_Gen, y = ratio, group=Group, color=Group))+
  geom_hline(yintercept = 0.5, color = "black", linetype = "solid", linewidth = 1) +
  geom_point(size = 3) +
  geom_line(linewidth = 2) + 
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Group),
            inherit.aes = FALSE, alpha = 0.1) +
  scale_fill_manual(values = c("GT3" = "#688ec8",  
                               "GT5" = "#deb180",  
                               "GT9" = "#8eba6b"))+
  scale_color_manual(values = c("GT3" = "#8faad9",  
                                "GT5" = "#fbd09c",  
                                "GT9" = "#abd18f")) + 
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(0, 2.4, by = 0.2))+
  scale_x_discrete(labels = x_labels)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "ts:tv ratio", title = "") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey", linetype = "solid"),  # lignes horizontales
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())



#------------------------------------------------------------------plot final figure ------------------------------------------------------------------------------

figure_global <- (p_moy_tstv_ancestral / p_tstv_ancestral) / (p_moy_tstv_denovo / p_tstv_denovo) 
figure_global


#------------------------------------------------------------------compute pi over time for each family (Supplementary figure)------------------------------------------------------------------------------


# Pivot all sample columns into two columns: 'Sample' and 'AltFreq'
# This converts the data from wide format (one column per sample) to long format
subst_long <- subst_table %>%
  pivot_longer(cols = starts_with("GT"), 
               names_to = "Sample", 
               values_to = "AltFreq")

# Extract Family, Passage, and Individual information from the Sample column
subst_long <- subst_long %>%
  mutate(
    Family = gsub("^(GT\\d+)_.*", "\\1", Sample),        # Extracts the family code: GT3, GT5, GT9
    Passage = gsub("^GT\\d+_(G\\d+)_.*", "\\1", Sample), # Extracts the passage: G0, G10, etc.
    Individual = Sample                                   # Keeps the original sample name as individual identifier
  )

# Calculate nucleotide diversity (π) for each site
# Formula: π = 2 * allele_frequency * (1 - allele_frequency)
subst_long <- subst_long %>%
  mutate(pi_site = 2 * AltFreq * (1 - AltFreq))

# Calculate mean π for each individual within each family
pi_per_ind <- subst_long %>%
  group_by(Family, Individual) %>%
  summarize(pi_ind = mean(pi_site, na.rm = TRUE))

# Calculate mean π per passage and individual
pi_time <- subst_long %>%
  group_by(Family, Individual, Passage) %>%
  summarize(pi_ind = mean(pi_site, na.rm = TRUE))

# Calculate mean π and standard deviation per family and passage
pi_summary <- pi_time %>%
  group_by(Family, Passage) %>%
  summarize(
    pi_mean = mean(pi_ind, na.rm = TRUE),  # Mean nucleotide diversity for the family at this passage
    pi_sd = sd(pi_ind, na.rm = TRUE)       # Standard deviation for error ribbon
  )

# Clean up the Passage column and rename families
pi_summary <- pi_summary %>%
  mutate(
    Passage = as.numeric(gsub("G", "", Passage)),        # Remove "G" and convert to numeric for proper ordering
    Family = recode(Family,                              
                    "GT3" = "FA",                        # Rename families for clarity
                    "GT5" = "C",
                    "GT9" = "NFA")
  ) %>%
  arrange(Passage)  # Order rows by passage

# Plot nucleotide diversity (π) per passage for each family
p_pi <- ggplot(pi_summary, aes(x = Passage, y = pi_mean, color = Family, group = Family)) +
  geom_line() +  # Line connecting mean π values per family
  geom_ribbon(aes(ymin = pi_mean - pi_sd, ymax = pi_mean + pi_sd, fill = Family), 
              alpha = 0.2, color = NA) +  # Ribbon shows ±1 SD
  ylab("Intra-individual mean nucleotide diversity (π)") +
  xlab("Passage") +
  # Set x-axis limits from 0 to 30 with breaks every 5 passages, custom labels
  scale_x_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, 5), 
                     labels = c("P0", "P5", "P10", "P15", "P20", "P25", "P30")) +
  theme_minimal() +
  # Manual color mapping for lines
  scale_color_manual(values = c("FA" = "#86a1da", "C" = "#f8c284", "NFA" = "#98d881", name = "Family")) +
  # Manual fill mapping for ribbons
  scale_fill_manual(values = c("FA" = "#86a1da", "C" = "#f8c284", "NFA" = "#98d881", name = "Family"))

