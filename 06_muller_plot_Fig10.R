#!/usr/bin/Rscript
##title: "Muller plot"
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
  c( "tidyverse", "dplyr", "knitr", "VennDiagram") # plot
) )

addline_format <- function(x,...){
  gsub('//s','/n',x)
}

#------------------------------------------------------------------import datas-------------------------------------------------------------------------------

ancestral_vcf<- read.csv2("./global_suspension_individuals_final.csv")
evol_expe_vcf<- read.csv2("./data_frame_allelic_frequencies_indiduals.csv")
attri <- read.csv2("./attributes_muller_plot.csv") ### table to define colors for plotting according to origin of the snps

#------------------------------------------------------------------prepare datas------------------------------------------------------------------------------
#------------------- 1 - create matrix of allelic frequencies for each individuals used to produce initial viral suspensions
mat_allFreq_ancestral <- matrix(data=0,nrow=nrow(ancestral_vcf),ncol=137)
colnames(mat_allFreq_ancestral) <- c("UniqueID", "POS", "REF", "ALT", "TYPE", "ORIGIN", "YEAR", colnames(ancestral_vcf[7:136]))

i=0

for (i in 1:nrow(ancestral_vcf)){
  TYPE_i <- ancestral_vcf[i,] %>% tidyr::separate(INFO, c("avant","après"), sep =c(";TYPE=")) %>%
    dplyr::select(-avant) %>%
    tidyr::separate(après, c("TYPE","après2"), sep =c(";technology")) %>% dplyr::select(TYPE)
  
  
  ####find informations to store in descriptive columns
  REF <- ancestral_vcf$REF[i]
  ALT <- ancestral_vcf$ALT[i]
  POS <- ancestral_vcf$POS[i]
  UNIQUEID <- ancestral_vcf[i,] %>% dplyr::mutate(UNIQUEID= glue::glue("{POS}_{REF}>{ALT}")) %>% dplyr::select(UNIQUEID)
  FORMAT <- ancestral_vcf$FORMAT[i]
  
  mat_allFreq_ancestral[i,1] <- as.character(UNIQUEID)
  mat_allFreq_ancestral[i,2] <- as.numeric(POS)
  mat_allFreq_ancestral[i,3] <- REF
  mat_allFreq_ancestral[i,4] <- ALT
  mat_allFreq_ancestral[i,5] <- as.character(TYPE_i)
  
  
  #####allelic frequencies for each individuals
  for (j in 7:136){
    FORMAT_j_tmp <- ancestral_vcf[i,j] %>% str_split(., ":")
    
    if (FORMAT == "GT:DP:AD:RO:QR:AO:QA:GL" ){
      
      AO_j <- as.numeric(FORMAT_j_tmp[[1]][6])
      AO_j[is.na(AO_j)] <- 0
      DP_j <- as.numeric(FORMAT_j_tmp[[1]][2])
      FREQ <- as.numeric(AO_j/DP_j)
      FREQ[is.na(FREQ)] <- 0
      mat_allFreq_ancestral[i,j+1] <- as.numeric(FREQ)
    }else if (FORMAT == "GT:DP:AD:RO:QR:AO:QA:GL:GQ:MIN_DP" ){
      AO_j <- as.numeric(FORMAT_j_tmp[[1]][6])
      AO_j[is.na(AO_j)] <- 0
      DP_j <- as.numeric(FORMAT_j_tmp[[1]][2])
      FREQ <- as.numeric(AO_j/DP_j)
      FREQ[is.na(FREQ)] <- 0
      mat_allFreq_ancestral[i,j+1] <- as.numeric(FREQ)
    }else if (FORMAT == "GT:GQ:DP:MIN_DP:QR:RO:QA:AO:AD:GL" ){
      AO_j <- as.numeric(FORMAT_j_tmp[[1]][8])
      AO_j[is.na(AO_j)] <- 0
      DP_j <- as.numeric(FORMAT_j_tmp[[1]][3])
      FREQ <- as.numeric(AO_j/DP_j)
      FREQ[is.na(FREQ)] <- 0
      mat_allFreq_ancestral[i,j+1] <- as.numeric(FREQ)
      
    }
  }
    ########add origin informations
    LINE_i <- mat_allFreq_ancestral[i,8:137]
    origin_i <- as.data.frame(LINE_i[LINE_i > 0]) %>% rownames(.) %>% str_split(., "_20") %>% unlist(.) %>% .[grep("[A-Z]", .)] %>% 
      sort() %>% str_replace(., "LT", "MO") %>% unique() %>% knitr::combine_words(., sep = "_", and="")
    
    year_i <- as.data.frame(LINE_i[LINE_i > 0]) %>% rownames(.) %>% str_split(., "_") %>% unlist(.) %>% .[grep("202[0-1]", .)] %>% 
      sort() %>% unique() %>% knitr::combine_words(., sep = "_", and="")
    
    mat_allFreq_ancestral[i,6] <- origin_i
    mat_allFreq_ancestral[i,7] <- year_i
  }

df_allFreq_ancestral <- as.data.frame(mat_allFreq_ancestral)

#------------------------------------------------------------------process datas------------------------------------------------------------------------------
#-------------1 - gneral statistics
nb_snps_MO <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "MO",])
nb_snps_LEU <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "LEU",])
nb_snps_BR <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "BR",])
nb_snps_ARC <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC",])
nb_snps_MO_LEU <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "LEU_MO",])
nb_snps_MO_BR <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "BR_MO",])
nb_snps_MO_ARC <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_MO",])
nb_snps_LEU_BR <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "BR_LEU",])
nb_snps_LEU_ARC <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_LEU",])
nb_snps_BR_ARC <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR",])
nb_snps_BR_ARC_LEU <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR_LEU",])
nb_snps_BR_LEU_MO <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "BR_LEU_MO",])
nb_snps_ARC_BR_MO <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR_MO",])
nb_snps_ARC_LEU_MO <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_LEU_MO",])
nb_snps_ALL <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR_LEU_MO",])

nb_snps_20_21 <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$YEAR == "2020_2021",])
nb_snps_20 <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$YEAR == "2020",])
nb_snps_21 <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$YEAR == "2021",])



#------------- 2 - find correspondence between ancestral population and evolved populations to assign an origin to snps in evol expe
######subsampling of alleles from only one origin
SNPS_ancestral_origin <-df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "MO" | df_allFreq_ancestral$ORIGIN == "LEU"| 
                   df_allFreq_ancestral$ORIGIN == "BR" | df_allFreq_ancestral$ORIGIN == "ARC",]

########find correspondence between ancestral population and evolved populations to assign an origin to snps in evol expe


i=0

for (i in 1:nrow(evol_expe_vcf)){
  LINE_i <- evol_expe_vcf[i,] 
  Unique_ID <- LINE_i$UniqueID
  
  ####find informations to store in descriptive columns
  origin_i <- df_allFreq_ancestral$ORIGIN[df_allFreq_ancestral$UniqueID == Unique_ID ]
  if (length(origin_i) > 0){
  evol_expe_vcf[i,103] <- origin_i
  }
}



#####subsample tables for each families 

GT3_tab <- cbind(evol_expe_vcf[,2:30], evol_expe_vcf[,103])
GT5_tab <- cbind(evol_expe_vcf[,2:6],evol_expe_vcf[,31:60] , evol_expe_vcf[,103])
GT9_tab <- cbind(evol_expe_vcf[,2:6],evol_expe_vcf[,61:103] )


#####average frequencies of SNPs of each generation

mat_meanFreq_GT3 <- matrix(data=0,nrow=nrow(GT3_tab),ncol=10)
colnames(mat_meanFreq_GT3) <- c("UniqueID", "POS", "REF", "ALT", "TYPE","ORIGIN", "G0", "G5", "G10", "G13")

i=0

for (i in 1:nrow(GT3_tab)){
  LINE_i <- GT3_tab[i,] 
  G0 <- mean(as.numeric(LINE_i[,6:11]))
  G10 <- mean(as.numeric(LINE_i[,12:17]))
  G5 <- mean(as.numeric(LINE_i[,24:29]))
  G13 <- mean(as.numeric(LINE_i[,18:23]))
  origin_i <- LINE_i$`evol_expe_vcf[, 103]`
  
  mat_meanFreq_GT3[i,1] <- LINE_i$UniqueID
  mat_meanFreq_GT3[i,2] <- LINE_i$POS
  mat_meanFreq_GT3[i,3] <- LINE_i$REF
  mat_meanFreq_GT3[i,4] <- LINE_i$ALT
  mat_meanFreq_GT3[i,5] <- LINE_i$TYPE
  mat_meanFreq_GT3[i,6] <- origin_i
  mat_meanFreq_GT3[i,7] <- G0
  mat_meanFreq_GT3[i,8] <- G5
  mat_meanFreq_GT3[i,9] <- G10
  mat_meanFreq_GT3[i,10] <- G13
}
mat_meanFreq_GT3[is.na(mat_meanFreq_GT3)] <- "?"
mat_meanFreq_GT3 <- as.data.frame(mat_meanFreq_GT3)


mat_meanFreq_GT5 <- matrix(data=0,nrow=nrow(GT5_tab),ncol=11)
colnames(mat_meanFreq_GT5) <- c("UniqueID", "POS", "REF", "ALT", "TYPE","ORIGIN", "G0", "G5", "G10", "G15", "G18")

i=0

for (i in 1:nrow(GT5_tab)){
  LINE_i <- GT5_tab[i,] 
  G0 <- mean(as.numeric(LINE_i[,6:11]))
  G10 <- mean(as.numeric(LINE_i[,12:17]))
  G5 <- mean(as.numeric(LINE_i[,30:35]))
  G15 <- mean(as.numeric(LINE_i[,18:23]))
  G18 <- mean(as.numeric(LINE_i[,24:29]))
  origin_i <- LINE_i$`evol_expe_vcf[, 103]`
  
  mat_meanFreq_GT5[i,1] <- LINE_i$UniqueID
  mat_meanFreq_GT5[i,2] <- LINE_i$POS
  mat_meanFreq_GT5[i,3] <- LINE_i$REF
  mat_meanFreq_GT5[i,4] <- LINE_i$ALT
  mat_meanFreq_GT5[i,5] <- LINE_i$TYPE
  mat_meanFreq_GT5[i,6] <- origin_i
  mat_meanFreq_GT5[i,7] <- G0
  mat_meanFreq_GT5[i,8] <- G5
  mat_meanFreq_GT5[i,9] <- G10
  mat_meanFreq_GT5[i,10] <- G15
  mat_meanFreq_GT5[i,11] <- G18
}
mat_meanFreq_GT5[is.na(mat_meanFreq_GT5)] <- "?"
mat_meanFreq_GT5 <- as.data.frame(mat_meanFreq_GT5)

mat_meanFreq_GT9 <- matrix(data=0,nrow=nrow(GT9_tab),ncol=13)
colnames(mat_meanFreq_GT9) <- c("UniqueID", "POS", "REF", "ALT", "TYPE","ORIGIN", "G0", "G5", "G10", "G15", "G20", "G24", "G28")

i=0

for (i in 1:nrow(GT9_tab)){
  LINE_i <- GT9_tab[i,] 
  G0 <- mean(as.numeric(LINE_i[,6:11]))
  G10 <- mean(as.numeric(LINE_i[,12:17]))
  G5 <- mean(as.numeric(LINE_i[,42:47]))
  G15 <- mean(as.numeric(LINE_i[,18:23]))
  G20 <- mean(as.numeric(LINE_i[,24:29]))
  G24 <- mean(as.numeric(LINE_i[,30:35]))
  G28 <- mean(as.numeric(LINE_i[,36:41]))
  origin_i <- LINE_i$V103
  
  mat_meanFreq_GT9[i,1] <- LINE_i$UniqueID
  mat_meanFreq_GT9[i,2] <- LINE_i$POS
  mat_meanFreq_GT9[i,3] <- LINE_i$REF
  mat_meanFreq_GT9[i,4] <- LINE_i$ALT
  mat_meanFreq_GT9[i,5] <- LINE_i$TYPE
  mat_meanFreq_GT9[i,6] <- origin_i
  mat_meanFreq_GT9[i,7] <- G0
  mat_meanFreq_GT9[i,8] <- G5
  mat_meanFreq_GT9[i,9] <- G10
  mat_meanFreq_GT9[i,10] <- G15
  mat_meanFreq_GT9[i,11] <- G20
  mat_meanFreq_GT9[i,12] <- G24
  mat_meanFreq_GT9[i,13] <- G28
}
mat_meanFreq_GT9[is.na(mat_meanFreq_GT9)] <- "?"
mat_meanFreq_GT9 <- as.data.frame(mat_meanFreq_GT9)


#------------- 3 - create dataframe of abundances for each origin  


mat_muller_GT3 <- matrix(data=0,nrow=5,ncol=4)
colnames(mat_muller_GT3) <- c( "0", "5", "10", "13")

row.names(mat_muller_GT3) <- c("MO", "BR", "ARC", "LEU","?")

i=0
for (i in 1:nrow(mat_muller_GT3)){

origin <- mat_meanFreq_GT3[mat_meanFreq_GT3$ORIGIN == rownames(mat_muller_GT3)[i], ] 

mat_muller_GT3[i,1] <- length(origin$G0[origin$G0 >0])
mat_muller_GT3[i,2] <- length(origin$G5[origin$G5 >0])
mat_muller_GT3[i,3] <- length(origin$G10[origin$G10 >0])
mat_muller_GT3[i,4] <- length(origin$G13[origin$G13 >0])
}
origin <- mat_meanFreq_GT3[mat_meanFreq_GT3$ORIGIN == rownames(mat_muller_GT3)[5], ] 
mat_muller_GT3[5,1] <- length(origin$G0[origin$G0 >0])
mat_muller_GT3[5,2] <- length(origin$G5[origin$G5 >0])
mat_muller_GT3[5,3] <- length(origin$G10[origin$G10 >0])
mat_muller_GT3[5,4] <- length(origin$G13[origin$G13 >0])



layout(matrix(c(1,2), nrow = 1), widths = c(0.8, 0.2))
par(mar=c(3, 3, 1, 0))
m.plot_GT3 <- Muller.plot(attributes = attri, population.data = mat_muller_GT3,
                      data.method = "table", time.interval.method = "linear",mgp=c(2,0.5,0))




mat_muller_GT5 <- matrix(data=0,nrow=5,ncol=6)
colnames(mat_muller_GT5) <- c( "0", "5", "10", "15","18")
row.names(mat_muller_GT5) <- c("MO", "BR", "ARC", "LEU", "?")

i=0
for (i in 1:nrow(mat_muller_GT5)){
  
  origin <- mat_meanFreq_GT5[mat_meanFreq_GT5$ORIGIN == rownames(mat_muller_GT5)[i], ] 
  #pattern <- rownames(mat_muller_GT5)[i]
  #origin <- mat_meanFreq_GT5[grep(pattern = pattern, mat_meanFreq_GT5$ORIGIN),]
  
  
  mat_muller_GT5[i,1] <- length(origin$G0[origin$G0 >0])
  mat_muller_GT5[i,2] <- length(origin$G5[origin$G5 >0])
  mat_muller_GT5[i,3] <- length(origin$G10[origin$G10 >0])
  mat_muller_GT5[i,4] <- length(origin$G15[origin$G15 >0])
  mat_muller_GT5[i,5] <- length(origin$G18[origin$G18 >0])
}
origin <- mat_meanFreq_GT5[mat_meanFreq_GT5$ORIGIN == rownames(mat_muller_GT5)[5], ] 
mat_muller_GT5[5,1] <- length(origin$G0[origin$G0 >0])
mat_muller_GT5[5,2] <- length(origin$G5[origin$G5 >0])
mat_muller_GT5[5,3] <- length(origin$G10[origin$G10 >0])
mat_muller_GT5[5,4] <- length(origin$G15[origin$G15 >0])
mat_muller_GT5[5,5] <- length(origin$G18[origin$G18 >0])



layout(matrix(c(1,2), nrow = 1), widths = c(0.8, 0.2))
par(mar=c(3, 3, 1, 0))
m.plot_GT5 <- Muller.plot(attributes = attri, population.data = mat_muller_GT5,
                          data.method = "table", time.interval.method = "linear",mgp=c(2,0.5,0))





mat_muller_GT9 <- matrix(data=0,nrow=5,ncol=7)
colnames(mat_muller_GT9) <- c( "0", "5", "10", "15","20", "24", "28")
row.names(mat_muller_GT9) <- c("MO", "BR", "ARC", "LEU", "?")

nb_snps_MO_LEU <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "LEU_MO",])
nb_snps_MO_BR <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "BR_MO",])
nb_snps_MO_ARC <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_MO",])
nb_snps_LEU_BR <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "BR_LEU",])
nb_snps_LEU_ARC <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_LEU",])
nb_snps_BR_ARC <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR",])
nb_snps_BR_ARC_LEU <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR_LEU",])
nb_snps_BR_LEU_MO <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "BR_LEU_MO",])
nb_snps_ARC_BR_MO <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR_MO",])
nb_snps_ARC_LEU_MO <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_LEU_MO",])
nb_snps_ALL <- nrow(df_allFreq_ancestral[df_allFreq_ancestral$ORIGIN == "ARC_BR_LEU_MO",])

i=0
for (i in 1:nrow(mat_muller_GT9)){
  
  origin <- mat_meanFreq_GT9[mat_meanFreq_GT9$ORIGIN == rownames(mat_muller_GT9)[i], ] 
  #pattern <- rownames(mat_muller_GT9)[i]
  #origin <- mat_meanFreq_GT9[grep(pattern = pattern, mat_meanFreq_GT9$ORIGIN),]
  
  
  mat_muller_GT9[i,1] <- length(origin$G0[origin$G0 >0])
  mat_muller_GT9[i,2] <- length(origin$G5[origin$G5 >0])
  mat_muller_GT9[i,3] <- length(origin$G10[origin$G10 >0])
  mat_muller_GT9[i,4] <- length(origin$G15[origin$G15 >0])
  mat_muller_GT9[i,5] <- length(origin$G20[origin$G20 >0])
  mat_muller_GT9[i,6] <- length(origin$G24[origin$G24 >0])
  mat_muller_GT9[i,7] <- length(origin$G28[origin$G28 >0])
}
origin <- mat_meanFreq_GT9[mat_meanFreq_GT9$ORIGIN == rownames(mat_muller_GT9)[5], ] 
mat_muller_GT9[5,1] <- length(origin$G0[origin$G0 >0])
mat_muller_GT9[5,2] <- length(origin$G5[origin$G5 >0])
mat_muller_GT9[5,3] <- length(origin$G10[origin$G10 >0])
mat_muller_GT9[5,4] <- length(origin$G15[origin$G15 >0])
mat_muller_GT9[5,5] <- length(origin$G20[origin$G20 >0])
mat_muller_GT9[5,6] <- length(origin$G24[origin$G24 >0])
mat_muller_GT9[5,7] <- length(origin$G28[origin$G28 >0])


layout(matrix(c(1,2), nrow = 1), widths = c(0.8, 0.2))
par(mar=c(3, 3, 1, 0))
m.plot_GT9 <- Muller.plot(attributes = attri, population.data = mat_muller_GT9,
                          data.method = "table", time.interval.method = "linear",mgp=c(2,0.5,0))



#------------------------------------------------------------------plot datas as Muller plot------------------------------------------------------------------------------

plot_grid(m.plot_GT3, m.plot_GT5, m.plot_GT9, labels = c('GT3', 'GT5', 'GT9'), label_size = 12)
