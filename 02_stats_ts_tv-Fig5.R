#!/usr/bin/Rscript
##title: "Plot transition and transversions experimental evolution -  Figure 5"
##author: "Pelletier Camille"
#date: 20/09/2022

rm(list=ls())
'%!in%' <- function(x,y)!('%in%'(x,y))
options( "digits"=3)

###################################### Function to install / load package if it's not installed ###############################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(unique(
  c( "tidyverse", "dplyr") # plot
) )

addline_format <- function(x,...){
  gsub('//s','/n',x)
}

#------------------------------------------------------------------import datas-------------------------------------------------------------------------------
setwd(/path/to/folder/output/04_vcffiles_bcftools)
##using the output of all files from bcftools stats, add a color per substitution type and a Transition (Ts) or Transversion (Tv) category for each SNP
subst_type<- read.csv2("./all_subst_type_filtered.csv")


#------------------------------------------------------------------prepare datas------------------------------------------------------------------------------

####add a column with details about family, generation and substitution type
subst_type <- subst_type %>% tidyr::separate(Individual, c("Family","Generation","Replicate","Ind"), sep =c("_") ,remove = FALSE) %>%
  dplyr::mutate(Family_generation= glue::glue("{Family}_{Generation}"))  %>%
  dplyr::mutate(Family_generation_subs= glue::glue("{Family}_{Generation}_{Substitution}")) %>%  select(-Replicate, -Ind, -Family, -Generation)

###Create a list of all possible combination of family, generation of infection and substitution type
list_FGS <- subst_type$Family_generation_subs %>% sort() %>% unique()

###create a matrix of mean number of Ts and Tv for each generations and family 
mat_moy_tstv <- matrix(data=NA,nrow=length(list_FGS),ncol=5)
colnames(mat_moy_tstv) <- c("Family_generation","Substitution", "Number","Color","TSTV")


i=0
for (i in 1:length(list_FGS)){
  FGS <- list_FGS[i]
   subsi<- subst_type[subst_type$Family_generation_subs == FGS, ]
  
   mat_moy_tstv[i,1] <- subsi$Family_generation[1]
   mat_moy_tstv[i,2] <- subsi$Substitution[1]
   mat_moy_tstv[i,4] <- subsi$Color[1]
   mat_moy_tstv[i,5] <- subsi$TSTV[1]
   
   mean_tstv <- as.numeric(mean(subsi$Number))
   mat_moy_tstv[i,3] <- as.numeric(mean_tstv)
}

df_moy_tstv <- as.data.frame(mat_moy_tstv)
df_moy_tstv$Number <- as.numeric(df_moy_tstv$Number)
write.csv2(df_moy_tstv, "./mean_subst_type.csv")


#------------------------------------------------------------------plot datas------------------------------------------------------------------------------

##plot all substituion type for all family and generation of infection
x11()
ggplot(data = df_moy_tstv, aes(x = Family_generation, y = Number)) +
  geom_col(aes(fill = Substitution, group=TSTV), width = 0.7)+
  scale_fill_manual(values=subst_type$Color) + theme_classic()



