#!/usr/bin/Rscript
##title: "Genomic selection detection, CMH - Figure6, 7 and 8"
##author: "Pelletier Camille"
#date: 20/09/2022

rm(list=ls())
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
  c( "tidyverse", "dplyr", "factoextra", "poolSeq", "ACER", "NbClust") # plot
) )

addline_format <- function(x,...){
  gsub('//s','/n',x)
}

#------------------------------------------------------------------import datas-------------------------------------------------------------------------------
setwd(/path/to/folder/output/04_vcffiles_bcftools)
vcf_file <- read.csv2("./merge_vcf_files_evol_expe_decompose_blocksub_sort_uniq_modif.csv")



#------------------------------------------------------------------prepare datas------------------------------------------------------------------------------
##create a matrix of allelic frequencies and allelic coverage for each polymorphims detected oyster population and infection generations
mat_allFreq <- matrix(data=0,nrow=nrow(vcf_file),ncol=101)
colnames(mat_allFreq) <- c("UniqueID", "POS", "REF", "ALT", "TYPE", colnames(vcf_file[10:105]))

mat_allCov <- matrix(data=0,nrow=nrow(vcf_file),ncol=101)
colnames(mat_allCov) <- c("UniqueID", "POS", "REF", "ALT", "TYPE", colnames(vcf_file[10:105]))

i=0
for (i in 1:nrow(vcf_file)){
  TYPE_i <- vcf_file[i,] %>% tidyr::separate(INFO, c("avant","après"), sep =c(";TYPE=")) %>%
    dplyr::select(-avant) %>%
    tidyr::separate(après, c("TYPE","après2"), sep =c(";technology")) %>% dplyr::select(TYPE)
  
	## first step add information about the polymorphism within matrix
  REF <- vcf_file$REF[i] ## reference allele
  ALT <- vcf_file$ALT[i] ## alternative allele
  POS <- vcf_file$POS[i] ## position of the polymrophism
  UNIQUEID <- vcf_file[i,] %>% dplyr::mutate(UNIQUEID= glue::glue("{POS}_{REF}>{ALT}")) %>% dplyr::select(UNIQUEID) ## create a unique ID using position, ref and alt alleles
  FORMAT <- vcf_file$FORMAT[i] ## format of the polymorphism
  
  mat_allFreq[i,1] <- as.character(UNIQUEID)
  mat_allFreq[i,2] <- as.numeric(POS)
  mat_allFreq[i,3] <- REF
  mat_allFreq[i,4] <- ALT
  mat_allFreq[i,5] <- as.character(TYPE_i)
  
  mat_allCov[i,1] <- as.character(UNIQUEID)
  mat_allCov[i,2] <- as.numeric(POS)
  mat_allCov[i,3] <- REF
  mat_allCov[i,4] <- ALT
  mat_allCov[i,5] <- as.character(TYPE_i)
  
  	## second step add allelic frequencies and coverage of the polymorphism within matrix
  for (j in 10:105){
    FORMAT_j_tmp <- vcf_file[i,j] %>% str_split(., ":")
    
    if (FORMAT == "GT:DP:AD:RO:QR:AO:QA:GL" ){ 
     
    AO_j <- as.numeric(FORMAT_j_tmp[[1]][6]) ## AO = Count of full observations of this alternate haplotype
    DP_j <- as.numeric(FORMAT_j_tmp[[1]][2]) ## DP = Total read depth at the locus
    FREQ <- as.numeric(AO_j/DP_j)
        mat_allFreq[i,j-4] <- as.numeric(FREQ)
        mat_allCov[i,j-4] <- as.numeric(AO_j)
    }else if (FORMAT == "GT:DP:AD:RO:QR:AO:QA:GL:GQ:MIN_DP" ){
      AO_j <- as.numeric(FORMAT_j_tmp[[1]][6])
      DP_j <- as.numeric(FORMAT_j_tmp[[1]][2])
      FREQ <- as.numeric(AO_j/DP_j)
      mat_allFreq[i,j-4] <- as.numeric(FREQ)
      mat_allCov[i,j-4] <- as.numeric(AO_j)
    }else if (FORMAT == "GT:GQ:DP:MIN_DP:QR:RO:QA:AO:AD:GL" ){
      AO_j <- as.numeric(FORMAT_j_tmp[[1]][8])
      DP_j <- as.numeric(FORMAT_j_tmp[[1]][3])
      FREQ <- as.numeric(AO_j/DP_j)
      mat_allFreq[i,j-4] <- as.numeric(FREQ)
      mat_allCov[i,j-4] <- as.numeric(AO_j)
      
    }
  
  
}}

mat_allFreq[is.na(mat_allFreq)] <- 0
mat_allCov[is.na(mat_allCov)] <- 0

df_allFreq <- as.data.frame(mat_allFreq)
df_allFreq<- df_allFreq %>% mutate_at(6:101, as.numeric)
df_allFreq <- as.data.frame(mat_allFreq)
df_allFreq<- df_allFreq %>% mutate_at(6:101, as.numeric)

df_allCov <- as.data.frame(mat_allCov)
df_allCov<- df_allCov %>% mutate_at(6:101, as.numeric)
df_ind_Cov <- df_allCov[6:101]
df_allCov<- df_allCov %>% mutate_at(6:101, as.numeric)

## write dataframes to files
write.csv2(df_allFreq, "./output_R/data_frame_allelic_frequencies_individuals.csv")
write.csv2(df_allCov, "./output_R/data_frame_allelic_coverage_indiduals.csv")



#------------------------------------------------------------------process datas------------------------------------------------------------------------------
#-------------1 - perform an ACP analysis for all samplesusing similarities between observed SNPs

df_allFreq <- df_allFreq[df_allFreq$TYPE == "snp",]
dat = na.omit(2*asin(sqrt(df_allFreq[6:101]))) 
dat <- dat[rowSums(dat) > 0,]


####add to the matrix columns to describe shared SNPs,etc

i=0
for (i in 1:nrow(dat)){
  GT3_i <- as.numeric(dat[i,1:24]) 
  GT5_i <- as.numeric(dat[i,25:54]) 
  GT9_i <- as.numeric(dat[i,55:96]) 
  
  sum_GT3_i <- as.numeric(sum(GT3_i >0))
  sum_GT5_i <- as.numeric(sum(GT5_i >0))
  sum_GT9_i <- as.numeric(sum(GT9_i >0))
  
  #####add information of shared alleles between oyster populations
  if (sum_GT3_i > 0 & sum_GT5_i > 0 & sum_GT9_i >0){
    dat[i,97] <- "Shared"
  }else if (sum_GT3_i > 0 & sum_GT5_i  ==  0 & sum_GT9_i  == 0 ){ ## SNP only in GT3
    dat[i,97] <- "GT3" 
  }else if (sum_GT3_i  ==  0 & sum_GT5_i > 0 & sum_GT9_i  == 0 ){ ## SNP only in GT5
    dat[i,97] <- "GT5" 
  }else if (sum_GT3_i  ==  0 & sum_GT5_i  ==  0 & sum_GT9_i >0 ){ ## SNP only in GT9
    dat[i,97] <- "GT9" 
      }else if (sum_GT3_i > 0 & sum_GT5_i > 0 & sum_GT9_i  == 0 ){ ## SNP only in GT3  and GT5
    dat[i,97] <- "GT3_GT5"
    }else if (sum_GT3_i > 0 & sum_GT5_i  ==  0 & sum_GT9_i >0 ){## SNP only in GT3  and GT9
    dat[i,97] <- "GT3_GT9"
  }else if (sum_GT3_i == 0 & sum_GT5_i > 0 & sum_GT9_i >0 ){ ## SNP only in GT5  and GT9
    dat[i,97] <- "GT5_GT9"
  }
  
  
  GEN0_i <- as.numeric(cbind(dat[i,1:6], dat[i,25:30],dat[i,55:60] ))
  GEN5_i <- as.numeric(cbind(dat[i,19:24],dat[i,49:54],dat[i,91:96] ) )
  GEN10_i <- as.numeric(cbind(dat[i,7:12],dat[i,31:36],dat[i,61:66]) )
  GEN13_i <- as.numeric(cbind(dat[i,13:18]))
  GEN15_i <- as.numeric(cbind(dat[i,37:42],dat[i,67:72]))
  GEN18_i <- as.numeric(cbind(dat[i,43:48]))
  GEN20_i <- as.numeric(cbind(dat[i,73:78]))
  GEN24_i <- as.numeric(cbind(dat[i,79:84]))
  GEN28_i <- as.numeric(cbind(dat[i,85:90]))
  
  sum_GEN0_i <- sum(GEN0_i >0)
  sum_GEN5_i <- sum(GEN5_i >0)
  sum_GEN10_i <- sum(GEN10_i >0)
  sum_GEN13_i <- sum(GEN13_i >0)
  sum_GEN15_i <- sum(GEN15_i >0)
  sum_GEN18_i <- sum(GEN18_i >0)
  sum_GEN20_i <- sum(GEN20_i >0)
  sum_GEN24_i <- sum(GEN24_i >0)
  sum_GEN28_i <- sum(GEN28_i >0)
  
	#####add information of shared alleles between infection generations
  if (sum_GEN0_i > 0 & sum_GEN5_i > 0 & sum_GEN10_i >0 & sum_GEN13_i >0 & sum_GEN15_i >0 & sum_GEN18_i>0 & sum_GEN20_i>0 & sum_GEN24_i>0 & sum_GEN28_i>0){
    dat[i,98] <- "Shared"
  }else if (sum_GEN0_i > 0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
    dat[i,98] <- "G0" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i > 0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
    dat[i,98] <- "G5" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i >0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
    dat[i,98] <- "G10" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i >0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
    dat[i,98] <- "G13" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i >0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
    dat[i,98] <- "G15" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i>0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
    dat[i,98] <- "G18" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i>0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
    dat[i,98] <- "G20" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i>0 & sum_GEN28_i == 0 ){
    dat[i,98] <-"G24" 
  }else if (sum_GEN0_i  ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i>0 ){
    dat[i,98] <- "G28" 
  } else {
    dat[i,98] <- "shared_not_all" 
  }
#   }  else if (sum_GEN0_i > 0 & sum_GEN5_i  >  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G0_G5" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i== 0 & sum_GEN10_i  > 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G0_G10" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i  ==  0 & sum_GEN10_i ==0 & sum_GEN13_i  > 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G0_G13" 
#   }else if (sum_GEN0_i>0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i ==0 & sum_GEN15_i  > 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G0_G15" 
#   }else if (sum_GEN0_i>0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i==0 & sum_GEN18_i>0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G0_G18" 
#   }else if (sum_GEN0_i>0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i>0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G0_G20" 
#   }else if (sum_GEN0_i>0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i>0 & sum_GEN28_i == 0 ){
#     dat[i,98] <-"G0_G24" 
#   }else if (sum_GEN0_i>0 & sum_GEN5_i  ==  0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i>0 ){
#     dat[i,98] <- "G0_G28" 
#   } else if (sum_GEN0_i == 0 & sum_GEN5_i  >  0 & sum_GEN10_i>0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G5_G10" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i >0 & sum_GEN10_i==0 & sum_GEN13_i  > 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G5_G13" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i >0 & sum_GEN10_i  == 0 & sum_GEN13_i ==0 & sum_GEN15_i  > 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G5_G15" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i >0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i==0 & sum_GEN18_i>0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G5_G18" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i >0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i>0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#     dat[i,98] <- "G5_G20" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i >0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i>0 & sum_GEN28_i == 0 ){
#     dat[i,98] <-"G5_G24" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i >0 & sum_GEN10_i  == 0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i>0 ){
#     dat[i,98] <- "G5_G28" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i>0 & sum_GEN13_i  > 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#   dat[i,98] <- "G10_G13" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i >0 & sum_GEN13_i ==0 & sum_GEN15_i  > 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#   dat[i,98] <- "G10_G15" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i >0 & sum_GEN13_i  == 0 & sum_GEN15_i==0 & sum_GEN18_i>0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#   dat[i,98] <- "G10_G18" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i >0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i>0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#   dat[i,98] <- "G10_G20" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i >0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i>0 & sum_GEN28_i == 0 ){
#   dat[i,98] <-"G10_G24" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i >0 & sum_GEN13_i  == 0 & sum_GEN15_i  == 0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i>0 ){
#   dat[i,98] <- "G10_G28" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i ==0 & sum_GEN13_i  == 0 & sum_GEN15_i==0 & sum_GEN18_i>0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#   dat[i,98] <- "G15_G18" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i ==0 & sum_GEN13_i  == 0 & sum_GEN15_i >0 & sum_GEN18_i == 0 & sum_GEN20_i>0 & sum_GEN24_i == 0 & sum_GEN28_i == 0 ){
#   dat[i,98] <- "G15_G20" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i ==0 & sum_GEN13_i  == 0 & sum_GEN15_i >0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i>0 & sum_GEN28_i == 0 ){
#   dat[i,98] <-"G15_G24" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i ==0 & sum_GEN13_i  == 0 & sum_GEN15_i >0 & sum_GEN18_i == 0 & sum_GEN20_i == 0 & sum_GEN24_i == 0 & sum_GEN28_i>0 ){
#   dat[i,98] <- "G15_G28" 
#   }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i ==0 & sum_GEN13_i  == 0 & sum_GEN15_i ==0 & sum_GEN18_i == 0 & sum_GEN20_i>0 & sum_GEN24_i>0 & sum_GEN28_i == 0 ){
#   dat[i,98] <-"G20_G24" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i ==0 & sum_GEN13_i  == 0 & sum_GEN15_i==0 & sum_GEN18_i == 0 & sum_GEN20_i>0 & sum_GEN24_i == 0 & sum_GEN28_i>0 ){
#   dat[i,98] <- "G20_G28" 
# }else if (sum_GEN0_i ==  0 & sum_GEN5_i ==0 & sum_GEN10_i ==0 & sum_GEN13_i  == 0 & sum_GEN15_i==0 & sum_GEN18_i == 0 & sum_GEN20_i==0 & sum_GEN24_i>0 & sum_GEN28_i>0 ){
#   dat[i,98] <- "G24_G28" }
  
  #####add information of shared alleles between replicates
  R1_i <- as.numeric(cbind(dat[i,1:3], dat[i,7:9],dat[i,13:15],dat[i,19:21],dat[i,25:27],dat[i,31:33],dat[i,37:39],dat[i,43:45],dat[i,49:51],
                           dat[i,55:57],dat[i,61:63],dat[i,67:69],dat[i,73:75],dat[i,79:81],dat[i,85:87],dat[i,91:93] ))
  R2_i <- as.numeric(cbind(dat[i,4:6], dat[i,10:12],dat[i,16:18],dat[i,22:24],dat[i,28:30],dat[i,34:36],dat[i,40:42],dat[i,46:48],dat[i,52:54],
                           dat[i,58:60],dat[i,64:66],dat[i,70:72],dat[i,76:78],dat[i,82:84],dat[i,88:90],dat[i,94:96] ) )
  
  sum_R1_i <- sum(R1_i >0)
  sum_R2_i <- sum(R2_i >0)
  
  if (sum_R1_i > 0 & sum_R2_i > 0){
    dat[i,99]<- "Shared"
  }else if (sum_R1_i > 0 & sum_R2_i==0 ){
    dat[i,99]<- "R1" 
  }else if (sum_R1_i==0 & sum_R2_i > 0 ){
    dat[i,99]<- "R2" 
 
  }
  
  
  }

sampleIDs<-colnames(dat[1:96])
splits <- unlist(strsplit(sampleIDs, "[_]"))
fam_oyster <- factor(splits[grep("GT",splits)])
splits <- splits[-grep("T",splits)]
generation <- as.factor(splits[grep("G",splits)])
generation <- as.numeric(splits[grep("G",splits)] %>% strsplit(., "G") %>% unlist(.) %>% .[grep("[0-9]",.)])
 
###process PCA
rep<- factor(splits[grep("R",splits)])
pcadata<-as.data.frame(t(dat[1:96]))
row.names(pcadata)<-sampleIDs

pcadata_explain<-as.data.frame(t(dat[97:99]))
row.names(pcadata_explain)<-c("Family", "Generation", "Replicate")



pca.res <- prcomp(pcadata, retx=TRUE, center = TRUE, scale. =TRUE)
pca.var<-get_eigenvalue(pca.res)

data_pca <- as.data.frame(pca.res$x)
NbClust(data_pca[,1:2], method = 'kmeans', index = 'all')$Best.nc
NbClust(data_pca[,1:2], method = 'complete', index = 'all')$Best.partition
clust <- kmeans(data_pca[,1:2], 5)$cluster %>%
  as.factor()
data_pca <- data_pca %>%
  mutate(groups = clust)

#------------------------------------------------------------------plot datas as PCA------------------------------------------------------------------------------

###global PCA
x11()
p1_gen_fam<-ggplot(as.data.frame(pca.res$x),aes(x=PC1,y=PC2, col=fam_oyster, shape=rep, alpha=generation))+
  geom_point(size=4)+ 
  scale_shape_manual(values = c(15,17,18))+
  scale_color_manual(values=c( "#000099", "#FF6600", "#009900"))+
  xlab(paste("PC1(",round(pca.var$variance.percent[1], 2),"%)",sep=""))+
  ylab(paste("PC2(",round(pca.var$variance.percent[2], 2),"%)",sep=""))+ theme_classic()

x11()
p1_gen_rep_clust<-ggplot(as.data.frame(pca.res$x),aes(x=PC1,y=PC2, col=clust, shape=rep, alpha=generation))+
  geom_point(size=4)+ 
  scale_shape_manual(values = c(15,17,18))+
  xlab(paste("PC1(",round(pca.var$variance.percent[1], 2),"%)",sep=""))+
  ylab(paste("PC2(",round(pca.var$variance.percent[2], 2),"%)",sep=""))+ theme_classic()

# 
# p2_gen_fam<-ggplot(as.data.frame(pca.res$x),aes(x=PC3,y=PC4, col=fam_oyster, shape=generation))+
#   geom_point(size=2)+
#   scale_shape_manual(values = c(0,15,17,18, 1,2,5,6,19))+
#   xlab(paste("PC3(",round(pca.var$variance.percent[3], 2),"%)",sep=""))+
#   ylab(paste("PC4(",round(pca.var$variance.percent[4], 2),"%)",sep=""))
# 
# p1_rep_fam<-ggplot(as.data.frame(pca.res$x),aes(x=PC1,y=PC2, col=fam_oyster, shape=rep))+
#   geom_point(size=2)+ 
#   scale_shape_manual(values = c(15,17))+
#   xlab(paste("PC1(",round(pca.var$variance.percent[1], 2),"%)",sep=""))+
#   ylab(paste("PC2(",round(pca.var$variance.percent[2], 2),"%)",sep=""))
# 
# 
# p2_rep_fam<-ggplot(as.data.frame(pca.res$x),aes(x=PC3,y=PC4, col=fam_oyster, shape=rep))+
#   geom_point(size=2)+
#   scale_shape_manual(values = c(15,17))+
#   xlab(paste("PC3(",round(pca.var$variance.percent[3], 2),"%)",sep=""))+
#   ylab(paste("PC4(",round(pca.var$variance.percent[4], 2),"%)",sep=""))
# 
# 
# p1_gen_rep<-ggplot(as.data.frame(pca.res$x),aes(x=PC1,y=PC2, col=rep, shape=generation))+
#   geom_point(size=2)+ 
#   scale_shape_manual(values = c(0,15,17,18, 1,2,5,6,19))+
#   xlab(paste("PC1(",round(pca.var$variance.percent[1], 2),"%)",sep=""))+
#   ylab(paste("PC2(",round(pca.var$variance.percent[2], 2),"%)",sep=""))
# 
# 
# p2_gen_rep<-ggplot(as.data.frame(pca.res$x),aes(x=PC3,y=PC4, col=rep, shape=generation))+
#   geom_point(size=2)+
#   scale_shape_manual(values = c(0,15,17,18, 1,2,5,6,19))+
#   xlab(paste("PC3(",round(pca.var$variance.percent[3], 2),"%)",sep=""))+
#   ylab(paste("PC4(",round(pca.var$variance.percent[4], 2),"%)",sep=""))

x11()
p1_gen_fam

x11()
p2_gen_fam


#####plot individuals as a PCA for GT3 samples 
df_GT3 <- df_allFreq[6:29]
df_GT3<- df_GT3 %>% mutate_at(1:24, as.numeric)

var <- apply(df_GT3, 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_GT3))>0)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf_GT3<-df_GT3[-ind,]
}

dat_GT3 = na.omit(2*asin(sqrt(my_freqf_GT3)))

sampleIDs<-colnames(dat_GT3)
splits <- unlist(strsplit(sampleIDs, "[_]"))
#fam_oyster <- factor(splits[grep("GT",splits)])
splits <- splits[-grep("T",splits)]
generation_GT3 <- as.factor(splits[grep("G",splits)])

rep_GT3<- factor(splits[grep("R",splits)])
pcadata_GT3<-as.data.frame(t(dat_GT3))
row.names(pcadata_GT3)<-sampleIDs


pca.res <- prcomp(pcadata_GT3, retx=TRUE, center = TRUE, scale. =TRUE)
pca.var<-get_eigenvalue(pca.res)

p1_GT3<-ggplot(as.data.frame(pca.res$x),aes(x=PC1,y=PC2, col=rep_GT3, shape=generation_GT3))+
  geom_point()+
  xlab(paste("PC1(",round(pca.var$variance.percent[1], 2),"%)",sep=""))+
  ylab(paste("PC2(",round(pca.var$variance.percent[2], 2),"%)",sep=""))


p2_GT3<-ggplot(as.data.frame(pca.res$x),aes(x=PC3,y=PC4, col=rep_GT3, shape=generation_GT3))+
  geom_point()+
  xlab(paste("PC3(",round(pca.var$variance.percent[3], 2),"%)",sep=""))+
  ylab(paste("PC4(",round(pca.var$variance.percent[4], 2),"%)",sep=""))


#####plot individuals as a PCA for GT5 samples 
df_GT5 <- df_allFreq[30:59]
df_GT5<- df_GT5 %>% mutate_at(1:30, as.numeric)

var <- apply(df_GT5, 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_GT5))>0)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf_GT5<-df_GT5[-ind,]
}

dat_GT5 = na.omit(2*asin(sqrt(my_freqf_GT5)))

sampleIDs<-colnames(dat_GT5)
splits <- unlist(strsplit(sampleIDs, "[_]"))
~#fam_oyster <- factor(splits[grep("GT",splits)])
splits <- splits[-grep("T",splits)]
generation_GT5 <- as.factor(splits[grep("G",splits)])

rep_GT5<- factor(splits[grep("R",splits)])
pcadata_GT5<-as.data.frame(t(dat_GT5))
row.names(pcadata_GT5)<-sampleIDs


pca.res <- prcomp(pcadata_GT5, retx=TRUE, center = TRUE, scale. =TRUE)
pca.var<-get_eigenvalue(pca.res)

p1_GT5<-ggplot(as.data.frame(pca.res$x),aes(x=PC1,y=PC2, col=rep_GT5, shape=generation_GT5))+
  geom_point()+
  xlab(paste("PC1(",round(pca.var$variance.percent[1], 2),"%)",sep=""))+
  ylab(paste("PC2(",round(pca.var$variance.percent[2], 2),"%)",sep=""))


p2_GT5<-ggplot(as.data.frame(pca.res$x),aes(x=PC3,y=PC4, col=rep_GT5, shape=generation_GT5))+
  geom_point()+
  xlab(paste("PC3(",round(pca.var$variance.percent[3], 2),"%)",sep=""))+
  ylab(paste("PC4(",round(pca.var$variance.percent[4], 2),"%)",sep=""))



#####plot individuals as a PCAfor GT9 samples 
df_GT9 <- df_allFreq[60:101]
df_GT9<- df_GT9 %>% mutate_at(1:42, as.numeric)

var <- apply(df_GT9, 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_GT9))>0)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf_GT9<-df_GT9[-ind,]
}

dat_GT9 = na.omit(2*asin(sqrt(my_freqf_GT9)))

sampleIDs<-colnames(dat_GT9)
splits <- unlist(strsplit(sampleIDs, "[_]"))
#fam_oyster <- factor(splits[grep("GT",splits)])
splits <- splits[-grep("T",splits)]
generation_GT9 <- as.factor(splits[grep("G",splits)])

rep_GT9<- factor(splits[grep("R",splits)])
pcadata_GT9<-as.data.frame(t(dat_GT9))
row.names(pcadata_GT9)<-sampleIDs


pca.res <- prcomp(pcadata_GT9, retx=TRUE, center = TRUE, scale. =TRUE)
pca.var<-get_eigenvalue(pca.res)

p1_GT9<-ggplot(as.data.frame(pca.res$x),aes(x=PC1,y=PC2, col=rep_GT9, shape=generation_GT9))+
  geom_point()+
  scale_shape_manual(values = c(15,16,17,18,19, 1,2,5,6))+
  xlab(paste("PC1(",round(pca.var$variance.percent[1], 2),"%)",sep=""))+
  ylab(paste("PC2(",round(pca.var$variance.percent[2], 2),"%)",sep=""))


p2_GT9<-ggplot(as.data.frame(pca.res$x),aes(x=PC3,y=PC4, col=rep_GT9, shape=generation_GT9))+
  geom_point()+
  scale_shape_manual(values = c(15,16,17,18,19, 1,2,5,6))+
  xlab(paste("PC3(",round(pca.var$variance.percent[3], 2),"%)",sep=""))+
  ylab(paste("PC4(",round(pca.var$variance.percent[4], 2),"%)",sep=""))


#------------------------------------------------------------------process datas------------------------------------------------------------------------------
#-------------2 - testing for effective size

##initial parameter
input_gen=c(0,10,20,30,40,50,60)
input_gen2 = c(0,60)
input_rep=2



var <- apply(df_allFreq[,6:101], 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_allFreq[,6:101]))>0)
ind3 <- which(apply(df_allCov[,6:101], 1, function(x) min(x))<5)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf<-df_allFreq[-ind,]
  my_cov<-df_allCov[-ind,]
}

df_allCov[df_allCov == 0] <- 1
my_freqf <- df_allFreq %>% mutate_at(6:101, as.numeric)
my_cov <- df_allCov %>% mutate_at(6:101, as.numeric)


sampleIDs<-colnames(my_freqf)
splits <- unlist(strsplit(sampleIDs, "_[1234567]$"))


med_Ne=NULL
for (i in 1:input_rep){
  ne=NULL
  index=splits
  af=my_freqf
  cov=my_cov
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(af)[1]
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 814)
    ne <- rbind(ne, 
                data.frame(trial = j,
                           ne = estimateNe(p0 = af[ind,colnames(af)[1]], pt = af[ind,colnames(af)[2]], 
                                           cov0 = cov[ind,colnames(af)[1]], covt = cov[ind,colnames(af)[2]],
                                           t = 5, 
                                           ploidy=2, 
                                           truncAF=0.05, 
                                           method=c("P.planI"), 
                                           poolSize=rep(1000, times=2), 
                                           Ncensus=1250)))
    off <- off+1
  }
  Ne <- median(ne$ne)
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  med_Ne<-c(med_Ne,round(median(ne$ne)))
}
med_Ne<-as.integer(med_Ne)




#-------------3 - process adapted CMH

df_allFreq <- read.csv2("./Mean_data_frame_allelic_frequencies_individuals.csv")
df_allCov <- read.csv2("./Mean_data_frame_allelic_coverage_individuals.csv")

df_allCov[df_allCov == 0] <- 1
my_freqf <- df_allFreq %>% mutate_at(6:37, as.numeric)
my_cov <- df_allCov %>% mutate_at(6:37, as.numeric)


sampleIDs<-colnames(my_freqf[6:37])
splits <- unlist(strsplit(sampleIDs, "M_"))
Generations <- splits[grep("G",splits)]

index=Generations
af=my_freqf
cov=my_cov


p.values<-adapted.cmh.test(freq = as.matrix(my_freqf[6:37]), 
                           coverage = as.matrix(my_cov[6:37]), 
                           Ne = as.integer(c(70,70)), 
                           gen = c(0,10,13,5, 0, 10, 15, 18, 5, 0, 10, 15, 20,24,28,5), 
                           repl = 1:2,
                           order = 1,
                           mincov = 1,
                           poolSize = c(3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), 
                           IntGen = TRUE)




#-------------4 - allelic trajectories

df_allFreq <- read.csv2("./Mean_data_frame_allelic_frequencies_individuals.csv")
df_allCov <- read.csv2("./Mean_data_frame_allelic_coverage_individuals.csv")

df_allFreq[,7:13]
var <- apply(df_allFreq[,7:13], 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_allFreq[,7:13]))>0)
ind3 <- which(apply(df_allCov[,7:13], 1, function(x) min(x))<5)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf_GT3<-df_allFreq[-ind,]
  my_cov_GT3<-df_allCov[-ind,]
}

my_cov_GT3[my_cov_GT3 == 0] <- 1
my_freqf_GT3 <- my_freqf_GT3 %>% mutate_at(6:13, as.numeric)
my_cov_GT3 <- my_cov_GT3 %>% mutate_at(6:13, as.numeric)



sampleIDs<-colnames(my_cov_GT3[,6:13])
splits <- unlist(strsplit(sampleIDs, "M_GT3_"))
Generations <- splits[grep("G",splits)]
input_rep=2
med_Ne=NULL
for (i in 1:input_rep){
  ne=NULL
  index=Generations
  af=my_freqf_GT3[6:13]
  cov=my_cov_GT3[6:13]
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(af)[1]
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 476)
    ne <- rbind(ne, 
                data.frame(trial = j,
                           ne = estimateNe(p0 = af[ind,colnames(af)[1]], pt = af[ind,colnames(af)[2]], 
                                           cov0 = cov[ind,colnames(af)[1]], covt = cov[ind,colnames(af)[2]],
                                           t = 5, 
                                           ploidy=1, 
                                           truncAF=0.05, 
                                           method=c("P.planI"), 
                                           poolSize=rep(1000, times=2), 
                                           Ncensus=1250)))
    off <- off+1
  }
  Ne <- median(ne$ne)
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  med_Ne<-c(med_Ne,round(median(ne$ne)))
}
med_Ne<-as.integer(med_Ne)



p.values<-adapted.cmh.test(freq = as.matrix(my_freqf_GT3[6:13]), 
                           coverage = as.matrix(my_cov_GT3[6:13]), 
                           Ne = as.integer(c(85,85)), 
                           gen = c(0,10,13,5), 
                           repl = 1:2,
                           order = 1,
                           mincov = 1,poolSize = c(3,3,3,3,3,3,3,3), IntGen = TRUE)


adj_p<-p.adjust(p.values,method = "fdr", n=length(p.values))
p.values<-cbind(my_freqf_GT3$POS, p.values, adj_p)
colnames(p.values) <- c("pos", "p.values", "adj_p")
p.values <- as.data.frame(p.values)



#------------------------------------------------------------------plot CMH and allelic trajectories datas------------------------------------------------------------------------------
p1_GT3<-ggplot(p.values, 
           aes(x=pos, y=-log10(adj_p)))+
  geom_point(cex=1)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "CMH", x="position", y ="-log10(FDR)") + theme_classic()





p.values_GT5<-adapted.cmh.test(freq = as.matrix(my_freqf[14:23]), 
                           coverage = as.matrix(my_cov[14:23]), 
                           Ne = as.integer(c(70,70)), 
                           gen = c(0,10,15,18,5), 
                           repl = 1:2,
                           order = 1,
                           mincov = 1,poolSize = c(3,3,3,3,3,3,3,3,3,3), IntGen = TRUE)


adj_p_GT5<-p.adjust(p.values_GT5,method = "fdr", n=length(p.values_GT5))
p.values_GT5 <-cbind(my_freqf$POS, p.values_GT5, adj_p_GT5)
colnames(p.values_GT5) <- c("pos", "p.values", "adj_p")
p.values_GT5 <- as.data.frame(p.values_GT5)


p1_GT5<-ggplot(p.values_GT5, 
           aes(x=pos, y=-log10(adj_p)))+
  geom_point(cex=1)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "CMH", x="position", y ="-log10(FDR)") + theme_classic() +
  xlim(1, 186683)



p.values_GT9<-adapted.cmh.test(freq = as.matrix(my_freqf[24:37]), 
                               coverage = as.matrix(my_cov[24:37]), 
                               Ne = as.integer(c(70,70)), 
                               gen = c(0,10,15,20,24,28,5), 
                               repl = 1:2,
                               order = 1,
                               mincov = 1,poolSize = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3), IntGen = TRUE)


adj_p_GT9<-p.adjust(p.values_GT9,method = "fdr", n=length(p.values_GT9))
p.values_GT9 <-cbind(my_freqf$POS, p.values_GT9, adj_p_GT9)
colnames(p.values_GT9) <- c("pos", "p.values", "adj_p")
p.values_GT9 <- as.data.frame(p.values_GT9)

p1_GT9<-ggplot(p.values_GT9, 
           aes(x=pos, y=-log10(adj_p)))+
  geom_point(cex=1)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "CMH", x="position", y ="-log10(FDR)") + theme_classic()



x11()
p1_GT3


x11()
p1_GT5

x11()
p1_GT9

figure <- ggarrange(p1_GT3, p1_GT5, p1_GT9 + font("x.text", size = 10),
                    ncol = 1, nrow = 3)



#### addd informations of genomic architecture of OsHV-1
Unique_long=164344
Repeat_long=7411
X_region=1513
Repeat_short=10042
Unique_short=3373
nr_genoe_size = Unique_long+Repeat_long+X_region+Repeat_short+Unique_short


UL_color = "gray79" ###"#999999"
RL_color = "olivedrab3"
X_color = "lightsalmon1"
US_color = "gray79"
RS_color = "red3"


# Genome drawing
OsHV_uVarA_scheme <- data.table::fread(glue::glue("./OsHV1_annotations2.csv") ) %>%
  dplyr::mutate(direction = if_else(strand == "forward",1,-1 ))


OsHV_structure <- OsHV_uVarA_scheme %>%
  filter(molecule == "OsHV genome") %>%
  ggplot(aes(xmin = start, xmax = end, y = yposition, fill=region , forward = direction, label = name )) +
  #facet_wrap(~ molecule, ncol = 1) +
  geom_gene_label(align = "left") +
  geom_gene_arrow(arrowhead_height = unit(16, "mm"),
                  arrowhead_width = unit(3, "mm"),
                  arrow_body_height= unit(16, "mm")) +
  scale_fill_manual(values=list(UL = UL_color, 
                                RL = RL_color,
                                X = X_color,
                                RS = RS_color,
                                US = US_color)) +
  labs(x = "",
       y = "") +
  theme_classic()+
  xlim(7585, 197665) +
  ylim(0,1) +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank(), 
        axis.ticks=element_blank())

x11()
OsHV_structure
