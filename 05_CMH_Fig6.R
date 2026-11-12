#!/usr/bin/Rscript
##title: "Adapted CMH, pvalues and plotting"
##author: "Pelletier Camille"
##Figures: Fig6
#date: 12/11/2026

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
  c( "tidyverse", "dplyr", "factoextra", "poolSeq", "ACER", "ggpubr", "gggenes") # plot
) )

addline_format <- function(x,...){
  gsub('//s','/n',x)
}

################################## data import#############################################################################

vcf_file <- read.csv2("./merge_vcf_files_evol_expe_decompose_blocksub_sort_uniq_modif.csv")

#create a matrix of mean number for each generations and family 
mat_allFreq <- matrix(data=0,nrow=nrow(vcf_file),ncol=101)
colnames(mat_allFreq) <- c("UniqueID", "POS", "REF", "ALT", "TYPE", colnames(vcf_file[10:105]))

mat_allCov <- matrix(data=0,nrow=nrow(vcf_file),ncol=101)
colnames(mat_allCov) <- c("UniqueID", "POS", "REF", "ALT", "TYPE", colnames(vcf_file[10:105]))

i=0
for (i in 1:nrow(vcf_file)){
  TYPE_i <- vcf_file[i,] %>% tidyr::separate(INFO, c("avant","après"), sep =c(";TYPE=")) %>%
    select(-avant) %>%
    tidyr::separate(après, c("TYPE","après2"), sep =c(";technology")) %>% select(TYPE)
  
  REF <- vcf_file$REF[i]
  ALT <- vcf_file$ALT[i]
  POS <- vcf_file$POS[i]
  UNIQUEID <- vcf_file[i,] %>% dplyr::mutate(UNIQUEID= glue::glue("{POS}_{REF}>{ALT}")) %>% select(UNIQUEID)
  FORMAT <- vcf_file$FORMAT[i]
  
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
  
  
  for (j in 10:105){
    FORMAT_j_tmp <- vcf_file[i,j] %>% str_split(., ":")
    print(FORMAT_j_tmp)
    
    if (FORMAT == "GT:DP:AD:RO:QR:AO:QA:GL" ){
      
      AO_j <- as.numeric(FORMAT_j_tmp[[1]][6])
      DP_j <- as.numeric(FORMAT_j_tmp[[1]][2])
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

df_allFreq_GT3 <- df_allFreq[6:29]
df_allFreq_GT5 <- df_allFreq[30:59]
df_allFreq_GT9 <- df_allFreq[60:101]

df_allCov <- as.data.frame(mat_allCov)
df_allCov<- df_allCov %>% mutate_at(6:101, as.numeric)
df_ind_Cov <- df_allCov[6:101]

df_allCov_GT3 <- df_allCov[6:29]
df_allCov_GT5 <- df_allCov[30:59]
df_allCov_GT9 <- df_allCov[60:101]

write.csv2(df_allFreq, "./allelic_frequencies.csv")
write.csv2(df_allCov, "./allelic_coverage.csv")
####################################allelic trajectories #######################################

df_allFreq <- read.csv2("./Mean_data_frame_allelic_frequencies_individuals.csv")
df_allCov <- read.csv2("./Mean_data_frame_allelic_coverage_individuals.csv")


#------------------------------------------------------------------process datas------------------------------------------------------------------------------
#-------------2 - testing for effective size

#############GT3
df_allFreq[,6:13]
var <- apply(df_allFreq[,6:13], 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_allFreq[,6:13]))>0)
ind3 <- which(apply(df_allCov[,6:13], 1, function(x) min(x))<5)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf<-df_allFreq[-ind,]
  my_cov<-df_allCov[-ind,]
}

my_cov[my_cov == 0] <- 1
my_freqf_GT3 <- my_freqf[,1:13] %>% mutate_at(6:13, as.numeric)
my_cov_GT3 <- my_cov[,1:13] %>% mutate_at(6:13, as.numeric)



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
                                           t = 13, 
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
med_Ne_GT3<-as.integer(med_Ne)



########################################GT5
df_allFreq[,14:23]
var <- apply(df_allFreq[,14:23], 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_allFreq[,14:23]))>0)
ind3 <- which(apply(df_allCov[,14:23], 1, function(x) min(x))<5)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf<-df_allFreq[-ind,]
  my_cov<-df_allCov[-ind,]
}

my_cov[my_cov == 0] <- 1
my_freqf_GT5 <- cbind(my_freqf[,1:5],my_freqf[,14:23]) %>% mutate_at(6:15, as.numeric)
my_cov_GT5 <- cbind(my_cov[,1:5],my_cov[,14:23]) %>% mutate_at(6:15, as.numeric)

sampleIDs<-colnames(my_cov_GT5[,6:15])
splits <- unlist(strsplit(sampleIDs, "M_GT5_"))
Generations <- splits[grep("G",splits)]
input_rep=2
med_Ne=NULL

for (i in 1:input_rep){
  ne=NULL
  index=Generations
  af=my_freqf_GT5[6:15]
  cov=my_cov_GT5[6:15]
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
                                           t = 18, 
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
med_Ne_GT5<-as.integer(med_Ne)

#########################GT9
df_allFreq[,24:37]
var <- apply(df_allFreq[,24:37], 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_allFreq[,24:37]))>0)
ind3 <- which(apply(df_allCov[,24:37], 1, function(x) min(x))<5)
ind<-sort(unique(c(ind1,ind2)))

if (length(ind) != 0) {
  my_freqf<-df_allFreq[-ind,]
  my_cov<-df_allCov[-ind,]
}

my_cov[my_cov == 0] <- 1
my_freqf_GT9 <- cbind(my_freqf[,1:5],my_freqf[,24:37]) %>% mutate_at(6:19, as.numeric)
my_cov_GT9 <- cbind(my_cov[,1:5],my_cov[,24:37]) %>% mutate_at(6:19, as.numeric)

sampleIDs<-colnames(my_cov_GT9[,6:19])
splits <- unlist(strsplit(sampleIDs, "M_GT9_"))
Generations <- splits[grep("G",splits)]
input_rep=2
med_Ne=NULL

for (i in 1:input_rep){
  ne=NULL
  index=Generations
  af=my_freqf_GT9[6:19]
  cov=my_cov_GT9[6:19]
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
                                           t = 24, 
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
med_Ne_GT9<-as.integer(med_Ne)


#-------------3 - compute pvalues from adapted CMH
####GT3
p.values<-adapted.cmh.test(freq = as.matrix(my_freqf_GT3[6:13]), 
                           coverage = as.matrix(my_cov_GT3[6:13]), 
                           Ne = med_Ne_GT3, 
                           gen = c(0,10,13,5), 
                           repl = 1:2,
                           order = 1,
                           mincov = 1,poolSize = c(3,3,3,3,3,3,3,3), IntGen = TRUE)


adj_p<-p.adjust(p.values,method = "fdr", n=length(p.values))
p.values<-cbind(my_freqf_GT3$POS, p.values, adj_p)
colnames(p.values) <- c("pos", "p.values", "adj_p")
p.values_GT3 <- as.data.frame(p.values)

#####GT5
p.values_GT5<-adapted.cmh.test(freq = as.matrix(my_freqf_GT5[6:15]), 
                               coverage = as.matrix(my_cov_GT5[6:15]), 
                               Ne = med_Ne, 
                               gen = c(0,10,15,18,5), 
                               repl = 1:2,
                               order = 1,
                               mincov = 1,poolSize = c(3,3,3,3,3,3,3,3,3,3), IntGen = TRUE)


adj_p_GT5<-p.adjust(p.values_GT5,method = "fdr", n=length(p.values_GT5))
p.values_GT5 <-cbind(my_freqf_GT5$POS, p.values_GT5, adj_p_GT5)
colnames(p.values_GT5) <- c("pos", "p.values", "adj_p")
p.values_GT5 <- as.data.frame(p.values_GT5)



#####GT9
p.values_GT9<-adapted.cmh.test(freq = as.matrix(my_freqf[24:37]), 
                               coverage = as.matrix(my_cov[24:37]), 
                               Ne = med_Ne, 
                               gen = c(0,10,15,20,24,28,5), 
                               repl = 1:2,
                               order = 1,
                               mincov = 1,poolSize = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3), IntGen = TRUE)


adj_p_GT9<-p.adjust(p.values_GT9,method = "fdr", n=length(p.values_GT9))
p.values_GT9 <-cbind(my_freqf$POS, p.values_GT9, adj_p_GT9)
colnames(p.values_GT9) <- c("pos", "p.values", "adj_p")
p.values_GT9 <- as.data.frame(p.values_GT9)


#-------------4 - plot adjusted pvalues
p1_GT3<-ggplot(p.values_GT3, 
               aes(x=pos, y=-log10(adj_p)))+
  geom_point(cex=1)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "FA", x="position", y ="-log10(FDR)") + theme_classic() + ylim(0,2)

p1_GT3


p1_GT5<-ggplot(p.values_GT5, 
               aes(x=pos, y=-log10(adj_p)))+
  geom_point(cex=1)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "C", x="position", y ="-log10(FDR)") + theme_classic() +
  xlim(1, 186683)+ ylim(0,2)

p1_GT5

p1_GT9<-ggplot(p.values_GT9, 
               aes(x=pos, y=-log10(adj_p)))+
  geom_point(cex=1, aes(colour= ifelse(-log10(adj_p)>=1.3, 'red', 'black')))+
  scale_colour_manual(labels = c("<1.3", ">1.3"), values=c('black', 'red')) +
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "NFA", x="position", y ="-log10(FDR)") + theme_classic() + ylim(0,2) +
  theme(legend.position = "none")


p1_GT9

figure <- ggarrange(p1_GT3, p1_GT5, p1_GT9 + font("x.text", size = 10),
                    ncol = 1, nrow = 3)


signif_CMH_GT9 <- p.values_GT9[-log10(p.values_GT9$adj_p) > -log10(0.05),]
not_signif_CMH_GT9 <- p.values_GT9[-log10(p.values_GT9$adj_p) < -log10(0.05),]
write.csv2(signif_CMH_GT9, "./significant_pvalues_CMH_GT9.csv")

#-------------5- plot adjusted pvalues colored according to ORF functions
orf_func <- read.csv2("./significant_pvalues_CMH_GT9_orf_function.csv")
signif_CMH_GT9 <- read.csv2( "./significant_pvalues_CMH_GT9.csv")
OsHV_uVarA_scheme <- data.table::fread(glue::glue("./OsHV1_annotations2.csv") ) %>%
  dplyr::mutate(direction = if_else(strand == "forward",1,-1 )) 
OsHV_uVarA_scheme <- OsHV_uVarA_scheme[OsHV_uVarA_scheme$molecule == "OsHV genome",]
OsHV_uVarA_scheme$start <- as.numeric(OsHV_uVarA_scheme$start)
OsHV_uVarA_scheme$end <- as.numeric(OsHV_uVarA_scheme$end)

orf_func$start <- as.numeric(orf_func$start)
orf_func$stop <- as.numeric(orf_func$stop)
i=0
for (i in 1:nrow(signif_CMH_GT9)){
  line_i <- signif_CMH_GT9[i,]
  pos_i <- line_i$pos
  if (line_i$Function != "intergenic" ){
    orf_func_i <- orf_func[orf_func$start < pos_i & orf_func$stop > pos_i,]
    signif_CMH_GT9$color[i] <- orf_func_i$Colors
  }else if (line_i$Function == "intergenic")
    signif_CMH_GT9$color[i] <- "#FFCC99"
}



p1_GT9_signif<- ggplot()+
  geom_point(data=not_signif_CMH_GT9,aes(x=pos, y=-log10(adj_p)),colour= 'black')+
  geom_point(data=signif_CMH_GT9, aes(x=pos, y=-log10(adj_p)), col=signif_CMH_GT9$color)+
  geom_rect(data=orf_func, aes(xmin =start, xmax=stop, ymin=rep(-0.1,nrow(orf_func)), ymax=rep(-0.05,nrow(orf_func)) ), fill=orf_func$Colors)+
  geom_rect(aes(xmin =1, xmax=163344, ymin=rep(-0.11,nrow(OsHV_uVarA_scheme)), ymax=rep(-0.14,nrow(OsHV_uVarA_scheme))), fill="gray79" )+
  geom_rect(aes(xmin =163345, xmax=171755, ymin=rep(-0.11,nrow(OsHV_uVarA_scheme)), ymax=rep(-0.14,nrow(OsHV_uVarA_scheme))), fill="olivedrab3" )+
  geom_rect(aes(xmin =171756, xmax=173268, ymin=rep(-0.11,nrow(OsHV_uVarA_scheme)), ymax=rep(-0.14,nrow(OsHV_uVarA_scheme))), fill="lightsalmon1" )+
  geom_rect(aes(xmin =173269, xmax=183310, ymin=rep(-0.11,nrow(OsHV_uVarA_scheme)), ymax=rep(-0.14,nrow(OsHV_uVarA_scheme))), fill="red3" )+
  geom_rect(aes(xmin =183311, xmax=186683, ymin=rep(-0.11,nrow(OsHV_uVarA_scheme)), ymax=rep(-0.14,nrow(OsHV_uVarA_scheme))), fill="gray79" )+
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "NFA", x="position", y ="-log10(FDR)") + theme_classic() + ylim(-0.2,2)

figure <- ggarrange(p1_GT3, p1_GT5, p1_GT9_signif + font("x.text", size = 10),
                    ncol = 1, nrow = 3)

figure

