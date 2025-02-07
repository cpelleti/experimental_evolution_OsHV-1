#!/usr/bin/Rscript
##title: "allelic trajectories - Figure 9"
##author: "Pelletier Camille"
#date: 20/09/2022
rm(list=ls())
'%!in%' <- function(x,y)!('%in%'(x,y))
options( "digits"=3)

#------------------------------------------------------------------import datas-------------------------------------------------------------------------------
setwd(/path/to/folder/output/02_allelic_frequencies)

df_allFreq <- read.csv2("./Mean_data_frame_allelic_frequencies_individuals.csv")
df_allCov <- read.csv2("./Mean_data_frame_allelic_coverage_individuals.csv")


#------------------------------------------------------------------prepare datas------------------------------------------------------------------------------
var <- apply(df_allFreq[,7:13], 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_allFreq[,7:13]))>0)
ind3 <- which(apply(df_allCov[,7:13], 1, function(x) min(x))<5)
ind<-sort(unique(c(ind1,ind2)))



##prepare tables for GT3 population
if (length(ind) != 0) {
  my_freqf_GT3<-df_allFreq[-ind,]
  my_cov_GT3<-df_allCov[-ind,]
}

my_cov_GT3[my_cov_GT3 == 0] <- 1
my_freqf_GT3 <- my_freqf_GT3 %>% mutate_at(6:13, as.numeric)

my_cov_GT3 <- my_cov_GT3 %>% mutate_at(6:13, as.numeric)



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
colnames(p.values) <- c("POS", "p.values", "adj_p")


freqs <- merge(my_freqf_GT3, p.values, by = c("POS"))
freqs$POS= as.numeric(freqs$POS)



##prepare tables for GT9 population
var <- apply(df_allFreq[,24:37], 1, function(x) var(na.omit(x)))
ind1 <- which(var == 0)
ind2 <- which(rowSums(is.na(df_allFreq[,24:37]))>0)
ind3 <- which(apply(df_allCov[,24:37], 1, function(x) min(x))<5)
ind<-sort(unique(c(ind1,ind2)))

df_allFreq_GT9 <- cbind(df_allFreq[,2],df_allFreq[,24:37])
df_allCov_GT9 <- cbind(df_allCov[,2],df_allCov[,24:37])
if (length(ind) != 0) {
  my_freqf_GT9<-df_allFreq_GT9[-ind,]
  my_cov_GT9<-df_allCov_GT9[-ind,]
}

my_cov_GT9[my_cov_GT9 == 0] <- 1
my_freqf_GT9 <- my_freqf_GT9 %>% mutate_at(2:15, as.numeric)
colnames(my_freqf_GT9) <- c("POS", colnames(my_freqf_GT9[2:15]))
my_cov_GT9 <- my_cov_GT9 %>% mutate_at(2:15, as.numeric)



p.values_GT9<-adapted.cmh.test(freq = as.matrix(my_freqf_GT9[2:15]), 
                               coverage = as.matrix(my_cov_GT9[2:15]), 
                               Ne = as.integer(c(70,70)), 
                               gen = c(0,10,15,20,24,28,5), 
                               repl = 1:2,
                               order = 1,
                               mincov = 1,poolSize = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3), IntGen = TRUE)


adj_p_GT9<-p.adjust(p.values_GT9,method = "fdr", n=length(p.values_GT9))
p.values_GT9 <-cbind(my_freqf_GT9$POS, p.values_GT9, adj_p_GT9)
colnames(p.values_GT9) <- c("POS", "p.values", "adj_p")
p.values_GT9 <- as.data.frame(p.values_GT9)
p.values_GT9 <- p.values_GT9 %>% mutate_at(1:3, as.numeric)


p1_GT9<-ggplot(p.values_GT9, 
               aes(x=POS, y=-log10(adj_p_GT9)))+
  geom_point(cex=1)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype=2)+
  labs(title = "GT9", x="position", y ="-log10(FDR)") + theme_classic()



freqs_GT9 <- merge(my_freqf_GT9, p.values_GT9, by = c("POS"))
freqs_GT9$POS= as.numeric(freqs_GT9$POS)

#------------------------------------------------------------------process datas, replicate 1------------------------------------------------------------------------------
#-------------1 - allelic trajectories on the global artifical NR-consensus (only for top snps)
region = freqs_GT9[ freqs_GT9$adj_p<0.05,]
topsnps = region[order(region$adj_p),][1:100,]


#median frequency start and end
topsnps$G0.median =apply(topsnps[,c(2:3)], 1, median)
topsnps$G28.median=apply(topsnps[,c(12:13)], 1, median)
#save afc
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#if negative, change allele frequency
topsnps[topsnps$AFC < 0,c(2:15,20,19)]=1-topsnps[topsnps$AFC < 0,c(2:15,20,19)]
#restimate AFC
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#topsnps <- topsnps[topsnps$AFC >0,]

gens=c("G0","G05","G10","G15","G20","G24","G28")
nrep=2
 

topsnps_R1 <- cbind(topsnps[1:2],topsnps[14], topsnps[4], topsnps[6], topsnps[8], topsnps[10],topsnps[12])


mat_R1 <- matrix(data=NA,nrow=nrow(topsnps_R1)*7,ncol=3)
colnames(mat_R1) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R1)){
  
  line <- topsnps_R1[i,]
    for (j in 1:7){
  mat_R1[k,3] <- line[,j+1]
  mat_R1[k,1] <- line[,1]
  gen <- gens[[j]]
  mat_R1[k,2] <- as.character(gen)
    k<- k+1
    
  }
}
df_R1 <- as.data.frame(mat_R1)
df_R1$pos <- as.character(df_R1$pos)
df_R1$gen <- as.character(df_R1$gen)
df_R1$freq <- as.numeric(df_R1$freq)

medianaf=apply(topsnps_R1[,2:8], 2, median)
medianaf <- as.data.frame(medianaf)
row.names(medianaf) <- gens
medianaf<- cbind(gens, medianaf, rep("medianaf", 7))
colnames(medianaf) <- c("gen", "medianaf", "pos")


#------------------------------------------------------------------plot top snps as allellic trajectories, replicate 1-----------------------------------------------------
p <- ggplot() +
  geom_line(data=df_R1, aes(x=gen, y=freq, group=pos),colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")+
  geom_line(data = medianaf,aes(x=gen, y=medianaf, group=NA), col = "red", linewidth=1)


x11()
p

#------------------------------------------------------------------process datas, replicate 2------------------------------------------------------------------------------
#-------------1 - allelic trajectories on the global artifical NR-consensus (only for top snps)

topsnps_R2 <- cbind(topsnps[1],topsnps[3],topsnps[15], topsnps[5], topsnps[7], topsnps[9], topsnps[11],topsnps[13])
mat_R2 <- matrix(data=NA,nrow=nrow(topsnps_R2)*7,ncol=3)
colnames(mat_R2) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R2)){
  
  line <- topsnps_R2[i,]
  
  for (j in 1:7){
    mat_R2[k,3] <- line[,j+1]
    mat_R2[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R2[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R2 <- as.data.frame(mat_R2)
df_R2$pos <- as.character(df_R2$pos)
df_R2$gen <- as.character(df_R2$gen)
df_R2$freq <- as.numeric(df_R2$freq)

medianaf=apply(topsnps_R2[,2:8], 2, median)
medianaf <- as.data.frame(medianaf)
row.names(medianaf) <- gens
medianaf<- cbind(gens, medianaf, rep("medianaf", 7))
colnames(medianaf) <- c("gen", "medianaf", "pos") 

#------------------------------------------------------------------plot top snps as allellic trajectories, replicate 2-----------------------------------------------------
p <- ggplot(data=df_R2, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")+
  geom_line(data = medianaf,aes(x=gen, y=medianaf, group=NA), col = "red", linewidth=1)


x11()
p


#------------------------------------------------------------------plot top snps as allellic trajectories, in specific regions with high snps density-----------------------------------------------------
######################allelic trajectories between 50000 and 80000 pb #################################
start=50000
end=80000
region = freqs_GT9[freqs_GT9$POS>start & freqs_GT9$POS<end & freqs_GT9$adj_p<0.05,]
topsnps = region[order(region$adj_p),][1:22,]


#median frequency start and end
topsnps$G0.median =apply(topsnps[,c(2:3)], 1, median)
topsnps$G28.median=apply(topsnps[,c(12:13)], 1, median)
#save afc
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#if negative, change allele frequency
topsnps[topsnps$AFC < 0,c(2:15,20,19)]=1-topsnps[topsnps$AFC < 0,c(2:15,20,19)]
#restimate AFC
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#topsnps <- topsnps[topsnps$AFC >0,]

gens=c("G0","G05","G10","G15","G20","G24","G28")
nrep=2


topsnps_R1 <- cbind(topsnps[1:2],topsnps[14], topsnps[4], topsnps[6], topsnps[8], topsnps[10],topsnps[12])

mat_R1 <- matrix(data=NA,nrow=nrow(topsnps_R1)*7,ncol=3)
colnames(mat_R1) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R1)){
  
  line <- topsnps_R1[i,]
  
  for (j in 1:7){
    mat_R1[k,3] <- line[,j+1]
    mat_R1[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R1[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R1 <- as.data.frame(mat_R1)
df_R1$pos <- as.character(df_R1$pos)
df_R1$gen <- as.character(df_R1$gen)
df_R1$freq <- as.numeric(df_R1$freq)

p <- ggplot(data=df_R1, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p


topsnps_R2 <- cbind(topsnps[1],topsnps[3],topsnps[15], topsnps[5], topsnps[7], topsnps[9], topsnps[11],topsnps[13])
mat_R2 <- matrix(data=NA,nrow=nrow(topsnps_R2)*7,ncol=3)
colnames(mat_R2) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R2)){
  
  line <- topsnps_R2[i,]
  
  for (j in 1:7){
    mat_R2[k,3] <- line[,j+1]
    mat_R2[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R2[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R2 <- as.data.frame(mat_R2)
df_R2$pos <- as.character(df_R2$pos)
df_R2$gen <- as.character(df_R2$gen)
df_R2$freq <- as.numeric(df_R2$freq)

p <- ggplot(data=df_R2, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p


######################allelic trajectories between 97000 and 110000 pb #################################
start=97000
end=110000
region = freqs_GT9[freqs_GT9$POS>start & freqs_GT9$POS<end & freqs_GT9$adj_p<0.05,]
topsnps = region[order(region$adj_p),][1:8,]


#median frequency start and end
topsnps$G0.median =apply(topsnps[,c(2:3)], 1, median)
topsnps$G28.median=apply(topsnps[,c(12:13)], 1, median)
#save afc
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#if negative, change allele frequency
topsnps[topsnps$AFC < 0,c(2:15,20,19)]=1-topsnps[topsnps$AFC < 0,c(2:15,20,19)]
#restimate AFC
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#topsnps <- topsnps[topsnps$AFC >0,]

gens=c("G0","G05","G10","G15","G20","G24","G28")
nrep=2


topsnps_R1 <- cbind(topsnps[1:2],topsnps[14], topsnps[4], topsnps[6], topsnps[8], topsnps[10],topsnps[12])

mat_R1 <- matrix(data=NA,nrow=nrow(topsnps_R1)*7,ncol=3)
colnames(mat_R1) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R1)){
  
  line <- topsnps_R1[i,]
  
  for (j in 1:7){
    mat_R1[k,3] <- line[,j+1]
    mat_R1[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R1[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R1 <- as.data.frame(mat_R1)
df_R1$pos <- as.character(df_R1$pos)
df_R1$gen <- as.character(df_R1$gen)
df_R1$freq <- as.numeric(df_R1$freq)

p <- ggplot(data=df_R1, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,0.10)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p


topsnps_R2 <- cbind(topsnps[1],topsnps[3],topsnps[15], topsnps[5], topsnps[7], topsnps[9], topsnps[11],topsnps[13])
mat_R2 <- matrix(data=NA,nrow=nrow(topsnps_R2)*7,ncol=3)
colnames(mat_R2) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R2)){
  
  line <- topsnps_R2[i,]
  
  for (j in 1:7){
    mat_R2[k,3] <- line[,j+1]
    mat_R2[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R2[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R2 <- as.data.frame(mat_R2)
df_R2$pos <- as.character(df_R2$pos)
df_R2$gen <- as.character(df_R2$gen)
df_R2$freq <- as.numeric(df_R2$freq)

p <- ggplot(data=df_R2, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,0.10)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p
######################allelic trajectories in the IRL  #################################
start=164345
end=171680
region = freqs_GT9[freqs_GT9$POS>start & freqs_GT9$POS<end & freqs_GT9$adj_p<0.05,]
topsnps = region[order(region$adj_p),][1:45,]


#median frequency start and end
topsnps$G0.median =apply(topsnps[,c(2:3)], 1, median)
topsnps$G28.median=apply(topsnps[,c(12:13)], 1, median)
#save afc
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#if negative, change allele frequency
topsnps[topsnps$AFC < 0,c(2:15,20,19)]=1-topsnps[topsnps$AFC < 0,c(2:15,20,19)]
#restimate AFC
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#topsnps <- topsnps[topsnps$AFC >0,]

gens=c("G0","G05","G10","G15","G20","G24","G28")
nrep=2


topsnps_R1 <- cbind(topsnps[1:2],topsnps[14], topsnps[4], topsnps[6], topsnps[8], topsnps[10],topsnps[12])

ind2 <- which((rowSums(topsnps_R1[,2:8])>0.1))
if (length(ind2) != 0) {
  topsnps_R1<-topsnps_R1[ind2,]
}

mat_R1 <- matrix(data=NA,nrow=nrow(topsnps_R1)*7,ncol=3)
colnames(mat_R1) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R1)){
  
  line <- topsnps_R1[i,]
  
  for (j in 1:7){
    mat_R1[k,3] <- line[,j+1]
    mat_R1[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R1[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R1 <- as.data.frame(mat_R1)
df_R1$pos <- as.character(df_R1$pos)
df_R1$gen <- as.character(df_R1$gen)
df_R1$freq <- as.numeric(df_R1$freq)

p <- ggplot(data=df_R1, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p


topsnps_R2 <- cbind(topsnps[1],topsnps[3],topsnps[15], topsnps[5], topsnps[7], topsnps[9], topsnps[11],topsnps[13])
ind2 <- which((rowSums(topsnps_R2[,2:8])>0.1))
if (length(ind2) != 0) {
  topsnps_R2<-topsnps_R2[ind2,]
}


mat_R2 <- matrix(data=NA,nrow=nrow(topsnps_R2)*7,ncol=3)
colnames(mat_R2) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R2)){
  
  line <- topsnps_R2[i,]
  
  for (j in 1:7){
    mat_R2[k,3] <- line[,j+1]
    mat_R2[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R2[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R2 <- as.data.frame(mat_R2)
df_R2$pos <- as.character(df_R2$pos)
df_R2$gen <- as.character(df_R2$gen)
df_R2$freq <- as.numeric(df_R2$freq)

p <- ggplot(data=df_R2, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p


######################allelic trajectories in the IRS  #################################
start=173399
end=183310
region = freqs_GT9[freqs_GT9$POS>start & freqs_GT9$POS<end & freqs_GT9$adj_p<0.05,]
topsnps = region[order(region$adj_p),][1:7,]


#median frequency start and end
topsnps$G0.median =apply(topsnps[,c(2:3)], 1, median)
topsnps$G28.median=apply(topsnps[,c(12:13)], 1, median)
#save afc
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#if negative, change allele frequency
topsnps[topsnps$AFC < 0,c(2:15,20,19)]=1-topsnps[topsnps$AFC < 0,c(2:15,20,19)]
#restimate AFC
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#topsnps <- topsnps[topsnps$AFC >0,]

gens=c("G0","G05","G10","G15","G20","G24","G28")
nrep=2


topsnps_R1 <- cbind(topsnps[1:2],topsnps[14], topsnps[4], topsnps[6], topsnps[8], topsnps[10],topsnps[12])

ind2 <- which((rowSums(topsnps_R1[,2:8])>0.05))
if (length(ind2) != 0) {
  topsnps_R1<-topsnps_R1[ind2,]
}

mat_R1 <- matrix(data=NA,nrow=nrow(topsnps_R1)*7,ncol=3)
colnames(mat_R1) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R1)){
  
  line <- topsnps_R1[i,]
  
  for (j in 1:7){
    mat_R1[k,3] <- line[,j+1]
    mat_R1[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R1[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R1 <- as.data.frame(mat_R1)
df_R1$pos <- as.character(df_R1$pos)
df_R1$gen <- as.character(df_R1$gen)
df_R1$freq <- as.numeric(df_R1$freq)

p <- ggplot(data=df_R1, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p


topsnps_R2 <- cbind(topsnps[1],topsnps[3],topsnps[15], topsnps[5], topsnps[7], topsnps[9], topsnps[11],topsnps[13])
ind2 <- which((rowSums(topsnps_R2[,2:8])>0.05))
if (length(ind2) != 0) {
  topsnps_R2<-topsnps_R2[ind2,]
}


mat_R2 <- matrix(data=NA,nrow=nrow(topsnps_R2)*7,ncol=3)
colnames(mat_R2) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R2)){
  
  line <- topsnps_R2[i,]
  
  for (j in 1:7){
    mat_R2[k,3] <- line[,j+1]
    mat_R2[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R2[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R2 <- as.data.frame(mat_R2)
df_R2$pos <- as.character(df_R2$pos)
df_R2$gen <- as.character(df_R2$gen)
df_R2$freq <- as.numeric(df_R2$freq)

p <- ggplot(data=df_R2, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p

######################allelic trajectories in the US  #################################
start=183311
end=186683
region = freqs_GT9[freqs_GT9$POS>start & freqs_GT9$POS<end & freqs_GT9$adj_p<0.05,]
topsnps = region[order(region$adj_p),][1:6,]


#median frequency start and end
topsnps$G0.median =apply(topsnps[,c(2:3)], 1, median)
topsnps$G28.median=apply(topsnps[,c(12:13)], 1, median)
#save afc
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#if negative, change allele frequency
topsnps[topsnps$AFC < 0,c(2:15,20,19)]=1-topsnps[topsnps$AFC < 0,c(2:15,20,19)]
#restimate AFC
topsnps$AFC =topsnps$G28.median-topsnps$G0.median
#topsnps <- topsnps[topsnps$AFC >0,]

gens=c("G0","G05","G10","G15","G20","G24","G28")
nrep=2


topsnps_R1 <- cbind(topsnps[1:2],topsnps[14], topsnps[4], topsnps[6], topsnps[8], topsnps[10],topsnps[12])

ind2 <- which((rowSums(topsnps_R1[,2:8])>0.05))
if (length(ind2) != 0) {
  topsnps_R1<-topsnps_R1[ind2,]
}

mat_R1 <- matrix(data=NA,nrow=nrow(topsnps_R1)*7,ncol=3)
colnames(mat_R1) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R1)){
  
  line <- topsnps_R1[i,]
  
  for (j in 1:7){
    mat_R1[k,3] <- line[,j+1]
    mat_R1[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R1[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R1 <- as.data.frame(mat_R1)
df_R1$pos <- as.character(df_R1$pos)
df_R1$gen <- as.character(df_R1$gen)
df_R1$freq <- as.numeric(df_R1$freq)

p <- ggplot(data=df_R1, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p


topsnps_R2 <- cbind(topsnps[1],topsnps[3],topsnps[15], topsnps[5], topsnps[7], topsnps[9], topsnps[11],topsnps[13])
ind2 <- which((rowSums(topsnps_R2[,2:8])>0.05))
if (length(ind2) != 0) {
  topsnps_R2<-topsnps_R2[ind2,]
}


mat_R2 <- matrix(data=NA,nrow=nrow(topsnps_R2)*7,ncol=3)
colnames(mat_R2) <- c("pos", "gen", "freq")

i=1
k=1
for (i in 1:nrow(topsnps_R2)){
  
  line <- topsnps_R2[i,]
  
  for (j in 1:7){
    mat_R2[k,3] <- line[,j+1]
    mat_R2[k,1] <- line[,1]
    gen <- gens[[j]]
    mat_R2[k,2] <- as.character(gen)
    k<- k+1
  }
}
df_R2 <- as.data.frame(mat_R2)
df_R2$pos <- as.character(df_R2$pos)
df_R2$gen <- as.character(df_R2$gen)
df_R2$freq <- as.numeric(df_R2$freq)

p <- ggplot(data=df_R2, aes(x=gen, y=freq, group=pos)) +
  geom_line(colour = "blue", alpha=0.3) + 
  theme_bw()+ ylim(0,1)+
  theme(legend.position="none") + 
  labs(y="Allele Frequency", x="Generation")


x11()
p
