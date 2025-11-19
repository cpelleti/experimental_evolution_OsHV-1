#!/usr/bin/Rscript
##title: "analysis of mortalities and load viral in mantle"
##author: "Pelletier Camille"
#date: 06/05/2022


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
  c( "ggplot2", "tidyr","ggpubr","tidyverse","rstatix") # plot
) )

addline_format <- function(x,...){
  gsub('//s','/n',x)
}

###################################### MORTALITIES ###############################################
#------------------------------------------------------------------import mortality datas---------------------------------------------------------------------

mortalities<- read.csv2("./data_mortality.csv", header=T)
viralLoad<- read.csv2("./load_viral_TR463.csv", header=F)
meanViralLoad<- read.csv2("./load_viral_moyen_TR463.csv", header=T) 


#------------------------------------------------------------------prepare datas------------------------------------------------------------------------------
mat_mortalities <- matrix(data=0,nrow=nrow(mortalities)*(ncol(mortalities)-1),ncol=4)
colnames(mat_mortalities) <- c("Fam_Gen", "Fam", "Gen", "Morta")

compteur=1
for (i in 1:nrow(mortalities)){
  Fam_i <- mortalities[i,1]
  for (j in 2:ncol(mortalities)){
    Gen_name <- colnames(mortalities)[j]
    Gen_i <- mortalities[i,j]
    
    mat_mortalities[compteur,1] <- paste0(Fam_i,"_",Gen_name)
    mat_mortalities[compteur,4] <- Gen_i
    mat_mortalities[compteur,2] <- Fam_i
    mat_mortalities[compteur,3] <-Gen_name
    compteur <- compteur+1
  }
  
}


mat_mortalities <- na.omit(as.data.frame(mat_mortalities))
mat_mortalities <- as.data.frame(mat_mortalities)
mat_mortalities$Morta <- as.numeric(mat_mortalities$Morta)
mat_mortalities <- mat_mortalities[mat_mortalities[,4]!=0,]


#------------------------------------------------------------------test for significance----------------------------------------------------------------------


kruskal.test(Morta ~ Gen, data = mat_mortalities)
kruskal.test(Morta ~ Fam, data = mat_mortalities)



# Construct linear model
model  <- lm(Morta ~ Gen, data = mat_mortalities)
# QQ plot of residus
ggqqplot(residuals(model))

# Normality test of Shapiro-Wilk
shapiro_test(residuals(model))


summary_stats_morta_gen<- mat_mortalities %>% group_by(Gen) %>% get_summary_stats(Morta, type="mean_sd")
summary_stats_morta_fam<- mat_mortalities %>% group_by(Fam) %>% get_summary_stats(Morta, type="mean_sd")


bxp <- ggboxplot(
  mat_mortalities, x = "Fam", y = "Morta",
  color = "Fam", palette = "jco"
)
bxp



res.aov1 <- mat_mortalities %>% anova_test(Morta ~ Fam)
res.aov1
res.aov2 <- mat_mortalities  %>% anova_test(Morta ~ Gen)
res.aov2



###################################### VIRAL LOAD ###############################################

#------------------------------------------------------------------prepare datas--------------------------------------------------------------------------

viralLoad_family <- viralLoad %>% 
  tidyr::separate(V1, c("Family","Replicate"), sep =c("_") ,remove = FALSE) %>%
  dplyr::mutate(param= glue::glue("{Family}_{V2}"))


#------------------------------------------------------------------plot datas--------------------------------------------------------------------------

p <- ggplot(viralLoad_family, aes(x=V2, y=V3, fill=Family)) + scale_y_continuous(trans='log10') +
  geom_boxplot(position=position_dodge(1))+  scale_x_discrete(limits=c("G0", "G5", "G10","G13","G15","G18","G20","G24","G28")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,
               position=position_dodge(1)) +scale_fill_manual(values=c("#0033FF", "#FF9933", "#99CC00"))+
  geom_line(data = meanViralLoad[meanViralLoad$Family == "GT3",], 
            mapping = aes(x = Generation, y = Load, group=1),color="#0033FF", linewidth=2)+
  geom_line(data = meanViralLoad[meanViralLoad$Family == "GT5",], 
            mapping = aes(x = Generation, y = Load, group=1),color="#FF9933", linewidth=2)+
  geom_line(data = meanViralLoad[meanViralLoad$Family == "GT9",], 
            mapping = aes(x = Generation, y = Load, group=1),color="#99CC00", linewidth=2) + 
  ggtitle("") +
  xlab("Infection generation") + ylab("log10(Viral copy number/ng of DNA") + theme_classic()


p

#------------------------------------------------------------------test for significance----------------------------------------------------------------------

# Construct linear model
model1 <- lm(V3 ~ V2, data = viralLoad_family)
model2  <- lm(V3 ~ Family, data = viralLoad_family)
# QQ plot of residus
ggqqplot(residuals(model1))
ggqqplot(residuals(model2))
# Normality test of Shapiro-Wilk
shapiro_test(residuals(model1))
shapiro_test(residuals(model2))

kruskal.test(V3 ~ V2, data = viralLoad_family)
kruskal.test(V3 ~ Family, data = viralLoad_family)


summary_stats_loadVir_gen<- viralLoad_family %>% group_by(V2) %>% get_summary_stats(V3, type="mean_sd")
summary_stats_loadVir_fam<- viralLoad_family %>% group_by(Family) %>% get_summary_stats(V3, type="mean_sd")
summary_stats_loadVir_fam_Gen<- viralLoad_family %>% group_by(Family,V2) %>% get_summary_stats(V3, type="mean_sd")

res.aov3 <- viralLoad_family %>% anova_test(V3 ~ Family)
res.aov3
res.aov4 <- viralLoad_family  %>% anova_test(V3 ~ V2)
res.aov4




anova_mat<- rbind(res.aov1, res.aov2, res.aov3, res.aov4)
anova_mat[1,1] <- "Mortalities ~ oyster populations"
anova_mat[2,1] <- "Mortalities ~ infection generations"
anova_mat[3,1] <- "Viral DNA copy number ~ oyster populations"
anova_mat[4,1] <- "Viral DNA copy number ~ infection generations"

anova_mat





