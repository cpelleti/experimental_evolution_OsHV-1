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



# Construire le modèle linéaire
model  <- lm(Morta ~ Gen, data = mat_mortalities)
# Créer un QQ plot des résidus
ggqqplot(residuals(model))

# Calculer le test de normalité de Shapiro-Wilk
shapiro_test(residuals(model))


mat_mortalities %>% group_by(Gen) %>% shapiro_test(Morta)


summary_stats_morta_gen<- mat_mortalities %>% group_by(Gen) %>% get_summary_stats(Morta, type="mean_sd")
summary_stats_morta_fam<- mat_mortalities %>% group_by(Fam) %>% get_summary_stats(Morta, type="mean_sd")


bxp <- ggboxplot(
  mat_mortalities, x = "Fam", y = "Morta",
  color = "Fam", palette = "jco"
)
bxp




res.aov <- mat_mortalities_ordered  %>% anova_test(Morta ~ Gen)
res.aov
res.aov <- mat_mortalities %>% anova_test(Morta ~ Fam)
res.aov

###significance test familly
GT3_GT5 <- mat_mortalities[mat_mortalities$Fam == "GT3"|mat_mortalities$Fam == "GT5_2",]
GT3_GT5 <- GT3_GT5[GT3_GT5[,4]!=0,]
t.test(Morta ~ Fam, data = GT3_GT5)
wilcox.test(Morta ~ Fam, data = GT3_GT5)


GT3_GT9 <- mat_mortalities[mat_mortalities$Fam == "GT3"|mat_mortalities$Fam == "GT9",]
GT3_GT9 <- GT3_GT9[GT3_GT9[,4]!=0,]
t.test(Morta ~ Fam, data = GT3_GT9)
wilcox.test(Morta ~ Fam, data = GT3_GT9)

GT5_GT9 <- mat_mortalities[mat_mortalities$Fam == "GT5_2"|mat_mortalities$Fam == "GT9",]
GT5_GT9 <- GT5_GT9[GT5_GT9[,4]!=0,]
t.test(Morta ~ Fam, data = GT5_GT9)
wilcox.test(Morta ~ Fam, data = GT5_GT9)

###significance test generation


results <- list()

# Nom des variables numériques pour lesquelles on veut faire le test
variables_to_test <- mat_mortalities$Fam_Gen %>% 
  
  colnames(mat_mortalities)[sapply(mat_mortalities, is.numeric) & colnames(mat_mortalities) != "Morta"]

for (var in variables_to_test) {
  
  # Sous-ensemble des données avec les familles GT3 et GT5_2
  GT3_GT5 <- mat_mortalities[mat_mortalities$Fam %in% c("GT3", "GT5_2"), c("Fam", var)]
  
  # Suppression des lignes avec valeur 0 ou NA dans la variable considérée
  GT3_GT5 <- GT3_GT5[GT3_GT5[[var]] != 0 & !is.na(GT3_GT5[[var]]), ]
  
  # Vérification qu'on a bien au moins deux familles présentes avec au moins une valeur
  if (length(unique(GT3_GT5$Fam)) == 2 && nrow(GT3_GT5) >= 3) {
    test_result <- t.test(as.formula(paste(var, "~ Fam")), data = GT3_GT5)
    results[[var]] <- list(p.value = test_result$p.value, statistic = test_result$statistic)
  } else {
    results[[var]] <- "Insufficient data"
  }
}



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
  ggtitle("Charge virale dans les échantillons au cours \n des générations d'infections") +
  xlab("Charge virale (Nombre de copies/ng ADN)") + ylab("Générations") + theme_classic()

x11()
p

#------------------------------------------------------------------test for significance----------------------------------------------------------------------

# Construire le modèle linéaire
model1 <- lm(V3 ~ V2, data = viralLoad_family)
model2  <- lm(V3 ~ Family, data = viralLoad_family)
# Créer un QQ plot des résidus
ggqqplot(residuals(model1))
ggqqplot(residuals(model2))
# Calculer le test de normalité de Shapiro-Wilk
shapiro_test(residuals(model1))
shapiro_test(residuals(model2))

kruskal.test(V3 ~ V2, data = viralLoad_family)
kruskal.test(V3 ~ Family, data = viralLoad_family)


summary_stats_loadVir_gen<- viralLoad_family %>% group_by(V2) %>% get_summary_stats(V3, type="mean_sd")
summary_stats_loadVir_fam<- viralLoad_family %>% group_by(Family) %>% get_summary_stats(V3, type="mean_sd")
summary_stats_loadVir_fam_Gen<- viralLoad_family %>% group_by(Family,V2) %>% get_summary_stats(V3, type="mean_sd")

res.aov1 <- viralLoad_family  %>% anova_test(V3 ~ V2)
res.aov1
res.aov2 <- viralLoad_family %>% anova_test(V3 ~ Family)
res.aov2

###significance test GT3
GT3_G0_G5 <- viralLoad_family[viralLoad_family$param == "GT3_G0"|viralLoad_family$param == "GT3_G5",]
t.test(V3 ~ param, data = GT3_G0_G5)
wilcox.test(V3 ~ param, data = GT3_G0_G5)

GT3_G5_G10 <- viralLoad_family[viralLoad_family$param == "GT3_G5"|viralLoad_family$param == "GT3_G10",]
t.test(V3 ~ param, data = GT3_G5_G10)
wilcox.test(V3 ~ param, data = GT3_G5_G10)

GT3_G10_G13 <- viralLoad_family[viralLoad_family$param == "GT3_G10"|viralLoad_family$param == "GT3_G13",]
t.test(V3 ~ param, data = GT3_G10_G13)
wilcox.test(V3 ~ param, data = GT3_G10_G13)

#G0-10
GT3_G0_G10 <- viralLoad_family[viralLoad_family$param == "GT3_G0"|viralLoad_family$param == "GT3_G10",]
t.test(V3 ~ param, data = GT3_G0_G10)
wilcox.test(V3 ~ param, data = GT3_G0_G10)

#G0-13
GT3_G0_G13 <- viralLoad_family[viralLoad_family$param == "GT3_G0"|viralLoad_family$param == "GT3_G13",]
t.test(V3 ~ param, data = GT3_G0_G13)
wilcox.test(V3 ~ param, data = GT3_G0_G13)

#G5-13
GT3_G5_G13 <- viralLoad_family[viralLoad_family$param == "GT3_G5"|viralLoad_family$param == "GT3_G13",]
t.test(V3 ~ param, data = GT3_G5_G13)
wilcox.test(V3 ~ param, data = GT3_G5_G13)

###significance test GT5

GT5_G0_G5 <- viralLoad_family[viralLoad_family$param == "GT5_G0"|viralLoad_family$param == "GT5_G5",]
t.test(V3 ~ param, data = GT5_G0_G5)
wilcox.test(V3 ~ param, data = GT5_G0_G5)

GT5_G5_G10 <- viralLoad_family[viralLoad_family$param == "GT5_G5"|viralLoad_family$param == "GT5_G10",]
t.test(V3 ~ param, data = GT5_G5_G10)
wilcox.test(V3 ~ param, data = GT5_G5_G10)

GT5_G10_G15 <- viralLoad_family[viralLoad_family$param == "GT5_G10"|viralLoad_family$param == "GT5_G15",]
t.test(V3 ~ param, data = GT5_G10_G15)
wilcox.test(V3 ~ param, data = GT5_G10_G15)

GT5_G15_G18 <- viralLoad_family[viralLoad_family$param == "GT5_G15"|viralLoad_family$param == "GT5_G18",]
t.test(V3 ~ param, data = GT5_G15_G18)
wilcox.test(V3 ~ param, data = GT5_G15_G18)

#G0-10
GT5_G0_G10 <- viralLoad_family[viralLoad_family$param == "GT5_G0"|viralLoad_family$param == "GT5_G10",]
t.test(V3 ~ param, data = GT5_G0_G10)
wilcox.test(V3 ~ param, data = GT5_G0_G10)

#G0-15
GT5_G0_G15 <- viralLoad_family[viralLoad_family$param == "GT5_G0"|viralLoad_family$param == "GT5_G15",]
t.test(V3 ~ param, data = GT5_G0_G15)
wilcox.test(V3 ~ param, data = GT5_G0_G15)

#G0-18
GT5_G0_G18 <- viralLoad_family[viralLoad_family$param == "GT5_G0"|viralLoad_family$param == "GT5_G18",]
t.test(V3 ~ param, data = GT5_G0_G18)
wilcox.test(V3 ~ param, data = GT5_G0_G18)


#G5-15
GT5_G5_G15 <- viralLoad_family[viralLoad_family$param == "GT5_G5"|viralLoad_family$param == "GT5_G15",]
t.test(V3 ~ param, data = GT5_G5_G15)
wilcox.test(V3 ~ param, data = GT5_G5_G15)

#G5-18
GT5_G5_G18 <- viralLoad_family[viralLoad_family$param == "GT5_G5"|viralLoad_family$param == "GT5_G18",]
t.test(V3 ~ param, data = GT5_G5_G18)
wilcox.test(V3 ~ param, data = GT5_G5_G18)

#G10-18
GT5_G10_G18 <- viralLoad_family[viralLoad_family$param == "GT5_G10"|viralLoad_family$param == "GT5_G18",]
t.test(V3 ~ param, data = GT5_G10_G18)
wilcox.test(V3 ~ param, data = GT5_G10_G18)

###significance test GT9

GT9_G0_G5 <- viralLoad_family[viralLoad_family$param == "GT9_G0"|viralLoad_family$param == "GT9_G5",]
t.test(V3 ~ param, data = GT9_G0_G5)
wilcox.test(V3 ~ param, data = GT9_G0_G5)

GT9_G5_G10 <- viralLoad_family[viralLoad_family$param == "GT9_G5"|viralLoad_family$param == "GT9_G10",]
t.test(V3 ~ param, data = GT9_G5_G10)
wilcox.test(V3 ~ param, data = GT9_G5_G10)

GT9_G10_G15 <- viralLoad_family[viralLoad_family$param == "GT9_G10"|viralLoad_family$param == "GT9_G15",]
t.test(V3 ~ param, data = GT9_G10_G15)
wilcox.test(V3 ~ param, data = GT9_G10_G15)

GT9_G15_G20 <- viralLoad_family[viralLoad_family$param == "GT9_G15"|viralLoad_family$param == "GT9_G20",]
t.test(V3 ~ param, data = GT9_G15_G20)
wilcox.test(V3 ~ param, data = GT9_G15_G20)

GT9_G20_G24 <- viralLoad_family[viralLoad_family$param == "GT9_G20"|viralLoad_family$param == "GT9_G24",]
t.test(V3 ~ param, data = GT9_G20_G24)
wilcox.test(V3 ~ param, data = GT9_G20_G24)

GT9_G24_G28 <- viralLoad_family[viralLoad_family$param == "GT9_G24"|viralLoad_family$param == "GT9_G28",]
t.test(V3 ~ param, data = GT9_G24_G28)
wilcox.test(V3 ~ param, data = GT9_G24_G28)


#G0-10
GT9_G0_G10 <- viralLoad_family[viralLoad_family$param == "GT9_G0"|viralLoad_family$param == "GT9_G10",]
t.test(V3 ~ param, data = GT9_G0_G10)
wilcox.test(V3 ~ param, data = GT9_G0_G10)

#G0-15
GT9_G0_G15 <- viralLoad_family[viralLoad_family$param == "GT9_G0"|viralLoad_family$param == "GT9_G15",]
t.test(V3 ~ param, data = GT9_G0_G15)
wilcox.test(V3 ~ param, data = GT9_G0_G15)

#G0-20
GT9_G0_G20 <- viralLoad_family[viralLoad_family$param == "GT9_G0"|viralLoad_family$param == "GT9_G20",]
t.test(V3 ~ param, data = GT9_G0_G20)
wilcox.test(V3 ~ param, data = GT9_G0_G20)

#G0-24
GT9_G0_G24 <- viralLoad_family[viralLoad_family$param == "GT9_G0"|viralLoad_family$param == "GT9_G24",]
t.test(V3 ~ param, data = GT9_G0_G24)
wilcox.test(V3 ~ param, data = GT9_G0_G24)

#G0-28
GT9_G0_G28 <- viralLoad_family[viralLoad_family$param == "GT9_G0"|viralLoad_family$param == "GT9_G28",]
t.test(V3 ~ param, data = GT9_G0_G28)
wilcox.test(V3 ~ param, data = GT9_G0_G28)

#G5-15
GT9_G5_G15 <- viralLoad_family[viralLoad_family$param == "GT9_G5"|viralLoad_family$param == "GT9_G15",]
t.test(V3 ~ param, data = GT9_G5_G15)
wilcox.test(V3 ~ param, data = GT9_G5_G15)

#G5-20
GT9_G5_G20 <- viralLoad_family[viralLoad_family$param == "GT9_G5"|viralLoad_family$param == "GT9_G20",]
t.test(V3 ~ param, data = GT9_G5_G20)
wilcox.test(V3 ~ param, data = GT9_G5_G20)

#G5-24
GT9_G5_G24 <- viralLoad_family[viralLoad_family$param == "GT9_G5"|viralLoad_family$param == "GT9_G24",]
t.test(V3 ~ param, data = GT9_G5_G24)
wilcox.test(V3 ~ param, data = GT9_G5_G24)

#G5-28
GT9_G5_G28 <- viralLoad_family[viralLoad_family$param == "GT9_G5"|viralLoad_family$param == "GT9_G28",]
t.test(V3 ~ param, data = GT9_G5_G28)
wilcox.test(V3 ~ param, data = GT9_G5_G28)

#G10-20
GT9_G10_G20 <- viralLoad_family[viralLoad_family$param == "GT9_G10"|viralLoad_family$param == "GT9_G20",]
t.test(V3 ~ param, data = GT9_G10_G20)
wilcox.test(V3 ~ param, data = GT9_G10_G20)

#G10-24
GT9_G10_G24 <- viralLoad_family[viralLoad_family$param == "GT9_G10"|viralLoad_family$param == "GT9_G24",]
t.test(V3 ~ param, data = GT9_G10_G24)
wilcox.test(V3 ~ param, data = GT9_G10_G24)

#G10-28
GT9_G10_G28 <- viralLoad_family[viralLoad_family$param == "GT9_G10"|viralLoad_family$param == "GT9_G28",]
t.test(V3 ~ param, data = GT9_G10_G28)
wilcox.test(V3 ~ param, data = GT9_G10_G28)

#G15-24
GT9_G15_G24 <- viralLoad_family[viralLoad_family$param == "GT9_G15"|viralLoad_family$param == "GT9_G24",]
t.test(V3 ~ param, data = GT9_G15_G24)
wilcox.test(V3 ~ param, data = GT9_G15_G24)

#G15-28
GT9_G15_G28 <- viralLoad_family[viralLoad_family$param == "GT9_G15"|viralLoad_family$param == "GT9_G28",]
t.test(V3 ~ param, data = GT9_G15_G28)
wilcox.test(V3 ~ param, data = GT9_G15_G28)

#G20-28
GT9_G20_G28 <- viralLoad_family[viralLoad_family$param == "GT9_G20"|viralLoad_family$param == "GT9_G28",]
t.test(V3 ~ param, data = GT9_G20_G28)
wilcox.test(V3 ~ param, data = GT9_G20_G28)

###significance test G0
GT3_GT5_G0 <- viralLoad_family[viralLoad_family$param == "GT3_G0"|viralLoad_family$param == "GT5_G0",]
t.test(V3 ~ param, data = GT3_GT5_G0)
wilcox.test(V3 ~ param, data = GT3_GT5_G0)

GT5_GT9_G0 <- viralLoad_family[viralLoad_family$param == "GT5_G0"|viralLoad_family$param == "GT9_G0",]
t.test(V3 ~ param, data = GT5_GT9_G0)
wilcox.test(V3 ~ param, data = GT5_GT9_G0)

GT3_GT9_G0 <- viralLoad_family[viralLoad_family$param == "GT3_G0"|viralLoad_family$param == "GT9_G0",]
t.test(V3 ~ param, data = GT3_GT9_G0)
wilcox.test(V3 ~ param, data = GT3_GT9_G0)


###significance test G5
GT3_GT5_G5 <- viralLoad_family[viralLoad_family$param == "GT3_G5"|viralLoad_family$param == "GT5_G5",]
t.test(V3 ~ param, data = GT3_GT5_G5)
wilcox.test(V3 ~ param, data = GT3_GT5_G5)

GT5_GT9_G5 <- viralLoad_family[viralLoad_family$param == "GT5_G5"|viralLoad_family$param == "GT9_G5",]
t.test(V3 ~ param, data = GT5_GT9_G5)
wilcox.test(V3 ~ param, data = GT5_GT9_G5)

GT3_GT9_G5 <- viralLoad_family[viralLoad_family$param == "GT3_G5"|viralLoad_family$param == "GT9_G5",]
t.test(V3 ~ param, data = GT3_GT9_G5)
wilcox.test(V3 ~ param, data = GT3_GT9_G5)


###significance test G10
GT3_GT5_G10 <- viralLoad_family[viralLoad_family$param == "GT3_G10"|viralLoad_family$param == "GT5_G10",]
t.test(V3 ~ param, data = GT3_GT5_G10)
wilcox.test(V3 ~ param, data = GT3_GT5_G10)

GT5_GT9_G10 <- viralLoad_family[viralLoad_family$param == "GT5_G10"|viralLoad_family$param == "GT9_G10",]
t.test(V3 ~ param, data = GT5_GT9_G10)
wilcox.test(V3 ~ param, data = GT5_GT9_G10)

GT3_GT9_G10 <- viralLoad_family[viralLoad_family$param == "GT3_G10"|viralLoad_family$param == "GT9_G10",]
t.test(V3 ~ param, data = GT3_GT9_G10)
wilcox.test(V3 ~ param, data = GT3_GT9_G10)


###significance test G15
GT3_GT5_G15 <- viralLoad_family[viralLoad_family$param == "GT3_G13"|viralLoad_family$param == "GT5_G15",]
t.test(V3 ~ param, data = GT3_GT5_G15)
wilcox.test(V3 ~ param, data = GT3_GT5_G15)

GT5_GT9_G15 <- viralLoad_family[viralLoad_family$param == "GT5_G15"|viralLoad_family$param == "GT9_G15",]
t.test(V3 ~ param, data = GT5_GT9_G15)
wilcox.test(V3 ~ param, data = GT5_GT9_G15)

GT3_GT9_G15 <- viralLoad_family[viralLoad_family$param == "GT3_G13"|viralLoad_family$param == "GT9_G15",]
t.test(V3 ~ param, data = GT3_GT9_G15)
wilcox.test(V3 ~ param, data = GT3_GT9_G15)

###significance test G20
GT5_GT9_G20 <- viralLoad_family[viralLoad_family$param == "GT5_G18"|viralLoad_family$param == "GT9_G20",]
t.test(V3 ~ param, data = GT5_GT9_G20)
wilcox.test(V3 ~ param, data = GT5_GT9_G20)



