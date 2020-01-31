#IMPORTATION DE DONNEES--------

donnee <- read.csv("TeaBag3_.csv",
                    header = TRUE, stringsAsFactors = TRUE,
                    dec=',', sep=';',  na.strings = "")

donnee_CNS <- read.csv("CNS_Fortier_45.csv",
                       header = TRUE, stringsAsFactors = TRUE,
                       dec=',', sep=';',  na.strings = "")
# Rapport C/N initial des sacs de the
RESIDUCAN_C_sur_Ni <- 50.033/0.90172
HIBISCUS_C_sur_Ni <- 46.04066667/0.872726667
ROOIBOS_C_sur_Ni <- 47.79166667/1.040022222
GREENT_C_sur_Ni <- 47.99466667/3.377
SENCHA_C_sur_Ni <- 47.08725/3.6811


# CHARGEMENT DES PACKAGES----------
library(lme4)
library(lattice)
library(ggplot2)
library(lmerTest)

# MODELISATION MIXTE----------
## DIFFERENCE POIDS-----

donnee$DifferencePoids <- donnee$Poids_Restant_g - donnee$PoidsContenu_g
donnee_CNS$C_sur_N <- donnee_CNS$Teneur_C / donnee_CNS$Teneur_N
donnee_CNS$TauxDecomposition_MO <- donnee_CNS$Difference_C * 100/ 90
donnee$DoseRep = paste0(donnee$Dose, '_', donnee$Repetition,'_', donnee$Champs)
str(donnee)
donnee$SorteThe = relevel(donnee$SorteThe, ref = 'ResidusCanneberge')
lm.out <- lmer(DifferencePoids ~ SorteThe +  Regie + TypeAzote +
                 DoseN + (1 | Producteur/Champs/DoseRep), data = donnee)
summary(lm.out)
anova(lm.out)

#lm.out1 <- lmer(DifferencePoids ~ SorteThe *  Regie * TypeAzote * 
#                  DoseN + (1 | Producteur/Champs) + (1 | Repetition), 
#                data = donnee)
#summary(lm.out1)
#anova(lm.out1)
## TENEUR EN CARBONE---------
donnee_CNS$SorteThe = relevel(donnee_CNS$SorteThe, ref = 'ResidusCanneberge')
donnee_CNS$DoseRep = paste0(donnee_CNS$Dose, '_', donnee_CNS$Repetition)

lm3 <- lmer(TauxDecomposition_MO ~ SorteThe + DoseN + (1 | DoseRep), data = donnee_CNS)
summary(lm3)
anova(lm3)
library(dplyr)
library(tidyverse)
donnee_CNS$DoseRep = paste0(donnee_CNS$DoseN, '_', donnee_CNS$Repetition)

data_CNS_mm = donnee_CNS %>%
  select(TauxDecomposition_MO, SorteThe, DoseN, C_sur_N, DoseRep) %>%
  drop_na(.) %>%
  droplevels(.)

lm3 <- lmer(TauxDecomposition_MO ~ DoseN + C_sur_N + (1 | DoseRep), data = data_CNS_mm)
summary(lm3)
anova(lm3)

summary(lm(residuals(lm3)~data_CNS_mm$SorteThe))

boxplot(C_sur_N~SorteThe, donnee_CNS)

lm3 <- lmer(donnee_CNS$Teneur_C ~ SorteThe + DoseN + C_sur_N + (1 | Repetition), data = donnee_CNS)
summary(lm3)
anova(lm3)


# GRAPHIQUES-----------

boxplot1 <- boxplot(TauxDecomposition_MO ~ SorteThe, donnee_CNS,
                    ylab = "Taux de Décomposition de la Matière Organique(%)",
                   xlab = "Sortes de thé", cex.lab = 0.7, cex.axis = 0.59)

boxplot2 <- boxplot(donnee$Poids_Restant_g ~ SorteThe, donnee,
                    ylab = "Poids restant(g)",
                    xlab = "Sortes de thé", cex.lab = 0.7, cex.axis = 0.59)


qplot(x = DoseN, y = TauxDecomposition_MO , data = donnee_CNS, geom = c("point", "smooth"),
      color = SorteThe,method = "lm", shape = SorteThe,
      xlab = "Dose d'azote(kg/ha)", ylab = "Taux de décomposition(%)",facets = . ~ SorteThe, par(mar = c(12 ,12 ,12 ,12)))

qplot(x = donnee_CNS$C_sur_N, y = TauxDecomposition_MO , data = donnee_CNS, geom = c("point", "smooth"),
      color = SorteThe,method = "lm", shape = SorteThe,
      xlab = "Rapport C/N", ylab = "Taux de décomposition(%)",facets = . ~ SorteThe, par(mar = c(12 ,12 ,12 ,12)))

#data <- aggregate(DifferencePoids ~ DoseN + SorteThe, donnee, FUN = mean)


#data1 <- aggregate(Poids_Restant_g ~ DoseN + SorteThe, donnee, FUN = mean)

data_C <- aggregate(Teneur_C ~ DoseN + SorteThe, donnee_CNS, FUN = mean)

#data_N <- aggregate(Teneur_N ~ DoseN + SorteThe, donnee_CNS, FUN = mean)

#data_S <- aggregate(Teneur_S ~ DoseN + SorteThe, donnee_CNS, FUN = mean)

data2 <- aggregate(TauxDecomposition_MO ~ DoseN + SorteThe, donnee_CNS, FUN = mean)

qplot(x= TauxDecomposition_MO , fill = SorteThe,data = donnee_CNS, 
      xlab = "Taux de décomposition", ylab = "")

ggplot(donnee_CNS, aes(SorteThe,TauxDecomposition_MO, fill = SorteThe)) +
  geom_bar(stat = "summary", fun.y = "mean") +
             labs(x = "Sorte de thé", y = "Taux de décomposition(%)")


qplot(x = DoseN, y = Teneur_C , data = data_C, geom = c("point", "smooth"),
      color = SorteThe, shape = SorteThe)


data8 <- aggregate(DifferencePoids ~ Champs + SorteThe + DoseN + Repetition, donnee, FUN = mean)

qplot(x = DoseN, y = DifferencePoids , data = data8, geom = c("point", "smooth"),
      color = SorteThe, shape = SorteThe, facets = . ~ SorteThe)

qplot(x = DoseN, y = DifferencePoids , data = data8, geom = c("point", "smooth"),
      color = SorteThe, shape = SorteThe, facets = . ~ SorteThe, method = "lm")

qplot(x = DoseN, y = Poids_Restant_g , data = donnee, geom = c("point", "smooth"),
      color = SorteThe, shape = SorteThe, facets = . ~ SorteThe, method = "lm")



#data7 <- aggregate(C_sur_N ~ DoseN + SorteThe, donnee_CNS, FUN = mean)

#data3 <- aggregate(C_sur_N ~ DoseN + SorteThe, donnee_CNS, FUN = mean)

#data4 <- aggregate(Difference_C ~ DoseN + SorteThe, donnee_CNS, FUN = mean)

ggplot(donnee_CNS, aes(SorteThe,TauxDecomposition_MO, fill = SorteThe)) +
  geom_bar(stat = "summary", fun.y = "mean") +
  labs(x = "Sorte de thé", y = "Taux de décomposition(%)")

