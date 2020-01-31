##########################
# IMPORTATION DE donnee_2018S #
##########################

donnee_2018 <- read.csv2("donnee_2018.csv", dec='.', sep=';')


###########################
# CHARGEMENT DES PACKAGES #
###########################
library(tidyverse) # imports ggplot2, dplyr, etc.
library(nlme)

# R²
pseudoR2 <- function(y, y_hat) {
  1 - (sum((y - y_hat)^2) / sum((y - mean(y))^2))
}

# Residuals
residues <- function(model, level = 0, computeResidues = TRUE, residues,
                     main = "") {
  library(e1071)
  if (computeResidues) {
    r = residuals(model, type="normalized", level)
  } else {
    r = residues
  }
  hist(r, freq=F, main = main, breaks=20)
  xfit<-seq(min(r),max(r),length=40)
  yfit<-dnorm(xfit)
  lines(xfit, yfit, col="red", lwd=2)
  print("Shapiro-Wilk:")
  print(shapiro.test(r))
  print(paste("Kurtosis:", kurtosis(r)))
  print(paste("Skewness:", skewness(r)))
}

##############
# PREPROCESS #
##############
# Weight difference
donnee_2018$Poids_final_g <- as.numeric(as.character(donnee_2018$Poids_final_g))
donnee_2018$Poids_Contenu..g. <- as.numeric(donnee_2018$Poids_Contenu..g.)

donnee_2018$PertePoids <- (donnee_2018$Poids_Contenu..g. - donnee_2018$Poids_final_g) / donnee_2018$Poids_Contenu..g.
glimpse(donnee_2018)
# Label repetitions in such a way they are unique per field and dose
donnee_2018$DoseRepPrel = paste0(donnee_2018$Dose, '_', donnee_2018$Repetition,'_', donnee_2018$Site, '_', 
                            donnee_2018$Prélèvement)
#donnee_2018_CNS$DoseRep = paste0(donnee_2018$Dose, '_', donnee_2018$Repetition)

# Cranberry residues are the reference to assess the effect
# of the kind of matter analysed
donnee_2018$TeaType = relevel(donnee_2018$TeaType, ref = 'Cranberry litter')
#donnee_2018$TeaType = relevel(donnee_2018$SorteThe, ref = 'ResidusCanneberge')
#donnee_2018_CNS$SorteThe = relevel(donnee_2018_CNS$SorteThe, ref = 'ResidusCanneberge')

######################
# MODELISATION MIXTE #
######################

# 1. What's affecting weight difference?
## prepare data
weightdiff_data_2018 = donnee_2018 %>%
  select(PertePoids,#C_over_N,
         TeaType, Farmers,Prélèvement, Site, DoseRepPrel) %>%
  drop_na(.) %>%
  droplevels(.)
weightdiff_data_factor_2018 <- weightdiff_data_2018
weightdiff_data_factor_2018$TeaType <- factor(weightdiff_data_2018$TeaType, 
                                         levels = levels(weightdiff_data_2018$TeaType)[c(2,3,1)]) # Reorder
ggplot(weightdiff_data_factor_2018, aes(x=Prélèvement, y=PertePoids)) +
  facet_grid(. ~ TeaType) +
  geom_smooth(method='lm') +
  labs(x = expression("Prélèvement (Semaines)"),
       y = "PertePoids (g)") +     #expression(frac(-ln(frac(Final~weight, Initial~weight)), 90))) +
  geom_point()
ggsave("Figure1.png", width = 10, height = 6, dpi = 600)# export plot high resolution

##############
library(ggplot2)
# Basic line plot with points
donnee_moy <- aggregate(PertePoids ~ Prélèvement +  TeaType, data = donnee_2018, FUN = mean)
str(donnee_moy)
ggplot(data=donnee_moy, aes(x=Prélèvement, y=PertePoids, group=TeaType)) +
  geom_line()+
  geom_point()
# Change the line type
ggplot(data=donnee_moy, aes(x=Prélèvement, y=PertePoids, group= TeaType)) +
  geom_line(linetype = "dashed")+
  geom_point()
# Change the color
ggplot(data=donnee_moy, aes(x=Prélèvement, y=PertePoids, group=TeaType)) +
  geom_line(aes(color=TeaType))+
  geom_point(aes(color=TeaType))+
geom_errorbar(aes(ymin=PertePoids-sd, ymax=PertePoids+sd), width=.1, 
              position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()
  






qplot(donnee_2018$Prélèvement,donnee_2018$PertePoids)+
geom_errorbar(aes(donnee_2018$Prélèvement, 
                    ymin=donnee_2018$PertePoids-sd, ymax=donnee_2018$PertePoids+sd), width=0.25)

greenTea <- subset(donnee_2018, TeaType == 'Green tea')
RooibosTea <- subset(donnee_2018, TeaType == 'Rooibos tea')
Canneberge <- subset(donnee_2018, TeaType == 'Cranberry litter')
boxplot(PertePoids ~ Prélèvement, greenTea )
boxplot(PertePoids ~ Prélèvement, RooibosTea )
boxplot(PertePoids ~ Prélèvement, Canneberge)

############
greenTea <- subset(donnee_2018, TeaType == 'Green tea')
greenTea_mean <- aggregate(PertePoids ~ Prélèvement, data = greenTea, FUN = mean)
qplot(greenTea_mean$Prélèvement,greenTea$PertePoids)+
  geom_errorbar(aes(greenTea_mean$Prélèvement, 
                    ymin=greenTea_mean$PertePoids - sd, ymax=greenTea_mean$PertePoids + sd), width=0.25)

#############
RooibosTea <- subset(donnee_2018, TeaType == 'Rooibos tea')
RooibosTea_mean <- aggregate(PertePoids ~ Prélèvement, data = RooibosTea, FUN = mean)
qplot(RooibosTea_mean$Prélèvement,RooibosTea_mean$PertePoids)+
  geom_errorbar(aes(RooibosTea_mean$Prélèvement, 
                    ymin=RooibosTea_mean$PertePoids - sd, ymax=RooibosTea_mean$PertePoids + sd), width=0.25)

###################
CranberryTea <- subset(donnee_2018, TeaType == 'Cranberry litter')
CranberryTea_mean <- aggregate(PertePoids ~ Prélèvement, data = CranberryTea, FUN = mean)
qplot(CranberryTea_mean$Prélèvement,CranberryTea_mean$PertePoids)
  geom_errorbar(aes(CranberryTea_mean$Prélèvement, 
                    ymin=CranberryTea_mean$PertePoids - sd, 
                    ymax=CranberryTea_mean$PertePoids + sd), width=0.25)



write.table(greenTea_mean, file = "greenTea.csv", 
            dec='.', sep=';')

# new plot (adjusts Yrange automatically)
with (
  data = d
  , expr = errbar(x, y, y+sd, y-sd, add=F, pch=1, cap=.015, log="x")
)
