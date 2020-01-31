##########################
# IMPORTATION DE DONNEES #
##########################
donnee <- read.csv("DONNEES.csv",
                    header = TRUE, stringsAsFactors = TRUE,
                    dec=',', sep=';',  na.strings = "")

#donnee_CNS <- read.csv("CNS_Fortier_45.csv",
 #                      header = TRUE, stringsAsFactors = TRUE,
  #                     dec=',', sep=';',  na.strings = "")

donnee$Teneur_Cf <- donnee$Pourcentage_C * donnee$Poids_Restant_g / 100
donnee$Teneur_Nf <- donnee$Pourcentage_N * donnee$Poids_Restant_g / 100
donnee$Teneur_Sf <- donnee$Pourcentage_S * donnee$Poids_Restant_g / 100  
  
#donnee$Teneur_Ci <- apply(donnee$PoidsContenu_g , MARGIN = 1, FUN = mean)
#donnee$Teneur_Ci <- subset(donnee, subset =SorteThe , select =PoidsContenu_g) * 47.99466667
#aggregate(donnee[, c("conc", "rate")], by = list(state = donnee$SorteThe), FUN = mean)

###########################
# CHARGEMENT DES PACKAGES #
###########################
library(tidyverse) # imports ggplot2, dplyr, etc.
library(nlme)

pseudoR2 = function(data, mm, y ,level=0) {
  R2 = cor(y, predict(mm))^2
  R2.1 = 1 - with(data, (sum((y-predict(mm, level=level))^2)/sum((y-mean(y))^2)))
  return(R2.1)
}

##############
# PREPROCESS #
##############
# Weight difference
donnee$DifferencePoids <- donnee$PoidsContenu_g - donnee$Poids_Restant_g

# C/N ratio
donnee$C_sur_N <- donnee$Teneur_Cf / donnee$Teneur_N


# C-decomposition rate, as carbon differencedivided by the time of the experiment
donnee$Difference_C_Total_g <- donnee$Teneur_Ci - donnee$Teneur_Cf

donnee$TauxDecomposition_MO <- donnee$Difference_C_Total_g * 100 / 90 # units: g/day

# Label repetitions in such a way they are unique per field and dose
donnee$DoseRep = paste0(donnee$Dose, '_', donnee$Repetition,'_', donnee$Champs)
#donnee_CNS$DoseRep = paste0(donnee$Dose, '_', donnee$Repetition)

# Cranberry residues are the reference to assess the effect
# of the kind of matter analysed
donnee$TeaType = relevel(donnee$TeaType, ref = 'Control / Cranberry Litter')
#donnee_CNS$SorteThe = relevel(donnee_CNS$SorteThe, ref = 'ResidusCanneberge')

######################
# MODELISATION MIXTE #
######################

# 1. What's affecting weight difference?
## prepare data
weightdiff_data = donnee %>%
  select(DifferencePoids, TeaType, System, NitrogenType, DoseN, Producteur, Champs, DoseRep) %>%
  drop_na(.) %>%
  droplevels(.)

ggplot(weightdiff_data, aes(x=DoseN, y=DifferencePoids)) +
  facet_grid(. ~ TeaType) +
  geom_smooth(method='lm') +
  labs(x = "Nitrogen doses (kg/ha)", y = "weight loss (g)") +
  geom_point()

## mixed model
weightdiff_data$DoseN_sc = scale(weightdiff_data$DoseN)

weightdiff_mm <- lme(fixed = DifferencePoids ~ TeaType +  System + NitrogenType + DoseN_sc*TeaType,
                     random = ~ 1 | Producteur/Champs/DoseRep,
                     data = weightdiff_data)
summary(weightdiff_mm)
anova(weightdiff_mm)
pseudoR2(data=weightdiff_data, mm=weightdiff_mm, y=weightdiff_data$DifferencePoids, level=0)

weightdiff_gg = data.frame(intervals(weightdiff_mm)$fixed)
weightdiff_gg$pvalue = summary(weightdiff_mm)$tTable[, 5]
weightdiff_gg$is_significant = ifelse(weightdiff_gg$pvalue <= 0.05,
                                      'Significant at 0.05 level',
                                      'Not significant at 0.05 level')
weightdiff_gg$variable = rownames(weightdiff_gg)
weightdiff_gg$facet = factor(c('Intercept',
                                'Tea', 'Tea', 'Tea', 'Tea',
                                'Culture',
                                'Nitrogen input type', 'Nitrogen input type', 'Nitrogen input type',
                                'Nitrogen input quantity',
                               'Interaction', 'Interaction', 'Interaction', 'Interaction'))
weightdiff_gg = weightdiff_gg[-1, ] # remove the intercept
ggplot(weightdiff_gg, aes(x = est., y=variable, colour=is_significant)) +
  facet_grid(facet ~ ., scales = 'free', space = 'free') +
  geom_point() +
  geom_segment(aes(x=lower, xend=upper, y=variable, yend=variable)) +
  geom_vline(xintercept = 0, colour='grey70') +
  xlab('Coefficient') +
  ylab('') +
  theme(strip.text.y = element_text(angle=0),
        legend.title = element_blank(),
        legend.position = "bottom")


# 2. What is affecting the organic matter decomposition rate?
omdecomp_data = donnee %>%
  select(TauxDecomposition_MO, TeaType, DoseN, Producteur,Champs, C_sur_N, DoseRep, System, NitrogenType) %>%
  drop_na(.) %>%
  droplevels(.)

ggplot(omdecomp_data, aes(x=DoseN, y=TauxDecomposition_MO)) +
  facet_grid(. ~ TeaType) +
  geom_smooth(method='lm') +
  labs(x = "Nitrogen doses (kg/ha)", y = "Decomposition rate (%)") +
  
  geom_point()

omdecomp_data$DoseN_sc = scale(omdecomp_data$DoseN)
omdecomp_data$C_sur_N_sc = scale(omdecomp_data$C_sur_N)
omdecomp_mm = lme(fixed = TauxDecomposition_MO ~ TeaType +  System + NitrogenType + DoseN_sc*TeaType,
                    random = ~ 1 | Producteur/Champs/DoseRep,
                    data = omdecomp_data)
summary(omdecomp_mm)
anova(omdecomp_mm)
pseudoR2(data=omdecomp_data, mm=omdecomp_mm, y=omdecomp_data$TauxDecomposition_MO, level=0)

omdecomp_gg = data.frame(intervals(omdecomp_mm)$fixed)
omdecomp_gg$pvalue = summary(omdecomp_mm)$tTable[, 5]
omdecomp_gg$is_significant = ifelse(omdecomp_gg$pvalue <= 0.05,
                                    'Significant at 0.05 level',
                                    'Not significant at 0.05 level')
omdecomp_gg$variable = rownames(omdecomp_gg)
omdecomp_gg$facet = factor(c('Intercept',
                             'Tea', 'Tea', 'Tea', 'Tea', 'Culture',
                             'Nitrogen input type', 'Nitrogen input type', 'Nitrogen input type',
                             'Nitrogen input quantity',
                             'Interaction', 'Interaction', 'Interaction', 'Interaction'))
omdecomp_gg = omdecomp_gg[-1, ] # remove the intercept

droplevels(omdecomp_gg)
ggplot(omdecomp_gg, aes(x = est., y=variable, colour=is_significant)) +
  facet_grid(facet ~ ., scales = 'free', space = 'free') +
  geom_point() +
  geom_segment(aes(x=lower, xend=upper, y=variable, yend=variable)) +
  geom_vline(xintercept = 0, colour='grey70') +
  xlab('Coefficient') +
  ylab('') +
  theme(strip.text.y = element_text(angle=0),
        legend.title = element_blank(),
        legend.position = "bottom")
------------------------------------------------------------------------------
## The problem with the former model is that C_sur_N and SorteThe are colinear
summary(lm(C_sur_N~TeaType, omdecomp_data))
ggplot(omdecomp_data, aes(y = C_sur_N, x=TeaType)) +
  geom_boxplot() +
  xlab('Organic matter type') +
  ylab(' C/N ratio')

## We can select one or the other. We select SorteThe for convenience
omdecomp_mm = lme(fixed = TauxDecomposition_MO ~ TeaType +  System + NitrogenType + DoseN_sc*TeaType,
                  random = ~ 1 | Producteur/Champs/DoseRep,
                  data = omdecomp_data)
summary(omdecomp_mm)
anova(omdecomp_mm)
pseudoR2(data=omdecomp_data, mm=omdecomp_mm, y=omdecomp_data$TauxDecomposition_MO, level=0)

omdecomp_gg = data.frame(intervals(omdecomp_mm, which = "fixed")$fixed)
omdecomp_gg$pvalue = summary(omdecomp_mm)$tTable[, 5]
omdecomp_gg$is_significant = ifelse(omdecomp_gg$pvalue <= 0.05,
                                    'Significant at 0.05 level',
                                    'Not significant at 0.05 level')
omdecomp_gg$variable = rownames(omdecomp_gg)
omdecomp_gg$facet = factor(c('Intercept',
                             'Tea', 'Tea', 'Tea', 'Tea',
                             'Nitrogen input quantity',
                             'Interaction', 'Interaction', 'Interaction', 'Interaction'))
omdecomp_gg = omdecomp_gg[-1, ] # remove the intercept
ggplot(omdecomp_gg, aes(x = est., y=variable, colour=is_significant)) +
  facet_grid(facet ~ ., scales = 'free', space = 'free') +
  geom_point() +
  geom_segment(aes(x=lower, xend=upper, y=variable, yend=variable)) +
  geom_vline(xintercept = 0, colour='grey70') +
  xlab('Coefficient') +
  ylab('') +
  theme(strip.text.y = element_text(angle=0),
        legend.title = element_blank(),
        legend.position = "bottom")

## we can check if residuals could still be predicted using the C/N ratio
summary(lm(residuals(omdecomp_mm)~omdecomp_data$C_sur_N)) # ... no, not really


## Now we can try to remove tea type and keep C/N ratio
omdecomp_mm = lme(fixed = TauxDecomposition_MO ~ TeaType +  System + NitrogenType + DoseN_sc*TeaType,
                  random = ~ 1 | DoseRep,
                  data = omdecomp_data)
summary(omdecomp_mm)
anova(omdecomp_mm)
pseudoR2(data=omdecomp_data, mm=omdecomp_mm, y=omdecomp_data$TauxDecomposition_MO, level=0)

omdecomp_gg = data.frame(intervals(omdecomp_mm, which = "fixed")$fixed)
omdecomp_gg$pvalue = summary(omdecomp_mm)$tTable[, 5]
omdecomp_gg$is_significant = ifelse(omdecomp_gg$pvalue <= 0.05,
                                    'Significant at 0.05 level',
                                    'Not significant at 0.05 level')
omdecomp_gg$variable = rownames(omdecomp_gg)
omdecomp_gg$facet = factor(c('Intercept',
                             'Tea',
                             'Nitrogen input quantity',
                             'Interaction'))
omdecomp_gg = omdecomp_gg[-1, ] # remove the intercept
ggplot(omdecomp_gg, aes(x = est., y=variable, colour=is_significant)) +
  facet_grid(facet ~ ., scales = 'free', space = 'free') +
  geom_point() +
  geom_segment(aes(x=lower, xend=upper, y=variable, yend=variable)) +
  geom_vline(xintercept = 0, colour='grey70') +
  xlab('Coefficient') +
  ylab('') +
  theme(strip.text.y = element_text(angle=0),
        legend.title = element_blank(),
        legend.position = "bottom")

donnee$pourcentage_pertePoids <- donnee$DifferencePoids * 100 / donnee$PoidsContenu_g
donnee$pourcentage_C_perdu <- donnee$Difference_C_Total_g * 100 / donnee$Teneur_Ci

aggregate(pourcentage_pertePoids ~ SorteThe + DoseN, data = donnee, FUN = min)
aggregate(pourcentage_pertePoids ~ SorteThe + DoseN, data = donnee, FUN = median)
aggregate(pourcentage_pertePoids ~ SorteThe + DoseN, data = donnee, FUN = max)
aggregate(pourcentage_pertePoids ~ SorteThe + DoseN, data = donnee, FUN = sd)



aggregate(pourcentage_C_perdu ~ SorteThe + DoseN, data = donnee, FUN = min)
aggregate(pourcentage_C_perdu ~ SorteThe + DoseN, data = donnee, FUN = median)
aggregate(pourcentage_C_perdu ~ SorteThe + DoseN, data = donnee, FUN = max)
aggregate(pourcentage_C_perdu ~ SorteThe + DoseN, data = donnee, FUN = sd)




aggregate(donnee$TauxDecomposition_MO ~ SorteThe + DoseN, data = donnee, FUN = mean)

aggregate(donnee[, c("SorteThe", "Champs")], by = list(donnee$Poids_Restant_g),
          FUN = mean)
