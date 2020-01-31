##########################
# IMPORTATION DE DONNEES #
##########################
donnee <- read.csv("DONNEES1.1.csv",
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
donnee$PertePoids <- -log((donnee$PoidsContenu_g - donnee$Poids_Restant_g)/donnee$PoidsContenu_g)/90

# C/N ratio
#donnee$C_over_N <- donnee$Teneur_Ci / donnee$Teneur_Ni


# C-decomposition rate, as carbon differencedivided by the time of the experiment
donnee$Difference_C_Total_g <- donnee$Teneur_Ci - donnee$Teneur_Cf

donnee$TauxDecomposition_C <- - log(donnee$Teneur_Cf/donnee$Teneur_Ci)/90

# Label repetitions in such a way they are unique per field and dose
donnee$DoseRep = paste0(donnee$Dose, '_', donnee$Repetition,'_', donnee$Site)
#donnee_CNS$DoseRep = paste0(donnee$Dose, '_', donnee$Repetition)

# Cranberry residues are the reference to assess the effect
# of the kind of matter analysed
donnee$TeaType = relevel(donnee$TeaType, ref = 'Cranberry litter')
#donnee$TeaType = relevel(donnee$SorteThe, ref = 'ResidusCanneberge')
#donnee_CNS$SorteThe = relevel(donnee_CNS$SorteThe, ref = 'ResidusCanneberge')

######################
# MODELISATION MIXTE #
######################

# 1. What's affecting weight difference?
## prepare data
weightdiff_data = donnee %>%
  select(PertePoids,#C_over_N,
         TeaType, NitrogenType, DoseN, Farmers, Site, DoseRep) %>%
  drop_na(.) %>%
  droplevels(.)
weightdiff_data_factor <- weightdiff_data
weightdiff_data_factor$TeaType <- factor(weightdiff_data$TeaType, 
                                  levels = levels(weightdiff_data$TeaType)[c(4,3,5,2,1)]) # Reorder
ggplot(weightdiff_data_factor, aes(x=DoseN, y=PertePoids)) +
  facet_grid(. ~ TeaType) +
  geom_smooth(method='lm') +
  labs(x = expression("Nitrogen dose (kg ha"^"-1"*")"),
       y = "Decomposition rate (k1)") +     #expression(frac(-ln(frac(Final~weight, Initial~weight)), 90))) +
  geom_point()
ggsave("Figure1.png", width = 10, height = 6, dpi = 600)# export plot high resolution


################
#note
data(mtcars)
png("Figure1.png", width = 3000, height=1000, res=300)
plot(mtcars$wt, mtcars$hp)
dev.off() #export plot high resolution
##################
## mixed model
weightdiff_data$DoseN_scaled = scale(weightdiff_data$DoseN)

weightdiff_mm <- lme(fixed = PertePoids ~ TeaType + #C_over_N + 
                       NitrogenType + DoseN_scaled*TeaType,
                     random = ~ 1 | Farmers/Site, #/DoseRep,
                     data = weightdiff_data)
summary(weightdiff_mm)
anova(weightdiff_mm)
pseudoR2(y = weightdiff_data$PertePoids,
         y_hat = predict(weightdiff_mm, level = 0))
residues(weightdiff_mm)

weightdiff_gg = data.frame(intervals(weightdiff_mm, which = "fixed")$fixed)
weightdiff_gg$pvalue = summary(weightdiff_mm)$tTable[, 5]
weightdiff_gg$is_significant = ifelse(weightdiff_gg$pvalue <= 0.05,
                                      'Significant at 0.05 level',
                                      'Not significant at 0.05 level')
weightdiff_gg$variable = rownames(weightdiff_gg)
weightdiff_gg$variable[c(2,3,4,5,6,7,8,9,11,12,13,14)] <- c('Green tea', 'Hibiscus tea', 'Rooibos tea',
                                   'Sencha tea', 'N-6-1-1','N-8-0-0', 'N-SCU-39-0-0',
                                   'N-growers', 'Green tea : DoseN_scaled',
                                   'Hibiscus tea : DoseN_scaled', 'Rooibos tea : DoseN_scaled',
                                   'Sencha tea : DoseN_scaled')
weightdiff_gg$facet = factor(c('Intercept',
                               'Tea type', 'Tea type', 'Tea type', 'Tea type',
                               'N sources', 'N sources', 'N sources',
                               'N sources','Nitrogen input quantity',
                               'Tea type x N dosage', 'Tea type x N dosage', 
                               'Tea type x N dosage', 'Tea type x N dosage'))
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
ggsave("Figure2.png", width = 8, height = 5, dpi = 600)# export plot high resolution



# 2. What is affecting the organic matter decomposition rate?
omdecomp_data = donnee %>%
  select(TauxDecomposition_C, TeaType, DoseN, Farmers,Site, #C_over_N, 
         DoseRep,NitrogenType) %>%
  drop_na(.) %>%
  droplevels(.)
omdecomp_data_factor <- omdecomp_data
omdecomp_data_factor$TeaType <- factor(omdecomp_data$TeaType, 
                                         levels = levels(omdecomp_data$TeaType)[c(4,3,5,2,1)]) # Reorder

ggplot(omdecomp_data_factor, aes(x=DoseN, y=TauxDecomposition_C)) +
  facet_grid(. ~ TeaType) +
  geom_smooth(method='lm') +
  labs(x = expression("Nitrogen dose (kg ha"^"-1"*")"),
       y = "Carbon decomposition rate (k2)")+      #expression( frac(Initial~carbon - Final~carbon, Initial~carbon) %*% 100)) +
  geom_point(aes())
ggsave("Figure3.png", width = 10, height = 6, dpi = 600)# export plot high resolution


omdecomp_data$DoseN_sc = scale(omdecomp_data$DoseN)
#omdecomp_data$C_over_N_sc = scale(omdecomp_data$C_over_N)

plot(density(omdecomp_data$TauxDecomposition_C))
omdecomp_mm = lme(fixed = TauxDecomposition_C ~ TeaType + #C_over_N+  
                    NitrogenType + DoseN_sc*TeaType,
                  random = ~ 1 | Farmers/Site, #/DoseRep,
                  data = omdecomp_data) # [omdecomp_data$TauxDecomposition_C > 0, ]
summary(omdecomp_mm)
ranef(omdecomp_mm)
summary(omdecomp_mm)$tTable
anova(omdecomp_mm)
pseudoR2(y=omdecomp_data$TauxDecomposition_C, # [omdecomp_data$TauxDecomposition_C > 0]
         y_hat = predict(omdecomp_mm, level=2))
residues(omdecomp_mm)



omdecomp_gg = data.frame(intervals(omdecomp_mm, which = "fixed")$fixed)
omdecomp_gg$pvalue = summary(omdecomp_mm)$tTable[, 5]
omdecomp_gg$is_significant = ifelse(omdecomp_gg$pvalue <= 0.05,
                                    'Significant at 0.05 level',
                                    'Not significant at 0.05 level')
omdecomp_gg$variable = rownames(omdecomp_gg)
omdecomp_gg$variable[c(2,3,4,5,6,7,8,9,10,11,12,13,14)] <- c('Green tea', 'Hibiscus tea', 'Rooibos tea',
                                                          'Sencha tea', 'N-6-1-1','N-8-0-0', 'N-SCU-39-0-0',
                                                          'N-growers', 'DoseN_scaled', 'Green tea : DoseN_scaled',
                                                          'Hibiscus tea : DoseN_scaled', 'Rooibos tea : DoseN_scaled',
                                                          'Sencha tea : DoseN_scaled')
omdecomp_gg$facet = factor(c('Intercept',
                             'Tea type', 'Tea type', 'Tea type', 'Tea type',
                             'N source', 'N source', 'N source',
                             'N source','Nitrogen input quantity',
                             'Tea type x N dosage', 'Tea type x N dosage', 'Tea type x N dosage', 'Tea type x N dosage'))
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
ggsave("Figure4.png", width = 8, height = 5, dpi = 600)# export plot high resolution

write.table(donnee, file = "donneeEXP.csv", 
            dec=',', sep=';')

#-----------------------------------------------------------------------
  #Tea bag index
head(donnee)

donnee$M_surM0 <- donnee$Poids_Restant_g / donnee$PoidsContenu_g
Hf_g = 0.805
Hf_r = 0.826


donnee$a = NA
donnee$k_tbi = NA

donnee$a[donnee$TeaType == "Green tea"] = (donnee$PoidsContenu_g[donnee$TeaType == "Green tea"] -
                                            donnee$Poids_Restant_g[donnee$TeaType == "Green tea"]) /
                                            donnee$PoidsContenu_g[donnee$TeaType == "Green tea"]

donnee$S = 1 - donnee$a/Hf_g
donnee$k_tbi[donnee$TeaType == "Rooibos tea"] = log((Hf_r/donnee$S[donnee$TeaType == "Green tea"]) / 
                                                      (donnee$M_surM0[donnee$TeaType == "Rooibos tea"]- 
                                                         (1-(Hf_r/donnee$S[donnee$TeaType == "Green tea"])))) / 90
sdr = paste(donnee$Site, donnee$DoseN, donnee$Repetition)
sdr_u = unique(sdr)

for (i in 1:length(sdr_u)) {
  Sg_i = mean(donnee$S[sdr == sdr_u[i] & donnee$TeaType == "Green tea"], na.rm = TRUE)
  
  # Rooibos
  donnee$a[sdr == sdr_u[i] & donnee$TeaType == "Rooibos tea"] = Hf_r * (1 - Sg_i)
  ar_i = donnee$a[sdr == sdr_u[i] & donnee$TeaType == "Rooibos tea"]
  M_surM0_i = donnee$M_surM0[sdr == sdr_u[i] & donnee$TeaType == "Rooibos tea"]
  
  donnee$k_tbi[sdr == sdr_u[i] & donnee$TeaType == "Rooibos tea"] = log(ar_i / (M_surM0_i - (1-ar_i))) / 90
  
  
}
donnee$k_tbi


PertePoids_rooibos <- subset(donnee, donnee$TeaType == "Rooibos tea")

plot(PertePoids_rooibos$PertePoids, PertePoids_rooibos$k_tbi)

TBI_data = donnee %>%
  select(k_tbi, TeaType, DoseN, Farmers,Site, #C_over_N, 
         NitrogenType, Regie) %>%
  drop_na(.) %>%
  droplevels(.)

TBI_rooibos <- subset(TBI_data, TeaType == "Rooibos tea")
lme_TBI <- lme(fixed = k_tbi ~  DoseN + NitrogenType,
    random = ~ 1 | Farmers/Site, #/DoseRep,
    data = TBI_rooibos) # [omdecomp_data$TauxDecomposition_C > 0, ]

summary(lme_TBI)
ranef(lme_TBI)
summary(lme_TBI)$tTable
attributes(summary(lme_TBI))
#intervals(lme_TBI)
anova(lme_TBI)
pseudoR2(y=TBI_rooibos$k_tbi, # [omdecomp_data$TauxDecomposition_C > 0]
         y_hat = predict(lme_TBI, level=2))
residues(lme_TBI)



donnee$M_surM0 <- donnee$Poids_Restant_g / donnee$PoidsContenu_g
Hf_g = 0.805
Hf_cr = 0.475


donnee$a = NA
donnee$k1_tbi = NA

donnee$a[donnee$TeaType == "Green tea"] = (donnee$PoidsContenu_g[donnee$TeaType == "Green tea"] -
                                             donnee$Poids_Restant_g[donnee$TeaType == "Green tea"]) /
  donnee$PoidsContenu_g[donnee$TeaType == "Green tea"]

donnee$S = 1 - donnee$a/Hf_g
donnee$k1_tbi[donnee$TeaType == "Cranberry litter"] = log((Hf_cr/donnee$S[donnee$TeaType == "Green tea"]) / 
                                                            (donnee$M_surM0[donnee$TeaType == "Cranberry litter"]- 
                                                               (1-(Hf_cr/donnee$S[donnee$TeaType == "Green tea"])))) / 90
sdr = paste(donnee$Site, donnee$DoseN, donnee$Repetition)
sdr_u = unique(sdr)

for (i in 1:length(sdr_u)) {
  Sg_i = mean(donnee$S[sdr == sdr_u[i] & donnee$TeaType == "Green tea"], na.rm = TRUE)
  
  # Rooibos
  donnee$a[sdr == sdr_u[i] & donnee$TeaType == "Cranberry litter"] = Hf_r * (1 - Sg_i)
  acr_i = donnee$a[sdr == sdr_u[i] & donnee$TeaType == "Cranberry litter"]
  M_surM0_i = donnee$M_surM0[sdr == sdr_u[i] & donnee$TeaType == "Cranberry litter"]
  
  donnee$k1_tbi[sdr == sdr_u[i] & donnee$TeaType == "Cranberry litter"] = log(acr_i / (M_surM0_i - (1-acr_i))) / 90
  
  
}
donnee$k1_tbi
TBI_cranberry <- subset(donnee, TeaType == "Cranberry litter")


mean_k1_tbi <- mean(TBI_cranberry$k1_tbi, na.rm = TRUE)






weightdiff_gg = data.frame(intervals(weightdiff_mm, which = "fixed")$fixed)
weightdiff_gg$pvalue = summary(weightdiff_mm)$tTable[, 5]
weightdiff_gg$is_significant = ifelse(weightdiff_gg$pvalue <= 0.05,
                                      'Significant at 0.05 level',
                                      'Not significant at 0.05 level')
weightdiff_gg$variable = rownames(weightdiff_gg)
weightdiff_gg$variable[c(2,3,4,5,6,7,8,9,11,12,13,14)] <- c('Green tea', 'Hibiscus tea', 'Rooibos tea',
                                                            'Sencha tea', 'N-6-1-1','N-8-0-0', 'N-SCU-39-0-0',
                                                            'N-growers', 'Green tea : DoseN_scaled',
                                                            'Hibiscus tea : DoseN_scaled', 'Rooibos tea : DoseN_scaled',
                                                            'Sencha tea : DoseN_scaled')
weightdiff_gg$facet = factor(c('Intercept',
                               'Tea type', 'Tea type', 'Tea type', 'Tea type',
                               'N sources', 'N sources', 'N sources',
                               'N sources','Nitrogen input quantity',
                               'Tea type x N dosage', 'Tea type x N dosage', 
                               'Tea type x N dosage', 'Tea type x N dosage'))
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
ggsave("Figure2.png", width = 8, height = 5, dpi = 600)# export plot high resolution



# 2. What is affecting the organic matter decomposition rate?
omdecomp_data = donnee %>%
  select(TauxDecomposition_C, TeaType, DoseN, Farmers,Site, #C_over_N, 
         DoseRep,NitrogenType) %>%
  drop_na(.) %>%
  droplevels(.)
omdecomp_data_factor <- omdecomp_data
omdecomp_data_factor$TeaType <- factor(omdecomp_data$TeaType, 
                                       levels = levels(omdecomp_data$TeaType)[c(4,3,5,2,1)]) # Reorder

ggplot(omdecomp_data_factor, aes(x=DoseN, y=TauxDecomposition_C)) +
  facet_grid(. ~ TeaType) +
  geom_smooth(method='lm') +
  labs(x = expression("Nitrogen dose (kg ha"^"-1"*")"),
       y = "Carbon decomposition rate (k2)")+      #expression( frac(Initial~carbon - Final~carbon, Initial~carbon) %*% 100)) +
  geom_point(aes())
ggsave("Figure3.png", width = 10, height = 6, dpi = 600)# export plot high resolution




S_data = donnee %>%
  select(S, TeaType, DoseN, Farmers,Site, #C_over_N, 
         NitrogenType) %>%
  drop_na(.) %>%
  droplevels(.)
#TBI_data_factor <- TBI_data
#TBI_data_factor$TeaType <- factor(TBI_data$TeaType, 
#                                      levels = levels(omdecomp_data$TeaType)[c(4,3,5,2,1)]) # Reorder
TBI_GreenTea <- subset(S_data, TeaType == "Green tea")
lme_S <- lme(fixed = S ~  DoseN + NitrogenType,
               random = ~ 1 | Farmers/Site, #/DoseRep,
               data = TBI_GreenTea) # [omdecomp_data$TauxDecomposition_C > 0, ]

summary(lme_S)
ranef(lme_S)
summary(lme_S)$tTable
attributes(summary(lme_S))
#intervals(lme_S)
anova(lme_S)
pseudoR2(y=TBI_GreenTea$S, # [omdecomp_data$TauxDecomposition_C > 0]
         y_hat = predict(lme_S, level=2))
residues(lme_S)





omdecomp_mm = lme(fixed = TauxDecomposition_C ~ TeaType + #C_over_N+  
                    NitrogenType + DoseN_sc*TeaType,
                  random = ~ 1 | Farmers/Site, #/DoseRep,
                  data = omdecomp_data) # [omdecomp_data$TauxDecomposition_C > 0, ]
summary(omdecomp_mm)
ranef(omdecomp_mm)
summary(omdecomp_mm)$tTable
anova(omdecomp_mm)
pseudoR2(y=omdecomp_data$TauxDecomposition_C, # [omdecomp_data$TauxDecomposition_C > 0]
         y_hat = predict(omdecomp_mm, level=2))
residues(omdecomp_mm)





donnee_extrait <- subset(donnee, DoseN == 45)
boxplot(TauxDecomposition_C ~ TeaType + DoseN, donnee, ylab="Carbon decomposition rate k (%)",
        las=1)
ggsave("Figure10.png", width = 10, height = 8, dpi = 600)# export plot high resolution

boxplot(donnee$PertePoids ~ TeaType + DoseN, donnee, ylab="Weight loss (%)",
        las=1)
ggsave("Figure6.png", width = 10, height = 8, dpi = 600)# export plot high resolution


boxplot(donnee_extrait$PertePoids ~ TeaType, donnee_extrait )

boxplot(TauxDecomposition_C ~ TeaType + DoseN, donnee_extrait )

--------------------------------------------------------
  #Carte
install.packages("tidyverse") 
library("tidyverse")
library("ggmap")
library("ggrepel")
library("cowplot")              #   c("Site 45", "Site 9", "Site A9", "Site 10")

sites <- tibble(sites = c("", "", "", ""),
                lat = c(46.2759555, 46.2759555, 46.2329064, 46.3277404),
                lon = c(-71.8584291, -71.8584291, -72.0461854, -71.7486645))
site <- tibble(sites = c("Site 45", "Site 9", "Site A9", "Site 10"),
                lat = c(46.2759555, 46.2759555, 46.2329064, 46.3277404),
                lon = c(-71.8584291, -71.8584291, -72.0461854, -71.7486645))

frame_map_zoom_out <- c(left = -90, bottom = 44, right = -53, top = 63)
frame_map_zoom_in <- c(left = -73, bottom = 46, right = -71, top = 47)

rect_zoom <- data.frame(t(data.frame(frame_map_zoom_in)))

# Zoom out
fetch_map_out <- get_map(frame_map_zoom_out, zoom = 4, source = "stamen",
                         maptype = "toner-2011")
map_out <- ggmap(fetch_map_out) +
  geom_rect(data = rect_zoom, aes(xmin = left, xmax = right,
                                  ymin = bottom, ymax = top,
                                  x = NULL, y = NULL),
            fill = rgb(0, 0, 0, 0.2), colour = "black") +
  geom_label_repel(data = sites, aes(x = lon, y = lat, label = sites)) +
  xlab("Longitude") +
  ylab("Latitude")
map_out
ggsave("site_map.png", width = 12, height = 5, dpi = 300)
# Zoom in
fetch_map_in <- get_map(frame_map_zoom_in, zoom = 10, source = "stamen",
                        maptype = "toner-2011")
map_in <- ggmap(fetch_map_in) +
  geom_label_repel(data = site, aes(x = lon, y = lat, label = sites)) +
  xlab("Longitude") +
  ylab("Latitude")

plot_grid(map_out, map_in, labels = c("A", "B"))
ggsave("site_map.png", width = 12, height = 5, dpi = 300)




--------------------------
  
  donnee$  
  
  
  
  
donnee$Pourcentage_pertePoids <- donnee$DifferencePoids * 100 / donnee$PoidsContenu_g
donnee$Pourcentage_C_perdu <- donnee$Difference_C_Total_g * 100 / donnee$Teneur_Ci

aggregate(PertePoids ~ TeaType + DoseN, data = donnee, FUN = min)
aggregate(PertePoids ~ TeaType + DoseN, data = donnee, FUN = median)
aggregate(PertePoids ~ TeaType + DoseN, data = donnee, FUN = max)
aggregate(PertePoids ~ TeaType + DoseN, data = donnee, FUN = sd)
aggregate(PertePoids ~ TeaType + DoseN, data = donnee, FUN = mean)
aggregate(PertePoids ~ TeaType, data = donnee, FUN = mean)
aggregate(PertePoids ~ TeaType, data = donnee, FUN = sd)

aggregate(donnee$PoidsContenu_g ~ TeaType, data = donnee, FUN = mean)

aggregate(donnee$Poids_Restant_g ~ TeaType, data = donnee, FUN = mean)

aggregate(TauxDecomposition_C ~ TeaType+DoseN, data = donnee, FUN = mean)

aggregate(TauxDecomposition_C ~ TeaType+DoseN, data = donnee, FUN = sd)



aggregate(TauxDecomposition_C ~ TeaType + DoseN, data = donnee, FUN = min)
aggregate(TauxDecomposition_C ~ TeaType + DoseN, data = donnee, FUN = median)
aggregate(TauxDecomposition_C ~ TeaType + DoseN, data = donnee, FUN = max)
aggregate(TauxDecomposition_C ~ TeaType + DoseN, data = donnee, FUN = sd)
aggregate(TauxDecomposition_C ~ TeaType + DoseN, data = donnee, FUN = mean)
aggregate(TauxDecomposition_C ~ TeaType, data = donnee, FUN = mean)
aggregate(TauxDecomposition_C ~ TeaType, data = donnee, FUN = sd)

library(tidyverse)
library(weathercan)

head(stations)
stations_search("Laurierville", interval = "month")
weather_dl(station_ids = 5392, start = "2007-05-15", end = "2017-08-15")

