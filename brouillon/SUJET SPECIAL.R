library(psych)
var(bfi[,1:4], na.rm=T)->ex.var

donnee <- read.csv("TeaBag3_.csv",
                     header = TRUE, stringsAsFactors = TRUE,
                     dec=',', sep=';',  na.strings = "")
str(donnee)
