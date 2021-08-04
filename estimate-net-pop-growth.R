#Analyse catch curve to predict abundance trend
# CJ Brown 2018-08-31
# not used in paper. 
#This looks to see if s simple lifehistory model 
#predictions are consistent
# with the observed rate of decline 

#TODO: Z mortality is pretty low so we don't actually predict any decline from fishing. 
#Could use a model like this to partition the causes of decline. 
#Z ~= M so F is estimated to be virtually zero here. What is wrong? 
# Perhaps life history model is wrong, M must be lower than it suggests. 
#It could be as low as 0.08 (looking at error bounds)
#If I can get reliable parameters I can actually estimate alpha
# from the other parameters and the magnitude of decline
# Then can I infer habitat loss somehow? 
#
#Thorson model probably underestimates m slightly., need to update it. 

library(dplyr)
library(ggplot2)
library(readr)
library(rethinking)
library(PlotTools)


dat <- read_csv("sites-abundance.csv")
maxn <- read_csv("max-catch-data.csv")
agefreq <- read_csv("topa-length-age-roviana.csv")
nyears <- 18
#
# Estimate pop growth prior 
#
library(nleqslv)
#I had to download fishlife package 
devtools::load_all("C:/Documents/R-packages-development/FishLife-master/Fishlife")


#Estimate Z 
af2 <- agefreq %>% group_by(age) %>%
  summarize(freq = n()) %>% ungroup()

imax <- which.max(af2$freq)
af3 <- af2[imax:nrow(af2),]
plot(af2$age, af2$freq)
lines(af3$age, af3$freq)
afm1 <- lm(log(freq) ~ age, data = af3)
spp_Z <- -coef(afm1)[2]
plot(af2$age, log(af2$freq))
abline(afm1)

# FSA::catchCurve(freq ~ age, data = af3) %>% summary()

#Estimate M
lhtraits <- Search_species(Genus="Bolbometopon",Species="muricatum",
                           add_ancestors=FALSE)$match_taxonomy
lhtraits <- Plot_taxa( lhtraits, mfrow=c(2,2) )

spp_m <-  exp(lhtraits[[1]]$Mean_pred[6]) # Instataneous natural mortality rate (ie not Z)
#Estimate R
#Functiont to estimate r (intrinsic popn growth) from age at maturity,
#instantaneous nat mortality and alpha (steepness). Use in nleqslv()
pop_form <- function(rm, amat, M, alpha){
  exp(rm * amat) - exp(rm*(amat - 1) - M) - alpha
}

#Estimate F 
spp_F <- spp_Z - spp_m

#Therefore based on life history invariants: 

#Estimate this from life history invariates M=~0.2 - just doing
#this visually
spp_amat <- 7 #Taylor median of 5-7 yrs for males and females 
ahat <- 4 #represents the number of spawners produced by each spawner per year, 
#if there was no fishing mortality, at very low spawner abundance Myers 2001
rval <- nleqslv(0.2, fn=pop_form, jac = NULL, 
                amat = spp_amat, M = spp_m, alpha = ahat)$x
rval

#Estimate Net growth 
#Now estimate population growth, accounting for F 
#dndt = rN - FN = (r-F)N
#Nt = N0 * exp((r-F)*t)
#deltaN = Nt/N0 = exp((r-F)*t)
#times decline = 
exp(-(rval - spp_F) * nyears)

#Picking an F to get a roviana decline over 18 years
exp(-(rval - 0.35) * nyears)

#BUT this doesn't account for hyperstability, because F will be increasing
# over time, so need to account for that here... 


