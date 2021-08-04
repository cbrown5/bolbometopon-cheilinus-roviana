#Analyse abundance for Roviana data for habili
# CJ Brown 2018-11-13
# v4 stan model 
#Adults only 

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(rethinking)
library(PlotTools)


datin <- read_csv("roviana-habili-2018-10-19.csv")
hab_mat <- 550


#add in zeros
#add in zeros
datin$number[is.na(datin$number)] <- 0
sitevars <- datin %>% select(site, date, month, year, lunar, area) %>%
  distinct()
dat <- datin %>% group_by(site, date) %>%
  filter(size > hab_mat) %>%
  summarize(number = sum(number)) %>%
  ungroup() %>%
  full_join(sitevars) 
dat$number[is.na(dat$number)] <- 0

#check no duplicates
which(duplicated(paste(dat$site, dat$date)))

#
# Main data only 
#
sites <- unique(dat$site)[5:8]
dat2 <- filter(dat, site %in% sites)
dat2$jsite <- as.integer(factor(dat2$site))
dat2$period <- factor(ifelse(dat2$year < 2002, "2000s", "2018"))
dat2$xyear <- as.integer(dat2$period)-1

#
# Data for models 
#


datlist_habili <- list(
  N = nrow(dat2), 
  N_site = length(unique(dat2$jsite)),
  abund = dat2$number, 
  jsite = dat2$jsite, 
  xyear = dat2$xyear
)

initfun <- function(x){list(scale = 2, scale_m = 2)}

#
# Fit models 
#
 
m3_adults_habili <- stan(file = "habili-model-4.stan", data = datlist_habili,
           iter=5000, chains=3, 
           init = initfun, 
           control = list(max_treedepth = 10))

save(datlist_habili, m3_adults_habili, file = "habili-models-adults.rda")

precis(m3_adults_habili)


#
# check against lme4 
#

library(lme4)
glm1 <- glmer(number ~ period + (1|site), data = dat2, family = "poisson")
summary(glm1)

