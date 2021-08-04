#Analyse abundance for Roviana data
# CJ Brown 2018-11-13
# v4 stan model 
#Adults onlyl model, for IUCN declines

rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(rethinking)
library(PlotTools)


datin <- read_csv("roviana-topa-2018-10-19.csv")
maxn <- read_csv("max-catch-data.csv")
topa_mat <- 650
nyears <- 18


#add in zeros
datin$number[is.na(datin$number)] <- 0
sitevars <- datin %>% select(site, date, month, year, lunar, area) %>%
  distinct()
dat <- datin %>% 
   filter(size > topa_mat) %>%
   group_by(site, date) %>%
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

# maxn2 <- filter(maxn, Year != "1980s")
maxn$xyear <- ifelse(maxn$Year == "2018", 1, 0)
maxn$xyearh <- ifelse(maxn$Year == "1980s", 1, 0)

maxn2 <- filter(maxn, Year != "1980s")
maxn2$xyear <- as.integer(factor(maxn2$Year))-1


#
# Data for models 
#


datlist <- list(
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
 
m3_adults_topa <- stan(file = "habili-model-4.stan", data = datlist,
           iter=5000, chains=3, 
           init = initfun, 
           control = list(max_treedepth = 10))
#using habili model because want poisson distribution for data

save(datlist, m3_adults_topa, file = "topa-models-adults.rda")

precis(m3_adults_topa)
post4 <- extract.samples(m3_adults_topa) %>% data.frame()
minx <- 0.5; maxx <- 10

dens(exp(-(post4$b_n)), xlab = "Times decline", xlim = c(minx, maxx), 
     ylab = "Probability density", col = "#63b8ba", lwd = 2)

#
# check against lme4 
#

dat2 %>% group_by(period) %>%
  summarize(mn = mean(number))

library(lme4)
glm1 <- glmer(number ~ period + (1|site), data = dat2, family ="poisson")
summary(glm1)

plot(glm1)

