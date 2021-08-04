#Analyse abundance for Roviana data - lunar periods pre 2017 
# CJ Brown 2018-10-23

library(dplyr)
library(ggplot2)
library(readr)
library(INLA)
library(PlotTools)


tdat <- read_csv("roviana-topa-2018-10-19.csv")
hdat <- read_csv("roviana-habili-2018-10-19.csv")

datalist <- c(list(tdat), list(hdat))
sprior <- list(theta = list(prior="loggamma", param=c(1, 1)))

models <- NULL
#
# Fit models 
#
for (i in 1:length(datalist)){

  #add in zeros
  datalist[[i]]$number[is.na(datalist[[i]]$number)] <- 0
  sitevars <- datalist[[i]] %>% select(site, date, month, year, lunar, area) %>%
    distinct()
  dat <- datalist[[i]] %>% group_by(site, date) %>%
    summarize(number = sum(number)) %>%
    ungroup() %>%
    inner_join(sitevars)
  sites <- unique(dat$site)[5:8]
  dat2 <- filter(dat, site %in% sites)
  dat2$jsite <- as.integer(factor(dat2$site))
  
  #
  # Lunar model
  #
  datl <- filter(dat2, year < 2017)
  datl$lunar <- factor(datl$lunar)
  
  lunar_m1 <- inla(number ~ lunar  + 
                     f(site, model = "iid", hyper = sprior), 
                  data = datl, family = "poisson", 
                  quantiles = c(0.055, 0.5, 0.945))
   
  summary(lunar_m1)  
  models <- c(models, list(lunar_m1))
  
}
save(models, file = "lunar-models.rda")
lapply(models, summary)
