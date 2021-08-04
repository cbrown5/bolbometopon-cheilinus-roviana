#Analyse trends in size distribution and maturity
#v2 uses stan
# CJ Brown 2018-10-23

rm(list = ls())

library(dplyr)
library(ggplot2)
library(readr)
library(rethinking)
library(PlotTools)
library(tidyr)
  
spp_mat_size <- c(650, 550)
tdat <- read_csv("roviana-topa-2018-10-19.csv")
hdat <- read_csv("roviana-habili-2018-10-19.csv")

datalist <- c(list(tdat), list(hdat))
sprior <- list(theta = list(prior="pc.prec", param=c(5, 0.025)))

size_models <- NULL
dat_sizes <- NULL
#
# Fit models 
#
for (i in 1:length(datalist)){

  #add in zeros
  dat <- datalist[[i]][!is.na(datalist[[i]]$number),]
  sites <- unique(dat$site)[1:4]
  dat2 <- filter(dat, site %in% sites)
  dat2$jsite <- as.integer(factor(dat2$site))
  dat2$site_survey <- paste0(dat2$site, dat2$date)
  dat2$jsurvey <- as.integer(factor(dat2$site_survey))
  dat2$period <- factor(ifelse(dat2$year < 2002, "2000s", "2018"))
  dat2$xyear <- as.integer(dat2$period)-1
  datsz <- uncount(dat2, number)
  datsz$lnsz <- log(datsz$size)
  # estimate prob mature 
   print(datsz %>% group_by(period) %>%
     summarize(sum(size>spp_mat_size[i])/n(), mean(size))) #650 for topa
  # boxplot(datsz$size ~ datsz$period)
   dat_sizes <- c(dat_sizes, list(datsz))
  datstan <- list(N = nrow(datsz),
                  N_survey = max(datsz$jsurvey), 
                  N_site = max(datsz$jsite),
                  jsurvey = datsz$jsurvey, 
                  jsite = datsz$jsite, 
                  lnsize = datsz$lnsz,
                  xyear = datsz$xyear)

  size_m1 <- stan(file = "size-model.stan", data = datstan,
                 iter=5000, chains=3, 
                 control = list(max_treedepth = 10))
  
  # shinystan::launch_shinystan(size_m1)
  # precis(size_m1)
  # exp(6.26)*exp(0.31)
  # exp(6.26 - 0.47)*exp(0.31)
  
  postdat <- extract.samples(size_m1) %>% data.frame()
  ipred <- grep('mupred', names(postdat))
  preds <- postdat[,ipred]
  predmed <- apply(preds, 2, median)
  mresids <- datsz$lnsz - predmed
  dev.new()
  plot(predmed, mresids)
  abline(h=0)
  dev.new()
  qqnorm(mresids)
  qqline(mresids, col = "red")
  size_models <- c(size_models, list(postdat))
  
}
save(size_models, file = "size-models.rda")
save(dat_sizes, file = "size-data.rda")
lapply(size_models, summary)
