rm(list=ls())
load("WorkingTests/Andrew_flanker.RData") # attached
library(EMC2)
library(dplyr)
names(dat)[1] <- "subjects"
dat <- dat[,c("subjects","S","CI","E","R","rt")]
dat <- dat[dat$rt>.2,]
design_WDM <- design(formula = list(a ~ E, v ~ 0 + S/CI, t0 ~ 1, Z ~ 1),data= dat, model = DDM,LT=.2)
pmean <- c(a=log(.17),a_Espeed=0,v_Sleft=-2.25,v_Sright=2.25,
 'v_Sleft:CIincongruent'=0,'v_Sright:CIincongruent'=0,t0=log(.45),Z=qnorm(.5))
psd <- c(.7,.5,2.5,2.5,1,1,.4,.4)
priorWDM <- prior(design_WDM,pmean=pmean,psd=psd,type="single")
sWDM <-  make_emc(dat,design_WDM,type="single",prior=priorWDM)
sWDM <- fit(sWDM)
ppWDM <- predict(sWDM)
correct_fun <-  function(data) factor(data$S == data$R)
plot_cdf(dat,ppWDM,layout=c(2,2), factors = c("S", "CI", "E"),
  functions = list(correct = correct_fun), defective_factor = "correct") 
plot_density(dat,ppWDM,layout=c(2,2), factors = c("S", "CI", "E"),
         functions = list(correct = correct_fun), defective_factor = "correct") 
plot_stat(dat, ppWDM, factors = c("S", "CI", "E"), stat_name = "MeanCorrect",
                        stat_fun = function(d) mean(as.numeric(d$correct)-1,na.rm=TRUE), functions = list(correct = correct_fun),
          support = c(0,1))
summary = ppWDM %>%
  mutate(Correct = as.numeric(as.logical(correct_fun(.)))) %>%
  group_by(S, CI, E) %>%
  summarise(MeanAcc = mean(Correct,na.rm=TRUE))
