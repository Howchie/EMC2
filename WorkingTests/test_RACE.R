# Test that variable RACE levels still work
rm(list=ls())
library(EMC2)

matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR) |
  (d$lR=="pm" & as.numeric(d$S)>2)
designLBA <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=LBA,constants=c(v_RACE3=0),
  formula=list(v~RACE*lM,B~1,t0~1,A~1),
)

## First test recovery against Old C
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(2), log(4), log(1), log(2), log(0.2),log(.5))

# Make square data so can remove pm in RACE = 2
template <- make_data(p_vector,designLBA,n_trials=1000)
template <- template[!(template$RACE==2 & (template$S %in% c("leftpm","rightpm"))),]
dat <- make_data(p_vector,designLBA,data=template)

emc <- make_emc(dat, designLBA, type = "single", compress = T)
emc1 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30)

print(recovery(emc1,p_vector,selection="alpha"))

## Also test recovery with RACE plus censoring and truncation
TC=list(LT=.5,UT=1.8,LC=.55,UC=1.6,verbose=TRUE)
datCT <- make_data(p_vector,designLBA,n_trials=10000,TC=TC)
emc <- make_emc(dat, designLBA, type = "single", compress = T)
emc2 <- fit(emc,stop_criteria = list(
  sample = list(
    iter = 1000,
    max_gd = 1.10,
    max_flat_loc = 0.5,
    flat_selection = c("alpha", "subj_ll"),
    flat_p1 = 1/3,
    flat_p2 = 1/3,
    max_sample_iter = 5000
  ),cores_per_chain=3, cores_for_chains = 3), max_tries=30)

print(recovery(emc2,p_vector,selection="alpha"))
save.image("testRACE.RData")