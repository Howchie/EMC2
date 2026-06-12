rm(list=ls())
library(EMC2)

check_dat <- function(dat) {
  datok <- dat[!is.na(dat$R),]
  if (length(levels(dat$R))>1) print(c(Accuracy=mean(datok$S==datok$R)))
  print(c(MRT=mean(datok$rt)))
  pStopSignal=mean(is.finite(dat$SSD))
  print(c(pStopSignal=pStopSignal,
          pStop=mean(is.finite(dat$SSD)&is.na(dat$R))/pStopSignal))
  print(c(SSD=sort(unique(dat$SSD))))
}


# These data sets all have 10 SSDs
load("WorkingTests/EXGintegrate.RData")


# "st" level is stop-triggered (lI=1), the rest are go (lI=2)
lIfun_st <- function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)
p_ssd     = c(0.25)


##### n = 1 ----

# For n = 1 used the following, n>1 case uses the defaults

rt_resolution <- stairstep <- .2

#### n = 1, untruncated ----

designSSEXG1 <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=-Inf,exgS_lb=-Inf)
)
p_vector <- sampled_pars(designSSEXG1, doMap = FALSE)
p_vector[] <- c(log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=400
dat1 <- make_data(p_vector, designSSEXG1, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=stairstep,stairmin=0,stairmax=Inf))
check_dat(dat1)

emcR1f <- make_emc(dat1,designSSEXG1,type="single",rt_resolution = rt_resolution)
tAf <- system.time({emcR1f <- fit(emcR1f)})[3]
recovery(emcR1f,p_vector,main="n=1, auto untruncated")


designSSEXG1 <- design(model = SSEXG("gl"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=-Inf,exgS_lb=-Inf)
)
emcR1gl <- make_emc(dat1,designSSEXG1gl,type="single",rt_resolution = rt_resolution)
tAfGL <- system.time({emcR1gl <- fit(emcR1gl)})[3]
recovery(emcR1gl,p_vector,main="n=1, GL truncated")

designSSEXG1 <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=-Inf,exgS_lb=-Inf)
)
emcR1i <- make_emc(dat1,designSSEXG1i,type="single",rt_resolution = rt_resolution)
tAfi <- system.time({emcR1i <- fit(emcR1i)})[3]
recovery(emcR1i,p_vector,main="n=1, integrate untruncated")


#### n = 1, truncated ----

designSSEXG1 <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=.2,exgS_lb=.2)
)
p_vector <- sampled_pars(designSSEXG1, doMap = FALSE)
p_vector[] <- c(log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

# n_trials=400
# dat1t <- make_data(p_vector, designSSEXG1, n_trials = n_trials,
#                  staircase = list(SSD0=.25,stairstep=stairstep,stairmin=0,stairmax=Inf))
# check_dat(dat1t)

emcR1t <- make_emc(dat1t,designSSEXG1,type="single",rt_resolution = rt_resolution)
tAt <- system.time({emcR1t <- fit(emcR1t)})[3]
recovery(emcR1t,p_vector,main="n=1, auto truncated")


designSSEXG1 <- design(model = SSEXG("gl"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=.2,exgS_lb=.2)
)
emcR1tGL <- make_emc(dat1t,designSSEXG1,type="single",rt_resolution = rt_resolution)
tAtGL <- system.time({emcR1tGL <- fit(emcR1tGL)})[3]
recovery(emcR1tGL,p_vector,main="n=1, GL truncated")

designSSEXG1 <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 1, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1,
                  exg_lb~1,exgS_lb~1),constants=c(exg_lb=.2,exgS_lb=.2)
)
emcR1ti <- make_emc(dat1t,designSSEXG1,type="single",rt_resolution = rt_resolution)
tAti <- system.time({emcR1ti <- fit(emcR1ti)})[3]
recovery(emcR1ti,p_vector,main="n=1, integrate truncated")


#### 2 choice EXG ----
designSSEXG2 <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSEXG2, doMap = FALSE)
p_vector[] <- c(log(.69),log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

# n_trials=200
# dat2 <- make_data(p_vector, designSSEXG2, n_trials = n_trials,
#                  staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
# check_dat(dat2)

emcR2 <- make_emc(dat2,designSSEXG2,type="single")
t2 <- system.time({emcR2 <- fit(emcR2)})[3]
recovery(emcR2,p_vector,main="n=2, auto")

designSSEXG2i <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR2i <- make_emc(dat2,designSSEXG2i,type="single")
t2i <- system.time({emcR2i <- fit(emcR2i)})[3]
recovery(emcR2i,p_vector,main="n=2, integrate")

#### 4 choice EXG ----
designSSEXG4 <- design(model = SSEXG,
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSEXG4, doMap = FALSE)
p_vector[] <- c(log(.8),log(.5),log(.1),log(0.25),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

# n_trials=100
# dat4 <- make_data(p_vector, designSSEXG4, n_trials = n_trials,
#                  staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
# check_dat(dat4)

emcR4 <- make_emc(dat4,designSSEXG4,type="single")
t4=system.time({emcR4 <- fit(emcR4)})[3]
recovery(emcR4,p_vector,main="n=4, auto")

designSSEXG4i <- design(model = SSEXG("integrate"),
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR4i <- make_emc(dat4,designSSEXG4i,type="single")
t4i=system.time({emcR4i <- fit(emcR4i)})[3]
recovery(emcR4i,p_vector,main="n=4, integrate")

# save(dat1,dat1t,dat2,dat4,file="WorkingTests/EXGintegrate.RData")

#### print timing results
print(c(tAf=tAf,tAfGL=tAfGL,tAfi=tAfi))
print(c(tAt=tAt,tAtGL=tAtGL,tAti=tAti))
print(c(t2=t2,t2i=t2i))
print(c(t4=t4,t4i=t4i))

