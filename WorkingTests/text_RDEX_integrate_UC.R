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
load("WorkingTests/RDEXintegrate.RData")

# "st" level is stop-triggered (lI=1), the rest are go (lI=2)
lIfun_st <- function(d) factor(ifelse(as.character(d$lR) == "st", 1, 2), levels = 1:2)
p_ssd     = c(0.25)



#### 1 choice RDEX ----
designSSRDEX1 <- design(model = SSRDEX,
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 1, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSRDEX1, doMap = FALSE)
p_vector[] <- c(log(1.5),log(1),log(0.2),log(0.2),
                log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

n_trials=400

# undebug(make_data)
# undebug(make_missing)
#
# designSSRDEX1 <- design(model = SSRDEX,
#   factors  = list(subjects = 1, S = 1),Rlevels  = 1,
#   TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
#   matchfun = function(d) as.character(d$S) == as.character(d$lR),
#   functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
#   formula  = list(v ~ 1, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
# )
#
# p_vector[] <- c(log(1.5),log(1),log(0.2),log(0.2),
#                 log(0.05),log(0.1),qnorm(0),qnorm(0))
# dat1 <- make_data(p_vector, designSSRDEX1, n_trials = 1000,
#                  staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))


dat1 <- make_data(p_vector, designSSRDEX1, n_trials = n_trials,
                 staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
check_dat(dat1)


emcR1 <- make_emc(dat1,designSSRDEX1,type="single")
t1=system.time({emcR1 <- fit(emcR1gl)})[3]
recovery(emcR1,p_vector,main="n=1, auto")

designSSRDEX1i <- design(model = SSRDEX("integrate"),
  factors  = list(subjects = 1, S = 1),Rlevels  = 1,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 1, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR1i <- make_emc(dat1,designSSRDEX1i,type="single")
t1i=system.time({emcR1i <- fit(emcR1i)})[3]
recovery(emcR1i,p_vector,main="n=1, integrate")

#### 2 choice RDEX ----
designSSRDEX2 <- design(model = SSRDEX,
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSRDEX2, doMap = FALSE)
p_vector[] <- c(log(0.25),log(1.5),log(1),log(0.2),
                log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

# n_trials=200
# dat2 <- make_data(p_vector, designSSRDEX2, n_trials = n_trials,
#                  staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
# check_dat(dat2)

emcR2 <- make_emc(dat2,designSSRDEX2,type="single")
t2 <- system.time({emcR2 <- fit(emcR2)})[3]
recovery(emcR2,p_vector,main="n=2, auto")

designSSRDEX2i <- design(model = SSRDEX("integrate"),
  factors  = list(subjects = 1, S = 1:2),Rlevels  = 1:2,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR2i <- make_emc(dat2,designSSRDEX2i,type="single")
t2i <- system.time({emcR2i <- fit(emcR2i)})[3]
recovery(emcR2i,p_vector,main="n=2, integrate")

#### 4 choice RDEX ----
designSSRDEX4 <- design(model = SSRDEX,
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
p_vector <- sampled_pars(designSSRDEX4, doMap = FALSE)
p_vector[] <- c(log(0.1),log(3),log(2),log(0.2),log(0.2),log(0.05),log(0.1),qnorm(0.08),qnorm(0.08))

# n_trials=100
# dat4 <- make_data(p_vector, designSSRDEX4, n_trials = n_trials,
#                  staircase = list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf))
# check_dat(dat4)

emcR4 <- make_emc(dat4,designSSRDEX4,type="single")
t4=system.time({emcR4 <- fit(emcR4)})[3]
# Time difference of 3.419095 mins
recovery(emcR4,p_vector,main="n=4, auto")

designSSRDEX4i <- design(model = SSRDEX("integrate"),
  factors  = list(subjects = 1, S = 1:4),Rlevels  = 1:4,
  TC=list(UC=\(d)quantile(d$rt,probs=.9,na.rm=T),verbose=T),
  matchfun = function(d) as.character(d$S) == as.character(d$lR),
  functions = list(lI = lIfun_st,SSD = function(d) SSD_function(d, SSD = NA, pSSD = p_ssd)),
  formula  = list(v ~ 0 + lM, B ~ 1, t0 ~ 1,muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
)
emcR4i <- make_emc(dat4,designSSRDEX4i,type="single")
t4i=system.time({emcR4i <- fit(emcR4i)})[3]
recovery(emcR4i,p_vector,main="n=4, integrate")

# save(dat1,dat2,dat4,file="WorkingTests/RDEXintegrate.RData")

#### print timing results
print(c(t1=t1,t1i=t1i))
print(c(t2=t2,t2i=t2i))
print(c(t4=t4,t4i=t4i))

