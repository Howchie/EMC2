rm(list = ls())
library(EMC2)

matchfun <- function(d) d$S == d$lR
UC <- UT <- function(d) quantile(d$rt,.9)
LC <- LT <- function(d) quantile(d$rt,.1)


design_LNR <- design(factors=list(subjects=1,S=c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = matchfun,TC=list(rt_resolution=1/60,LT=LT,UT=UT),
                     formula =list(m~lM,s~1, t0~1),
                     model = LNR)

prior_LNR <- prior(design_LNR, type = "single",
                  pmean = c(-.7, -.5, log(1), log(.2)),
                  psd = c(.2, .1, .1, .05))

run_sbc(design_LNR, prior_LNR, replicates = 500, trials = 100,
                          fileName = "SBC_LNR_T.RData", cores_per_chain = 3)


design_RDM <- design(factors=list(subjects=1,S=c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = matchfun,TC=list(rt_resolution=1/60,LT=LT,UT=UT),
                     formula =list(v~lM,B~1, t0~1, A~1, s ~ 1),
                     constants=c(s=log(1)),
                     model = RDM)

prior_RDM <- prior(design_RDM, type = "single",
                  pmean = c(1.4, .3, log(1.5), log(.2), log(.3)),
                  psd = c(.05, .1, .1, .05, .05))

run_sbc(design_RDM, prior_RDM, replicates = 500, trials = 100,
                        n_subjects = n_subjects, fileName = "SBC_RDM_T.RData", cores_per_chain = 3)


design_LBA <- design(factors=list(subjects=1,S=c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = matchfun,TC=list(rt_resolution=1/60,LC=LC,UC=UC),
                     formula =list(v~lM,B~1, t0~1, sv~1, A~1),
                     constants=c(sv=log(1)),
                     model = LBA)

prior_LBA <- prior(design_LBA, type = "single",
                  pmean = c(1.3, .7, log(.8), log(.2), log(.3)),
                  psd = c(.2, .1, .1, .05, .05))

run_sbc(design_LBA, prior_LBA, replicates = 500, trials = 100, plot_data = FALSE,
                  iter = 1000, n_post = 1000, fileName = "SBC_LBA_C.RData",
                  cores_per_chain = 3)


design_WDM <- design(factors=list(subjects=1,S=c("left", "right")),
                       Rlevels = c("left", "right"),TC=list(rt_resolution=1/60,LT=LT,UT=UT),
                       formula =list(v~1,a~1, t0~1, s~1, Z~1),
                       constants=c(s=log(1)),
                       model = DDM)

prior_WDM <- prior(design_WDM, type = "single",
                  pmean = c(1, log(.8), log(.3), qnorm(.5)),
                  psd = c(.15, .15, .1, .05))

run_sbc(design_WDM, prior_WDM, replicates = 500, trials = 100,
        fileName = "SBC_WDM_T.RData", cores_per_chain = 3)
