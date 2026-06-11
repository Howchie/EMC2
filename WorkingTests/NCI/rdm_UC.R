rm(list = ls())
# remotes::install_github("https://github.com/Howchie/EMC2/", ref="dev-oo",dependencies=TRUE)
library(EMC2)

matchfun <- function(d) d$S == d$lR
n_subjects <- 30
UC <- function(d) quantile(d$rt,.9)


design_RDM <- design(factors=list(subjects=1:n_subjects,S=c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = matchfun,
                     formula =list(v~lM,B~1, t0~1, A~1, s ~ 1),
                     constants=c(s=log(1)),
                     model = RDM,UC=UC)

prior_RDM <- prior(design_RDM, type = "diagonal-gamma",
                  mu_mean = c(1.4, .3, log(1.5), log(.2), log(.3)),
                  mu_sd = c(.05, .1, .1, .05, .05),
                  shape = 10,
                  rate = c(.1, .2, .2, .2, .1))

SBC_RDM <- run_sbc(design_RDM, prior_RDM, replicates = 500, trials = 100,
      n_subjects = n_subjects, fileName = "SBC_RDM_UC.RData", cores_per_chain = 30)

save(SBC_RDM,file="SBC_RDM_UC1.RData")

pdf("RDM_UC.pdf")
plot_sbc_ecdf(SBC_RDM)
dev.off()

