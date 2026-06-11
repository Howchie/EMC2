rm(list = ls())
# remotes::install_github("https://github.com/Howchie/EMC2/", ref="dev-oo",dependencies=TRUE)
library(EMC2)

matchfun <- function(d) d$S == d$lR
n_subjects <- 30
UC <- function(d) quantile(d$rt,.9)


design_LNR <- design(factors=list(subjects=1:n_subjects,S=c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = matchfun,
                     formula =list(m~lM,s~1, t0~1),
                     model = LNR,UC=UC)

prior_LNR <- prior(design_LNR, type = "diagonal-gamma",
                  mu_mean = c(-.7, -.5, log(1), log(.2)),
                  mu_sd = c(.2, .1, .1, .05),
                  shape = 10,
                  rate = c(.2, .2, .2, .1))

SBC_LNR <- run_sbc(design_LNR, prior_LNR, replicates = 500, trials = 100,
      n_subjects = n_subjects, fileName = "SBC_LNR_UC.RData", cores_per_chain = 30)

save(SBC_LNR,file="SBC_LNR_UC1.RData")

pdf("LNR_UC.pdf")
plot_sbc_ecdf(SBC_LNR, layout = c(2,2))
dev.off()

