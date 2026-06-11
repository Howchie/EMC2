rm(list = ls())
# remotes::install_github("https://github.com/Howchie/EMC2/", ref="dev-oo",dependencies=TRUE)
library(EMC2)

matchfun <- function(d) d$S == d$lR
n_subjects <- 30
UC <- function(d) quantile(d$rt,.9)


design_LBA <- design(factors=list(subjects=1:n_subjects,S=c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = matchfun,
                     formula =list(v~lM,B~1, t0~1, sv~1, A~1),
                     constants=c(sv=log(1)),
                     model = LBA,UC=UC)

prior_LBA <- prior(design_LBA, type = "diagonal-gamma",
                  mu_mean = c(1.3, .7, log(.8), log(.2), log(.3)),
                  mu_sd = c(.2, .1, .1, .05, .05),
                  shape = 10,
                  rate = c(.2, .2, .2, .1, .1))

SBC_LBA <- run_sbc(design_LBA, prior_LBA, replicates = 500, trials = 100,
      n_subjects = n_subjects, fileName = "SBC_LBA_UC.RData", cores_per_chain = 30)

save(SBC_LBA,file="SBC_LBA_UC1.RData")

pdf("LBA_UC.pdf")
plot_sbc_ecdf(SBC_LBA)
dev.off()



