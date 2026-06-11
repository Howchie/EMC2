rm(list = ls())
# remotes::install_github("https://github.com/Howchie/EMC2/", ref="dev-oo",dependencies=TRUE)
library(EMC2)

matchfun <- function(d) d$S == d$lR
n_subjects <- 30


design_DDM <- design(factors=list(subjects=1:n_subjects,S=c("go","nogo")),
                       Rlevels = c("go","nogo"),
                       formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1, st0 ~ 1),
                       constants=c(s=log(1)),
                       model = DDM)


prior_DDM <- prior(design_DDM, type = "diagonal-gamma",
                  pmean = c(1.2, log(.8), log(.3), qnorm(.5), log(.1), qnorm(.05), log(.05)),
                  psd = c(.15, .15, .1, .05, .1, .1, .15),
                  shape = 10,
                  rate = c(.2, .2, .2, .1, .1, .1, .1))

SBC_DDM <- run_sbc(design_DDM, prior_DDM, replicates = 250, trials = 200,
      n_subjects = n_subjects, fileName = "SBC_DDMgng.RData", cores_per_chain = 30)

save(SBC_DDM,file="SBC_DDMgng1.RData")

pdf("DDMgng.pdf")
plot_sbc_ecdf(SBC_DDM)
dev.off()
