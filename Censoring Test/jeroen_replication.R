rm(list = ls())
# remotes::install_github("timmerj1/EMC2-advanced",ref="dev_rename_jeroen")
library(EMC2)

ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

n_trials_small = 200
n_samples_small = 10
n_trials_large = 10000
n_samples_large = 1

# Standard models ----
# * LBA ====
designLBA_small <- design(
  factors = list(S = c("left", "right"),
                 subjects = 1:n_samples_small),
  Rlevels = c("left", "right"),
  model = LBA,
  matchfun = function(d) d$S == d$lR,
  formula = list(B ~ 1, v ~ lM, A ~ 1, sv ~ lM, t0 ~ 1),
  constants = c(sv = 0),
  contrasts=list(lM=ADmat)
)

designLBA_large <- design(
  factors = list(S = c("left", "right"),
                 subjects = 1:n_samples_large),
  Rlevels = c("left", "right"),
  model = LBA,
  matchfun = function(d) d$S == d$lR,
  formula = list(B ~ 1, v ~ lM, A ~ 1, sv ~ lM, t0 ~ 1),
  constants = c(sv = 0),
  contrasts=list(lM=ADmat)
)

parsLBA <- sampled_pars(designLBA_small,doMap = F)
parsLBA[1:6] <- c(log(2),3,1,log(2),log(.75),log(.2))
datLBA_small <- make_data(parsLBA, designLBA_small, n_trials = n_trials_small)
datLBA_large <- make_data(parsLBA, designLBA_large, n_trials = n_trials_large)
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)
tapply(Cfun(datLBA_large),datLBA_large[,c("S")],mean)
# * RDM ====

designRDM_small <- design(
  data=datLBA_small, model = RDM,
  matchfun = function(d) d$S == d$lR,
  formula = list(B ~ 1, v ~ lM, s ~ lM, t0 ~ 1),
  constants = c(s = 0),
  contrasts=list(lM=ADmat)
)

designRDM_large <- design(
  data=datLBA_large, model = RDM,
  matchfun = function(d) d$S == d$lR,
  formula = list(B ~ 1, v ~ lM, s ~ lM, t0 ~ 1),
  constants = c(s = 0),
  contrasts=list(lM=ADmat)
)


parsRDM <- sampled_pars(designRDM_small,doMap = F)
parsRDM[1:5] <- c(log(3),1,.4,log(.75),log(.2))
datRDM_large <- make_data(parsRDM, designRDM_large, data = datLBA_large)
datRDM_small <- make_data(parsRDM, designRDM_small, data = datLBA_small)

# * LNR ====
designLNR_small <- design(
  data=datLBA_small, model = LNR,
  matchfun = function(d) d$S == d$lR,
  formula = list(m ~ lM, s ~ lM, t0 ~ 1),
  contrasts=list(lM=ADmat)
)

designLNR_large <- design(
  data=datLBA_large, model = LNR,
  matchfun = function(d) d$S == d$lR,
  formula = list(m ~ lM, s ~ lM, t0 ~ 1),
  contrasts=list(lM=ADmat)
)

parsLNR <- sampled_pars(designLNR_small,doMap = F)
parsLNR[1:5] <- c(log(.75),log(.65),log(.5),log(.8),log(.4))
datLNR_large <- make_data(parsLNR, designLNR_large, data = datLBA_large)
datLNR_small <- make_data(parsLNR, designLNR_small, data = datLBA_small)

# Simulation ----

run_txt <- "library(EMC2)
load(\"./01_simulation/baseline/EMCs/s%1$s_%2$s.Rdata\")
s <- fit(s, fileName = \"./01_simulation/baseline/tmp_files/tmps%1$s_%2$s.Rdata\")
save(s, file = \"./01_simulation/baseline/EMCs/s%1$s_%2$s.Rdata\")"


cases <- expand.grid(
  model = c("LBA", "LNR", "RDM"),
  size = c("large", "small")
)

run_files <- sprintf(run_txt, cases[,1], cases[,2])

file_names <- sprintf(
  "./01_simulation/baseline/run_files/run%1$s_%2$s.R",
  cases[,1], cases[,2]
)

s_file_names <- sprintf(
  "./01_simulation/baseline/EMCs/s%1$s_%2$s.Rdata",
  cases[,1], cases[,2]
)

for (i in 1:nrow(cases)) {
  run_dir <- dirname(file_names[i])
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)
  write(run_files[i], file = file_names[i])

  s_dir <- dirname(s_file_names[i])
  if (!dir.exists(s_dir)) dir.create(s_dir, recursive = TRUE)
  s <- make_emc(get(paste0("dat", as.character(cases[i,1]), "_", cases[i,2])),
           get(paste0("design", as.character(cases[i,1]), "_", cases[i,2])),
           type = "single")
  save(s, file = s_file_names[i])
}
