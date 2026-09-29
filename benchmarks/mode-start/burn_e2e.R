# End-to-end burn with the mode start on forstmann LBA + 10-level dummy factor
# (p = 71, 19 subjects x 3 chains): prints the mode-start summary and burn time.
#   Rscript benchmarks/mode-start/burn_e2e.R
library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])); cat("EMC2 from", find.package("EMC2"), "\n")
set.seed(5)
dat <- forstmann
dat$F <- factor(sample(1:10, nrow(dat), TRUE))
invisible(capture.output(des <- suppressMessages(design(data = dat, model = LBA,
  formula = list(v ~ lM * E * F, B ~ E * lR, A ~ 1, t0 ~ E, sv ~ lM),
  constants = c(sv = log(1)), matchfun = function(d) d$S == d$lR))))
emc <- make_emc(dat, des, n_chains = 3, type = "diagonal-gamma")
cat("p =", emc[[1]]$n_pars, "\n")
emc <- run_emc(emc, stage = "preburn", cores_per_chain = 2, cores_for_chains = 3, verbose = FALSE, stop_criteria = list(iter = 150))
t_burn <- system.time(emc <- run_emc(emc, stage = "burn", cores_per_chain = 2, cores_for_chains = 3,
  verbose = TRUE, stop_criteria = list(iter = 100, max_gd = Inf, min_unique = 0, min_es = 0)))
print(t_burn)
print(table(unlist(lapply(emc, function(x) x$burn_mode_stats$reason))))
