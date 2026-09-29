# Post-preburn states for the mode-start benchmarks (2026-09-29).
#   CASE=lba71|joint|corr|freeA|wide Rscript benchmarks/mode-start/setup.R
# lba71: forstmann LBA + 10-level dummy factor (p = 71, 19 subjects).
# joint: forstmann LBA + LNR joint model (p = 74, 19 subjects).
# corr/freeA/wide: RDM fixtures from benchmarks/group-move/fixture_fits.R (run
# that first; the fixture .rds files are local, not in git).
# Writes <out>/preburn_<case>.rds after 150 preburn iterations, 3 chains.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
cat("EMC2 from", find.package("EMC2"), "\n")
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

case <- Sys.getenv("CASE"); set.seed(11)
q <- function(e) { invisible(capture.output(r <- suppressMessages(e))); r }
if (case == "lba71") {
  # forstmann with a 10-level dummy factor in the drift formula: p = 71.
  set.seed(5)
  dat <- forstmann; dat$F <- factor(sample(1:10, nrow(dat), TRUE))
  des <- q(design(data = dat, model = LBA,
    formula = list(v ~ lM * E * F, B ~ E * lR, A ~ 1, t0 ~ E, sv ~ lM),
    constants = c(sv = log(1)), matchfun = function(d) d$S == d$lR))
  emc <- q(make_emc(dat, des, n_chains = 3, type = "diagonal-gamma"))
} else if (case == "joint") {
  dat <- forstmann; dat$F <- factor(sample(1:5, nrow(dat), TRUE))
  d1 <- q(design(data = dat, model = LBA, formula = list(v ~ lM * E * F, B ~ E * lR, A ~ 1, t0 ~ E, sv ~ lM),
    constants = c(sv = log(1)), matchfun = function(d) d$S == d$lR))
  d2 <- q(design(data = dat, model = LNR, formula = list(m ~ lM * E * F, s ~ lM, t0 ~ 1), matchfun = function(d) d$S == d$lR))
  emc <- q(make_emc(list(dat, dat), list(a = d1, b = d2), type = "diagonal-gamma", n_chains = 3))
} else {
  fx <- readRDS(file.path("benchmarks/group-move/results/fixture_A0v15", paste0("fixture-", case, ".rds")))
  des <- q(design(data = fx$data, model = RDM, matchfun = function(d) d$S == d$lR,
    formula = fx$settings$formula, constants = fx$settings$constants))
  emc <- q(make_emc(fx$data, des, n_chains = 3, type = "standard"))
}
cat(case, "p =", emc[[1]]$n_pars, " subjects =", emc[[1]]$n_subjects, "\n")
emc <- q(run_emc(emc, stage = "preburn", cores_per_chain = as.integer(Sys.getenv("CORES_PER_CHAIN", "1")), cores_for_chains = 3, verbose = FALSE, stop_criteria = list(iter = 150)))
saveRDS(emc, file.path(out, paste0("preburn_", case, ".rds")))
