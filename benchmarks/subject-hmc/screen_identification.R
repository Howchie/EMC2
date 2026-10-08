# Identification screen: every subject's likelihood at its own multi-start
# MLE must have a positive definite, numerically stable observed information.
#   EMC_LIB=<lib> Rscript benchmarks/subject-hmc/screen_identification.R <fixture.rds> [n_starts]
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB", .libPaths()[1])))
ns <- asNamespace("EMC2"); library(parallel)
args <- commandArgs(TRUE); fx <- readRDS(args[1]); n_starts <- if (length(args) > 1) as.integer(args[2]) else 6L
emc <- suppressMessages(make_emc(fx$data, fx$design, n_chains = 1))
ch <- emc[[1]]; pn <- ch$par_names; p <- length(pn); truth <- fx$subject_effects[, pn, drop = FALSE]
res <- mclapply(seq_len(ch$n_subjects), function(s) {
  set.seed(s)
  f <- function(x) { v <- ns$.emc_burn_ll_batch(matrix(x, 1), ch$data[[s]], ch$model, pn); if (!is.finite(v)) -1e10 else v }
  starts <- c(list(truth[s, ]), lapply(seq_len(n_starts - 1), function(k) fx$group_means[pn] + rnorm(p, 0, sqrt(diag(fx$intended_covariance[pn, pn])) * 1.5)))
  fits <- lapply(starts, function(x0) optim(x0, function(x) -f(x), method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)))
  vals <- sapply(fits, function(o) -o$value); best <- fits[[which.max(vals)]]
  # Distinct optima: within 0.5 nats of the best but more than 0.1 apart in parameters.
  modes <- sum(vals > max(vals) - 0.5 & sapply(fits, function(o) max(abs(o$par - best$par))) > 0.1)
  H <- -optimHess(best$par, f)
  H2 <- -optimHess(best$par, f, control = list(ndeps = rep(3e-4, p)))
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  ev2 <- eigen((H2 + t(H2)) / 2, symmetric = TRUE, only.values = TRUE)$values
  se <- if (min(ev) > 0) sqrt(diag(solve(H))) else rep(NA, p)
  list(s = s, ll = max(vals), spread = max(vals) - min(vals), other_modes = modes, min_eig = min(ev), min_eig_alt = min(ev2),
       cond = max(ev) / min(ev), se = se, err = best$par - truth[s, ])
}, mc.cores = min(20, ch$n_subjects))
tab <- do.call(rbind, lapply(res, function(r) data.frame(s = r$s, ll = round(r$ll, 1), other_modes = r$other_modes,
  min_eig = signif(r$min_eig, 3), min_eig_alt = signif(r$min_eig_alt, 3), cond = signif(r$cond, 3),
  se_A = round(r$se[match("A", pn)], 3), se_B = round(r$se[match("B", pn)], 3), se_t0 = round(r$se[match("t0", pn)], 3),
  z_A = round(r$err[match("A", pn)] / r$se[match("A", pn)], 2))))
print(tab)
sd <- sqrt(diag(fx$intended_covariance[pn, pn]))
cat(sprintf("PD at MLE: %d/%d (alt step %d/%d) | max condition %.3g | median SE/group SD: A %.2f B %.2f t0 %.2f\n",
  sum(tab$min_eig > 0), nrow(tab), sum(tab$min_eig_alt > 0), nrow(tab), max(tab$cond),
  median(tab$se_A) / sd["A"], median(tab$se_B) / sd["B"], median(tab$se_t0) / sd["t0"]))
