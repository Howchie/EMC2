# optimHess vs the batched gradient-difference Hessian on one forstmann subject,
# LBA widened with a dummy factor (KS = its level counts; p = 23, 41, 71, 101).
#   KS=2,5,10,15 Rscript benchmarks/mode-start/hessian_timing.R
lib <- Sys.getenv("EMC_LIB"); library(EMC2, lib.loc = lib)
out <- Sys.getenv("MODE_START_OUT", "benchmarks/mode-start/results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2"); for (f in c(".emc_burn_ll_batch", ".emc_optimize_subject", ".emc_gradient_difference_hessian")) assign(f, get(f, ns))
options(emc2.burn_mode_start = FALSE, emc2.laplace_refresh = FALSE)
set.seed(3)
dat <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1], ]
dat$subjects <- droplevels(dat$subjects)
for (k in as.integer(strsplit(Sys.getenv("KS", "2,5,10"), ",")[[1]])) {
  dat$F <- factor(sample(seq_len(k), nrow(dat), TRUE))
  des <- design(data = dat, model = LBA,
    formula = list(v ~ lM * E * F, B ~ E * lR, A ~ 1, t0 ~ E, sv ~ lM),
    constants = c(sv = log(1)), matchfun = function(d) d$S == d$lR)
  emc <- suppressMessages(make_emc(dat, des, n_chains = 1, type = "single"))
  sm <- emc[[1]]
  p <- sm$n_pars
  ll_batch <- function(x) .emc_burn_ll_batch(x, sm$data[[1]], sm$model, sm$par_names)
  prior <- list(mu = sm$prior$theta_mu_mean, var = sm$prior$theta_mu_var)
  x0 <- as.numeric(prior$mu); x0[sm$par_names == "t0"] <- log(0.1)
  opt <- system.time(res <- .emc_optimize_subject(x0, prior, ll_batch))
  if (is.null(res$alpha)) { print(res$reason); str(prior); print(ll_batch(x0)); next }
  xm <- res$alpha
  t_old <- system.time(h_old <- stats::optimHess(xm, function(z) -ll_batch(z)))
  t_new <- system.time(h_new <- .emc_gradient_difference_hessian(xm, ll_batch))
  t_one <- system.time(for (r in 1:20) ll_batch(xm))[["user.self"]] / 20
  cat(sprintf(paste0("p=%3d  BFGS+Hessian %.1fs (%s, conv=%s)  optimHess %.2fs  batched %.2fs  ",
    "speedup %.0fx  rel diff %.1e  single ll %.1fms\n"), p, opt[["user.self"]], res$reason,
    res$converged, t_old[["user.self"]], t_new[["user.self"]],
    t_old[["user.self"]] / t_new[["user.self"]],
    max(abs(h_new - h_old)) / max(abs(h_old)), 1000 * t_one))
}
