# FRQ heavy validation -- NOT part of R CMD check.
#
# tests/testthat/test-frq.R keeps the fast oracles.  The two things it cannot
# afford are here:
#
#   1. Large-sample Monte Carlo of the literal finite-reservoir process across
#      a full parameter grid, at a tolerance tight enough to detect a small
#      systematic error rather than only a gross one.
#   2. A sampler-based recovery study.
#
# Run with:  Rscript WorkingTests/test_frq_validation.R

library(EMC2)
library(testthat)

# ---------------------------------------------------------------------------
# 1. Large-sample Monte Carlo against the literal generative process
# ---------------------------------------------------------------------------

cat("== literal finite-reservoir Monte Carlo ==\n")
set.seed(20260816)
nsim <- 4e5
cfgs <- list(c(N = 7, K = 3, p = 0.80, lambda = 2.5),
             c(N = 4, K = 4, p = 0.60, lambda = 1.2),
             c(N = 5, K = 1, p = 0.45, lambda = 4.0),
             c(N = 12, K = 5, p = 0.95, lambda = 3.0),
             c(N = 3, K = 2, p = 0.30, lambda = 0.8))
for (cfg in cfgs) {
  N <- cfg[["N"]]; K <- cfg[["K"]]; p <- cfg[["p"]]; lam <- cfg[["lambda"]]
  a <- K; b <- N - K + 1
  avail <- matrix(runif(nsim * N) < p, nrow = nsim)
  lat <- matrix(rexp(nsim * N, lam), nrow = nsim)
  lat[!avail] <- Inf
  Tk <- apply(lat, 1, function(r) sort(r)[K])

  h <- pbeta(p, a, b)
  tau <- -log1p(-qbeta(0.5 * h, a, b) / p) / lam
  probe <- quantile(Tk[is.finite(Tk)], c(.05, .1, .25, .5, .75, .9, .95))
  emp <- vapply(probe, function(x) mean(Tk <= x), numeric(1))
  thy <- EMC2:::pfrq(as.numeric(probe), a, b, h, tau)
  # 4e5 draws gives a binomial SE of at most 8e-4, so 3e-3 is ~4 SE.
  cat(sprintf("  N=%2d K=%d p=%.2f lam=%.1f | max|F_emp - F_thy| = %.2e | "
              , N, K, p, lam, max(abs(emp - thy))))
  cat(sprintf("omit emp %.4f thy %.4f\n", mean(is.infinite(Tk)), 1 - h))
  expect_lt(max(abs(emp - thy)), 3e-3)
  expect_lt(abs(mean(is.infinite(Tk)) - (1 - h)), 3e-3)
}

# ---------------------------------------------------------------------------
# 2. Sampler-based recovery
# ---------------------------------------------------------------------------

cat("== single-subject recovery ==\n")
set.seed(20260816)
matchfun <- function(d) d$S == d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
dat$subjects <- droplevels(dat$subjects)
ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
des <- design(data = dat, model = FRQ, matchfun = matchfun,
              formula = list(alpha ~ 1, beta ~ 1, h ~ lM, tau ~ lM, t0 ~ 1),
              contrasts = list(h = list(lM = ADmat), tau = list(lM = ADmat)))
p <- c(alpha = log(2), beta = log(3), h = qnorm(0.93), h_lMd = 0.9,
       tau = log(0.35), tau_lMd = -0.5, t0 = log(0.15))[names(sampled_pars(des))]

sim <- make_data(p, design = des, n_trials = 2000)
cat(sprintf("  simulated %d rows, omission rate %.4f\n",
            nrow(sim), mean(!is.finite(sim$rt))))

emc <- make_emc(sim, des, type = "single", n_chains = 3, rt_resolution = NULL)
emc <- fit(emc, cores_per_chain = 3, verbose = TRUE, fileName = NULL)

est <- summary(emc)
print(est)
cat("\n  truth:\n"); print(round(p, 3))

# The shapes (alpha, beta) are the weakly-informed coordinates -- see the
# identifiability note in ?FRQ -- so they get a looser band than h/tau/t0.
cat("\nRecovery check (credible intervals should cover the truth):\n")
ci <- credint(emc)
print(ci)
