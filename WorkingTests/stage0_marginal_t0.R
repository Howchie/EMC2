# ===========================================================================
# Stage 0: does marginalizing the shared subject-level t0 (option 3) remove the
# t0 <-> race-speed trade-off that breaks chain mixing in go/no-go race models?
#
# No sampler surgery. We build the model-agnostic marginalized likelihood
# (WorkingTests/marginal_t0_lib.R) and, on simulated GNG data with a finite
# deadline D, we:
#   (A) validate the t0 quadrature against a fine-grid reference;
#   (B) exhibit the pathology: the naive joint posterior has a sharp
#       t0 <-> B(go) ridge (near-collinear), the classic slow-mixing geometry;
#   (C) compare mixing head-to-head with a controlled adaptive RWMH, giving the
#       explicit-t0 and marginalized samplers the SAME t0 hyperprior p(t0|eta)
#       so the target over the SHARED parameters is identical -- only the nuisance
#       t0 dimension differs. We report ESS on the shared parameters;
#   (D) recover parameters (and a post-hoc t0 from the quadrature weights).
#
# Model is LBA standard-race GNG here; the machinery is model-agnostic (swap
# `model = LBA` for RDM etc.), which we spot-check at the end.
# ===========================================================================

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("WorkingTests/marginal_t0_lib.R")
set.seed(2026)

OUT <- "WorkingTests"
pdf(file.path(OUT, "stage0_marginal_t0.pdf"), width = 9, height = 7)
on.exit(dev.off(), add = TRUE)

## Effective sample size via Geyer's initial positive sequence estimator.
ess1 <- function(z) {
  z <- z[is.finite(z)]; n <- length(z); z <- z - mean(z)
  v0 <- sum(z * z) / n
  if (v0 == 0) return(n)
  ac <- as.numeric(acf(z, lag.max = min(n - 1L, 2000L), plot = FALSE, demean = FALSE)$acf)
  # pair sums Gamma_k = rho_{2k}+rho_{2k+1}; truncate at first non-positive
  g <- ac[seq(2, length(ac) - 1, by = 2)] + ac[seq(3, length(ac), by = 2)]
  k <- which(g <= 0)[1]; if (is.na(k)) k <- length(g)
  tau <- 1 + 2 * sum(ac[2:(2 * k)])
  n / max(tau, 1)
}

## --------------------------------------------------------------------------
## Build a GNG data set + fast likelihood closure.
## --------------------------------------------------------------------------
build <- function(model = LBA, n_trials = 250L, UC = 0.7, seed = 1, noise = "sv") {
  set.seed(seed)
  mf <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  # LBA's between-trial drift sd is "sv"; RDM's diffusion sd is "s".
  fml <- c(list(B ~ 0 + lR, v ~ 0 + mismatch + S:match, A ~ 1, t0 ~ 1),
           list(stats::as.formula(paste(noise, "~ 1"))))
  consts <- stats::setNames(c(log(1), log(0.4)), c(noise, "A"))
  des <- design(
    factors = list(subjects = 1, S = c("go","nogo","nogo","nogo")), Rlevels = c("go","nogo"),
    matchfun = mf,
    functions = list(match = function(d) ifelse(d$lM == TRUE, 1, 0),
                     mismatch = function(d) ifelse(d$lM == TRUE, 0, 1)),
    model = model, formula = fml, constants = consts)
  p <- sampled_pars(des, doMap = FALSE)
  p[["B_lRgo"]] <- log(.75); p[["B_lRnogo"]] <- log(.6); p[["t0"]] <- log(.2)
  p[["v_mismatch"]] <- .5; p[["v_Sgo:match"]] <- 2.1; p[["v_Snogo:match"]] <- 2.5
  dat <- make_data(p, des, n_trials = n_trials, TC = list(UC = UC))
  emc <- suppressMessages(make_emc(dat, des, type = "single", compress = TRUE, n_chains = 1))
  list(p = p, dat = dat, dadm = emc[[1]]$data[[1]],
       ll_fun = make_ll_fun(emc[[1]]$data[[1]], emc[[1]]$model),
       prop = matrix(p, nrow = 1, dimnames = list(NULL, names(p))))
}

cat("======================================================================\n")
cat("Stage 0: marginalizing subject-level t0 in a GNG race (LBA)\n")
cat("======================================================================\n")
B <- build(model = LBA, n_trials = 600L)
p <- B$p; ll_fun <- B$ll_fun; prop <- B$prop
free_shared <- c("B_lRgo","B_lRnogo","v_Sgo:match","v_Snogo:match","v_mismatch")
eta_mu <- log(0.20); eta_sd <- 0.30                      # t0 hyperprior p(t0|eta)
cat(sprintf("data: %d rows compressed; omission rate %.2f; min rt %.3f\n",
            nrow(B$dadm), mean(!is.finite(B$dat$rt) | is.na(B$dat$rt)),
            min(B$dat$rt[is.finite(B$dat$rt)])))

## (A) validation ------------------------------------------------------------
gl64 <- marginal_ll_t0(prop, ll_fun, eta_mu, eta_sd, n_nodes = 64L)
ref  <- grid_marginal_ll_t0(prop, ll_fun, eta_mu, eta_sd, n_grid = 2000L)
cat(sprintf("\n(A) quadrature check: GL(64)=%.5f  grid=%.5f  |diff|=%.2e\n", gl64, ref, abs(gl64 - ref)))

## (B) ridge geometry --------------------------------------------------------
t0g <- log(seq(0.12, 0.32, length.out = 41))
bg  <- p[["B_lRgo"]] + seq(-0.6, 0.6, length.out = 41)
M <- outer(t0g, bg, Vectorize(function(a, b) { q <- prop; q[,"t0"] <- a; q[,"B_lRgo"] <- b; ll_fun(q) }))
w <- exp(M - max(M)); hd <- which(w > 1e-4, arr.ind = TRUE)
ridge_cor <- cor(t0g[hd[,1]], bg[hd[,2]])
cat(sprintf("(B) naive joint (t0, B_go) high-density ridge correlation = %+.3f\n", ridge_cor))
image(exp(t0g), bg, w, xlab = "t0 (s)", ylab = "B_lRgo (log)",
      main = sprintf("Naive joint posterior slice: t0<->B_go ridge (r=%+.2f)", ridge_cor),
      col = hcl.colors(40, "YlOrRd", rev = TRUE))
contour(exp(t0g), bg, w, add = TRUE, drawlabels = FALSE, nlevels = 6)
points(exp(p[["t0"]]), p[["B_lRgo"]], pch = 3, lwd = 2, col = "blue")

## (C) ISOLATED t0<->B_go ridge: does carrying t0 slow B_go, and does
## integrating it out fix it? Others held at truth so the ONLY difference is
## whether the sampler carries the collinear t0 dimension. A block *isotropic*
## RWMH is used deliberately -- a sampler that does not know the ridge in
## advance, exactly the kind that suffers on it.
lp_prior_shared <- function(x) sum(dnorm(x[free_shared], 0, 3, log = TRUE))

# naive: sample (t0, B_go) jointly; marginal: sample B_go with t0 integrated.
lp_pair <- function(x) {                                 # x = c(t0, B_lRgo)
  q <- prop; q[, "t0"] <- x[1]; q[, "B_lRgo"] <- x[2]
  as.numeric(ll_fun(q)) + dnorm(x[1], eta_mu, eta_sd, log = TRUE) + dnorm(x[2], 0, 3, log = TRUE)
}
lp_bg_marg <- function(x) {                              # x = B_lRgo (t0 integrated)
  q <- prop; q[, "B_lRgo"] <- x
  as.numeric(marginal_ll_t0(q, ll_fun, eta_mu, eta_sd, n_nodes = 32L)) + dnorm(x, 0, 3, log = TRUE)
}

# (C1) deterministic argument: condition number of the naive (t0,B_go) posterior.
# RW-Metropolis mixing time scales with kappa = lambda_max/lambda_min of the
# precision, so a large kappa *is* the slow-mixing pathology, MCMC-noise-free.
mode2 <- c(p[["t0"]], p[["B_lRgo"]])
H <- matrix(0, 2, 2); h <- 1e-3
for (a in 1:2) for (b in 1:2) {
  e_a <- replace(numeric(2), a, h); e_b <- replace(numeric(2), b, h)
  H[a, b] <- (lp_pair(mode2 + e_a + e_b) - lp_pair(mode2 + e_a - e_b) -
              lp_pair(mode2 - e_a + e_b) + lp_pair(mode2 - e_a - e_b)) / (4 * h^2)
}
prec <- -H; ev <- sort(eigen(prec, symmetric = TRUE)$values)
kappa_naive <- ev[2] / ev[1]
# marginal is 1-D in B_go: curvature = -d2/dB^2 of the marginal log-post; kappa=1.
cm <- -(lp_bg_marg(p[["B_lRgo"]] + h) - 2 * lp_bg_marg(p[["B_lRgo"]]) + lp_bg_marg(p[["B_lRgo"]] - h)) / h^2
cat(sprintf("\n(C1) naive (t0,B_go) posterior condition number kappa = %.1f  (marginal B_go: kappa = 1)\n", kappa_naive))

# (C2) confirm with isolated block-isotropic RWMH; compare ESS(B_go).
block_rwmh <- function(lp, init, n_iter = 40000L, burn = 5000L, s0 = 0.08, target = 0.30) {
  d <- length(init); x <- init; l <- lp(x); s <- s0; acc <- 0
  D <- matrix(NA_real_, n_iter, d); t0 <- proc.time()[["elapsed"]]
  for (it in seq_len(n_iter)) {
    xp <- x + rnorm(d, 0, s); lpp <- lp(xp)
    if (is.finite(lpp) && log(runif(1)) < lpp - l) { x <- xp; l <- lpp; acc <- acc + 1 }
    D[it, ] <- x
    if (it <= burn && it %% 200 == 0) { s <- s * exp((acc/200 - target) * 0.7); acc <- 0 }
  }
  list(draws = D[(burn + 1L):n_iter, , drop = FALSE], sec = proc.time()[["elapsed"]] - t0)
}
cat("(C2) running isolated block-isotropic RWMH (naive 2-D vs marginal 1-D)...\n")
Rn <- block_rwmh(lp_pair,    c(p[["t0"]], p[["B_lRgo"]]))
Rm <- block_rwmh(lp_bg_marg, p[["B_lRgo"]])
ess_bg_naive <- ess1(Rn$draws[, 2]); ess_bg_marg <- ess1(Rm$draws[, 1])
cat(sprintf("   ESS(B_go): naive-2D = %.0f   marginal-1D = %.0f   (%.1fx)  | %d draws each\n",
            ess_bg_naive, ess_bg_marg, ess_bg_marg / ess_bg_naive, nrow(Rn$draws)))
cat(sprintf("   ESS/sec  : naive = %.1f    marginal = %.1f    (%.1fx)  | naive %.1fs, marg %.1fs\n",
            ess_bg_naive / Rn$sec, ess_bg_marg / Rm$sec,
            (ess_bg_marg / Rm$sec) / (ess_bg_naive / Rn$sec), Rn$sec, Rm$sec))
cat(sprintf("   ESS(t0) in naive chain (the sticky nuisance dimension): %.0f\n", ess1(Rn$draws[, 1])))

par(mfrow = c(2, 2))
plot(Rn$draws[, 2], type = "l", col = "firebrick", main = sprintf("naive B_go (ESS=%.0f)", ess_bg_naive),
     ylab = "B_lRgo", xlab = "iter"); abline(h = p[["B_lRgo"]], col = "blue", lwd = 2)
plot(Rm$draws[, 1], type = "l", col = "seagreen", main = sprintf("marginal B_go (ESS=%.0f)", ess_bg_marg),
     ylab = "B_lRgo", xlab = "iter"); abline(h = p[["B_lRgo"]], col = "blue", lwd = 2)
plot(exp(Rn$draws[, 1]), Rn$draws[, 2], pch = ".", col = "firebrick",
     main = sprintf("naive ridge (kappa=%.0f)", kappa_naive), xlab = "t0 (s)", ylab = "B_lRgo")
points(exp(p[["t0"]]), p[["B_lRgo"]], pch = 3, lwd = 2, col = "blue")
image(exp(t0g), bg, w, xlab = "t0 (s)", ylab = "B_lRgo", main = "posterior slice",
      col = hcl.colors(40, "YlOrRd", rev = TRUE)); contour(exp(t0g), bg, w, add = TRUE, drawlabels = FALSE)
par(mfrow = c(1, 1))

## (D) recovery via direct optimization of the marginal objective + post-hoc t0
neg_marg <- function(x) {
  v <- setNames(x, free_shared); q <- prop; q[, free_shared] <- v
  -(as.numeric(marginal_ll_t0(q, ll_fun, eta_mu, eta_sd, n_nodes = 32L)) + lp_prior_shared(v))
}
set.seed(11)
opt <- optim(unlist(p[free_shared]) + rnorm(length(free_shared), 0, 0.3), neg_marg,
             method = "Nelder-Mead", control = list(maxit = 1200, reltol = 1e-8))
cat("\n(D) recovery of shared params by maximizing the marginal likelihood:\n")
rec <- data.frame(param = free_shared, true = round(unlist(p[free_shared]), 3),
                  est = round(opt$par, 3))
print(rec, row.names = FALSE)
# post-hoc t0 (writeup's step): quadrature-weighted E[t0|data] at the recovered mode
xm <- prop; xm[, free_shared] <- opt$par
gl <- gl_rule(64L); lo <- eta_mu - 6 * eta_sd; hi <- min(eta_mu + 6 * eta_sd, attr(ll_fun, "log_ms"))
xs <- (hi + lo)/2 + (hi - lo)/2 * gl$x
lw <- log(gl$w) + dnorm(xs, eta_mu, eta_sd, log = TRUE) +
      sapply(xs, function(z) { q <- xm; q[, "t0"] <- z; ll_fun(q) })
wq <- exp(lw - max(lw)); wq <- wq / sum(wq)
cat(sprintf("post-hoc E[t0 | data] = %.3f s  (true %.3f)\n", sum(exp(xs) * wq), exp(p[["t0"]])))

## model-agnostic spot check (RDM) ------------------------------------------
cat("\n(model-agnostic check) same wrapper on an RDM GNG data set:\n")
Br <- build(model = RDM, seed = 7, noise = "s")
glr <- marginal_ll_t0(Br$prop, Br$ll_fun, eta_mu, eta_sd, n_nodes = 64L)
refr <- grid_marginal_ll_t0(Br$prop, Br$ll_fun, eta_mu, eta_sd, n_grid = 2000L)
cat(sprintf("  RDM: GL(64)=%.5f  grid=%.5f  |diff|=%.2e\n", glr, refr, abs(glr - refr)))

cat("\nWrote figure: ", file.path(OUT, "stage0_marginal_t0.pdf"), "\n")
cat("STAGE0 DONE\n")
