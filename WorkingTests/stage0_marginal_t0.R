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
# The primary benchmark below is an LBA GNG/no-go model (deadline-censored
# data).  The machinery is model-agnostic; the spot checks at the end include
# both a GNG/no-go RDM and a standard two-choice RDMSWTN with sv > 0.
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
build <- function(model = LBA, n_trials = 250L, UC = 0.7, seed = 1,
                  noise = "sv", sample_sv = FALSE) {
  set.seed(seed)
  mf <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  # LBA's between-trial drift sd is "sv"; RDM's diffusion sd is "s".
  fml <- c(list(B ~ 0 + lR, v ~ 0 + mismatch + S:match, A ~ 1, t0 ~ 1),
           list(stats::as.formula(paste(noise, "~ 1"))))
  consts <- stats::setNames(c(log(1), log(0.4)), c(noise, "A"))
  if (isTRUE(sample_sv) && identical(noise, "sv")) {
    consts <- consts[names(consts) != "sv"]
  }
  des <- design(
    factors = list(subjects = 1, S = c("go","nogo","nogo","nogo")), Rlevels = c("go","nogo"),
    matchfun = mf,
    functions = list(match = function(d) ifelse(d$lM == TRUE, 1, 0),
                     mismatch = function(d) ifelse(d$lM == TRUE, 0, 1)),
    model = model, formula = fml, constants = consts)
  p <- sampled_pars(des, doMap = FALSE)
  p[["B_lRgo"]] <- log(.75); p[["B_lRnogo"]] <- log(.6); p[["t0"]] <- log(.2)
  p[["v_mismatch"]] <- .5; p[["v_Sgo:match"]] <- 2.1; p[["v_Snogo:match"]] <- 2.5
  # RDMSWTN exposes sv as a sampled parameter; keep this strictly positive in
  # that check.  For LBA, sv is fixed by `constants` and is not in `p`.
  if ("sv" %in% names(p)) p[["sv"]] <- log(.35)
  dat <- make_data(p, des, n_trials = n_trials, TC = list(UC = UC))
  emc <- suppressMessages(make_emc(dat, des, type = "single", compress = TRUE, n_chains = 1))
  list(p = p, dat = dat, dadm = emc[[1]]$data[[1]], model = emc[[1]]$model,
       ll_fun = make_ll_fun(emc[[1]]$data[[1]], emc[[1]]$model),
       prop = matrix(p, nrow = 1, dimnames = list(NULL, names(p))))
}

# A standard two-choice race (no `nogo` level and no deadline) used to ensure
# the node-count check is not accidentally relying on GNG censoring semantics.
build_standard_rdmswtn <- function(n_trials = 600L, seed = 17) {
  set.seed(seed)
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  des <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"), matchfun = matchfun,
    functions = list(match = function(d) ifelse(d$lM == TRUE, 1, 0),
                     mismatch = function(d) ifelse(d$lM == TRUE, 0, 1)),
    model = RDMSWTN,
    formula = list(B ~ 0 + lR, v ~ 0 + lM, A ~ 1, t0 ~ 1,
                   s ~ 1, sv ~ 1),
    constants = c(A = log(.2), s = log(1)))
  p <- sampled_pars(des, doMap = FALSE)
  p[["B_lRleft"]] <- log(.8); p[["B_lRright"]] <- log(1.0)
  p[["v_lMFALSE"]] <- log(1.2); p[["v_lMTRUE"]] <- log(2.0)
  p[["t0"]] <- log(.2); p[["sv"]] <- log(.35)
  dat <- make_data(p, des, n_trials = n_trials)
  emc <- suppressMessages(make_emc(dat, des, type = "single",
                                   compress = TRUE, n_chains = 1))
  list(p = p, dat = dat, dadm = emc[[1]]$data[[1]],
       model = emc[[1]]$model,
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

## (A) validation and node-count benchmark ----------------------------------
# Use the registered C++ marginal likelihood for the node-count comparison.
# Its adaptive route first uses a 7-node pilot to locate the common
# response-informed t0 mass, then spends n_nodes in a composite GL rule with
# a dense central panel plus explicit tail panels.  The older R grid remains as
# an independent check, but is explicitly clipped at the model's t0 lower bound
# so it cannot spend nodes in the min_ll dead zone.
make_cpp_marginal <- function(B) {
  model <- B$model()
  dadm <- EMC2:::.cache_ll_data_attrs(B$dadm)
  constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
  designs <- EMC2:::.oo_expanded_designs(dadm)
  marginal <- function(proposals, n_nodes = NULL, mu = eta_mu, sigma = eta_sd,
                       adaptive = TRUE, n_scan = 7L) {
    if (is.null(dim(proposals))) {
      proposals <- matrix(proposals, nrow = 1L,
                          dimnames = list(NULL, names(proposals)))
    }
    marginalise <- list(param = "t0", mu = mu, sigma = sigma,
                         adaptive = adaptive, n_scan = as.integer(n_scan))
    if (!is.null(n_nodes)) marginalise$n_nodes <- as.integer(n_nodes)
    EMC2:::calc_ll_oo(
      proposals, dadm, constants = constants, designs = designs,
      type = model$c_name, bounds = model$bound, transforms = model$transform,
      pretransforms = model$pre_transform, p_types = names(model$p_types),
      min_ll = log(1e-10), trend = model$trend,
      marginalise = marginalise
    )
  }
  nodes <- function(proposals, n_nodes = 40L, mu = eta_mu, sigma = eta_sd,
                    adaptive = TRUE, n_scan = 7L) {
    if (is.null(dim(proposals))) {
      proposals <- matrix(proposals, nrow = 1L,
                          dimnames = list(NULL, names(proposals)))
    }
    EMC2:::calc_ll_oo_marginal_nodes(
      proposals, dadm, constants = constants, designs = designs,
      type = model$c_name, bounds = model$bound, transforms = model$transform,
      pretransforms = model$pre_transform, p_types = names(model$p_types),
      min_ll = log(1e-10), trend = model$trend,
      marginalise = list(param = "t0", mu = mu, sigma = sigma,
                         n_nodes = as.integer(n_nodes), adaptive = adaptive,
                         n_scan = as.integer(n_scan))
    )
  }
  list(model = model, marginal = marginal, nodes = nodes)
}
cpp <- make_cpp_marginal(B)
t0_lower <- cpp$model$bound$minmax[1, "t0"]
production_nodes <- 40L
node_sizes <- c(6L, 8L, 10L, 12L, 16L, 20L, 24L, 28L, 32L, 40L)
reference_nodes <- 256L
ll_reference <- as.numeric(cpp$marginal(prop, reference_nodes, adaptive = FALSE))
ll_by_nodes <- vapply(node_sizes, function(k) as.numeric(cpp$marginal(prop, k)), numeric(1L))
ll_default <- as.numeric(cpp$marginal(prop))

# Reproduce the old interval only to quantify the bug this benchmark guards:
# with t0 ~ N(log(.2), .3^2), 7/25 (28%) of its nodes are below the .05 bound.
legacy_lo <- eta_mu - 6 * eta_sd
log_ms <- log_ms_from_dadm(B$dadm)
legacy_nodes <- gl_rule(25L)
legacy_nodes <- (min(eta_mu + 6 * eta_sd, log_ms) + legacy_lo) / 2 +
  (min(eta_mu + 6 * eta_sd, log_ms) - legacy_lo) / 2 * legacy_nodes$x
cat(sprintf("\n(A) legacy interval: %d/%d nodes (%.1f%%) below t0 lower bound %.3f\n",
            sum(exp(legacy_nodes) < t0_lower), length(legacy_nodes),
            100 * mean(exp(legacy_nodes) < t0_lower), t0_lower))
stopifnot(sum(exp(legacy_nodes) < t0_lower) == 7L)

cpp_nodes <- cpp$nodes(prop, production_nodes)
stopifnot(all(exp(cpp_nodes$nodes) >= t0_lower))
node_table <- data.frame(
  final_nodes = node_sizes,
  kernel_calls = node_sizes + 7L,
  log_likelihood = ll_by_nodes,
  abs_error_vs_256 = abs(ll_by_nodes - ll_reference)
)
print(node_table, row.names = FALSE, digits = 8)
ref_stability <- abs(as.numeric(cpp$marginal(prop, 128L, adaptive = FALSE)) - ll_reference)
cat(sprintf("Fixed C++ GL(128) vs GL(256) stability: %.2e log units\n", ref_stability))
cat(sprintf("C++ omitted-n_nodes default agrees with adaptive GL(%d): %.2e log units\n",
            production_nodes, abs(ll_default - ll_by_nodes[node_sizes == production_nodes])))
stopifnot(ref_stability < 1e-6,
          abs(ll_default - ll_by_nodes[node_sizes == production_nodes]) < 1e-12,
          node_table$abs_error_vs_256[node_table$final_nodes == production_nodes] < 1e-3)

# Timing is measured on a modest particle batch to expose the linear node cost
# without making this scratch test depend on wall-clock precision for one call.
bench_particles <- prop[rep(1L, 32L), , drop = FALSE]
bench_cpp <- function(n_nodes, adaptive = TRUE, reps = 20L) {
  elapsed <- system.time(for (i in seq_len(reps)) {
    cpp$marginal(bench_particles, n_nodes, adaptive = adaptive)
  })[["elapsed"]]
  1000 * elapsed / reps
}
timing_ms <- vapply(node_sizes, bench_cpp, numeric(1L))
node_table$milliseconds_per_32_particles <- timing_ms
fixed80_ms <- bench_cpp(80L, adaptive = FALSE)
cat("\nC++ timing (milliseconds per 32-particle marginal call):\n")
print(node_table[, c("final_nodes", "kernel_calls", "abs_error_vs_256",
                     "milliseconds_per_32_particles")],
      row.names = FALSE, digits = 6)
cat(sprintf("Fixed GL(80) baseline: %.2f ms per 32-particle call\n", fixed80_ms))

gl40 <- as.numeric(cpp$marginal(prop, production_nodes))
ref_grid <- grid_marginal_ll_t0(
  prop, ll_fun, eta_mu, eta_sd, n_grid = 2000L,
  log_lo = log(t0_lower), log_ms = log_ms
)
cat(sprintf("\nC++ adaptive GL(%d)=%.5f  clipped R grid=%.5f  |diff|=%.2e\n",
            production_nodes, gl40, ref_grid, abs(gl40 - ref_grid)))

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
  as.numeric(cpp$marginal(q, n_nodes = production_nodes)) + dnorm(x, 0, 3, log = TRUE)
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
  -(as.numeric(cpp$marginal(q, n_nodes = production_nodes)) + lp_prior_shared(v))
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
nd <- cpp$nodes(xm, production_nodes)
lw <- nd$log_terms[1L, ]
wq <- exp(lw - max(lw)); wq <- wq / sum(wq)
cat(sprintf("post-hoc E[t0 | data] = %.3f s  (true %.3f)\n",
            sum(exp(nd$nodes) * wq), exp(p[["t0"]])))

## model-agnostic spot checks -------------------------------------------------
cat("\n(model-agnostic check) RDM GNG/no-go data set:\n")
Br <- build(model = RDM, seed = 7, noise = "s")
cppr <- make_cpp_marginal(Br)
rdm_nodes <- c(8L, 12L, 16L, 20L, 24L, 28L, 32L, 40L)
rdm_reference <- as.numeric(cppr$marginal(Br$prop, n_nodes = 256L, adaptive = FALSE))
rdm_ll <- vapply(rdm_nodes, function(k) as.numeric(cppr$marginal(Br$prop, n_nodes = k)), numeric(1L))
cat("  RDM node convergence (C++ log-likelihood error vs GL(256)):\n")
print(data.frame(nodes = rdm_nodes, abs_error_vs_256 = abs(rdm_ll - rdm_reference)),
      row.names = FALSE, digits = 8)
stopifnot(abs(rdm_ll[rdm_nodes == production_nodes] - rdm_reference) < 1e-3)
glr <- as.numeric(cppr$marginal(Br$prop, n_nodes = production_nodes))
refr <- grid_marginal_ll_t0(
  Br$prop, Br$ll_fun, eta_mu, eta_sd, n_grid = 2000L,
  log_lo = log(cppr$model$bound$minmax[1, "t0"]),
  log_ms = log_ms_from_dadm(Br$dadm)
)
cat(sprintf("  RDM: C++ adaptive GL(%d)=%.5f  clipped grid=%.5f  |diff|=%.2e\n",
            production_nodes,
            glr, refr, abs(glr - refr)))

cat("\n(model-agnostic check) RDMSWTN GNG/no-go with sv > 0:\n")
Bswtn_gng <- build(model = RDMSWTN, seed = 9, noise = "sv", sample_sv = TRUE)
stopifnot(exp(Bswtn_gng$p[["sv"]]) > 0,
          any(as.character(Bswtn_gng$dat$S) == "nogo"))
cppswtn_gng <- make_cpp_marginal(Bswtn_gng)
swtn_gng_nodes <- c(16L, 20L, 24L, 28L, 32L, 40L)
swtn_gng_reference <- as.numeric(cppswtn_gng$marginal(
  Bswtn_gng$prop, n_nodes = 256L, adaptive = FALSE))
swtn_gng_ll <- vapply(
  swtn_gng_nodes,
  function(k) as.numeric(cppswtn_gng$marginal(Bswtn_gng$prop, n_nodes = k)),
  numeric(1L))
cat(sprintf("  RDMSWTN GNG sv=%.2f node convergence (error vs fixed GL(256)):\n",
            exp(Bswtn_gng$p[["sv"]])))
print(data.frame(nodes = swtn_gng_nodes,
                 abs_error_vs_256 = abs(swtn_gng_ll - swtn_gng_reference)),
      row.names = FALSE, digits = 8)
stopifnot(abs(swtn_gng_ll[swtn_gng_nodes == production_nodes] - swtn_gng_reference) < 1e-3)

cat("\n(model-agnostic check) standard race: RDMSWTN with sv > 0:\n")
Bswtn <- build_standard_rdmswtn()
stopifnot(exp(Bswtn$p[["sv"]]) > 0,
          !any(as.character(Bswtn$dat$S) == "nogo"))
cppswtn <- make_cpp_marginal(Bswtn)
swtn_nodes <- c(12L, 16L, 20L, 24L, 28L, 32L, 40L)
swtn_reference <- as.numeric(cppswtn$marginal(Bswtn$prop, n_nodes = 256L,
                                               adaptive = FALSE))
swtn_ll <- vapply(swtn_nodes,
                  function(k) as.numeric(cppswtn$marginal(Bswtn$prop,
                                                          n_nodes = k)),
                  numeric(1L))
cat(sprintf("  RDMSWTN sv=%.2f node convergence (C++ log-likelihood error vs fixed GL(256)):\n",
            exp(Bswtn$p[["sv"]])))
print(data.frame(nodes = swtn_nodes,
                 abs_error_vs_256 = abs(swtn_ll - swtn_reference)),
      row.names = FALSE, digits = 8)
stopifnot(abs(swtn_ll[swtn_nodes == production_nodes] - swtn_reference) < 1e-3)
glswtn <- as.numeric(cppswtn$marginal(Bswtn$prop, n_nodes = production_nodes))
refswtn <- grid_marginal_ll_t0(
  Bswtn$prop, Bswtn$ll_fun, eta_mu, eta_sd, n_grid = 2000L,
  log_lo = log(cppswtn$model$bound$minmax[1, "t0"]),
  log_ms = log_ms_from_dadm(Bswtn$dadm)
)
cat(sprintf("  RDMSWTN: C++ adaptive GL(%d)=%.5f  clipped grid=%.5f  |diff|=%.2e\n",
            production_nodes,
            glswtn, refswtn, abs(glswtn - refswtn)))

cat("\nWrote figure: ", file.path(OUT, "stage0_marginal_t0.pdf"), "\n")
cat("STAGE0 DONE\n")
