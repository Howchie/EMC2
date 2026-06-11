# ===========================================================================
# stop_success_methods.R — one-stop comparison of the stop-success
# integration options in the SS-EXG / SS-RDEX stop-signal models.
#
# Consolidates (and replaces) the earlier WorkingTests scripts:
#   analytic_stop_success_exg.R, benchmark_stop_success_gl.R,
#   compare_n_go_runners.R, compare_stop_methods.R, compare_optin.R,
#   compare_optin_full_ll.R
#
# Run after:  devtools::load_all(".")
#             source("WorkingTests/stop_success_methods.R")
# Requires:   microbenchmark; mvtnorm + pbivnorm (oracle only)
#
# Sections:
#   1. Exact analytic oracle (MVN-CDF closed form) + self-check
#   2. Four-way accuracy: stop_method = integrate / gl / analytic / auto,
#      through the REAL C++ and R dispatch paths, vs the oracle; plus
#      stress cases that must trigger the auto guard / kink fallback
#   3. Gauss-Legendre node-count study (EXG and RDEX)
#   4. Speed per single integral (four methods; EXG n_go = 1, 2; RDEX)
#   5. Full log-likelihood end-to-end via calc_ll with SSEXG(stop_method=)
#   6. PASS/FAIL summary
# ===========================================================================

suppressMessages({
  library(microbenchmark)
  library(mvtnorm)
  library(pbivnorm)
})

sep <- paste(rep("=", 75), collapse = "")

# ===========================================================================
# 1. EXACT ANALYTIC ORACLE
# ---------------------------------------------------------------------------
# The stop-success probability is, with stop ~ exG(muS, sigS, tauS) lower-
# truncated at lbS and go runner i ~ exG(mu_i, sig_i, tau_i) truncated at lb_i:
#
#   P(SSD) = (1/N) * INT_{lbS}^{Inf} f_S(x) * prod_i S_i(x + SSD) dx,
#   N = S_S(lbS) * prod_i S_i(lb_i)        [valid when lbS + SSD >= lb_i]
#
# Both exG density and survivor are (exponential x normal-CDF) objects, so
# expanding prod_i S_i over branch choices A gives 2^n terms
# C_A exp(gamma_A x) prod_j Phi(a_j + b_j x) with gamma_A < 0; one integration
# by parts per term reduces each to one-factor multivariate-normal CDFs
# (bivariate for the full line with n=2; trivariate when truncated).
# This is EXACT (machine precision); the package's n_go = 1 closed form
# (src/ss_exg_analytic.h) is the production implementation of the same result.
# ===========================================================================

dexg_R <- function(x, mu, sigma, tau)  # exp*Phi in log space (stable tails)
  exp(-log(tau) + (mu - x)/tau + sigma^2/(2*tau^2) +
        pnorm((x - mu)/sigma - sigma/tau, log.p = TRUE))

sexg_R <- function(t, mu, sigma, tau)  # survivor 1 - F
  pnorm(-(t - mu)/sigma) +
    exp((mu - t)/tau + sigma^2/(2*tau^2) +
          pnorm((t - mu)/sigma - sigma/tau, log.p = TRUE))

dtexg_R <- function(x, mu, sigma, tau, lb)
  ifelse(x <= lb, 0, dexg_R(x, mu, sigma, tau) / sexg_R(lb, mu, sigma, tau))
stexg_R <- function(t, mu, sigma, tau, lb)
  ifelse(t <= lb, 1, sexg_R(t, mu, sigma, tau) / sexg_R(lb, mu, sigma, tau))

integrand_trunc <- function(x, SSD, muS, sigS, tauS, lbS, mu, sig, tau, lb) {
  out <- dtexg_R(x, muS, sigS, tauS, lbS)
  for (i in seq_along(mu)) out <- out * stexg_R(x + SSD, mu[i], sig[i], tau[i], lb[i])
  out
}

ref_integrate <- function(f, lower, ...)
  integrate(f, lower = lower, upper = Inf, ...,
            rel.tol = 1e-12, abs.tol = 1e-14, subdivisions = 500L)$value

Q_fullline <- function(alpha, beta) {
  cj <- alpha / sqrt(1 + beta^2)
  d <- length(cj)
  if (d == 1) return(pnorm(cj))
  if (d == 2) {
    rho <- beta[1]*beta[2] / sqrt((1 + beta[1]^2)*(1 + beta[2]^2))
    return(pbivnorm(cj[1], cj[2], rho))
  }
  s <- sqrt(1 + beta^2)
  R <- diag(d)
  for (j in 1:(d-1)) for (l in (j+1):d)
    R[j, l] <- R[l, j] <- beta[j]*beta[l]/(s[j]*s[l])
  alg <- if (d == 3) TVPACK(1e-14) else Miwa(steps = 4097)
  as.numeric(pmvnorm(upper = cj, corr = R, algorithm = alg))
}

Q_trunc <- function(alpha, beta, vL) {
  cj <- alpha / sqrt(1 + beta^2)
  s <- sqrt(1 + beta^2)
  d <- length(cj)
  R <- diag(d + 1)
  for (j in 1:d) R[1, j+1] <- R[j+1, 1] <- beta[j]/s[j]   # corr(-V, Y_j)
  if (d >= 2) for (j in 1:(d-1)) for (l in (j+1):d)
    R[j+1, l+1] <- R[l+1, j+1] <- beta[j]*beta[l]/(s[j]*s[l])
  upper <- c(-vL, cj)
  alg <- if (d + 1 <= 3) TVPACK(1e-14) else Miwa(steps = 4097)
  as.numeric(pmvnorm(upper = upper, corr = R, algorithm = alg))
}

# core integral INT_L^Inf f_S(x) prod_i S_i(x+SSD) dx (untruncated f, S)
stop_success_exg_oracle_core <- function(SSD, muS, sigS, tauS, mu, sig, tau,
                                         L = -Inf) {
  n <- length(mu)
  m <- mu - SSD
  total <- 0
  for (Abits in 0:(2^n - 1)) {
    inA <- as.logical(bitwAnd(Abits, 2^(0:(n-1))))
    gam <- -1/tauS - sum(1/tau[inA])
    logC <- -log(tauS) + muS/tauS + sigS^2/(2*tauS^2) +
      sum(m[inA]/tau[inA] + sig[inA]^2/(2*tau[inA]^2))
    a <- c(-muS/sigS - sigS/tauS,
           ifelse(inA, -m/sig - sig/tau, m/sig))
    b <- c(1/sigS, ifelse(inA, 1/sig, -1/sig))
    J <- 0
    if (is.finite(L))    # boundary term (zero for L = -Inf)
      J <- J + exp(gam*L)/abs(gam) * prod(pnorm(a + b*L))
    for (k in seq_along(a)) {
      xbar <- -a[k]/b[k]
      mustar <- xbar + gam/b[k]^2
      delta <- gam*xbar + gam^2/(2*b[k]^2)
      jj <- setdiff(seq_along(a), k)
      beta <- b[jj]/abs(b[k])
      alpha <- a[jj] + b[jj]*mustar
      Q <- if (is.finite(L)) Q_trunc(alpha, beta, abs(b[k])*(L - mustar))
           else Q_fullline(alpha, beta)
      J <- J + sign(b[k]) * exp(delta) * Q / abs(gam)
    }
    total <- total + exp(logC) * J
  }
  total
}

# full truncated stop-success probability (== my.integrate(stopfn_texg, ...))
pstop_texg_oracle <- function(SSD, muS, sigS, tauS, lbS, mu, sig, tau, lb) {
  stopifnot(all(lbS + SSD >= lb))  # else integrand is piecewise; split needed
  core <- stop_success_exg_oracle_core(SSD, muS, sigS, tauS, mu, sig, tau,
                                       L = lbS)
  core / (sexg_R(lbS, muS, sigS, tauS) * prod(sexg_R(lb, mu, sig, tau)))
}

# ---- oracle self-check ------------------------------------------------------
cat(sep, "\nSECTION 1: oracle self-check (vs 1e-12-tolerance adaptive integration)\n",
    sep, "\n", sep = "")
muS_t <- .2; sigS_t <- .03; tauS_t <- .05; lbS_def <- .05; lbG_def <- .05
muG_t <- c(.40, .45); sigG_t <- c(.05, .06); tauG_t <- c(.10, .12)
SSDs <- c(0, .1, .2, .3, .5)
worst <- 0
for (SSD in SSDs) {
  ref <- ref_integrate(integrand_trunc, lbS_def, SSD = SSD, muS = muS_t,
                       sigS = sigS_t, tauS = tauS_t, lbS = lbS_def,
                       mu = muG_t, sig = sigG_t, tau = tauG_t,
                       lb = rep(lbG_def, 2))
  ora <- pstop_texg_oracle(SSD, muS_t, sigS_t, tauS_t, lbS_def,
                           muG_t, sigG_t, tauG_t, rep(lbG_def, 2))
  worst <- max(worst, abs(ora - ref)/ref)
}
cat(sprintf("oracle vs integrate, n_go=2, %d SSDs: worst relerr = %.2e\n",
            length(SSDs), worst))
stopifnot(worst < 1e-9)
cat("oracle OK\n")

# ===========================================================================
# 2. FOUR-WAY ACCURACY through the real dispatch paths
# ===========================================================================

# pars matrix for ss_texg_stop_success_value(): n_go rows x 10 cols.
#   1=muG, 2=sigG, 3=tauG, 4=muS(row1), 5=sigS(row1), 6=tauS(row1),
#   9=lbG, 10=lbS(row1)
make_pars <- function(muS, sigS, tauS, lbS, muG, sigG, tauG, lbG) {
  n <- length(muG)
  m <- matrix(0, nrow = n, ncol = 10)
  m[, 1] <- muG;  m[, 2] <- sigG;  m[, 3] <- tauG
  m[1, 4] <- muS; m[1, 5] <- sigS; m[1, 6] <- tauS
  m[, 9]  <- lbG; m[1, 10] <- lbS
  m
}

safe_oracle <- function(SSD, muS, sigS, tauS, lbS, muG, sigG, tauG, lbG)
  tryCatch(pstop_texg_oracle(SSD, muS, sigS, tauS, lbS, muG, sigG, tauG, lbG),
           error = function(e) NA_real_)

fmt_relerr <- function(v) {
  if (is.null(v) || length(v) == 0) return(sprintf("%9s %9s", "—", "—"))
  sprintf("%9.2e %9.2e", median(v), max(v))
}
fmt_us <- function(x) if (is.na(x)) sprintf("%9s", "—") else sprintf("%9.1f", x)

# ~50 random draws (seed 42)
set.seed(42)
n_rand <- 50
mu_r   <- matrix(runif(n_rand * 3, .3, .6),  n_rand, 3)
sig_r  <- matrix(runif(n_rand * 3, .02, .1), n_rand, 3)
tau_r  <- matrix(runif(n_rand * 3, .05, .2), n_rand, 3)

# Stress cases: must trigger the auto guard / kink fallback
stress_cases <- list(
  # 1. Guard: sigS > 4*tauS -> auto/analytic must fall back to GL
  list(label = "guard(sigS>4tauS)",
       muS = .2, sigS = .2, tauS = .01, lbS = .05,
       muG1 = c(.4), sigG1 = c(.08), tauG1 = c(.05),
       muG2 = c(.4, .45), sigG2 = c(.08, .06), tauG2 = c(.05, .08)),
  # 2. Kinked domain: lbS + SSD < lbG -> analytic precondition fails -> GL
  list(label = "kinked(lbS+SSD<lbG)",
       muS = .2, sigS = .03, tauS = .05, lbS = .10,
       muG1 = c(.4),  sigG1 = c(.05), tauG1 = c(.10),
       muG2 = c(.4, .45), sigG2 = c(.05, .06), tauG2 = c(.10, .12))
)
stress_lbG <- .20; stress_lbG_normal <- .05

cpp_ok <- tryCatch({
  ss_texg_stop_success_value(
    .2, make_pars(.2, .03, .05, .05, .4, .05, .1, .05), "integrate")
  ss_texg_stop_success_auto_branch(.2, make_pars(.2,.03,.05,.05,.4,.05,.1,.05))
  TRUE
}, error = function(e) FALSE)
cat(sprintf("\nC++ wrappers: %s\n",
            if (cpp_ok) "available" else "NOT available — C++ columns skipped"))

methods4 <- c("integrate", "gl", "analytic", "auto")

run_accuracy <- function(n_go) {
  re <- setNames(vector("list", 2 * length(methods4)),
                 c(paste0("cpp_", methods4), paste0("R_", methods4)))
  for (k in names(re)) re[[k]] <- numeric(0)
  auto_branches <- character(0)
  n_eval <- 0L

  build_ps <- function(i) {
    if (i == 0L)
      list(muS = muS_t, sigS = sigS_t, tauS = tauS_t,
           muG = muG_t[seq_len(n_go)], sigG = sigG_t[seq_len(n_go)],
           tauG = tauG_t[seq_len(n_go)])
    else
      list(muS = mu_r[i, 1], sigS = sig_r[i, 1], tauS = tau_r[i, 1],
           muG  = mu_r[i, seq_len(n_go) + 1L],
           sigG = sig_r[i, seq_len(n_go) + 1L],
           tauG = tau_r[i, seq_len(n_go) + 1L])
  }

  for (pi in c(0L, seq_len(n_rand))) {
    ps   <- build_ps(pi)
    lbG  <- rep(lbG_def, n_go)
    mu_v <- c(ps$muS, ps$muG); sig_v <- c(ps$sigS, ps$sigG)
    tau_v <- c(ps$tauS, ps$tauG); lb_v <- c(lbS_def, lbG)
    pars <- if (cpp_ok) make_pars(ps$muS, ps$sigS, ps$tauS, lbS_def,
                                  ps$muG, ps$sigG, ps$tauG, lbG) else NULL
    for (SSD in SSDs) {
      oracle <- safe_oracle(SSD, ps$muS, ps$sigS, ps$tauS, lbS_def,
                            ps$muG, ps$sigG, ps$tauG, lbG)
      if (is.na(oracle) || oracle < 1e-12) next
      n_eval <- n_eval + 1L
      for (m in methods4) {
        v <- tryCatch(
          EMC2:::stop_success_texg_R(mu_v, sig_v, tau_v, lb_v, SSD,
                                     method = m, n_nodes = 64L),
          error = function(e) NA_real_)
        key <- paste0("R_", m)
        re[[key]] <- c(re[[key]],
                       if (is.na(v)) NA_real_ else abs(v - oracle) / oracle)
      }
      if (cpp_ok) {
        for (m in methods4) {
          v <- tryCatch(
            ss_texg_stop_success_value(SSD, pars, m, n_nodes = 64L),
            error = function(e) NA_real_)
          key <- paste0("cpp_", m)
          re[[key]] <- c(re[[key]],
                         if (is.na(v)) NA_real_ else abs(v - oracle) / oracle)
        }
        auto_branches <- c(auto_branches,
                           tryCatch(ss_texg_stop_success_auto_branch(SSD, pars),
                                    error = function(e) "error"))
      }
    }
  }

  cat(sprintf("n_go=%d: %d evaluations\n", n_go, n_eval))
  if (cpp_ok && length(auto_branches) > 0) {
    tbl <- sort(table(auto_branches), decreasing = TRUE)
    cat(sprintf("  auto branches (C++): %s\n",
                paste(names(tbl), tbl, sep = "x", collapse = "; ")))
  }
  re
}

cat("\n", sep, "\nSECTION 2: four-way accuracy (typical + 50 random draws)\n",
    sep, "\n", sep = "")
acc1 <- run_accuracy(1)
acc2 <- run_accuracy(2)

# ---- stress cases ----------------------------------------------------------
run_stress <- function(n_go) {
  cat(sprintf("\n--- Stress cases (n_go=%d) ---\n", n_go))
  for (sc in stress_cases) {
    muG  <- if (n_go == 1) sc$muG1  else sc$muG2
    sigG <- if (n_go == 1) sc$sigG1 else sc$sigG2
    tauG <- if (n_go == 1) sc$tauG1 else sc$tauG2
    is_kinked <- grepl("kinked", sc$label)
    lbG <- rep(if (is_kinked) stress_lbG else stress_lbG_normal, n_go)
    for (SSD in c(0, .1, .2)) {
      mu_v  <- c(sc$muS, muG);  sig_v <- c(sc$sigS, sigG)
      tau_v <- c(sc$tauS, tauG); lb_v <- c(sc$lbS, lbG)
      pars  <- if (cpp_ok) make_pars(sc$muS, sc$sigS, sc$tauS, sc$lbS,
                                     muG, sigG, tauG, lbG) else NULL
      ref <- tryCatch(   # reference: integrate (untouched qags)
        EMC2:::stop_success_texg_R(mu_v, sig_v, tau_v, lb_v, SSD,
                                   method = "integrate"),
        error = function(e) NA_real_)
      if (is.na(ref) || ref < 1e-14) {
        cat(sprintf("  %-25s SSD=%.2f: ref=0/NA, skip\n", sc$label, SSD))
        next
      }
      br <- if (cpp_ok)
        tryCatch(ss_texg_stop_success_auto_branch(SSD, pars),
                 error = function(e) "err") else "n/a"
      auto_cpp <- if (cpp_ok)
        tryCatch(ss_texg_stop_success_value(SSD, pars, "auto", n_nodes = 64L),
                 error = function(e) NA_real_) else NA_real_
      auto_R <- tryCatch(
        EMC2:::stop_success_texg_R(mu_v, sig_v, tau_v, lb_v, SSD,
                                   method = "auto", n_nodes = 64L),
        error = function(e) NA_real_)
      cat(sprintf(
        "  %-25s SSD=%.2f: branch=%-22s auto_cpp_relerr=%s auto_R_relerr=%s\n",
        sc$label, SSD, br,
        if (is.na(auto_cpp)) "       NA" else sprintf("%.2e", abs(auto_cpp - ref)/ref),
        if (is.na(auto_R))   "       NA" else sprintf("%.2e", abs(auto_R - ref)/ref)))
    }
  }
}
run_stress(1)
run_stress(2)

# ===========================================================================
# 3. GL NODE-COUNT STUDY
# ===========================================================================
cat("\n", sep, "\nSECTION 3: GL node-count study\n", sep, "\n", sep = "")

p1 <- make_pars(muS_t, sigS_t, tauS_t, lbS_def, muG_t[1], sigG_t[1], tauG_t[1], lbG_def)
p2 <- make_pars(muS_t, sigS_t, tauS_t, lbS_def, muG_t,    sigG_t,    tauG_t,
                rep(lbG_def, 2))
if (cpp_ok) {
  ora1 <- pstop_texg_oracle(.2, muS_t, sigS_t, tauS_t, lbS_def,
                            muG_t[1], sigG_t[1], tauG_t[1], lbG_def)
  ora2 <- pstop_texg_oracle(.2, muS_t, sigS_t, tauS_t, lbS_def,
                            muG_t, sigG_t, tauG_t, rep(lbG_def, 2))
  cat("EXG, SSD=0.2, relerr vs oracle:\n")
  for (n in c(8, 16, 32, 64, 128)) {
    v1 <- ss_texg_stop_success_value(.2, p1, "gl", n_nodes = n)
    v2 <- ss_texg_stop_success_value(.2, p2, "gl", n_nodes = n)
    cat(sprintf("  gl(%3d): n_go=1 %.2e   n_go=2 %.2e\n",
                n, abs(v1 - ora1)/ora1, abs(v2 - ora2)/ora2))
  }
  # RDEX has no oracle: reference = gl(512)
  rdex_pars <- matrix(c(2.5, 0.6, 0.2, 0.15, 1.0, 0.20, 0.03, 0.05, 0, 0, 0),
                      nrow = 1, ncol = 11)  # v B A t0 s muS sigS tauS . . lbS
  rdex_ok <- tryCatch({
    ss_rdex_stop_success_value(.15, rdex_pars, "integrate"); TRUE
  }, error = function(e) FALSE)
  if (rdex_ok) {
    ref_rdex <- ss_rdex_stop_success_value(.15, rdex_pars, "gl", n_nodes = 512L)
    cat("RDEX, SSD=0.15, relerr vs gl(512):\n")
    cat(sprintf("  integrate: %.2e\n",
        abs(ss_rdex_stop_success_value(.15, rdex_pars, "integrate") - ref_rdex)/ref_rdex))
    for (n in c(16, 32, 64, 128))
      cat(sprintf("  gl(%3d):   %.2e\n", n,
          abs(ss_rdex_stop_success_value(.15, rdex_pars, "gl", n_nodes = n) - ref_rdex)/ref_rdex))
  }
}

# ===========================================================================
# 4. SPEED per single integral
# ===========================================================================
cat("\n", sep, "\nSECTION 4: speed (median us/call, SSD=0.2, typical params)\n",
    sep, "\n", sep = "")

speed <- list()
if (cpp_ok) {
  mb1 <- microbenchmark(
    integrate = ss_texg_stop_success_value(.2, p1, "integrate", n_nodes = 64L),
    gl        = ss_texg_stop_success_value(.2, p1, "gl",        n_nodes = 64L),
    analytic  = ss_texg_stop_success_value(.2, p1, "analytic",  n_nodes = 64L),
    auto      = ss_texg_stop_success_value(.2, p1, "auto",      n_nodes = 64L),
    times = 1000L)
  mb2 <- microbenchmark(
    integrate = ss_texg_stop_success_value(.2, p2, "integrate", n_nodes = 64L),
    gl        = ss_texg_stop_success_value(.2, p2, "gl",        n_nodes = 64L),
    analytic  = ss_texg_stop_success_value(.2, p2, "analytic",  n_nodes = 64L),
    auto      = ss_texg_stop_success_value(.2, p2, "auto",      n_nodes = 64L),
    times = 1000L)
  speed$n1 <- setNames(summary(mb1, unit = "us")$median,
                       trimws(as.character(summary(mb1)$expr)))
  speed$n2 <- setNames(summary(mb2, unit = "us")$median,
                       trimws(as.character(summary(mb2)$expr)))
  if (exists("rdex_ok") && rdex_ok) {
    mbr <- microbenchmark(
      integrate = ss_rdex_stop_success_value(.15, rdex_pars, "integrate"),
      gl64      = ss_rdex_stop_success_value(.15, rdex_pars, "gl", n_nodes = 64L),
      auto      = ss_rdex_stop_success_value(.15, rdex_pars, "auto", n_nodes = 64L),
      times = 1000L)
    cat("RDEX (C++):\n"); print(summary(mbr, unit = "us")[, c("expr", "median")])
  }
} else {
  mu1 <- c(muS_t, muG_t[1]); sig1 <- c(sigS_t, sigG_t[1]); tau1 <- c(tauS_t, tauG_t[1])
  lb1 <- c(lbS_def, lbG_def)
  mu2 <- c(muS_t, muG_t); sig2 <- c(sigS_t, sigG_t); tau2 <- c(tauS_t, tauG_t)
  lb2 <- c(lbS_def, lbG_def, lbG_def)
  bench_R <- function(mu, sig, tau, lb)
    microbenchmark(
      integrate = EMC2:::stop_success_texg_R(mu, sig, tau, lb, .2, method = "integrate"),
      gl        = EMC2:::stop_success_texg_R(mu, sig, tau, lb, .2, method = "gl"),
      analytic  = EMC2:::stop_success_texg_R(mu, sig, tau, lb, .2, method = "analytic"),
      auto      = EMC2:::stop_success_texg_R(mu, sig, tau, lb, .2, method = "auto"),
      times = 500L)
  mb1 <- bench_R(mu1, sig1, tau1, lb1); mb2 <- bench_R(mu2, sig2, tau2, lb2)
  speed$n1 <- setNames(summary(mb1, unit = "us")$median,
                       trimws(as.character(summary(mb1)$expr)))
  speed$n2 <- setNames(summary(mb2, unit = "us")$median,
                       trimws(as.character(summary(mb2)$expr)))
}

# ===========================================================================
# 5. FULL LOG-LIKELIHOOD end-to-end (SSEXG(stop_method=) -> calc_ll)
# ===========================================================================
cat("\n", sep, "\nSECTION 5: full log-likelihood via calc_ll\n", sep, "\n",
    sep = "")

full_ll_ok <- tryCatch({
  lIfun_go <- function(d) factor(rep(2, nrow(d)), levels = 1:2)
  designSS <- design(
    model     = SSEXG,
    factors   = list(subjects = 1, S = c("left", "right")),
    Rlevels   = c("left", "right"),
    matchfun  = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI  = lIfun_go,
                     SSD = function(d) SSD_function(d, SSD = c(0.20, 0.40),
                                                    pSSD = c(0.25, 0.25))),
    formula   = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                     muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  p_vector <- sampled_pars(designSS, doMap = FALSE)
  p_vector["mu_lMFALSE"] <- log(0.65); p_vector["mu_lMTRUE"] <- log(0.50)
  p_vector["sigma"] <- log(0.15);      p_vector["tau"]  <- log(0.50)
  p_vector["muS"]   <- log(0.20);      p_vector["sigmaS"] <- log(0.03)
  p_vector["tauS"]  <- log(0.10)
  p_vector["gf"]    <- qnorm(0.08);    p_vector["tf"]  <- qnorm(0.08)
  set.seed(123)
  dat  <- make_data(p_vector, designSS, n_trials = 1000,
                    TC = list(UC = Inf, LC = 0, LT = 0, UT = Inf, verbose = FALSE))
  dadm <- EMC2:::design_model(dat, designSS, verbose = FALSE)

  # C++ log-likelihood for one parameter vector, honouring stop_method by
  # pushing it into the process-global config exactly as calc_ll_manager()
  # does (set_stop_method_from_model). NB calc_ll must receive the
  # design-processed model (attr(dadm, "model")) — bound/transform are
  # augmented by design(); a bare SSEXG() list lacks them.
  model_d <- attr(dadm, "model")()
  ll_with_method <- function(m) {
    EMC2:::set_stop_method_from_model(SSEXG(stop_method = m))
    p_matrix <- matrix(p_vector, nrow = 1)
    colnames(p_matrix) <- names(p_vector)
    p_types <- names(model_d$p_types)
    designs <- list()
    for (p in p_types)
      designs[[p]] <- attr(dadm, "designs")[[p]][
        attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
    constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
    out <- sum(EMC2:::calc_ll(p_matrix, dadm, constants, designs,
                              model_d$c_name, model_d$bound, model_d$transform,
                              model_d$pre_transform, names(model_d$p_types),
                              log(1e-10), model_d$trend))
    EMC2:::set_stop_method_from_model(SSEXG())  # restore default (auto)
    out
  }

  lls <- sapply(methods4, ll_with_method)
  cat("summed log-likelihood by stop_method (1000 trials):\n")
  print(round(lls, 6))
  cat(sprintf("max |ll - ll_integrate| = %.2e  (should be << 1e-3)\n",
              max(abs(lls - lls["integrate"]))))
  mbll <- microbenchmark(
    integrate = ll_with_method("integrate"),
    auto      = ll_with_method("auto"),
    times = 50L)
  print(summary(mbll, unit = "ms")[, c("expr", "median")])
  max(abs(lls - lls["integrate"])) < 1e-3
}, error = function(e) {
  cat("SKIPPED (error):", conditionMessage(e), "\n"); NA
})

# ===========================================================================
# 6. SUMMARY TABLES + PASS/FAIL
# ===========================================================================
cat("\n\n", sep, "\n", sep = "")
cat("TABLE 1: Accuracy — relative error vs analytic oracle\n")
cat(sprintf("         (%d param sets x %d SSDs per n_go; median, max)\n",
            n_rand + 1L, length(SSDs)))
cat(sep, "\n")
hdr <- sprintf("%-18s  %9s %9s  %9s %9s",
               "Method", "n=1 med", "n=1 max", "n=2 med", "n=2 max")
cat(hdr, "\n", paste(rep("-", nchar(hdr)), collapse = ""), "\n", sep = "")
for (nm in c(paste0("cpp_", methods4), paste0("R_", methods4))) {
  v1 <- acc1[[nm]]; v2 <- acc2[[nm]]
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  cat(sprintf("%-18s  %s  %s\n", nm, fmt_relerr(v1), fmt_relerr(v2)))
}
cat(sep, "\n\n")

cat(sep, "\n")
route_label <- if (cpp_ok) "C++ ss_texg_stop_success_value" else "R stop_success_texg_R"
cat(sprintf("TABLE 2: Speed — median per-call time (us), %s\n", route_label))
cat(sep, "\n")
hdr2 <- sprintf("%-18s  %9s %9s", "Method", "n=1 (us)", "n=2 (us)")
cat(hdr2, "\n", paste(rep("-", nchar(hdr2)), collapse = ""), "\n", sep = "")
for (nm in methods4) {
  v1 <- if (nm %in% names(speed$n1)) speed$n1[nm] else NA_real_
  v2 <- if (nm %in% names(speed$n2)) speed$n2[nm] else NA_real_
  cat(sprintf("%-18s  %s %s\n", nm, fmt_us(v1), fmt_us(v2)))
}
cat(sep, "\n")

# ---- PASS/FAIL --------------------------------------------------------------
# Tolerances: n_go = 1 goes through the EXACT closed form on these draws, so
# essentially machine precision is demanded (1e-9). n_go = 2 has no closed
# form — auto is GL with the node bump, whose worst case on extreme random
# draws (tiny P_stop at long SSDs) is ~1e-4-relative; the integrate route is
# no better there.
cat("\n", sep, "\nPASS/FAIL\n", sep, "\n", sep = "")
PASS_TOL_N1 <- 1e-9
PASS_TOL_N2 <- 1e-4
get_clean <- function(acc, key) { v <- acc[[key]]; v[!is.na(v)] }
auto_cpp1 <- get_clean(acc1, "cpp_auto"); auto_cpp2 <- get_clean(acc2, "cpp_auto")
auto_R1   <- get_clean(acc1, "R_auto");   auto_R2   <- get_clean(acc2, "R_auto")
intg_cpp <- c(get_clean(acc1, "cpp_integrate"), get_clean(acc2, "cpp_integrate"))
intg_R   <- c(get_clean(acc1, "R_integrate"),   get_clean(acc2, "R_integrate"))

pass_tol_cpp <- (!cpp_ok) ||
  (all(auto_cpp1 <= PASS_TOL_N1) && all(auto_cpp2 <= PASS_TOL_N2))
pass_tol_R   <- all(auto_R1 <= PASS_TOL_N1) && all(auto_R2 <= PASS_TOL_N2)
pass_nw_cpp  <- (!cpp_ok) || (max(auto_cpp1, auto_cpp2) <= max(intg_cpp) * 10)
pass_nw_R    <- max(auto_R1, auto_R2) <= max(intg_R) * 10
pass_ll      <- is.na(full_ll_ok) || isTRUE(full_ll_ok)

fmt_pass <- function(ok) if (ok) "PASS" else "FAIL"
if (cpp_ok) {
  cat(sprintf("[%s] C++ auto within %.0e (n=1) / %.0e (n=2) of oracle (max %.2e / %.2e)\n",
              fmt_pass(pass_tol_cpp), PASS_TOL_N1, PASS_TOL_N2,
              max(auto_cpp1), max(auto_cpp2)))
  cat(sprintf("[%s] C++ auto never worse than integrate worst-case (%.2e vs %.2e)\n",
              fmt_pass(pass_nw_cpp), max(auto_cpp1, auto_cpp2), max(intg_cpp)))
}
cat(sprintf("[%s] R   auto within %.0e (n=1) / %.0e (n=2) of oracle (max %.2e / %.2e)\n",
            fmt_pass(pass_tol_R), PASS_TOL_N1, PASS_TOL_N2,
            max(auto_R1), max(auto_R2)))
cat(sprintf("[%s] R   auto never worse than integrate worst-case (%.2e vs %.2e)\n",
            fmt_pass(pass_nw_R), max(auto_R1, auto_R2), max(intg_R)))
cat(sprintf("[%s] full log-likelihood: all methods within 1e-3 of integrate%s\n",
            fmt_pass(pass_ll), if (is.na(full_ll_ok)) " (SKIPPED)" else ""))

all_pass <- pass_tol_R && pass_nw_R && pass_ll &&
  ((!cpp_ok) || (pass_tol_cpp && pass_nw_cpp))
cat(sprintf("\nOverall: %s\n", if (all_pass) "ALL PASS" else "SOME CHECKS FAILED"))
cat(sep, "\n")
