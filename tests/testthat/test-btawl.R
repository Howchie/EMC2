skip_model_validation()

# BTAwL: ballistic Smith transient/sustained local race with state leak.
#
# dX/du = V (u/tau) exp(-u/tau) - k X, b = B + A, z ~ U(0, A).
#
# Checked against references written from the model definition (numerical
# minimisation of the required-launch curve and quadrature over start point),
# against Monte Carlo, against internal consistency (integrate(f) == F),
# and against the compiled race likelihood.

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

btawl_tmax <- function(k, tau) {
  EMC2:::btawl_tmax_vec(k, tau)
}

ref_H_btawl <- function(t, k, tau) {
  EMC2:::.btawl_h(t, k, tau)
}

ref_Hp_btawl <- function(t, k, tau) {
  g <- if (tau > 0 && t >= 0) (t / tau) * exp(-t / tau) else 0
  g - k * ref_H_btawl(t, k, tau)
}

ref_vstar_btawl <- function(u, z, b, k, tau) {
  h <- ref_H_btawl(u, k, tau)
  if (!isTRUE(h > 0)) return(Inf)
  E <- if (k <= 1e-12) 1 else exp(-k * u)
  (b - z * E) / h
}

ref_gbar_btawl <- function(w, p1, p2, launch, posdrift = TRUE) {
  if (launch == 1L) return(plnorm(w, p1, p2, lower.tail = FALSE))
  den <- if (posdrift) pnorm(0, p1, p2, lower.tail = FALSE) else 1
  pnorm(w, p1, p2, lower.tail = FALSE) / den
}

ref_g_btawl <- function(w, p1, p2, launch, posdrift = TRUE) {
  if (launch == 1L) return(dlnorm(w, p1, p2))
  den <- if (posdrift) pnorm(0, p1, p2, lower.tail = FALSE) else 1
  dnorm(w, p1, p2) / den
}

ref_F_btawl <- function(u, p1, p2, b, A, k, tau, launch = 1L, posdrift = TRUE) {
  if (!isTRUE(u > 0)) return(0)
  tm <- btawl_tmax(k, tau)
  vmin <- function(z) {
    hi <- min(u, if (is.finite(tm)) tm else u)
    o <- stats::optimize(function(s) ref_vstar_btawl(s, z, b, k, tau),
                         c(1e-12, hi), tol = 1e-12)$objective
    min(o, ref_vstar_btawl(hi, z, b, k, tau))
  }
  if (A <= 0) return(ref_gbar_btawl(vmin(0), p1, p2, launch, posdrift))
  stats::integrate(function(zz) {
    vapply(zz, function(z) ref_gbar_btawl(vmin(z), p1, p2, launch, posdrift),
           numeric(1))
  }, 0, A, rel.tol = 1e-10)$value / A
}

logdiff_btawl <- function(a, b) {
  if (!is.finite(a) || !is.finite(b) || a <= b) return(-Inf)
  a + log1p(-exp(b - a))
}

ref_log_launch_cdf_btawl <- function(w, p1, p2, launch, posdrift = TRUE) {
  if (!(w > 0)) return(-Inf)
  if (launch == 1L) return(pnorm((log(w) - p1) / p2, log.p = TRUE))
  upper <- pnorm((w - p1) / p2, log.p = TRUE)
  if (!posdrift) return(upper)
  lower <- pnorm(-p1 / p2, log.p = TRUE)
  logdiff_btawl(upper, lower) - pnorm(p1 / p2, log.p = TRUE)
}

ref_log_surv_btawl <- function(u, p1, p2, b, A, k, tau,
                               launch = 1L, posdrift = TRUE) {
  if (!isTRUE(u > 0)) return(0)
  one <- function(z) {
    objective <- function(s) {
      value <- ref_vstar_btawl(s, z, b, k, tau)
      if (is.finite(value)) value else .Machine$double.xmax
    }
    vmin <- optimize(objective, c(1e-8, min(u, btawl_tmax(k, tau))),
                     tol = 1e-12)$objective
    exp(ref_log_launch_cdf_btawl(vmin, p1, p2, launch, posdrift))
  }
  if (A <= 0) return(log(one(0)))
  integrand <- function(z) vapply(z, one, numeric(1))
  log(stats::integrate(integrand, 0, A, rel.tol = 1e-9,
                       subdivisions = 200L)$value / A)
}

ref_f_btawl <- function(u, p1, p2, b, A, k, tau, launch = 1L, posdrift = TRUE) {
  if (!isTRUE(u > 0)) return(0)
  tm <- btawl_tmax(k, tau)
  if (is.finite(tm) && u >= tm) return(0)
  h <- ref_H_btawl(u, k, tau)
  hp <- ref_Hp_btawl(u, k, tau)
  if (!isTRUE(h > 0) || !isTRUE(hp > 0)) return(0)
  E <- exp(-k * u)
  dens <- function(z) {
    w <- ref_vstar_btawl(u, z, b, k, tau)
    jac <- ((b - z * E) * hp - z * k * E * h) / (h * h)
    if (!isTRUE(jac > 0)) return(0)
    ref_g_btawl(w, p1, p2, launch, posdrift) * jac
  }
  if (A <= 0) return(dens(0))
  stats::integrate(function(zz) vapply(zz, dens, numeric(1)), 0, A,
                   rel.tol = 1e-10)$value / A
}

ref_fp_btawl <- function(V, z, b, k, tau) {
  if (!isTRUE(b > z)) return(0)
  if (!isTRUE(V > 0)) return(Inf)
  tm <- btawl_tmax(k, tau)
  H <- function(t) ref_H_btawl(t, k, tau)
  X <- function(t) z * exp(-k * t) + V * H(t)
  if (is.finite(tm)) {
    if (X(tm) < b) return(Inf)
    return(stats::uniroot(function(t) X(t) - b, c(0, tm), tol = 1e-12)$root)
  }
  hi <- 2 * max(tau, 1 / max(k, 1e-6))
  while (X(hi) < b && hi < 1e8) hi <- 2 * hi
  if (X(hi) < b) Inf else stats::uniroot(function(t) X(t) - b, c(0, hi), tol = 1e-12)$root
}

d_btawl <- function(t, A, b, p1, p2, k, tau, launch = 1L, posdrift = TRUE,
                    log_out = FALSE) {
  EMC2:::dbtawl_transient(t, A, b, p1, p2, k, tau, launch = as.integer(launch),
                posdrift = posdrift, log_out = log_out)
}
p_btawl <- function(t, A, b, p1, p2, k, tau, launch = 1L, posdrift = TRUE,
                    log_out = FALSE) {
  EMC2:::pbtawl_transient(t, A, b, p1, p2, k, tau, launch = as.integer(launch),
                posdrift = posdrift, log_out = log_out)
}

btawl_rows <- function(n, launch = 0L, A = 0.3, k = 0.5) {
  if (launch == 1L) {
    data.frame(mu = rep(0.2, n), sigma = rep(0.5, n),
               b = rep(1.3, n), A = rep(A, n), t0 = rep(0, n),
               k = rep(k, n), tau = rep(1, n))
  } else {
    data.frame(v = rep(1, n), sv = rep(1, n),
               b = rep(1.3, n), A = rep(A, n), t0 = rep(0, n),
               k = rep(k, n), tau = rep(1, n))
  }
}

btawl_local_race_rows <- function(n, launch = 0L, pi = 0.5, k = 0.5) {
  p <- btawl_rows(n, launch = launch, k = k)
  p$tau_s <- 2
  p$tau_t <- 1
  p$pi <- pi
  p
}

# ---------------------------------------------------------------------------
# 1. Geometry & Charts
# ---------------------------------------------------------------------------

test_that("T_max = 2*tau at k*tau = 1 identity and roundtrips with tau solver", {
  tau <- c(0.05, 0.2, 0.5, 1, 2)
  k <- 1 / tau
  tm <- EMC2:::btawl_tmax_vec(k, tau)
  expect_equal(tm, 2 * tau, tolerance = 2e-12)
  expect_equal(EMC2:::btawl_tau_vec(k, tm), tau, tolerance = 2e-10)
  expect_true(is.infinite(EMC2:::btawl_tmax_vec(0, 1)[1L]))
})

test_that("V_c(0) is b / H(T_max) and scales with b", {
  tau <- 1; k <- 1; b <- c(1.0, 2.5)
  vc <- EMC2:::btawl_vcrit_vec(rep(k, 2), rep(tau, 2), b)
  h_max <- EMC2:::.btawl_h(2 * tau, k, tau)
  expect_equal(vc, b / h_max, tolerance = 1e-12)
})

test_that("the stiff transient window remains finite and monotone", {
  k <- 20; tau <- 1; tm <- EMC2:::btawl_tmax_vec(k, tau)
  expect_equal(tm, 1.05263, tolerance = 2e-4)
  p <- btawl_rows(5, k = k, A = 0.3)
  tt <- tm * c(0.95, 0.99, 0.999, 1, 1.1)
  cdf <- EMC2:::pBTAwLTransient(tt, p)
  expect_true(all(is.finite(cdf)))
  expect_true(all(diff(cdf) >= -1e-10))
  expect_equal(cdf[4], cdf[5], tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 2. Both launches against independent start-point quadrature
# ---------------------------------------------------------------------------

btawl_cases <- list(
  list(A = 0.4, b = 1.3, p1 = 0.5, p2 = 0.5, k = 1.0, tau = 1.0, launch = 1L),
  list(A = 0.0, b = 1.0, p1 = 0.8, p2 = 0.4, k = 0.8, tau = 0.5, launch = 1L),
  list(A = 0.5, b = 1.5, p1 = 3.0, p2 = 1.0, k = 1.2, tau = 0.8, launch = 0L),
  list(A = 0.3, b = 1.2, p1 = 3.5, p2 = 1.2, k = 2.0, tau = 0.5, launch = 0L),
  list(A = 0.4, b = 1.4, p1 = 0.2, p2 = 0.6, k = 0.5, tau = 2.0, launch = 1L)
)

test_that("the CDF matches independent start-point quadrature", {
  skip_on_cran()
  for (cs in btawl_cases) {
    Tm <- btawl_tmax(cs$k, cs$tau)
    us <- Tm * c(0.05, 0.2, 0.5, 0.75, 0.9, 0.98)
    got <- p_btawl(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
    ref <- vapply(us, ref_F_btawl, numeric(1), cs$p1, cs$p2, cs$b, cs$A, cs$k,
                  cs$tau, cs$launch)
    expect_equal(got, ref, tolerance = 1e-6)
  }
})

test_that("the density matches independent start-point quadrature", {
  skip_on_cran()
  for (cs in btawl_cases) {
    Tm <- btawl_tmax(cs$k, cs$tau)
    us <- Tm * c(0.05, 0.2, 0.5, 0.75, 0.9, 0.98)
    got <- d_btawl(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
    ref <- vapply(us, ref_f_btawl, numeric(1), cs$p1, cs$p2, cs$b, cs$A, cs$k,
                  cs$tau, cs$launch)
    expect_equal(got, ref, tolerance = 1e-6)
  }
})

test_that("integrate(f) reproduces F through the live/frozen seam", {
  skip_on_cran()
  for (cs in btawl_cases) {
    Tm <- btawl_tmax(cs$k, cs$tau)
    for (frac in c(0.5, 0.95)) {
      u1 <- Tm * frac
      I <- stats::integrate(function(x) {
        d_btawl(x, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
      }, 1e-10, u1, subdivisions = 2000, rel.tol = 1e-10)$value
      expect_lt(abs(I - p_btawl(u1, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau,
                                launch = cs$launch)), 1e-7)
    }
  }
})

# ---------------------------------------------------------------------------
# 3. The endpoint & shutdown
# ---------------------------------------------------------------------------

test_that("the CDF is flat past T_max and F(Inf) is the finishing mass", {
  for (cs in btawl_cases) {
    Tm <- btawl_tmax(cs$k, cs$tau)
    at <- p_btawl(Tm, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
    beyond <- p_btawl(c(Tm * 1.2, Tm * 10, Inf), cs$A, cs$b, cs$p1, cs$p2,
                      cs$k, cs$tau, launch = cs$launch)
    expect_equal(beyond, rep(at, 3), tolerance = 1e-12)
    expect_lt(at, 1)  # defective tail
    expect_equal(d_btawl(c(Tm, Tm * 1.5), cs$A, cs$b, cs$p1, cs$p2, cs$k,
                         cs$tau, launch = cs$launch), c(0, 0))
  }
})

test_that("density shuts down quadratically for A > 0, linearly for A = 0", {
  k <- 1.0; tau <- 1.0; b <- 1.3; Tm <- 2.0
  eps <- c(1e-2, 5e-3, 2.5e-3, 1.25e-3)
  ord <- function(A) {
    fv <- d_btawl(Tm - eps, A, b, 0.5, 0.5, k, tau, launch = 1L)
    log(fv[-1] / fv[-length(fv)]) / log(eps[-1] / eps[-length(eps)])
  }
  expect_equal(ord(0), rep(1, 3), tolerance = 0.05)
  expect_equal(ord(0.4), rep(2, 3), tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# 4. Monte Carlo
# ---------------------------------------------------------------------------

test_that("Monte Carlo agrees for both launch distributions", {
  skip_on_cran()
  set.seed(42)
  N <- 1000
  for (cs in btawl_cases[c(1, 3)]) {
    V <- if (cs$launch == 1L) rlnorm(N, cs$p1, cs$p2)
         else msm::rtnorm(N, cs$p1, cs$p2, lower = 0)
    z <- runif(N, 0, cs$A)
    hit <- vapply(seq_len(N), function(i)
      ref_fp_btawl(V[i], z[i], cs$b, cs$k, cs$tau), numeric(1))
    Tm <- btawl_tmax(cs$k, cs$tau)
    us <- Tm * c(0.25, 0.5, 0.75, 0.95)
    emp <- vapply(us, function(u) mean(hit <= u), numeric(1))
    got <- p_btawl(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
    expect_lt(max(abs(got - emp)), 4 * sqrt(0.25 / N) + 1e-3)
    expect_lt(abs(mean(is.finite(hit)) -
                    p_btawl(Inf, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau,
                            launch = cs$launch)),
              4 * sqrt(0.25 / N) + 1e-3)
  }
})

# ---------------------------------------------------------------------------
# 5. Numerical hygiene
# ---------------------------------------------------------------------------

test_that("log output equals log of natural output", {
  for (cs in btawl_cases) {
    Tm <- btawl_tmax(cs$k, cs$tau)
    us <- Tm * c(0.01, 0.3, 0.6, 0.95)
    for (fn in list(d_btawl, p_btawl)) {
      nat <- fn(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
      lg <- fn(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch,
               log_out = TRUE)
      keep <- nat > 0
      expect_equal(lg[keep], log(nat[keep]), tolerance = 1e-8)
    }
  }
})

test_that("the CDF is monotone and density non-negative across seams", {
  for (cs in btawl_cases) {
    Tm <- btawl_tmax(cs$k, cs$tau)
    # Representative support/interior probes are sufficient for this package
    # contract; dense CDF sweeps are development-time validation.
    us <- Tm * c(1e-4, 0.01, 0.1, 0.5, 0.9, 0.99, 0.9999)
    Fv <- p_btawl(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
    fv <- d_btawl(us, cs$A, cs$b, cs$p1, cs$p2, cs$k, cs$tau, launch = cs$launch)
    expect_false(anyNA(c(Fv, fv)))
    expect_true(all(diff(Fv) >= -1e-12))
    expect_true(all(fv >= 0))
  }
})

test_that("log survivors stay finite, monotone, and below the raw floor", {
  tt <- c(0.25, 0.5, 1, 1.5, 1.9, 1.99, 2)
  ls_n <- EMC2:::btawl_transient_log_surv_vec(
    tt, A = 0.4, b = 1.3, p1 = 10, p2 = 0.5, k = 1, tau = 1,
    launch = 0L
  )
  ls_l <- EMC2:::btawl_transient_log_surv_vec(
    tt, A = 0.4, b = 1.3, p1 = 4, p2 = 0.2, k = 1, tau = 1,
    launch = 1L
  )
  for (ls in list(ls_n, ls_l)) {
    expect_true(all(is.finite(ls)))
    expect_true(all(diff(ls) <= 1e-8))
    expect_lt(min(ls), log(1e-10))
  }
  expect_equal(
    ls_n[1:3],
    log1p(-p_btawl(tt[1:3], 0.4, 1.3, 10, 0.5, 1, 1, launch = 0L)),
    tolerance = 2e-6
  )
})

test_that("log survivors match the direct start-point launch-CDF integral", {
  cases <- list(
    list(p1 = 10, p2 = 0.5, launch = 0L),
    list(p1 = 20, p2 = 0.5, launch = 0L),
    list(p1 = 4, p2 = 0.2, launch = 1L),
    list(p1 = 7, p2 = 0.2, launch = 1L)
  )
  tt <- c(1.5, 1.9, 1.99, 2)
  for (cs in cases) {
    got <- EMC2:::btawl_transient_log_surv_vec(
      tt, A = 0.4, b = 1.3, p1 = cs$p1, p2 = cs$p2, k = 1, tau = 1,
      launch = cs$launch
    )
    ref <- vapply(tt, ref_log_surv_btawl, numeric(1), cs$p1, cs$p2,
                  1.3, 0.4, 1, 1, cs$launch)
    expect_equal(got, ref, tolerance = 2e-5)
    expect_true(all(got < log(1e-10)))
  }
})

test_that("local-race all-live survivors survive interval collapse", {
  for (launch in 0:1) {
    p1 <- if (launch == 0L) 20 else 7
    p2 <- 0.2
    got <- EMC2:::btawl_local_race_log_surv_vec(
      c(50, 100), A = 0.4, b = 1.3, p1 = p1, p2 = p2, k = 1,
      tau_s = 1, tau_t = 1, pi = 1, launch = launch
    )
    ref <- if (launch == 0L)
      pnorm((1.3 - p1) / p2, log.p = TRUE) else
      pnorm((log(1.3) - p1) / p2, log.p = TRUE)
    expect_true(all(is.finite(got)))
    expect_equal(got, rep(ref, 2), tolerance = 1e-6)
  }
})

test_that("boundary parameters produce no NaN", {
  grid <- expand.grid(A = c(0, 1e-10, 0.5), k = c(1e-4, 1, 50),
                      tau = c(0.05, 1, 10), u = c(1e-6, 0.3, 10, Inf))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    for (launch in c(0L, 1L)) {
      expect_false(is.na(d_btawl(g$u, g$A, 1.2, 1, 0.5, g$k, g$tau, launch = launch)))
      expect_false(is.na(p_btawl(g$u, g$A, 1.2, 1, 0.5, g$k, g$tau, launch = launch)))
    }
  }
})

# ---------------------------------------------------------------------------
# 6. Constructors & Charts
# ---------------------------------------------------------------------------

test_that("constructors expose full local-race and pure-process charts", {
  endpoint <- BTAwL(drift_distribution = "normal", chart = "endpoint")
  rate <- BTAwL(drift_distribution = "normal", chart = "rate")
  expect_equal(endpoint$c_name, "BTAwL")
  expect_equal(rate$c_name, "BTAwL_RATE")
  expect_true("Ttrans" %in% names(endpoint$p_types))
  expect_true(all(c("tau_s", "Ttrans", "pi") %in% endpoint$p_types_canonical))
  expect_true(all(c("tau_s", "tau_t", "pi") %in% rate$p_types_canonical))
  expect_equal(BTAwLTransient()$c_name, "BTAwL_TRANSIENT")
  expect_equal(BTAwLSustained()$c_name, "BTAwL_SUSTAINED")
})

test_that("Ttransform reports b, tau, Tmax, rt_max and Vcrit", {
  m <- BTAwLTransient(drift_distribution = "normal", chart = "endpoint")
  pars <- cbind(v = c(2, 2), sv = c(0.5, 0.5), B = c(0.8, 1.6),
                A = c(0.3, 0.3), t0 = c(0.15, 0.15), k = c(1, 1),
                Ttrans = c(2, 2))
  out <- m$Ttransform(pars, NULL)
  expect_equal(out[, "b"], pars[, "B"] + pars[, "A"])
  expect_equal(out[, "tau"], EMC2:::btawl_tau_vec(pars[, "k"], pars[, "Ttrans"]), tolerance = 1e-10)
  expect_equal(out[, "Tmax"], pars[, "Ttrans"])
  expect_equal(out[, "rt_max"], pars[, "t0"] + out[, "Tmax"])
  expect_equal(out[, "Vcrit"], EMC2:::btawl_vcrit_vec(pars[, "k"], out[, "tau"], out[, "b"]))

  rate <- BTAwL(drift_distribution = "normal", chart = "rate")
  rate_pars <- cbind(v = 2, sv = 0.5, B = 0.8, A = 0.3, t0 = 0.15,
                     k = 1, tau_s = 1.2, tau_t = 0.7, pi = 0.4)
  rate_out <- rate$Ttransform(rate_pars, NULL)
  expect_equal(sum(colnames(rate_out) == "tau_t"), 1L)
  expect_equal(unname(rate_out[, "tau_t"]), unname(rate_pars[, "tau_t"]))
})

test_that("dfun/pfun use the same launch distribution as c_name", {
  pars <- cbind(mu = 1, sigma = 0.5, v = 1, sv = 0.5, b = 1.3, A = 0.3,
                t0 = 0.15, k = 1.0, tau = 1.0)
  rt <- 0.5
  expect_equal(BTAwLTransient(drift_distribution = "lognormal", chart = "rate")$dfun(rt, pars),
               d_btawl(rt - 0.15, 0.3, 1.3, 1, 0.5, 1.0, 1.0, launch = 1L))
  expect_equal(BTAwLTransient(drift_distribution = "normal", chart = "rate")$dfun(rt, pars),
               d_btawl(rt - 0.15, 0.3, 1.3, 1, 0.5, 1.0, 1.0, launch = 0L))
})

# ---------------------------------------------------------------------------
# 7. Compiled Race Likelihood
# ---------------------------------------------------------------------------

ref_race_ll_btawl <- function(dadm, pars, launch, min_ll = log(1e-10)) {
  nm <- if (launch == 1L) c("mu", "sigma") else c("v", "sv")
  clear_col <- if ("Ttrans" %in% colnames(pars)) "Ttrans" else "tau"
  tau_vals <- if ("Ttrans" %in% colnames(pars))
    EMC2:::btawl_tau_vec(pars[, "k"], pars[, "Ttrans"]) else pars[, "tau"]
  Fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    ref_F_btawl(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                pars[i, "k"], tau_vals[i], launch)
  }
  fu <- function(i, u) {
    if (!isTRUE(u > 0)) return(0)
    ref_f_btawl(u, pars[i, nm[1]], pars[i, nm[2]], pars[i, "b"], pars[i, "A"],
                pars[i, "k"], tau_vals[i], launch)
  }
  n_lR <- length(levels(dadm$lR))
  n_tr <- nrow(dadm) / n_lR
  lt <- numeric(n_tr)
  for (j in seq_len(n_tr)) {
    idx <- ((j - 1) * n_lR + 1):(j * n_lR)
    rt <- dadm$rt[idx[1]]
    if (is.infinite(rt) && rt > 0) {
      term <- prod(vapply(idx, function(i) 1 - Fu(i, Inf), numeric(1)))
    } else {
      w <- idx[which(dadm$winner[idx])]
      term <- fu(w, rt - pars[w, "t0"])
      for (i in setdiff(idx, w))
        term <- term * (1 - Fu(i, rt - pars[i, "t0"]))
    }
    lt[j] <- max(log(max(term, 0)), min_ll)
  }
  sum(lt[attr(dadm, "expand")])
}

btawl_mk <- function(model, form, consts, dat) {
  des <- design(data = dat, model = model, matchfun = function(d) d$S == d$lR,
                formula = form, constants = consts)
  list(des = des, emc = suppressMessages(
    make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
             rt_resolution = NULL)))
}

btawl_ll <- function(fx, p, p_types_override = NULL) {
  model <- fx$emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- fx$emc[[1]]$data[[1]]
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  if (!is.null(p_types_override)) p_types <- p_types_override
  EMC2:::calc_ll_oo(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
                    constants = attr(dadm, "constants"), designs = designs,
                    type = model$c_name, bounds = model$bound,
                    transforms = model$transform,
                    pretransforms = model$pre_transform, p_types = p_types,
                    min_ll = log(1e-10), trend = model$trend)
}

btawl_dat <- function(n = 60) {
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
  dat$subjects <- droplevels(dat$subjects)
  dat[seq_len(n), ]
}

test_that("the compiled race likelihood matches the R reference", {
  skip_on_cran()
  dat <- btawl_dat()
  p <- c(mu = 0.5, sigma = log(0.5), v = 2.5, sv = log(1), B = log(0.8),
         A = log(0.3), t0 = log(0.15), k = log(1), Ttrans = log(2))
  fx <- list(
    ln = suppressMessages(btawl_mk(function() BTAwLTransient(drift_distribution = "lognormal"),
      list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, Ttrans ~ 1), NULL, dat)),
    no = suppressMessages(btawl_mk(function() BTAwLTransient(drift_distribution = "normal"),
      list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, Ttrans ~ 1), c(sv = log(1)), dat)))
  for (which in c("ln", "no")) {
    launch <- if (which == "ln") 1L else 0L
    pp <- p[names(sampled_pars(fx[[which]]$des))]
    got <- btawl_ll(fx[[which]], pp)
    expect_true(is.finite(got))
    dadm <- fx[[which]]$emc[[1]]$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(pp, dadm, fx[[which]]$emc[[1]]$model())
    expect_equal(got, ref_race_ll_btawl(dadm, pars, launch), tolerance = 1e-5)
  }
})

test_that("a p_types reordering is caught by the column contract", {
  skip_on_cran()
  dat <- btawl_dat()
  fx <- suppressMessages(btawl_mk(function() BTAwLTransient(drift_distribution = "lognormal"),
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, Ttrans ~ 1), NULL, dat))
  p <- c(mu = 0.5, sigma = log(0.5), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(1), Ttrans = log(2))
  swapped <- names(BTAwLTransient(drift_distribution = "lognormal")$p_types)
  swapped[1:2] <- swapped[2:1]
  expect_error(btawl_ll(fx, p[names(sampled_pars(fx$des))],
                        p_types_override = swapped),
               "BTAwL_LOGN kernels expect parameter column")
})

test_that("omissions past t0 + T_max stay well posed", {
  skip_on_cran()
  dat <- btawl_dat(80)
  dat$rt[1:8] <- Inf
  dat$R[1:8] <- NA
  fx <- suppressMessages(btawl_mk(function() BTAwLTransient(drift_distribution = "lognormal"),
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, Ttrans ~ 1), NULL, dat))
  p <- c(mu = 0.5, sigma = log(0.5), B = log(0.8), A = log(0.3),
         t0 = log(0.15), k = log(1), Ttrans = log(2))
  got <- btawl_ll(fx, p[names(sampled_pars(fx$des))])
  expect_true(is.finite(got))
  dadm <- fx$emc[[1]]$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(p[names(sampled_pars(fx$des))], dadm,
                                    fx$emc[[1]]$model())
  expect_equal(got, ref_race_ll_btawl(dadm, pars, 1L), tolerance = 1e-5)
  # Beyond endpoint floored
  dat2 <- dat
  dat2$rt[1:8] <- 0.15 + 2 + 5
  fx2 <- suppressMessages(btawl_mk(function() BTAwLTransient(drift_distribution = "lognormal"),
    list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, Ttrans ~ 1), NULL, dat2))
  expect_true(is.finite(btawl_ll(fx2, p[names(sampled_pars(fx2$des))])))
})

# ---------------------------------------------------------------------------
# 8. Local-race invariants
# ---------------------------------------------------------------------------

test_that("the local race is exactly transient BTAwL at pi = 0", {
  tt <- c(0.5, 1, 1.5, 2, 2.5, 3, Inf)
  for (launch in 0:1) {
    p0 <- btawl_rows(length(tt), launch = launch)
    pm <- btawl_local_race_rows(length(tt), launch = launch, pi = 0)
    expect_equal(
      EMC2:::pBTAwL(tt, pm, launch = launch),
      EMC2:::pBTAwLTransient(tt, p0, launch = launch), tolerance = 1e-12
    )
    expect_equal(
      EMC2:::dBTAwL(tt, pm, launch = launch),
      EMC2:::dBTAwLTransient(tt, p0, launch = launch), tolerance = 1e-12
    )
  }
})

test_that("pi = 1 handles the sustained-only member", {
  tt <- c(1, 2, 3, 5, 10, 20, Inf)
  p <- btawl_local_race_rows(length(tt), pi = 1)
  cdf <- EMC2:::pBTAwL(tt, p)
  expect_true(all(diff(cdf) >= -1e-10))
  expect_gt(cdf[6], cdf[5])
  expect_gt(cdf[7], cdf[6])
  expect_gt(EMC2:::dBTAwL(10, p[5, ]), 0)
})

test_that("the no-leak sustained limit is retained", {
  tt <- c(1, 2, 5, Inf)
  p <- btawl_local_race_rows(length(tt), pi = 1, k = 0)
  cdf <- EMC2:::pBTAwL(tt, p)
  expect_true(all(diff(cdf) >= -1e-12))
  expect_equal(cdf[length(cdf)], 1, tolerance = 1e-12)
})

test_that("within-accumulator race CDF and PDF match pure member composition", {
  tt <- c(0.2, 0.5, 1.0, 1.8, 2.5)
  for (launch in 0:1) {
    for (pi_val in c(0.2, 0.5, 0.8)) {
      p_local <- btawl_local_race_rows(length(tt), launch = launch, pi = pi_val, k = 0.8)
      p_local$tau_s <- 1.5; p_local$tau_t <- 0.6
      
      # Transient member parameters
      p_T <- p_local
      p_T$tau <- p_local$tau_t
      if (launch == 1L) {
        p_T$mu <- p_local$mu + log(1 - pi_val)
      } else {
        p_T$v <- p_local$v * (1 - pi_val)
        p_T$sv <- p_local$sv * (1 - pi_val)
      }
      
      # Sustained member parameters
      p_S <- p_local
      p_S$pi <- 1
      if (launch == 1L) {
        p_S$mu <- p_local$mu + log(pi_val)
      } else {
        p_S$v <- p_local$v * pi_val
        p_S$sv <- p_local$sv * pi_val
      }
      
      F_T <- EMC2:::pBTAwLTransient(tt, p_T, launch = launch)
      f_T <- EMC2:::dBTAwLTransient(tt, p_T, launch = launch)
      F_S <- EMC2:::pBTAwL(tt, p_S, launch = launch)
      f_S <- EMC2:::dBTAwL(tt, p_S, launch = launch)
      
      F_expected <- 1 - (1 - F_T) * (1 - F_S)
      f_expected <- f_T * (1 - F_S) + f_S * (1 - F_T)
      
      F_got <- EMC2:::pBTAwL(tt, p_local, launch = launch)
      f_got <- EMC2:::dBTAwL(tt, p_local, launch = launch)
      
      expect_equal(F_got, F_expected, tolerance = 1e-12)
      expect_equal(f_got, f_expected, tolerance = 1e-12)
    }
  }
})

test_that("integrate(dBTAwL) matches pBTAwL for within-accumulator race", {
  for (launch in 0:1) {
    p_local <- btawl_local_race_rows(1, launch = launch, pi = 0.4, k = 0.5)
    p_local$tau_s <- 1.2; p_local$tau_t <- 0.7
    for (t_eval in c(0.5, 1.2, 2.0)) {
      I <- stats::integrate(function(x) {
        p_eval <- p_local[rep(1, length(x)), ]
        EMC2:::dBTAwL(x, p_eval, launch = launch)
      }, 1e-8, t_eval, rel.tol = 1e-8)$value
      expect_equal(I, EMC2:::pBTAwL(t_eval, p_local, launch = launch), tolerance = 1e-6)
    }
  }
})

# ---------------------------------------------------------------------------
# 9. Simulation
# ---------------------------------------------------------------------------

test_that("BTAwL simulators select a winner per trial", {
  set.seed(42)
  lR <- factor(rep(c("left", "right"), 4))
  pars <- data.frame(v = rep(3, 8), sv = rep(0.05, 8),
                     B = rep(1, 8), A = rep(0.2, 8), t0 = rep(0.1, 8),
                     k = rep(0.5, 8), tau = rep(1, 8), b = rep(1.2, 8))
  out <- EMC2:::rBTAwLTransient(lR, pars)
  expect_equal(nrow(out), 4)
  expect_true(all(is.finite(out$rt)))

  local_pars <- data.frame(
    v = rep(3, 8), sv = rep(0.2, 8), B = rep(1, 8), A = rep(0.2, 8),
    t0 = rep(0.1, 8), k = rep(0.5, 8), tau_s = rep(1.2, 8),
    tau_t = rep(0.7, 8), pi = rep(0.4, 8), b = rep(1.2, 8)
  )
  local_out <- EMC2:::rBTAwL(lR, local_pars)
  expect_equal(nrow(local_out), 4)
  expect_true(all(is.finite(local_out$rt)))
  expect_true(all(!is.na(local_out$R)))
})

test_that("the transient design produces omissions that can be fit back through", {
  skip_on_cran()
  dat <- btawl_dat(40)
  des <- suppressMessages(design(
    data = dat, model = BTAwLTransient, matchfun = function(d) d$S == d$lR,
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, Ttrans ~ 1),
    constants = c(sv = log(1))))
  p <- c(v = 2.5, B = log(0.8), A = log(0.3), t0 = log(0.15),
         k = log(1), Ttrans = log(2))
  set.seed(3)
  sim <- make_data(p, design = des, n_trials = 60)
  expect_true(any(is.infinite(sim$rt)))  # intrinsic omissions
  fin <- is.finite(sim$rt)
  expect_true(all(sim$rt[fin] <= 0.15 + 2 + 1e-8))
  emc <- suppressMessages(make_emc(sim, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  expect_true(is.finite(btawl_ll(list(emc = emc, des = des), p)))
})

# ---------------------------------------------------------------------------
# Nested submodels reachable from the full kernel
# ---------------------------------------------------------------------------

test_that("pi is permitted on both endpoints", {
  for (chart in c("rate", "endpoint")) {
    exc <- BTAwL(chart = chart)$bound$exception
    expect_equal(sort(unname(exc[names(exc) == "pi"])), c(0, 1))
    # The bound machinery must accept both endpoints, and nothing between them
    # that would otherwise be excluded by the strict [0,1] test.
    mm <- BTAwL(chart = chart)$bound$minmax
    expect_equal(unname(mm[, "pi"]), c(0, 1))
  }
})

test_that("bound exceptions admit several permitted values per parameter", {
  bound <- list(minmax = cbind(x = c(0, 1), y = c(0, 1)),
                exception = c(x = 0, x = 1, y = 0))
  pars <- cbind(x = c(0, 1, 0.5, 2), y = c(0, 0, 0, 0))
  expect_equal(unname(EMC2:::do_bound(pars, bound)),
               c(TRUE, TRUE, TRUE, FALSE))
  # y permits only 0, so 1 must still be rejected on the strict upper test.
  pars2 <- cbind(x = c(0.5, 0.5), y = c(0, 1))
  expect_equal(unname(EMC2:::do_bound(pars2, bound)), c(TRUE, FALSE))
})

test_that("pi = 0 and pi = 1 reproduce the dedicated wrappers exactly", {
  skip_on_cran()
  dat <- btawl_dat()
  vals <- c(mu = 0.5, sigma = log(0.5), B = log(0.8), A = log(0.3),
            t0 = log(0.15), k = log(1), tau = log(0.6), tau_t = log(0.6),
            tau_s = log(0.6), pi = qnorm(0.5))
  ll_of <- function(model, form, consts) {
    fx <- btawl_mk(model, form, consts, dat)
    p <- vals[names(sampled_pars(fx$des))]
    btawl_ll(fx, p)
  }
  full <- list(mu ~ 1, B ~ 1, A ~ 1, sigma ~ 1, t0 ~ 1, k ~ 1, tau_t ~ 1,
               tau_s ~ 1, pi ~ 1)

  # Pure transient: pi = 0, tau_s inert.
  expect_equal(
    ll_of(BTAwL(chart = "rate"), full, c(pi = qnorm(0), tau_s = log(1))),
    ll_of(BTAwLTransient(chart = "rate"),
          list(mu ~ 1, B ~ 1, A ~ 1, sigma ~ 1, t0 ~ 1, k ~ 1, tau ~ 1),
          NULL))

  # Pure sustained: pi = 1, tau_t inert.
  expect_equal(
    ll_of(BTAwL(chart = "rate"), full, c(pi = qnorm(1), tau_t = log(1))),
    ll_of(BTAwLSustained(chart = "rate"),
          list(mu ~ 1, B ~ 1, A ~ 1, sigma ~ 1, t0 ~ 1, k ~ 1, tau_s ~ 1),
          NULL))
})

test_that("an inert time constant really is inert at a pinned pi", {
  # This is what makes the nesting safe: at pi = 0 the sustained clearance
  # cannot affect the likelihood, so pinning it at any legal value is a free
  # choice rather than an assumption.
  for (ts in c(log(.2), log(1), log(5))) {
    m <- BTAwL(chart = "rate")
    pars <- cbind(mu = log(6), sigma = .5, b = 1.4, A = .4, t0 = .2, k = 2,
                  tau_t = .6, tau_s = exp(ts), pi = 0)
    expect_equal(m$dfun(0.5, pars),
                 dBTAwLTransient(0.5, cbind(mu = log(6), sigma = .5, b = 1.4,
                                            A = .4, t0 = .2, k = 2, tau = .6),
                                 launch = 1L))
  }
})
