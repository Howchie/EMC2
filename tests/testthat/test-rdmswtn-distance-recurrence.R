# RDMSWTN no-clock CDF, survivor and density with drift and start-point
# variability, computed by the distance recurrence.  References are nested
# adaptive quadrature of the fixed-(v, r) Wald CDF / survivor / density over
# the drift law and the uniform distance (no bivariate-normal algebra),
# evaluated at rel.tol 1e-12 on a log-shifted integrand.

library(EMC2)

rdmswtn_ref_cases <- data.frame(
  t   = c(3, 3, 40, 30, 0.002, 2, 10, 1, 0.8),
  mu  = c(0.01, 0.01, 2, 1.5, 1, -1.5, 0.1, -2, 2.5),
  b   = c(4, 4, 1.5, 2, 2, 1.2, 2, 1, 1.1),
  A   = c(0, 1.5, 1, 0.8, 1, 1, 1.5, 0.5, 0),
  s   = c(1, 1, 1, 0.9, 1, 1, 1, 1, 1.2),
  sv  = c(2, 2, 0.3, 0.4, 1, 1, 8, 1, 0.7),
  pos = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE),
  logF = c(-0.5772787813908, -0.428988098838, -2.941781799623e-10,
           -0.000128597375148, -257.9049122239, -1.422355357266,
           -0.007617605370874, -0.5463569460736, -0.1483399297273),
  logS = c(-0.8242222719837, -1.053163927212, -21.94683532411,
           -8.958888455346, NA, -0.2759450403299, -4.88109959932,
           -0.8652542898563, -1.981502085386)
)

# Absolute agreement: on the log scale this is relative agreement of F or S.
expect_abs <- function(object, expected, tol, label = NULL) {
  expect_lte(max(abs(object - expected)), tol, label = label)
}

swtn_log_p <- function(x, lower_tail = TRUE, ...) {
  EMC2:::pSWTNspv(x$t, x$mu, x$b, x$A, x$s, 0, x$sv, log_out = TRUE,
                  posdrift = x$pos, lower_tail = lower_tail, ...)
}

test_that("log CDF and log survivor match nested quadrature", {
  for (i in seq_len(nrow(rdmswtn_ref_cases))) {
    x <- rdmswtn_ref_cases[i, ]
    expect_abs(swtn_log_p(x), x$logF, 1e-10, label = paste("log F, case", i))
    if (!is.na(x$logS))
      expect_abs(swtn_log_p(x, FALSE), x$logS, 1e-10, label = paste("log S, case", i))
  }
})

test_that("large tilt at a fixed threshold is not lost to the BVN tail", {
  # ell = 2 b mu / s^2 + 2 b^2 sv^2 / s^4 = 128: exp(ell) amplifies any
  # absolute error in the bivariate-normal seed.
  F <- EMC2:::prdmswtn(3, 0.01, 4, 0, s = 1, sv = 2, posdrift = TRUE)
  expect_abs(F, 0.5614240473, 1e-10)
  F_A <- EMC2:::prdmswtn(3, 0.01, 4, 1.5, s = 1, sv = 2, posdrift = TRUE)
  expect_abs(log(F_A), -0.428988098838, 1e-10)
})

test_that("CDF and survivor are complementary and continuous in A", {
  set.seed(11)
  n <- 60
  P <- data.frame(t = exp(runif(n, log(0.05), log(5))), mu = runif(n, -1, 3),
                  A = runif(n, 0, 1.5), s = runif(n, 0.6, 1.3),
                  sv = exp(runif(n, log(0.1), log(3))), pos = runif(n) < 0.5)
  P$b <- runif(n, 0.2, 2) + P$A
  P$mu[P$pos] <- abs(P$mu[P$pos]) + 0.01
  for (i in seq_len(n)) {
    x <- P[i, ]
    expect_equal(exp(swtn_log_p(x)) + exp(swtn_log_p(x, FALSE)), 1, tolerance = 1e-12)
  }
  x <- P[1, ]
  x$A <- 0
  lF0 <- swtn_log_p(x)
  x$A <- 1e-9
  expect_abs(swtn_log_p(x), lF0, 1e-8)
})

test_that("deep survivors keep relative accuracy where 1 - F underflows", {
  x <- data.frame(t = 100, mu = 4, b = 1.5, A = 1, s = 0.8, sv = 0.05, pos = TRUE)
  expect_abs(swtn_log_p(x, FALSE), -903.0439711922, 1e-9)
  x$pos <- FALSE
  expect_abs(swtn_log_p(x, FALSE), -903.0439711922, 1e-9)
  expect_equal(swtn_log_p(x), 0)
  x <- data.frame(t = 60, mu = 3, b = 1.5, A = 0, s = 1, sv = 0.05, pos = TRUE)
  expect_abs(swtn_log_p(x, FALSE), -238.844028687, 1e-9)
})

test_that("distance-averaged density matches nested quadrature in bulk and tails", {
  cases <- data.frame(
    t   = c(0.8, 0.0367, 6, 1.2, 0.05, 9),
    mu  = c(1.5, 0.476, 1.2, -0.5, 1, 2),
    b   = c(1.6, 2.4, 1.5, 1.4, 2.2, 1),
    A   = c(0.6, 0.638, 0.9, 0.7, 0.5, 0.8),
    s   = c(1, 1, 0.9, 1.1, 1, 1),
    sv  = c(0.8, 1.423, 0.5, 0.9, 1, 0.3),
    pos = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    logf = c(-0.5249231526785, -36.90857231669, -5.33456756881,
             -2.077934589548, -24.60278591293, -14.19695834606)
  )
  for (i in seq_len(nrow(cases))) {
    x <- cases[i, ]
    lf <- EMC2:::dSWTNspv(x$t, x$mu, x$b, x$A, x$s, 0, x$sv, log_out = TRUE,
                          posdrift = x$pos)
    expect_abs(lf, x$logf, 1e-10, label = paste("log f, case", i))
  }
})

test_that("full-Gaussian drift at t = Inf splits into hit and miss masses", {
  set.seed(12)
  for (k in 1:12) {
    mu <- runif(1, -2, 2); sv <- exp(runif(1, log(0.1), log(2)))
    s <- runif(1, 0.6, 1.3); A <- runif(1, 0, 1.5); b <- runif(1, 0.2, 2) + A
    lH <- EMC2:::pSWTNspv(Inf, mu, b, A, s, 0, sv, log_out = TRUE, posdrift = FALSE)
    lM <- EMC2:::pSWTNspv(Inf, mu, b, A, s, 0, sv, log_out = TRUE, posdrift = FALSE,
                          lower_tail = FALSE)
    expect_equal(exp(lH) + exp(lM), 1, tolerance = 1e-13)
    hit_r <- function(r) vapply(r, function(rr) pnorm(mu / sv) + integrate(
      function(v) exp(2 * rr * v / s^2) * dnorm(v, mu, sv), -Inf, 0,
      rel.tol = 1e-12)$value, 0)
    ref <- if (A > 0) integrate(hit_r, b - A, b, rel.tol = 1e-12)$value / A else hit_r(b)
    expect_abs(lH, log(ref), 1e-9)
  }
  expect_equal(EMC2:::pSWTNspv(Inf, 1, 1, 0.5, 1, 0, 1, posdrift = TRUE), 1)
})

test_that("clocks compose with the no-clock CDF, survivor and density", {
  set.seed(13)
  for (k in 1:24) {
    t <- runif(1, 0.3, 3); t0 <- runif(1, 0.05, 0.2); mu <- runif(1, 0.3, 3)
    A <- if (k %% 3 == 0) 0 else runif(1, 0, 1); b <- runif(1, 0.3, 2) + A
    s <- runif(1, 0.7, 1.2); sv <- runif(1, 0.1, 1.5); lam <- rexp(1, 1)
    shape <- 1 + (k %% 2); pos <- k %% 4 < 2
    if (!pos) mu <- mu - 1
    lSG <- if (shape == 1) -lam * t else -lam * t + log1p(lam * t)
    lfG <- if (shape == 1) log(lam) - lam * t else 2 * log(lam) + log(t) - lam * t
    lSD <- EMC2:::pSWTNspv(t, mu, b, A, s, t0, sv, log_out = TRUE, posdrift = pos,
                           lower_tail = FALSE)
    lfD <- EMC2:::dSWTNspv(t, mu, b, A, s, t0, sv, log_out = TRUE, posdrift = pos)
    # Guess clock: S_R = S_D S_G, f_R = f_D S_G + f_G S_D.
    lF_g <- EMC2:::pSWTNspv(t, mu, b, A, s, t0, sv, lambda_g = lam, log_out = TRUE,
                            kill_shape = shape, posdrift = pos)
    expect_abs(lF_g, log1p(-exp(lSD + lSG)), 1e-12)
    lS_g <- EMC2:::pSWTNspv(t, mu, b, A, s, t0, sv, lambda_g = lam, log_out = TRUE,
                            kill_shape = shape, posdrift = pos, lower_tail = FALSE)
    expect_abs(lS_g, lSD + lSG, 1e-12)
    lf_g <- EMC2:::dSWTNspv(t, mu, b, A, s, t0, sv, lambda_g = lam, log_out = TRUE,
                            kill_shape = shape, posdrift = pos)
    expect_equal(exp(lf_g), exp(lfD + lSG) + exp(lfG + lSD), tolerance = 1e-12)
    # Kill clock: f_R = f_D S_K.
    lf_k <- EMC2:::dSWTNspv(t, mu, b, A, s, t0, sv, lambda_k = lam, log_out = TRUE,
                            kill_shape = shape, posdrift = pos)
    expect_abs(lf_k, lfD + lSG, 1e-12)
  }
})
