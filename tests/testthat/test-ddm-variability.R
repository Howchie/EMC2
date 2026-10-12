# Wiener DDM with start-point (SZ) and onset (st0) variability. The density and
# CDF are fixed-node integrals (src/ddm_variability.h), so they are smooth in
# the parameters; with st0 > 0 a binned response integrates its bin.

ddm_var_frame <- function(rt, R) {
  data.frame(subjects = factor(1), S = factor("s"),
             R = factor(R, levels = c("lower", "upper")), rt = rt)
}

ddm_var_data <- function(rt, R, rt_resolution = NULL, dat = ddm_var_frame(rt, R)) {
  des <- design(data = dat, model = DDM,
                formula = list(v ~ 1, a ~ 1, t0 ~ 1, s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1),
                constants = c(s = 0), verbose = FALSE)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   rt_resolution = rt_resolution, compress = FALSE))
  list(dadm = emc[[1]]$data[[1]], model = emc[[1]]$model())
}

# Trial-wise log likelihoods, one row per parameter vector (natural scale).
ddm_var_pointwise <- function(x, pars) {
  dadm <- x$dadm
  model <- x$model
  p_types <- names(model$p_types)
  designs <- lapply(setNames(p_types, p_types), function(p)
    attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE])
  p_mat <- with(pars, cbind(v = v, a = log(a), t0 = log(t0), Z = qnorm(Z),
                            sv = log(sv), SZ = qnorm(SZ), st0 = log(st0),
                            s = 0, pContaminant = -Inf, pGuess = -Inf))[, p_types, drop = FALSE]
  out <- EMC2:::calc_ll_oo_pw(p_mat, dadm, constants = attr(dadm, "constants"),
                              designs = designs, type = model$c_name, model$bound,
                              model$transform, model$pre_transform, p_types = p_types,
                              min_ll = log(1e-10), model$trend)
  matrix(out, nrow = nrow(p_mat))
}

test_that("DDM density with SZ and st0 is smooth in the parameters", {
  # One bin above onset, where the adaptive integrator the density used before
  # accepted a first-pass estimate over a band of `a` and jumped by ~5e-3.
  t0 <- 0.3
  x <- ddm_var_data(t0 + 1 / 60, "lower")
  la <- seq(-0.99, -0.94, length.out = 201)
  pars <- data.frame(v = 1, a = exp(la), t0 = t0, Z = 0.5, sv = 0.5, SZ = 0.1, st0 = 0.05)
  ll <- ddm_var_pointwise(x, pars)[, 1]
  expect_true(all(is.finite(ll)))
  expect_lt(max(abs(diff(ll, differences = 2))), 1e-6)
})

# Integral of the unbinned density of response R over [lo, hi], on Gauss-Legendre
# panels in sqrt(t - t0) split at the end of the onset window.
ddm_var_integral <- function(lo, hi, R, pars, n = 32) {
  g <- statmod::gauss.quad(n)
  cuts <- sort(unique(pmin(pmax(c(lo, hi, pars$t0 + pars$st0), max(lo, pars$t0)), hi)))
  total <- 0
  for (i in seq_len(length(cuts) - 1)) {
    zl <- sqrt(cuts[i] - pars$t0); zh <- sqrt(cuts[i + 1] - pars$t0)
    if (!(zh > zl)) next
    z <- (zl + zh) / 2 + (zh - zl) / 2 * g$nodes
    t <- pars$t0 + z^2
    dens <- exp(ddm_var_pointwise(ddm_var_data(t, rep(R, length(t))), pars)[1, ])
    total <- total + sum((zh - zl) / 2 * g$weights * 2 * z * dens)
  }
  total
}

test_that("DDM binned likelihood with st0 integrates the unbinned density over the bin", {
  h <- 1 / 60
  pars <- data.frame(v = 1.2, a = 1.1, t0 = 0.255, Z = 0.45, sv = 0.6, SZ = 0.5, st0 = 0.12)
  # Grid bins straddling onset, just above it, across the end of the onset
  # window (t0 + st0), and in the body.
  edges <- h * (floor(pars$t0 / h) + c(0, 1, 3, 7, 8, 10, 30))
  R <- rep(c("lower", "upper"), length.out = length(edges))
  binned <- ddm_var_pointwise(ddm_var_data(edges + h / 2, R, rt_resolution = h), pars)[1, ]
  reference <- vapply(seq_along(edges), function(i)
    log(ddm_var_integral(edges[i], edges[i] + h, R[i], pars) / h), numeric(1))
  expect_lt(max(abs(binned - reference)), 5e-4)
})

test_that("DDM binned likelihood with st0 is smooth in t0", {
  # The sweep takes both sides of the bin through their exact, blended and
  # corner-corrected forms, and onset into the bin.
  x <- ddm_var_data(0.42, "upper", rt_resolution = 1 / 60)
  pars <- data.frame(v = 1.2, a = 1.1, t0 = seq(0.2, 0.41, length.out = 2001),
                     Z = 0.45, sv = 0.6, SZ = 0.5, st0 = 0.12)
  ll <- ddm_var_pointwise(x, pars)[, 1]
  expect_true(all(is.finite(ll)))
  expect_lt(max(abs(diff(ll, differences = 3))), 1e-4)
})

test_that("DDM censored likelihood matches the integrated density", {
  pars <- data.frame(v = 1.2, a = 1.1, t0 = 0.255, Z = 0.45, sv = 0.6, SZ = 0.5, st0 = 0.12)
  LC <- 0.45; UC <- 0.8
  dat <- make_missing(ddm_var_frame(c(0.4, 0.6, 1.0), c("lower", "upper", "lower")),
                      LC = LC, UC = UC, verbose = FALSE)
  ll <- ddm_var_pointwise(ddm_var_data(dat = dat), pars)[1, ]
  below <- ddm_var_integral(0, LC, "lower", pars) + ddm_var_integral(0, LC, "upper", pars)
  inside <- ddm_var_integral(LC, UC, "lower", pars) + ddm_var_integral(LC, UC, "upper", pars)
  expect_lt(abs(ll[1] - log(below)), 1e-4)
  expect_lt(abs(ll[3] - log(1 - below - inside)), 1e-4)
})
