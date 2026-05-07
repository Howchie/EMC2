rdm_gl <- function(n = 80) {
  skip_if_not_installed("statmod")
  statmod::gauss.quad(n, "legendre")
}

rdm_q_integrate <- function(mu, sv, f, lower = -Inf, n = 80) {
  gl <- rdm_gl(n)
  lo <- if (is.finite(lower)) pnorm((lower - mu) / sv) else 0
  lo <- max(0, min(lo, 1 - 1e-15))
  0.5 * sum(gl$weights * vapply(gl$nodes, function(node) {
    u <- 0.5 * (node + 1)
    p <- lo + (1 - lo) * u
    drift <- mu + sv * qnorm(p)
    f(drift)
  }, numeric(1)))
}

expect_local_derivative <- function(pdf, cdf, t, tol = 4e-4) {
  h <- 1e-5
  slope <- (cdf(t + h) - cdf(t - h)) / (2 * h)
  expect_equal(slope, pdf(t), tolerance = tol)
}

test_that("positive-drift RDMSWTN A>0 matches q-transform drift integral", {
  cases <- list(
    list(label = "no_erlang", lambda_g = 0, lambda_k = 0, guess = FALSE),
    list(label = "kill",      lambda_g = 0, lambda_k = 0.35, guess = FALSE),
    list(label = "guess",     lambda_g = 0.45, lambda_k = 0, guess = TRUE)
  )

  pars <- list(t = 0.95, mu = -0.25, b = 1.45, A = 0.35,
               s = 0.9, t0 = 0.12, sv = 0.55)

  for (case in cases) {
    ref_p <- rdm_q_integrate(pars$mu, pars$sv, function(drift) {
      EMC2:::pwald(pars$t, drift, pars$b, A = pars$A, sigma = pars$s,
                   t0 = pars$t0, lambda_g = case$lambda_g,
                   lambda_k = case$lambda_k, guess = case$guess,
                   posdrift = FALSE)
    }, lower = 0, n = 80)
    ref_d <- rdm_q_integrate(pars$mu, pars$sv, function(drift) {
      EMC2:::dwald(pars$t, drift, pars$b, A = pars$A, sigma = pars$s,
                   t0 = pars$t0, lambda_g = case$lambda_g,
                   lambda_k = case$lambda_k, guess = case$guess,
                   posdrift = FALSE)
    }, lower = 0, n = 80)

    got_p <- EMC2:::prdmswtn(pars$t, pars$mu, pars$b, pars$A, s = pars$s,
                             t0 = pars$t0, sv = pars$sv,
                             lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                             n_gauss_nodes = 80, guess = case$guess,
                             posdrift = TRUE)
    got_d <- EMC2:::drdmswtn(pars$t, pars$mu, pars$b, pars$A, s = pars$s,
                             t0 = pars$t0, sv = pars$sv,
                             lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                             n_gauss_nodes = 80, guess = case$guess,
                             posdrift = TRUE)

    expect_equal(got_p, ref_p, tolerance = 1e-5, info = case$label)
    expect_equal(got_d, ref_d, tolerance = 1e-5, info = case$label)
  }
})

test_that("defective RDMSWTN A>0 matches independent threshold integration", {
  gl <- rdm_gl(100)
  pars <- list(t = 1.05, mu = -0.35, b = 1.4, A = 0.3,
               s = 1.0, t0 = 0.15, sv = 0.5, lambda_k = 0.25)
  center <- pars$b - 0.5 * pars$A
  half_width <- 0.5 * pars$A

  ref_p <- 0.5 * sum(gl$weights * vapply(gl$nodes, function(node) {
    threshold <- center + half_width * node
    EMC2:::pswtn(pars$t, pars$mu, threshold, s = pars$s, t0 = pars$t0,
                 sv = pars$sv, lambda_k = pars$lambda_k, posdrift = FALSE)
  }, numeric(1)))
  ref_d <- 0.5 * sum(gl$weights * vapply(gl$nodes, function(node) {
    threshold <- center + half_width * node
    EMC2:::dswtn(pars$t, pars$mu, threshold, s = pars$s, t0 = pars$t0,
                 sv = pars$sv, lambda_k = pars$lambda_k, posdrift = FALSE)
  }, numeric(1)))

  got_p <- EMC2:::prdmswtn(pars$t, pars$mu, pars$b, pars$A, s = pars$s, t0 = pars$t0,
                           sv = pars$sv, lambda_k = pars$lambda_k,
                           n_gauss_nodes = 80, posdrift = FALSE)
  expect_lt(abs(got_p - ref_p), 1e-4)
  expect_equal(
    EMC2:::drdmswtn(pars$t, pars$mu, pars$b, pars$A, s = pars$s, t0 = pars$t0,
                    sv = pars$sv, lambda_k = pars$lambda_k,
                    n_gauss_nodes = 80, posdrift = FALSE),
    ref_d,
    tolerance = 1e-4
  )
})

test_that("SWTN and RDMSWTN CDF derivatives match densities across variants", {
  variants <- list(
    list(posdrift = FALSE, lambda_g = 0,    lambda_k = 0,    guess = FALSE, tol = 4e-5),
    list(posdrift = TRUE,  lambda_g = 0,    lambda_k = 0,    guess = FALSE, tol = 2e-4),
    list(posdrift = FALSE, lambda_g = 0,    lambda_k = 0.45, guess = FALSE, tol = 3e-4),
    list(posdrift = TRUE,  lambda_g = 0,    lambda_k = 0.45, guess = FALSE, tol = 5e-4),
    list(posdrift = TRUE,  lambda_g = 0.35, lambda_k = 0,    guess = TRUE,  tol = 5e-4)
  )

  for (case in variants) {
    expect_local_derivative(
      pdf = function(t) EMC2:::dswtn(t, -0.15, 1.25, s = 1.0, t0 = 0.1, sv = 0.45,
                                     lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                                     guess = case$guess, posdrift = case$posdrift),
      cdf = function(t) EMC2:::pswtn(t, -0.15, 1.25, s = 1.0, t0 = 0.1, sv = 0.45,
                                     lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                                     guess = case$guess, posdrift = case$posdrift),
      t = 0.9,
      tol = case$tol
    )

    expect_local_derivative(
      pdf = function(t) EMC2:::drdmswtn(t, -0.15, 1.35, 0.3, s = 1.0, t0 = 0.1, sv = 0.45,
                                        lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                                        guess = case$guess, posdrift = case$posdrift),
      cdf = function(t) EMC2:::prdmswtn(t, -0.15, 1.35, 0.3, s = 1.0, t0 = 0.1, sv = 0.45,
                                        lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                                        guess = case$guess, posdrift = case$posdrift),
      t = 0.9,
      tol = case$tol
    )
  }
})

test_that("density integrals recover finite CDF mass and infinite tail ordering", {
  cases <- list(
    list(label = "defective", posdrift = FALSE, lambda_g = 0, lambda_k = 0, guess = FALSE, upper = 80, tol = 5e-3),
    list(label = "positive_trunc", posdrift = TRUE, lambda_g = 0, lambda_k = 0, guess = FALSE, upper = 80, tol = 6e-3),
    list(label = "kill", posdrift = TRUE, lambda_g = 0, lambda_k = 0.5, guess = FALSE, upper = 20, tol = 5e-4),
    list(label = "guess", posdrift = TRUE, lambda_g = 0.45, lambda_k = 0, guess = TRUE, upper = 30, tol = 1e-4)
  )

  for (case in cases) {
    swtn_int <- integrate(function(x) {
      vapply(x, function(t) EMC2:::dswtn(t, -0.2, 1.2, s = 1.0, t0 = 0.1, sv = 0.4,
                                         lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                                         guess = case$guess, posdrift = case$posdrift), numeric(1))
    }, lower = 0, upper = case$upper, subdivisions = 500L, rel.tol = 1e-6)$value
    swtn_mass <- EMC2:::pswtn(case$upper, -0.2, 1.2, s = 1.0, t0 = 0.1, sv = 0.4,
                              lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                              guess = case$guess, posdrift = case$posdrift)
    swtn_inf <- EMC2:::pswtn(Inf, -0.2, 1.2, s = 1.0, t0 = 0.1, sv = 0.4,
                             lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                             guess = case$guess, posdrift = case$posdrift)

    rdm_int <- integrate(function(x) {
      vapply(x, function(t) EMC2:::drdmswtn(t, -0.2, 1.3, 0.25, s = 1.0, t0 = 0.1, sv = 0.4,
                                           lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                                           guess = case$guess, posdrift = case$posdrift), numeric(1))
    }, lower = 0, upper = case$upper, subdivisions = 500L, rel.tol = 1e-6)$value
    rdm_mass <- EMC2:::prdmswtn(case$upper, -0.2, 1.3, 0.25, s = 1.0, t0 = 0.1, sv = 0.4,
                                lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                                guess = case$guess, posdrift = case$posdrift)
    rdm_inf <- EMC2:::prdmswtn(Inf, -0.2, 1.3, 0.25, s = 1.0, t0 = 0.1, sv = 0.4,
                               lambda_g = case$lambda_g, lambda_k = case$lambda_k,
                               guess = case$guess, posdrift = case$posdrift)

    expect_equal(swtn_int, swtn_mass, tolerance = case$tol, info = paste(case$label, "SWTN"))
    expect_equal(rdm_int, rdm_mass, tolerance = case$tol, info = paste(case$label, "RDMSWTN"))
    expect_lte(swtn_mass, swtn_inf + 1e-12)
    expect_lte(rdm_mass, rdm_inf + 1e-12)
  }
})

test_that("combined guess makes defective negative-drift response distribution proper", {
  p_decision <- EMC2:::prdmswtn(Inf, -0.8, 1.3, 0.25, s = 1.0, t0 = 0.1,
                                sv = 0.0, lambda_g = 0, lambda_k = 0,
                                guess = FALSE, posdrift = FALSE)
  p_guess <- EMC2:::prdmswtn(Inf, -0.8, 1.3, 0.25, s = 1.0, t0 = 0.1,
                             sv = 0.0, lambda_g = 0.5, lambda_k = 0,
                             guess = TRUE, posdrift = FALSE)
  p_kill <- EMC2:::prdmswtn(Inf, -0.8, 1.3, 0.25, s = 1.0, t0 = 0.1,
                            sv = 0.0, lambda_g = 0, lambda_k = 0.5,
                            guess = FALSE, posdrift = FALSE)

  expect_gt(p_decision, 0)
  expect_lt(p_decision, 1)
  expect_equal(p_guess, 1, tolerance = 1e-12)
  expect_gt(p_kill, 0)
  expect_lt(p_kill, p_decision)
})

test_that("posdrift fixed zero drift is invalid before Erlang clocks", {
  scalar_calls <- list(
    dswtn_guess = function(log_out) EMC2:::dswtn(0.25, 0, 1.2, sv = 0, s = 1, t0 = 0.4,
                                                lambda_g = 0.7, guess = TRUE,
                                                posdrift = TRUE, log_out = log_out),
    pswtn_guess = function(log_out) EMC2:::pswtn(0.25, 0, 1.2, sv = 0, s = 1, t0 = 0.4,
                                                lambda_g = 0.7, guess = TRUE,
                                                posdrift = TRUE, log_out = log_out),
    drdmswtn_kill = function(log_out) EMC2:::drdmswtn(0.8, 0, 1.3, 0.2, sv = 0, s = 1, t0 = 0.1,
                                                      lambda_k = 0.6, guess = FALSE,
                                                      posdrift = TRUE, log_out = log_out),
    prdmswtn_guess = function(log_out) EMC2:::prdmswtn(Inf, 0, 1.3, 0.2, sv = 0, s = 1, t0 = 0.1,
                                                       lambda_g = 0.6, guess = TRUE,
                                                       posdrift = TRUE, log_out = log_out)
  )

  for (f in scalar_calls) {
    expect_identical(f(FALSE), 0)
    expect_equal(f(TRUE), -Inf)
  }
})
