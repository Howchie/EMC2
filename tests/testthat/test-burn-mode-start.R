local_rng_guard()

test_that("Gaussian subject optimisation recovers the exact mode and covariance", {
  prior <- list(mu = c(-1, 0.5),
                var = matrix(c(2, 0.3, 0.3, 1), 2))
  observation <- c(1.2, -0.7)
  likelihood_var <- matrix(c(0.6, 0.1, 0.1, 0.4), 2)
  precision <- solve(prior$var) + solve(likelihood_var)
  exact_covariance <- solve(precision)
  exact_mode <- drop(exact_covariance %*%
                       (solve(prior$var, prior$mu) +
                          solve(likelihood_var, observation)))
  ll_batch <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    d <- sweep(x, 2L, observation)
    -0.5 * rowSums((d %*% solve(likelihood_var)) * d)
  }
  result <- EMC2:::.emc_optimize_subject(c(4, 3), prior, ll_batch)
  expect_true(result$converged)
  expect_identical(result$reason, "ok")
  expect_equal(result$alpha, exact_mode, tolerance = 1e-4)
  expect_equal(result$covariance, exact_covariance, tolerance = 1e-5)
})

test_that("unstable mixed differences retain a positive proposal covariance", {
  ll_batch <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    -0.5 * rowSums(x^2) +
      1e8 * x[, 1L] * x[, 2L] * rowSums(x^2)
  }
  covariance <- EMC2:::.emc_burn_hessian(
    c(0, 0), diag(2), diag(2), ll_batch)
  expect_true(EMC2:::.emc_empirical_covariance_ok(covariance))
  expect_equal(diag(covariance), c(0.5, 0.5), tolerance = 1e-5)
})

test_that("the second burn proposal window excludes every pre-mode draw", {
  chains <- lapply(seq_len(3L), function(i)
    list(samples = list(idx = 151L), burn_mode_index = 51L))
  expect_identical(EMC2:::.emc_chain_proposal_window(
    chains, stage = "burn"), 51:151)
  expect_identical(EMC2:::.emc_chain_proposal_window(
    chains, samples_idx = 45:75, stage = "burn"), 51:75)
  expect_identical(EMC2:::.emc_chain_proposal_window(
    chains, samples_idx = 45:75, stage = "sample"), 45:75)
})

test_that("proposal window indices select those exact stored draws", {
  chain <- list(type = "single", samples = list(
    stage = c("init", rep("preburn", 50L), rep("burn", 100L)),
    subj_ll = matrix(seq_len(151L), nrow = 1L), idx = 151L))
  filter <- EMC2:::.emc_proposal_filter_indices(list(chain), 51:151)
  expect_identical(filter, 50:150)
  actual <- EMC2::get_pars(list(chain), selection = "LL",
                           stage = c("preburn", "burn"),
                           filter = filter, merge_chains = TRUE,
                           return_mcmc = FALSE)
  expect_identical(as.numeric(actual), as.numeric(51:151))
  trimmed <- chain
  trimmed$samples$stage <- trimmed$samples$stage[-1L]
  trimmed$samples$idx <- 150L
  expect_identical(EMC2:::.emc_proposal_filter_indices(
    list(trimmed), 50:150), 50:150)
})

test_that("production entry preserves the adapted epsilon for every subject", {
  pm <- lapply(seq_len(3L), function(s)
    list(list(epsilon = c(0.003 * s, 0.04 * s),
              mix = c(0.2, 0.3, 0.5), iter = 100L)))
  chain <- list(samples = structure(list(stage = "adapt"), pm_settings = pm))
  before <- lapply(pm, function(x) x[[1L]]$epsilon)
  after_reset <- EMC2:::reset_pm_settings(list(chain), "sample")
  expect_equal(lapply(attr(after_reset[[1L]]$samples, "pm_settings"),
                      function(x) x[[1L]]$epsilon), before)
  expect_equal(lapply(attr(after_reset[[1L]]$samples, "pm_settings"),
                      function(x) x[[1L]]$mix),
               lapply(pm, function(x) x[[1L]]$mix))
})

test_that("a t0 pinned just below the fastest RT keeps its Laplace covariance", {
  # Wald-like leading edge: log density falls as -c / (rt - t0) with
  # t0 = exp(x1 + x2). A prior-scaled four-point stencil is indefinite on
  # this ridge, and shrinking its mixed terms collapsed the covariance to
  # ~3e-4 of the true variance along it (wide RDM fixture, 9 of 30 subjects).
  rt <- 0.2505
  ll_batch <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    dt <- rt - exp(x[, 1L] + x[, 2L])
    v <- -0.5 * rowSums(x[, 3:4, drop = FALSE]^2) - 0.02 / pmax(dt, 1e-300) +
      400 * (x[, 1L] + x[, 2L])
    v[dt <= 0] <- -Inf
    v
  }
  prior <- list(mu = c(rep(log(0.25) / 2, 2L), 0, 0), var = diag(0.25, 4L))
  result <- EMC2:::.emc_optimize_subject(
    c(rep(log(0.2) / 2, 2L), 0.1, 0.1), prior, ll_batch)
  expect_identical(result$reason, "ok")
  hessian <- solve(prior$var) +
    stats::optimHess(result$alpha, function(z) -ll_batch(z))
  e <- eigen(hessian, symmetric = TRUE)
  root <- e$vectors %*% diag(sqrt(e$values))
  whitened <- eigen(t(root) %*% result$covariance %*% root,
                    symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(abs(whitened - 1) < 0.05))
})

test_that("proposal covariances can use a wider group prior than the mode", {
  # Mode under the actual prior; covariance under a 4x wider one, so weakly
  # informed directions are not tied to a shrinking group variance.
  prior <- list(mu = c(0, 0), var = diag(c(0.01, 1)))
  likelihood_var <- diag(c(100, 0.1))
  ll_batch <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    d <- sweep(x, 2L, c(1, 1))
    -0.5 * rowSums((d %*% solve(likelihood_var)) * d)
  }
  base <- EMC2:::.emc_optimize_subject(c(0.5, 0.5), prior, ll_batch)
  wide <- EMC2:::.emc_optimize_subject(c(0.5, 0.5), prior, ll_batch,
                                       proposal_prior_scale = 4)
  expect_equal(wide$alpha, base$alpha, tolerance = 1e-5)
  expect_equal(wide$covariance,
               solve(solve(prior$var) / 4 + solve(likelihood_var)),
               tolerance = 1e-4)
})
