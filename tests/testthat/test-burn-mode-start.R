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
  # Its mean is the mode of the same widened density (exact for a Gaussian
  # likelihood), and the conditional covariance is kept for starts.
  P4 <- solve(prior$var) / 4; Pl <- solve(likelihood_var)
  expect_equal(wide$proposal_mode,
               drop(solve(P4 + Pl, P4 %*% prior$mu + Pl %*% c(1, 1))),
               tolerance = 1e-4)
  with_start <- EMC2:::.emc_optimize_subject(c(0.5, 0.5), prior, ll_batch,
                                             proposal_prior_scale = 4,
                                             start_covariance = TRUE)
  expect_equal(with_start$start_covariance,
               solve(solve(prior$var) + solve(likelihood_var)), tolerance = 1e-4)
  expect_equal(base$proposal_mode, base$alpha)
})

test_that("the batched gradient-difference Hessian matches stats::optimHess", {
  set.seed(11)
  for (p in c(1L, 3L, 12L)) {
    A <- crossprod(matrix(rnorm(p * p), p)) / p + diag(p)
    ll_batch <- function(x) {
      if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
      -0.5 * rowSums((x %*% A) * x) - rowSums(sin(x)^2) + rowSums(x^3) / 10
    }
    x <- rnorm(p) * 0.3
    reference <- stats::optimHess(x, function(z) -ll_batch(z))
    # Several likelihood batches, one split inside a pair's four points.
    for (max_rows in c(4096L, 7L)) {
      expect_equal(EMC2:::.emc_gradient_difference_hessian(
        x, ll_batch, max_rows = max_rows), reference,
        tolerance = 1e-8, ignore_attr = TRUE)
    }
  }
  outside <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    v <- -rowSums(x^2)
    v[x[, 1L] > 0.0015] <- -Inf
    v
  }
  expect_null(EMC2:::.emc_gradient_difference_hessian(c(0.001, 0), outside))
})

test_that("subjects with more than 40 parameters still get a mode covariance", {
  p <- 45L
  prior <- list(mu = rep(0, p), var = diag(p))
  likelihood_var <- diag(seq(0.2, 1, length.out = p))
  observation <- seq(-1, 1, length.out = p)
  ll_batch <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    d <- sweep(x, 2L, observation)
    -0.5 * rowSums(sweep(d^2, 2L, diag(likelihood_var), "/"))
  }
  result <- EMC2:::.emc_optimize_subject(rep(0.5, p), prior, ll_batch)
  expect_identical(result$reason, "ok")
  expect_equal(result$covariance,
               solve(solve(prior$var) + solve(likelihood_var)),
               tolerance = 1e-4)
})

test_that("the optimiser never re-evaluates the point it just evaluated", {
  # optim() calls the gradient at the point it has just evaluated; the smooth
  # Gaussian also never needs a one-sided difference, so every single-point
  # likelihood call is a new point.
  prior <- list(mu = c(-1, 0.5), var = matrix(c(2, 0.3, 0.3, 1), 2))
  seen <- list()
  ll_batch <- function(x) {
    if (is.null(dim(x))) {
      seen[[length(seen) + 1L]] <<- x
      x <- matrix(x, nrow = 1L)
    }
    -0.5 * rowSums(sweep(x, 2L, c(1.2, -0.7))^2 / 0.5)
  }
  result <- EMC2:::.emc_optimize_subject(c(4, 3), prior, ll_batch)
  expect_identical(result$reason, "ok")
  repeats <- vapply(seq_along(seen)[-1L], function(i)
    identical(seen[[i]], seen[[i - 1L]]), logical(1))
  expect_false(any(repeats))
})

test_that("joint models take per-model gradients and Hessians", {
  design_in <- get_design(samples_LNR)[[1]]
  data_in <- get_data(samples_LNR)
  joint <- suppressMessages(make_emc(
    list(data_in, data_in), list(a = design_in, b = design_in),
    prior_list = prior(list(a = design_in, b = design_in)), n_chains = 1))
  sampler <- joint[[1]]
  blocks <- EMC2:::.emc_joint_ll_blocks(sampler, 1L)
  expect_length(blocks, 2L)
  p <- sampler$n_pars
  expect_setequal(unlist(lapply(blocks, `[[`, "idx")), seq_len(p))
  par_names <- rownames(sampler$samples$alpha)
  ll_batch <- function(x) EMC2:::.emc_burn_ll_batch(
    x, sampler$data[[1L]], sampler$model, par_names)
  # A different point in each model, so a swapped block would show.
  x <- as.numeric(get_pars(samples_LNR, selection = "alpha", stage = "sample",
                           return_mcmc = FALSE, merge_chains = TRUE)[, 1L, 1L])
  x <- c(x, x + c(0.05, -0.05))[seq_len(p)]
  points <- rbind(x, x + 0.01, x - 0.01)
  expect_equal(ll_batch(points), Reduce(`+`, lapply(blocks, function(b)
    b$ll(points[, b$idx, drop = FALSE]))), tolerance = 1e-12)
  prior_var <- diag(p)
  full <- EMC2:::.emc_burn_hessian(x, prior_var, diag(p), ll_batch)
  blocked <- EMC2:::.emc_burn_hessian(x, prior_var, diag(p), ll_batch,
                                      blocks = blocks)
  expect_true(is.matrix(blocked))
  expect_equal(blocked, full, tolerance = 1e-6)
  expect_null(EMC2:::.emc_joint_ll_blocks(
    suppressMessages(make_emc(data_in, design_in, n_chains = 1))[[1]], 1L))
})

test_that("burn ensemble normal retains the wide Laplace metric", {
  emc <- samples_LNR
  covariance <- diag(c(0.2, 0.4, 0.6, 0.8))
  for (j in seq_along(emc)) {
    emc[[j]]$samples$stage[] <- "burn"
    emc[[j]]$burn_mode_index <- 1L
    n_subjects <- dim(emc[[j]]$samples$alpha)[2L]
    n_pars <- dim(emc[[j]]$samples$alpha)[1L]
    emc[[j]]$burn_mode_var <- rep(list(covariance), n_subjects)
    emc[[j]]$burn_mode_mu <- rep(list(rep(0, n_pars)), n_subjects)
  }
  out <- EMC2:::create_chain_proposals(emc, samples_idx = 1:50,
                                        stage = "burn")
  expect_equal(out[[1L]]$chains_var[[1L]], covariance)
  expect_equal(out[[1L]]$ensemble_q$var[[1L]], 3 * covariance)
})

test_that("Laplace importance starts are posterior draws with a Pareto diagnostic", {
  set.seed(21)
  # Gaussian posterior: the Laplace is exact, so every weight is equal (k low)
  # and the resampled starts have the posterior's moments.
  prior <- list(mu = c(0, 0), var = diag(2))
  likelihood_var <- diag(c(0.5, 0.2))
  ll_batch <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    -0.5 * rowSums(sweep(sweep(x, 2L, c(1, -1))^2, 2L, diag(likelihood_var), "/"))
  }
  post_var <- solve(diag(2) + solve(likelihood_var))
  post_mode <- drop(post_var %*% solve(likelihood_var, c(1, -1)))
  starts <- t(replicate(400, {
    out <- EMC2:::.emc_laplace_is_start(post_mode, post_var, prior, ll_batch)
    c(out$start, out$pareto_k)
  }))
  expect_true(all(starts[, 3L] < 0.5))
  expect_equal(colMeans(starts[, 1:2]), post_mode, tolerance = 0.05)
  expect_equal(cov(starts[, 1:2]), post_var, tolerance = 0.2)
  # A posterior cut off just above the mode (support boundary): the Gaussian
  # puts half its mass outside and every start stays inside the support.
  edge <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    v <- ll_batch(x)
    v[x[, 1L] > post_mode[1L] + 0.01] <- -Inf
    v
  }
  out <- replicate(50, unlist(EMC2:::.emc_laplace_is_start(
    post_mode, post_var, prior, edge)[c("start", "outside")]))
  expect_true(all(out[1L, ] <= post_mode[1L] + 0.01))
  expect_equal(mean(out[3L, ]), 0.5 - 0.01 / sqrt(post_var[1L, 1L]) * dnorm(0),
               tolerance = 0.05)
  # Nothing inside the support: no start, and a support-failure status.
  none <- EMC2:::.emc_laplace_is_start(
    post_mode, post_var, prior, function(x) rep(-Inf, NROW(x)))
  expect_null(none$start)
  expect_identical(none$status, "support")
})

test_that("Laplace gate treats infinite, missing and outside-support checks as poor", {
  st <- list(pareto_k = c(0.2, Inf, NA, 0.2, 0.9),
             outside = c(0, 0, 0.9, 0.9, NA))
  expect_identical(EMC2:::.emc_laplace_poor(st), c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_identical(EMC2:::.emc_laplace_poor(list()), logical(0))
  # A diagnostic with no finite weights still reports, as a support failure.
  set.seed(1)
  prior <- list(mu = c(0, 0), var = diag(2))
  out <- EMC2:::.emc_laplace_is_start(c(0, 0), diag(2), prior,
                                      function(x) rep(-Inf, nrow(x)), n_draws = 32L)
  expect_identical(out$status, "support")
  expect_true(is.na(out$pareto_k))
  expect_equal(out$outside, 1)
  expect_true(EMC2:::.emc_laplace_poor(out))
  # An errored check is recorded as unknown, never as the previous result.
  chk <- EMC2:::.emc_laplace_check(NULL)
  expect_identical(chk$status, "error")
  expect_true(EMC2:::.emc_laplace_poor(chk))
})

test_that("a poor Laplace fit swaps the ensemble pool mean and covariance together", {
  emc <- samples_LNR
  covariance <- diag(c(0.2, 0.4, 0.6, 0.8))
  for (j in seq_along(emc)) {
    emc[[j]]$samples$stage[] <- "burn"
    emc[[j]]$burn_mode_index <- 1L
    n_subjects <- dim(emc[[j]]$samples$alpha)[2L]
    n_pars <- dim(emc[[j]]$samples$alpha)[1L]
    emc[[j]]$burn_mode_var <- rep(list(covariance), n_subjects)
    emc[[j]]$burn_mode_mu <- rep(list(rep(-5, n_pars)), n_subjects)
    emc[[j]]$burn_mode_stats <- list(pareto_k = c(Inf, rep(0.1, n_subjects - 1L)),
                                     outside = rep(0, n_subjects))
  }
  out <- EMC2:::create_chain_proposals(emc, samples_idx = 1:50, stage = "burn")
  q <- out[[1L]]$ensemble_q
  # Subject 1 (k = Inf): window moments for both.
  expect_false(isTRUE(all.equal(q$var[[1L]], 3 * covariance)))
  expect_true(all(abs(q$mu[[1L]] + 5) > 1))
  # Subject 2 (good fit): the Laplace pair.
  expect_equal(q$var[[2L]], 3 * covariance)
  expect_equal(q$mu[[2L]], rep(-5, n_pars))
})

test_that("poor Laplace covariances use the window in every sampling stage", {
  covariance <- diag(c(0.2, 0.4, 0.6, 0.8))
  for (stage in c("burn", "adapt", "sample")) {
    emc <- samples_LNR
    for (j in seq_along(emc)) {
      emc[[j]]$samples$stage[] <- stage
      n_subjects <- dim(emc[[j]]$samples$alpha)[2L]
      emc[[j]]$burn_mode_var <- rep(list(covariance), n_subjects)
      emc[[j]]$burn_mode_mu <- rep(list(rep(-5, 4)), n_subjects)
      emc[[j]]$burn_mode_stats <- list(pareto_k = c(Inf, rep(0.1, n_subjects - 1L)),
                                      outside = rep(0, n_subjects))
    }
    reference <- emc
    for (j in seq_along(reference)) reference[[j]]$burn_mode_var <- NULL
    out <- EMC2:::create_chain_proposals(emc, samples_idx = 1:50, stage = stage)
    expected <- EMC2:::create_chain_proposals(reference, samples_idx = 1:50, stage = stage)
    expect_equal(out[[1]]$chains_var[[1]], expected[[1]]$chains_var[[1]])
    expect_equal(out[[1]]$chains_mu[[1]], expected[[1]]$chains_mu[[1]])
    expect_equal(out[[1]]$chains_var[[2]], covariance)
  }
})
