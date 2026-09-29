# The conditional importance-sampling particle update must leave the subject
# posterior invariant. Two past defects broke this while every chain still
# looked converged: (1) numbers_from_proportion() forced at least one particle
# per mixture component, so rare components were over-sampled relative to the
# mixture density used in the weights; (2) proposal components centred on the
# current state made the reference particle's proposal density depend on
# itself. Both inflated subject-level spread, and so the group variance, by an
# amount that depended on each fit's adapted proposal settings.

local_rng_guard()

test_that("particle allocation is plain multinomial (no forced minimum)", {
  set.seed(1)
  counts <- replicate(200, EMC2:::numbers_from_proportion(c(.98, .01, .01), 5))
  expect_true(any(counts[2, ] == 0))
  expect_true(all(colSums(counts) == 5))
})

test_that("burn runs the independent and local proposals as two IS steps", {
  calls <- 0L
  local_mocked_bindings(
    calc_ll_pooled = function(proposals, ...) {
      calls <<- calls + 1L
      rep(0, nrow(proposals))
    },
    .package = "EMC2")
  set.seed(124L)
  p <- 24L
  alpha <- stats::setNames(rep(0, p), paste0("p", seq_len(p)))
  settings <- list(list(epsilon = c(0.1, 0.1), mix = c(0.2, 0.4, 0.4),
                        n_particles = 100L, iter = 0L,
                        proposal_counts = integer(3L),
                        acc_counts = integer(3L), gd_good = FALSE))
  tune <- list(components = rep(1L, p), shared_ll_idx = rep(1L, p), n0 = Inf)
  result <- EMC2:::new_particle(
    1L, NULL, settings, chains_var = diag(p), prev_ll = 0,
    parameters = NULL, stage = "burn", type = "standard", tune = tune,
    current_alpha = alpha, population_mu = alpha, population_var = diag(p))
  expect_equal(calls, 2L)
  expect_true(all(is.finite(result$proposal)))
})

test_that("independence proposals are never shrunk by the local step size", {
  # A scale below one on a fitted covariance removes an independence
  # proposal's coverage in high dimension; only local steps are tuned.
  seen <- NULL
  local_mocked_bindings(
    calc_ll_pooled = function(proposals, ...) {
      if (is.null(seen)) seen <<- proposals
      rep(0, nrow(proposals))
    },
    .package = "EMC2")
  set.seed(7L)
  p <- 6L
  alpha <- stats::setNames(rep(0, p), paste0("p", seq_len(p)))
  settings <- list(list(epsilon = c(0.1, 0.1, 0.1), mix = c(0, 0, 1, 0),
                        n_particles = 800L, iter = 0L,
                        proposal_counts = integer(4L),
                        acc_counts = integer(4L), gd_good = FALSE))
  tune <- list(components = rep(1L, p), shared_ll_idx = rep(1L, p), n0 = Inf)
  result <- EMC2:::new_particle(
    1L, NULL, settings, chains_mu = alpha, chains_var = diag(p),
    eff_mu = alpha, eff_var = diag(p), prev_ll = 0, parameters = NULL,
    stage = "sample", type = "standard", tune = tune,
    current_alpha = alpha, population_mu = alpha, population_var = diag(p))
  sds <- apply(seen[-1L, , drop = FALSE], 2L, stats::sd)
  expect_true(all(sds > 0.85 & sds < 1.15))
  # Selection diagnostics cover every component of the production kernel.
  expect_length(result$pm_settings[[1L]]$sel$sample$n, 4L)
})

test_that("adapt tunes the full production mixture once efficient proposals exist", {
  local_mocked_bindings(
    calc_ll_pooled = function(proposals, ...) rep(0, nrow(proposals)),
    .package = "EMC2")
  set.seed(8L)
  p <- 4L
  alpha <- stats::setNames(rep(0, p), paste0("p", seq_len(p)))
  settings <- list(list(epsilon = c(0.5, 0.5, 0.5),
                        mix = EMC2:::get_default_mix("adapt"),
                        n_particles = 50L, iter = 0L,
                        proposal_counts = integer(4L),
                        acc_counts = integer(4L), gd_good = FALSE))
  tune <- list(components = rep(1L, p), shared_ll_idx = rep(1L, p), n0 = Inf)
  result <- EMC2:::new_particle(
    1L, NULL, settings, chains_mu = alpha, chains_var = diag(p),
    eff_mu = alpha, eff_var = diag(p), prev_ll = 0, parameters = NULL,
    stage = "adapt", type = "standard", tune = tune,
    current_alpha = alpha, population_mu = alpha, population_var = diag(p))
  expect_length(EMC2:::get_default_mix("adapt"), 4L)
  expect_length(result$pm_settings[[1L]]$sel$adapt$n, 4L)
})

test_that("new_particle leaves a conjugate Gaussian posterior invariant", {
  skip_model_validation()
  ybar <- c(1.5, -0.5); Li <- diag(1 / c(0.3, 0.5))
  mu <- c(0, 0); Sig <- matrix(c(1, .3, .3, 1), 2)
  post_var <- solve(solve(Sig) + Li)
  post_mu <- as.vector(post_var %*% (solve(Sig) %*% mu + Li %*% ybar))
  local_mocked_bindings(
    calc_ll_pooled = function(proposals, ...) {
      d <- sweep(proposals, 2, ybar)
      -0.5 * rowSums((d %*% Li) * d)
    },
    .package = "EMC2")
  bm_se <- function(x, nb = 40) {
    b <- split(x, cut(seq_along(x), nb)); sd(vapply(b, mean, 0)) / sqrt(nb)
  }
  run <- function(stage, mix, eps, iters = 20000) {
    set.seed(3)
    pm <- list(list(epsilon = eps, mix = mix, n_particles = 5, iter = 0,
                    proposal_counts = rep(0, length(mix)),
                    acc_counts = rep(0, length(mix)), gd_good = FALSE))
    tune <- list(components = c(1, 1), shared_ll_idx = c(1, 1), n0 = Inf)
    a <- c(a = 0, b = 0)
    ll <- -0.5 * sum((a - ybar)^2 * diag(Li))
    out <- matrix(NA_real_, iters, 2)
    for (t in seq_len(iters)) {
      r <- EMC2:::new_particle(
        1, NULL, pm, eff_mu = c(.8, -.2), eff_var = diag(2) * .3,
        chains_mu = c(.8, -.2), chains_var = diag(2) * .2, prev_ll = ll,
        parameters = NULL, stage = stage, type = "standard", tune = tune,
        current_alpha = a, population_mu = mu, population_var = Sig)
      a <- r$proposal; ll <- r$ll; out[t, ] <- a
    }
    out[-(1:1000), 1]
  }
  # Mixes with a state-centred component and a rare component: the old
  # kernel missed the mean here by z = -18 (sample) and z = -22 (adapt).
  for (cfg in list(list("sample", c(.2, .4, .2, .2), c(1, 1, 1)),
                   list("adapt", c(.3, .4, .3), c(1, 1)),
                   list("adapt", c(.2, .4, .2, .2), c(.3, .3, .3)),
                   list("sample", c(.3, .02, .3, .38), c(1, 1, 1)))) {
    x <- run(cfg[[1]], cfg[[2]], cfg[[3]])
    expect_lt(abs(mean(x) - post_mu[1]) / bm_se(x), 4.5)
    expect_lt(abs(var(x) - post_var[1, 1]) / bm_se((x - mean(x))^2), 4.5)
  }
})
