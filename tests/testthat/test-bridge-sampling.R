# Log inverse-Wishart density from log determinants and an explicit inverse.
iwish_reference <- function(W, v, S) {
  p <- nrow(W)
  0.5 * v * as.numeric(determinant(S, logarithm = TRUE)$modulus) -
    0.5 * (v + p + 1) * as.numeric(determinant(W, logarithm = TRUE)$modulus) -
    0.5 * sum(diag(S %*% solve(W))) -
    sum(lgamma((v + 1 - seq_len(p)) / 2)) - 0.5 * v * p * log(2) -
    0.25 * p * (p - 1) * log(pi)
}

log_chol <- function(W) {
  L <- t(chol(W))
  diag(L) <- log(diag(L))
  L[lower.tri(L, diag = TRUE)]
}

test_that("inverse-Wishart log density is exact on well- and ill-conditioned matrices", {
  local_test_rng(1)
  for (p in c(1L, 2L, 5L, 9L, 13L)) {
    for (v in c(p - 0.5, p + 1, p + 10)) {
      for (condition in c(1, 1e6, 1e9, 1e12)) {
        correlation <- cov2cor(crossprod(matrix(rnorm(p * (p + 3)), p + 3)))
        sds <- exp(seq(0, log(condition) / 2, length.out = p))[sample(p)]
        W <- outer(sds, sds) * correlation
        scale <- exp(rnorm(p, 0, log(condition + 1) / 4))
        expected <- iwish_reference(W, v, diag(scale, p))
        expect_equal(EMC2:::log_diwish(W, v, diag(scale, p)), expected, tolerance = 1e-8)
        expect_equal(EMC2:::log_diwish_log_chol(log_chol(W), v, log(scale)), expected,
                     tolerance = 1e-8)
        expect_equal(EMC2:::log_diwish_log_chol(log_chol(solve(W)), v, log(scale),
                                                covariance = FALSE),
                     expected, tolerance = 1e-6)
      }
    }
  }
  # A full scale matrix.
  W <- matrix(c(2, 0.5, 0.5, 3), 2)
  S <- matrix(c(3, -1, -1, 2), 2)
  expect_equal(EMC2:::log_diwish(W, 4, S), iwish_reference(W, 4, S), tolerance = 1e-12)
})

test_that("inverse-Wishart log density is not repaired off the positive-definite cone", {
  expect_identical(EMC2:::log_diwish(matrix(c(1, 2, 2, 1), 2), 3, diag(2)), -Inf)
  expect_identical(EMC2:::log_diwish(diag(c(1, 0)), 3, diag(2)), -Inf)
  expect_identical(EMC2:::log_diwish(diag(2), 3, diag(c(1, -1))), -Inf)
  expect_identical(EMC2:::log_diwish_log_chol(c(0, NaN, 0), 3, c(0, 0)), -Inf)
  expect_identical(EMC2:::log_diwish_log_chol(c(-800, 0.5, 0), 3, c(0, 0)), -Inf)
  # Variances this far apart stay on the log scale.
  W <- diag(c(1, 1e-8))
  S <- diag(c(3, 1e-8))
  expect_equal(EMC2:::log_diwish(W, 4, S), iwish_reference(W, 4, S), tolerance = 1e-12)
  expect_equal(EMC2:::log_diwish_log_chol(c(300, 0, -300), 3, c(600, -600)),
               -1 - (lgamma(1.5) + lgamma(1) + 3 * log(2) + 0.5 * log(pi)),
               tolerance = 1e-12)
  expect_error(EMC2:::log_diwish(diag(2), 1, diag(2)), "degrees of freedom")
})

# Draws from the standard hierarchy's prior in the layout bridge sampling
# uses: subject parameters, then the group block.
standard_prior_draws <- function(n, par_group, n_subjects, v, A, mu_sd = 1) {
  p <- length(par_group)
  theta_mu <- matrix(rnorm(p * n, 0, mu_sd), p, n)
  a_half <- matrix(1 / rgamma(p * n, shape = 0.5, rate = 1 / A^2), p, n)
  theta_var <- array(0, c(p, p, n))
  alpha <- array(0, c(p, n_subjects, n))
  for (i in seq_len(n)) {
    for (group in unique(par_group)) {
      idx <- which(par_group == group)
      d <- length(idx)
      theta_var[idx, idx, i] <- if (d > 1) {
        solve(rWishart(1, v + d - 1, diag(a_half[idx, i] / (2 * v), d))[, , 1])
      } else {
        1 / rgamma(1, shape = v / 2, rate = v / a_half[idx, i])
      }
    }
    alpha[, , i] <- theta_mu[, i] +
      t(chol(theta_var[, , i])) %*% matrix(rnorm(p * n_subjects), p)
  }
  samples <- list(
    subjects = seq_len(n_subjects), par_names = paste0("p", seq_len(p)), n_pars = p,
    n_subjects = n_subjects, par_group = par_group, group_designs = NULL,
    is_blocked = duplicated(par_group) | duplicated(par_group, fromLast = TRUE),
    prior = list(theta_mu_mean = rep(0, p), theta_mu_var = diag(mu_sd^2, p), v = v,
                 A = rep(A, p)),
    samples = list(theta_mu = theta_mu, a_half = a_half, theta_var = theta_var)
  )
  group <- EMC2:::bridge_add_group(matrix(numeric(), n, 0L), samples, seq_len(n), "standard")
  list(draws = cbind(t(matrix(alpha, p * n_subjects, n)), group),
       info = EMC2:::make_info(samples, "standard"),
       data = vector("list", n_subjects))
}

test_that("prior draws under a constant likelihood have log marginal likelihood zero", {
  local_test_rng(11)
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, model, ...) numeric(nrow(proposals)),
    .package = "EMC2"
  )
  # One covariance block with a wide spread of scales, blocks of three and
  # two parameters beside an uncorrelated one, and group means whose scales
  # differ by nine orders of magnitude.
  cases <- list(list(par_group = c(1, 1, 1, 1), v = 4, A = 0.3, mu_sd = 1),
                list(par_group = c(1, 1, 1, 2, 2, 3), v = 8, A = 1, mu_sd = 1),
                list(par_group = c(1, 1, 1, 1), v = 8, A = 1, mu_sd = c(1e6, 1, 1e-3, 0.3)))
  for (settings in cases) {
    case <- standard_prior_draws(30000L, settings$par_group, n_subjects = 5L,
                                 v = settings$v, A = settings$A, mu_sd = settings$mu_sd)
    n <- nrow(case$draws)
    for (split in list(seq_len(n / 2), seq.int(n / 2 + 1, n))) {
      fit <- EMC2:::bridge_fit(case$draws, split)
      log_ml <- EMC2:::bridge_sampling(case$draws, case$data, case$info, n / 2, split, fit,
                                       1L, 1L, 5000L, 1e-5, 1e-10, 1e-6, FALSE)
      expect_lt(abs(log_ml), 0.5)
    }
  }
})

test_that("factor prior draws under a constant likelihood have log marginal likelihood zero", {
  local_test_rng(12)
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, model, ...) numeric(nrow(proposals)),
    .package = "EMC2"
  )
  n <- 30000L
  p <- 4L
  n_subjects <- 5L
  prior <- list(theta_mu_mean = rep(0, p), theta_mu_var = rep(1, p),
                theta_lambda_var = rep(0.2, p), ap = 20, bp = 2, as = rep(20, p), bs = rep(2, p))
  Lambda_mat <- matrix(c(1, Inf, Inf, Inf), p, 1L)
  free <- Lambda_mat == Inf
  theta_mu <- matrix(rnorm(p * n), p, n)
  epsilon_inv <- matrix(rgamma(p * n, prior$as, prior$bs), p, n)
  psi_inv <- matrix(rgamma(n, prior$ap, prior$bp), 1L, n)
  lambda <- array(0, c(p, 1L, n))
  alpha <- array(0, c(p, n_subjects, n))
  for (i in seq_len(n)) {
    L <- Lambda_mat
    L[free] <- rnorm(sum(free), 0, sqrt(prior$theta_lambda_var[1]))
    lambda[, , i] <- L
    eta <- matrix(rnorm(n_subjects, 0, sqrt(1 / psi_inv[, i])), 1L)
    alpha[, , i] <- theta_mu[, i] + L %*% eta +
      matrix(rnorm(p * n_subjects, 0, sqrt(1 / epsilon_inv[, i])), p)
  }
  samples <- structure(list(
    subjects = seq_len(n_subjects), par_names = paste0("p", seq_len(p)), n_pars = p,
    n_subjects = n_subjects, n_factors = 1L, prior = prior,
    samples = list(theta_mu = theta_mu, lambda_untransf = lambda,
                   epsilon_inv = epsilon_inv, psi_inv = psi_inv)
  ), Lambda_mat = Lambda_mat)
  draws <- cbind(t(matrix(alpha, p * n_subjects, n)),
                 EMC2:::bridge_add_group(matrix(numeric(), n, 0L), samples, seq_len(n), "factor"))
  info <- EMC2:::make_info(samples, "factor")
  for (split in list(seq_len(n / 2), seq.int(n / 2 + 1, n))) {
    fit <- EMC2:::bridge_fit(draws, split)
    log_ml <- EMC2:::bridge_sampling(draws, vector("list", n_subjects), info, n / 2, split, fit,
                                     1L, 1L, 5000L, 1e-5, 1e-10, 1e-6, FALSE)
    expect_lt(abs(log_ml), 0.5)
  }
})

test_that("SEM prior draws under a constant likelihood have log marginal likelihood zero", {
  local_test_rng(13)
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, model, ...) numeric(nrow(proposals)),
    .package = "EMC2"
  )
  n <- 30000L
  p <- 4L
  k <- 2L
  n_subjects <- 5L
  prior <- list(theta_mu_mean = rep(0, p), theta_mu_var = rep(1, p), lambda_var = rep(0.5, k),
                B_var = rep(0.5, k), K_var = numeric(0), G_var = numeric(0),
                a_d = 20, b_d = 2, a_e = rep(20, p), b_e = rep(2, p))
  Lambda_mat <- matrix(c(1, Inf, 0, 0, 0, 0, 1, Inf), p, k)
  free <- Lambda_mat == Inf
  # Both factors share one covariance block, then each has its own precision.
  for (factor_groups in list(c(1, 1), c(1, 2))) {
    theta_mu <- matrix(rnorm(p * n), p, n)
    epsilon_inv <- matrix(rgamma(p * n, prior$a_e, prior$b_e), p, n)
    lambda <- array(0, c(p, k, n))
    delta_inv <- array(0, c(k, k, n))
    alpha <- array(0, c(p, n_subjects, n))
    for (i in seq_len(n)) {
      L <- Lambda_mat
      L[free] <- rnorm(sum(free), 0, sqrt(prior$lambda_var[col(Lambda_mat)[free]]))
      lambda[, , i] <- L
      for (group in unique(factor_groups)) {
        idx <- which(factor_groups == group)
        d <- length(idx)
        delta_inv[idx, idx, i] <- if (d > 1) {
          rWishart(1, prior$a_d, diag(1 / prior$b_d, d))[, , 1]
        } else {
          rgamma(1, prior$a_d, prior$b_d)
        }
      }
      V <- L %*% solve(delta_inv[, , i]) %*% t(L) + diag(1 / epsilon_inv[, i], p)
      alpha[, , i] <- theta_mu[, i] + t(chol(V)) %*% matrix(rnorm(p * n_subjects), p)
    }
    samples <- list(
      subjects = seq_len(n_subjects), par_names = paste0("p", seq_len(p)), n_pars = p,
      n_subjects = n_subjects, n_factors = k, prior = prior,
      sem_settings = list(Lambda_mat = Lambda_mat, B_mat = matrix(0, k, k),
                          K_mat = matrix(0, p, 0), G_mat = matrix(0, k, 0),
                          factor_groups = factor_groups, covariates = NULL),
      samples = list(theta_mu = theta_mu, lambda = lambda, B = array(0, c(k, k, n)),
                     K = array(0, c(p, 0, n)), G = array(0, c(k, 0, n)),
                     epsilon_inv = epsilon_inv, delta_inv = delta_inv)
    )
    draws <- cbind(t(matrix(alpha, p * n_subjects, n)),
                   EMC2:::bridge_add_group(matrix(numeric(), n, 0L), samples, seq_len(n), "SEM"))
    info <- EMC2:::make_info(samples, "SEM")
    for (split in list(seq_len(n / 2), seq.int(n / 2 + 1, n))) {
      fit <- EMC2:::bridge_fit(draws, split)
      log_ml <- EMC2:::bridge_sampling(draws, vector("list", n_subjects), info, n / 2, split,
                                       fit, 1L, 1L, 5000L, 1e-5, 1e-10, 1e-6, FALSE)
      expect_lt(abs(log_ml), 0.5)
    }
  }
})

test_that("bridge updates retain mass when both sums underflow on the natural scale", {
  expect_equal(EMC2:::bridge_update(c(1000, 1001), c(-1002, -1003), log(0.5), log(0.5)),
               -2, tolerance = 1e-12)
  expect_equal(EMC2:::bridge_log_add(c(-Inf, -1000, NA_real_), c(-Inf, -1000, -1000)),
               c(-Inf, -1000 + log(2), NA_real_))
})

test_that("Bayes factor comparison reproduces the seeded LNR marginal deviance", {
  local_test_rng(12, kind = "Mersenne-Twister")
  result <- compare(list(m = samples_LNR), BayesFactor = TRUE, WAIC = FALSE,
                    both_splits = FALSE, cores_for_props = 1L, print_summary = FALSE)
  expect_equal(result$MD, -571.823715913, tolerance = 1e-8)
  expect_equal(result$wMD, 1)
})

test_that("proposal blocks preserve draws and densities at singleton boundaries", {
  local_test_rng(5)
  samples <- EMC2:::merge_chains(samples_LNR)
  info <- EMC2:::make_info(samples, samples$type)
  draws <- EMC2:::bridge_draws(samples_LNR, "sample", info)
  fit <- EMC2:::bridge_fit(draws, seq_len(70))
  kinds <- c("posterior", "reflected", "negative", "positive")
  evaluate <- function(idx, gen, block_size, cores = 1L) {
    EMC2:::bridge_eval_kinds(kinds, draws, idx, gen, samples$data, fit$mean, fit$chol,
                             info, cores, FALSE, block_size)
  }
  # Only the last block contains a single draw.
  idx <- 71:77
  gen <- matrix(rnorm(length(idx) * ncol(draws)), length(idx))
  whole <- evaluate(idx, gen, length(idx))
  expect_named(whole, kinds)
  expect_equal(evaluate(idx, gen, 3L), whole, tolerance = 1e-10)

  # Row chunks shared between processes return the rows in order.
  skip_on_os("windows")
  idx <- seq_len(nrow(draws))
  gen <- matrix(rnorm(length(idx) * ncol(draws)), length(idx))
  expect_equal(evaluate(idx, gen, length(idx), cores = 2L), evaluate(idx, gen, length(idx)),
               tolerance = 1e-10)
})

# The standard group term from covariance matrices: multivariate normal
# densities of the subjects and of the group mean coefficients, inverse-Wishart
# and inverse-gamma priors, and the Jacobian of the sampled coordinates.
standard_group_reference <- function(proposals_group, proposals_list, info) {
  prior <- info$prior
  p <- info$n_pars
  designs <- EMC2:::add_group_design(info$par_names, info$group_designs, info$n_subjects)
  n_mu <- length(prior$theta_mu_mean)
  log_inv_gamma <- function(x, shape, rate) dgamma(1 / x, shape, rate, log = TRUE) - 2 * log(x)
  vapply(seq_len(nrow(proposals_group)), function(i) {
    row <- proposals_group[i, ]
    theta_mu <- row[seq_len(n_mu)]
    a_half <- exp(row[n_mu + seq_len(p)])
    at <- n_mu + p
    Sigma <- matrix(0, p, p)
    out <- sum(log_inv_gamma(a_half, 0.5, 1 / prior$A^2)) + sum(log(a_half)) +
      mvtnorm::dmvnorm(theta_mu, prior$theta_mu_mean, prior$theta_mu_var, log = TRUE)
    for (k in which(!info$is_blocked)) {
      at <- at + 1L
      Sigma[k, k] <- exp(row[at])
      out <- out + log_inv_gamma(Sigma[k, k], prior$v / 2, prior$v / a_half[k]) + row[at]
    }
    for (group in unique(info$par_group[info$is_blocked])) {
      idx <- which(info$par_group == group)
      d <- length(idx)
      L <- matrix(0, d, d)
      L[lower.tri(L, diag = TRUE)] <- row[at + seq_len(d * (d + 1) / 2)]
      at <- at + d * (d + 1) / 2
      log_diag <- diag(L)
      diag(L) <- exp(log_diag)
      Sigma[idx, idx] <- tcrossprod(L)
      out <- out + d * log(2) + sum((d - seq_len(d) + 2) * log_diag) +
        iwish_reference(Sigma[idx, idx], prior$v + d - 1, diag(2 * prior$v / a_half[idx], d))
    }
    mu <- EMC2:::calculate_subject_means(designs, theta_mu)
    alpha <- matrix(vapply(proposals_list, function(x) x[i, ], numeric(p)), p)
    out + sum(mvtnorm::dmvnorm(t(alpha - mu), sigma = Sigma, log = TRUE))
  }, numeric(1))
}

test_that("native standard group term matches densities built from covariance matrices", {
  local_test_rng(21)
  n_subjects <- 7L
  cases <- list(
    list(par_group = c(1, 1, 1, 2, 2, 3), design = TRUE),
    list(par_group = c(1, 1, 1, 1), design = FALSE),
    list(par_group = c(1, 2, 3), design = FALSE),
    list(par_group = 1, design = FALSE)
  )
  for (case in cases) {
    p <- length(case$par_group)
    par_names <- paste0("p", seq_len(p))
    group_designs <- if (case$design)
      list(p2 = cbind(1, rnorm(n_subjects)), p5 = cbind(1, rnorm(n_subjects), rnorm(n_subjects)))
    n_mu <- p + sum(vapply(group_designs, ncol, integer(1)) - 1L)
    blocked <- duplicated(case$par_group) | duplicated(case$par_group, fromLast = TRUE)
    info <- list(
      par_names = par_names, n_pars = p, n_subjects = n_subjects, par_group = case$par_group,
      is_blocked = blocked, group_designs = group_designs,
      prior = list(theta_mu_mean = rnorm(n_mu),
                   theta_mu_var = crossprod(matrix(rnorm(n_mu * (n_mu + 2)), n_mu + 2)) / n_mu,
                   v = 2.5, A = runif(p, 0.2, 2))
    )
    n_group <- n_mu + p + sum(!blocked) +
      sum(vapply(unique(case$par_group[blocked]), function(group) {
        d <- sum(case$par_group == group)
        d * (d + 1) / 2
      }, numeric(1)))
    n <- 25L
    proposals_group <- matrix(rnorm(n * n_group, sd = 0.7), n)
    proposals_list <- lapply(seq_len(n_subjects), function(s) matrix(rnorm(n * p), n))
    expect_equal(
      EMC2:::bridge_group_and_prior_and_jac(proposals_group, proposals_list, info, "standard"),
      standard_group_reference(proposals_group, proposals_list, info),
      tolerance = 1e-10
    )
  }

  # A proposal off the support has zero density.
  proposals_group[1L, 1L] <- NaN
  expect_identical(
    EMC2:::bridge_group_and_prior_and_jac(proposals_group, proposals_list, info, "standard")[1L],
    -Inf
  )
})

test_that("stored subject log-likelihoods stand in for the posterior likelihood pass only when they match", {
  emc <- samples_LNR
  samples <- emc[[1L]]
  info <- EMC2:::make_info(samples, samples$type)
  draws <- EMC2:::bridge_draws(emc, "sample", info)
  recomputed <- EMC2:::bridge_ll(EMC2:::bridge_subject_proposals(draws, info), samples$data, info)
  expect_equal(EMC2:::bridge_stored_ll(emc, "sample", draws, samples$data, info), recomputed,
               tolerance = 1e-10)

  shifted <- emc
  shifted[[2L]]$samples$subj_ll <- shifted[[2L]]$samples$subj_ll + 1e-3
  expect_null(EMC2:::bridge_stored_ll(shifted, "sample", draws, samples$data, info))
  missing <- emc
  missing[[3L]]$samples$subj_ll <- NULL
  expect_null(EMC2:::bridge_stored_ll(missing, "sample", draws, samples$data, info))

  # The estimate does not depend on which route supplied the likelihood.
  local_test_rng(7)
  with_stored <- suppressWarnings(run_bridge_sampling(emc, cores_for_props = 1L))
  set.seed(7)
  recomputing <- suppressWarnings(run_bridge_sampling(shifted, cores_for_props = 1L))
  expect_equal(with_stored, recomputing, tolerance = 1e-10)
})

test_that("repetitions reuse the fits and redraw only the proposals", {
  local_test_rng(3)
  together <- suppressWarnings(run_bridge_sampling(samples_LNR, repetitions = 3L,
                                                   cores_for_props = 1L))
  set.seed(3)
  apart <- suppressWarnings(vapply(1:3, function(i)
    run_bridge_sampling(samples_LNR, cores_for_props = 1L), numeric(1)))
  expect_equal(together, apart, tolerance = 1e-12)
  expect_gt(sd(together), 0)
})

test_that("the filter is applied once", {
  for (filter in list(20L, 11:40)) {
    local_test_rng(4)
    filtered <- suppressWarnings(run_bridge_sampling(samples_LNR, filter = filter,
                                                     cores_for_props = 1L))
    set.seed(4)
    expected <- suppressWarnings(run_bridge_sampling(subset(samples_LNR, filter = filter),
                                                     cores_for_props = 1L))
    expect_equal(filtered, expected, tolerance = 1e-12)
  }
})

test_that("the proposal moments are those of the fitting draws", {
  samples <- samples_LNR[[1L]]
  info <- EMC2:::make_info(samples, samples$type)
  draws <- EMC2:::bridge_draws(samples_LNR, "sample", info)
  idx <- seq.int(2L, nrow(draws), by = 2L)
  fit <- EMC2:::bridge_fit(draws, idx)
  expect_equal(fit$mean, colMeans(draws[idx, ]), tolerance = 1e-12)
  expect_equal(crossprod(fit$chol), var(draws[idx, ]), tolerance = 1e-10)

  # Columns on very different scales keep their own variances.
  scales <- 10^seq(-4, 4, length.out = ncol(draws))
  wide <- sweep(draws, 2L, scales, "*")
  fit <- EMC2:::bridge_fit(wide, idx)
  expect_equal(diag(crossprod(fit$chol)) / scales^2, diag(var(draws[idx, ])),
               tolerance = 1e-8)

  # Fewer fitting draws than columns still give a usable factor.
  fit <- EMC2:::bridge_fit(draws, seq_len(5L))
  expect_true(all(is.finite(fit$chol)) && all(diag(fit$chol) > 0))
})

test_that("two bridge estimates swap complementary halves and use all fitting draws", {
  local_test_rng()
  for (odd in c(FALSE, TRUE)) {
    emc <- samples_LNR
    if (odd) emc <- subset(emc, length.out = 49L)
    samples <- emc[[1L]]
    info <- EMC2:::make_info(samples, samples$type)
    draws <- EMC2:::bridge_draws(emc, "sample", info)
    n_eff <- round(ess_summary(emc, selection = "alpha", stat = "median", stat_only = TRUE) / 2)
    halves <- split(seq_len(nrow(draws)), seq_len(nrow(draws)) > round(nrow(draws) / 2))

    set.seed(12)
    expected <- vapply(halves, function(idx) {
      fitted <- draws[idx, , drop = FALSE]
      fit <- list(mean = colMeans(fitted), chol = chol(var(fitted)))
      EMC2:::bridge_sampling(draws, samples$data, info, n_eff, idx, fit,
                             1L, 1L, 5000L, 1e-5, 1e-10, 1e-6, FALSE)
    }, numeric(1))

    set.seed(12)
    actual <- suppressWarnings(run_bridge_sampling(emc, both_splits = TRUE, cores_for_props = 1L))
    expect_equal(actual, mean(expected), tolerance = 1e-10)
  }
})
