local_rng_guard()

cache_chain <- function() {
  list(data = list(data.frame(rt = c(0.2, 0.4), subjects = factor(c("a", "a")))),
       model = function() 1, par_names = "x", n_pars = 1L,
       type = "single", n_subjects = 1L)
}

test_that("likelihood configuration fingerprints survive serialization", {
  chain <- cache_chain()
  key <- EMC2:::.emc_likelihood_config_fingerprint(chain)
  expect_true(is.character(key) && nchar(key) == 32L)
  restored <- unserialize(serialize(chain, NULL, version = 3L))
  expect_identical(EMC2:::.emc_likelihood_config_fingerprint(restored), key)
  attr(restored$data[[1L]], "pooled_values") <- c(group = 0.1)
  expect_identical(EMC2:::.emc_likelihood_config_fingerprint(restored), key)
  restored$data[[1L]]$rt[1L] <- 0.3
  expect_false(identical(EMC2:::.emc_likelihood_config_fingerprint(restored), key))
})

test_that("configuration changes invalidate fit and diagnostic caches", {
  emc <- EMC2:::.emc_prepare_subject_fit_cache(list(cache_chain()),
                                                validate = TRUE)
  emc[[1L]]$subject_fit_cache[[1L]] <- list(sentinel = TRUE)
  emc[[1L]]$laplace_check_cache[[1L]] <- list(list(sentinel = TRUE))
  emc[[1L]]$lik_laplace <- list(list(prec = matrix(1), lin = 1))
  emc[[1L]]$burn_mode_mu <- list(1)
  emc[[1L]]$burn_mode_var <- list(matrix(1))
  emc[[1L]]$data[[1L]]$rt[1L] <- 0.3
  refreshed <- EMC2:::.emc_prepare_subject_fit_cache(emc, validate = TRUE)
  expect_null(refreshed[[1L]]$subject_fit_cache[[1L]])
  expect_null(refreshed[[1L]]$laplace_check_cache[[1L]])
  expect_null(refreshed[[1L]]$lik_laplace[[1L]])
  expect_null(refreshed[[1L]]$burn_mode_mu)
  expect_null(refreshed[[1L]]$burn_mode_var)
})

test_that("only benchmarked REXG uses fit reuse by default", {
  rexg <- cache_chain()
  rexg$model <- function() list(c_name = "REXG")
  lnr <- cache_chain()
  lnr$model <- function() list(c_name = "LNR")
  expect_true(is.infinite(EMC2:::.emc_subject_fit_cache_max_reuses(rexg)))
  expect_identical(EMC2:::.emc_subject_fit_cache_max_reuses(lnr), 0L)
  withr::local_options(emc2.subject_fit_cache_max_reuses = Inf)
  expect_true(is.infinite(EMC2:::.emc_subject_fit_cache_max_reuses(lnr)))
})

test_that("fit cache hits require exact prior and proposal scale", {
  emc <- EMC2:::.emc_prepare_subject_fit_cache(list(cache_chain()),
                                                validate = TRUE)
  chain <- emc[[1L]]
  prior <- list(mu = c(0), var = matrix(1))
  key <- EMC2:::.emc_subject_fit_key(chain, 1L, prior, 4)
  fit <- list(alpha = c(0.2), covariance = matrix(0.5),
              proposal_mode = c(0.2), start_covariance = matrix(1),
              lik = list(prec = matrix(1), lin = c(0)))
  chain <- EMC2:::.emc_subject_fit_store(chain, 1L, key, fit)
  hit <- EMC2:::.emc_subject_fit_lookup(chain, 1L, key)
  expect_true(hit$hit)
  expect_identical(hit$result$alpha, fit$alpha)
  expect_true(is.na(hit$result$gain))
  expect_identical(hit$warm_start, fit$alpha)
  changed_prior <- prior
  changed_prior$mu <- c(1e-12)
  miss <- EMC2:::.emc_subject_fit_lookup(
    chain, 1L, EMC2:::.emc_subject_fit_key(chain, 1L, changed_prior, 4))
  expect_false(miss$hit)
  expect_identical(miss$warm_start, fit$alpha)
  changed_scale <- EMC2:::.emc_subject_fit_lookup(
    chain, 1L, EMC2:::.emc_subject_fit_key(chain, 1L, prior, 2))
  expect_false(changed_scale$hit)
  fresh <- EMC2:::.emc_prepare_subject_fit_cache(list(cache_chain()),
                                                 validate = TRUE)[[1L]]
  failed <- EMC2:::.emc_subject_fit_store(
    fresh, 1L, key, list(alpha = c(0), covariance = NULL))
  expect_null(failed$subject_fit_cache[[1L]])
})

test_that("a configured cache age forces a fresh current-state fit", {
  emc <- EMC2:::.emc_prepare_subject_fit_cache(list(cache_chain()),
                                                validate = TRUE)
  chain <- emc[[1L]]
  prior <- list(mu = c(0), var = matrix(1))
  key <- EMC2:::.emc_subject_fit_key(chain, 1L, prior, 4)
  fit <- list(alpha = c(0.2), covariance = matrix(0.5), proposal_mode = c(0.2))
  chain <- EMC2:::.emc_subject_fit_store(chain, 1L, key, fit,
                                          reuse_count = 2L)
  expect_true(EMC2:::.emc_subject_fit_refresh_due(chain, 1L, key, 2L))
  expect_false(EMC2:::.emc_subject_fit_refresh_due(chain, 1L, key, 3L))
  forced <- EMC2:::.emc_subject_fit_lookup(
    chain, 1L, key, force_refresh = TRUE)
  expect_false(forced$hit)
  expect_null(forced$warm_start)
  expect_identical(forced$reuse_count, 2L)
})

test_that("diagnostic cache matches the assessed Gaussian and target prior exactly", {
  emc <- EMC2:::.emc_prepare_subject_fit_cache(list(cache_chain()),
                                                validate = TRUE)
  chain <- emc[[1L]]
  prior <- list(mu = c(0), var = matrix(1))
  gaussian_mu <- c(0.2)
  gaussian_var <- matrix(0.4)
  key <- EMC2:::.emc_laplace_check_key(
    chain, 1L, gaussian_mu, gaussian_var, prior)
  check <- list(pareto_k = 0.1, outside = 0, status = "ok")
  chain <- EMC2:::.emc_laplace_check_store(chain, 1L, key, check)
  expect_identical(EMC2:::.emc_laplace_check_lookup(chain, 1L, key), check)
  expect_null(EMC2:::.emc_laplace_check_lookup(
    chain, 1L, EMC2:::.emc_laplace_check_key(
      chain, 1L, gaussian_mu + 1e-12, gaussian_var, prior)))
  changed_prior <- prior
  changed_prior$var <- matrix(1 + 1e-12)
  expect_null(EMC2:::.emc_laplace_check_lookup(
    chain, 1L, EMC2:::.emc_laplace_check_key(
      chain, 1L, gaussian_mu, gaussian_var, changed_prior)))
  pooled_chain <- chain
  pooled_chain$pooled <- list(value = c(g = 0.2))
  pooled_key <- EMC2:::.emc_laplace_check_key(
    pooled_chain, 1L, gaussian_mu, gaussian_var, prior)
  pooled_chain <- EMC2:::.emc_laplace_check_store(
    pooled_chain, 1L, pooled_key, check)
  pooled_chain$pooled$value <- c(g = 0.3)
  expect_null(EMC2:::.emc_laplace_check_lookup(
    pooled_chain, 1L, EMC2:::.emc_laplace_check_key(
      pooled_chain, 1L, gaussian_mu, gaussian_var, prior)))
  chain <- EMC2:::.emc_laplace_check_store(
    chain, 1L, key, list(pareto_k = NA_real_, outside = NA_real_, status = "error"))
  expect_length(chain$laplace_check_cache[[1L]], 1L)
})

test_that("the subject cache is durable in saved fit state", {
  emc <- EMC2:::.emc_prepare_subject_fit_cache(list(cache_chain()),
                                                validate = TRUE)
  chain <- emc[[1L]]
  prior <- list(mu = c(0), var = matrix(1))
  key <- EMC2:::.emc_subject_fit_key(chain, 1L, prior, 4)
  fit <- list(alpha = c(0.2), covariance = matrix(0.5), proposal_mode = c(0.2))
  chain <- EMC2:::.emc_subject_fit_store(chain, 1L, key, fit)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(chain), path)
  restored <- readRDS(path)
  expect_identical(restored[[1L]]$likelihood_cache_version,
                   chain$likelihood_cache_version)
  expect_true(EMC2:::.emc_subject_fit_lookup(restored[[1L]], 1L, key)$hit)
})

test_that("duplicate stripping leaves chain-local fit caches intact", {
  emc <- EMC2:::.emc_prepare_subject_fit_cache(
    list(cache_chain(), cache_chain()), validate = TRUE)
  for (j in seq_along(emc)) {
    prior <- list(mu = c(0), var = matrix(1))
    key <- EMC2:::.emc_subject_fit_key(emc[[j]], 1L, prior, 4)
    fit <- list(alpha = c(j), covariance = matrix(0.5),
                proposal_mode = c(j))
    emc[[j]] <- EMC2:::.emc_subject_fit_store(emc[[j]], 1L, key, fit)
  }
  saved <- EMC2:::strip_duplicates(emc)
  expect_true(is.list(saved[[2L]]$subject_fit_cache[[1L]]))
  restored <- EMC2:::restore_duplicates(saved)
  expect_identical(restored[[1L]]$subject_fit_cache[[1L]]$result$alpha, c(1L))
  expect_identical(restored[[2L]]$subject_fit_cache[[1L]]$result$alpha, c(2L))
})

test_that("an exact single-subject refresh reuses its fit and diagnostic", {
  withr::local_options(emc2.subject_fit_cache_max_reuses = Inf)
  withr::local_seed(812)
  calls <- new.env(parent = emptyenv())
  calls$mode <- calls$check <- 0L
  local_mocked_bindings(
    .emc_burn_subject_mode = function(sampler, s, ...) {
      calls$mode <- calls$mode + 1L
      list(alpha = c(0.2, -0.1), ll = -0.5,
           covariance = diag(c(0.4, 0.7)),
           proposal_mode = c(0.2, -0.1),
           start_covariance = diag(c(0.5, 0.8)),
           lik = list(prec = diag(2), lin = c(0, 0)))
    },
    .emc_laplace_is_start = function(mode, ...) {
      calls$check <- calls$check + 1L
      rnorm(5L)
      list(pareto_k = 0.1, outside = 0, status = "ok")
    },
    .package = "EMC2")
  chain <- list(type = "single", n_pars = 2L, n_subjects = 1L,
    par_names = c("a", "b"), nuisance = FALSE,
    samples = list(idx = 1L, alpha = array(0, c(2L, 1L, 1L))),
    data = list(c(0, 0)), model = NULL,
    prior = list(theta_mu_mean = c(0, 0), theta_mu_var = diag(2)))
  emc <- EMC2:::.emc_prepare_chain_rng(list(chain))
  before <- .Random.seed
  first <- EMC2:::.emc_laplace_refresh(emc, 1L)
  expect_identical(.Random.seed, before)
  expect_identical(c(calls$mode, calls$check), c(1L, 1L))
  second <- EMC2:::.emc_laplace_refresh(first, 1L)
  expect_identical(c(calls$mode, calls$check), c(1L, 1L))
  expect_identical(second[[1L]]$rng$subjects, first[[1L]]$rng$subjects)
  expect_identical(second[[1L]]$subject_fit_cache[[1L]]$result$alpha,
                   first[[1L]]$subject_fit_cache[[1L]]$result$alpha)
  expect_identical(second[[1L]]$subject_fit_cache[[1L]]$result$covariance,
                   first[[1L]]$subject_fit_cache[[1L]]$result$covariance)
})

test_that("an exact likelihood diagnostic is reused but changed inputs are checked", {
  withr::local_seed(813)
  calls <- 0L
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, ...) {
      -0.5 * rowSums(sweep(proposals, 2L, dadm)^2)
    },
    .emc_laplace_is_start = function(mode, ...) {
      calls <<- calls + 1L
      rnorm(5L)
      list(pareto_k = 0.1, outside = 0, status = "ok")
    },
    .package = "EMC2")
  chain <- list(type = "single", n_pars = 2L, n_subjects = 1L,
    par_names = c("a", "b"), samples = list(idx = 1L),
    data = list(c(1, -1)), model = NULL,
    prior = list(theta_mu_mean = c(0, 0), theta_mu_var = diag(2)),
    lik_laplace = list(list(prec = diag(2), lin = c(1, -1))))
  emc <- EMC2:::.emc_prepare_chain_rng(list(chain))
  caller_seed <- .Random.seed
  first <- EMC2:::.emc_check_likelihood_fits(emc, 1L)
  subject_seed <- first[[1L]]$rng$subjects[[1L]]
  second <- EMC2:::.emc_check_likelihood_fits(first, 1L)
  expect_identical(calls, 1L)
  expect_identical(second[[1L]]$rng$subjects[[1L]], subject_seed)
  expect_identical(.Random.seed, caller_seed)
  second[[1L]]$lik_laplace[[1L]]$lin <- c(2, -2)
  third <- EMC2:::.emc_check_likelihood_fits(second, 1L)
  expect_identical(calls, 2L)
  expect_length(third[[1L]]$laplace_check_cache[[1L]], 2L)
})

test_that("a hierarchical prior change refits from the cached conditional mode", {
  calls <- list()
  local_mocked_bindings(
    .emc_burn_subject_mode = function(sampler, s, x0 = NULL, prior = NULL, ...) {
      calls[[length(calls) + 1L]] <<- if (is.null(x0)) "cold" else x0
      list(alpha = prior$mu, ll = 0, covariance = diag(2),
           proposal_mode = prior$mu, start_covariance = diag(2),
           lik = list(prec = diag(2), lin = prior$mu))
    },
    .emc_laplace_is_start = function(...) {
      list(pareto_k = 0.1, outside = 0, status = "ok")
    },
    .package = "EMC2")
  chain <- list(type = "factor", n_pars = 2L, n_subjects = 1L,
    par_names = c("a", "b"), nuisance = FALSE,
    samples = list(idx = 1L, alpha = array(0, c(2L, 1L, 1L)),
                   theta_mu = matrix(c(0, 0), ncol = 1L),
                   theta_var = array(diag(2), c(2L, 2L, 1L))),
    data = list(c(0, 0)), model = NULL)
  emc <- EMC2:::.emc_prepare_chain_rng(list(chain))
  first <- EMC2:::.emc_laplace_refresh(emc, 1L)
  first_mode <- first[[1L]]$subject_fit_cache[[1L]]$result$alpha
  first[[1L]]$samples$theta_mu[, 1L] <- c(0.5, -0.25)
  second <- EMC2:::.emc_laplace_refresh(first, 1L)
  expect_length(calls, 2L)
  expect_identical(calls[[1L]], "cold")
  expect_identical(calls[[2L]], first_mode)
  expect_identical(second[[1L]]$subject_fit_cache[[1L]]$key$prior_mu,
                   c(0.5, -0.25))
})

test_that("burn mode starts reuse cached fits but still draw a fresh start", {
  withr::local_options(emc2.subject_fit_cache_max_reuses = Inf)
  withr::local_seed(814)
  calls <- new.env(parent = emptyenv())
  calls$mode <- calls$check <- 0L
  local_mocked_bindings(
    .emc_burn_subject_mode = function(sampler, s, ...) {
      calls$mode <- calls$mode + 1L
      list(alpha = c(0.2, -0.1), ll = 0,
           covariance = diag(2), start_covariance = diag(2),
           proposal_mode = c(0.2, -0.1),
           lik = list(prec = diag(2), lin = c(0, 0)))
    },
    .emc_laplace_is_start = function(...) {
      calls$check <- calls$check + 1L
      rnorm(5L)
      list(pareto_k = 0.1, outside = 0, status = "ok")
    },
    .emc_burn_ll_batch = function(x, ...) rep(0, NROW(x)),
    .package = "EMC2")
  chain <- list(type = "single", n_pars = 2L, n_subjects = 1L,
    par_names = c("a", "b"), nuisance = FALSE,
    samples = list(idx = 1L, alpha = array(0, c(2L, 1L, 1L)),
                   subj_ll = matrix(0, 1L, 1L)),
    data = list(c(0, 0)), model = NULL,
    prior = list(theta_mu_mean = c(0, 0), theta_mu_var = diag(2)))
  emc <- EMC2:::.emc_prepare_chain_rng(list(chain))
  first <- EMC2:::.emc_burn_mode_start(emc, 1L)
  second <- EMC2:::.emc_burn_mode_start(first, 1L)
  expect_identical(c(calls$mode, calls$check), c(1L, 1L))
  expect_false(identical(second[[1L]]$rng$subjects, first[[1L]]$rng$subjects))
  expect_false("start" %in% names(second[[1L]]$subject_fit_cache[[1L]]$result))
})
