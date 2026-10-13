test_that("the likelihood approximation recombines exactly with a new group prior", {
  set.seed(3)
  p <- 4
  A <- matrix(rnorm(p * p), p); S <- crossprod(A) + diag(p)   # likelihood covariance
  y <- rnorm(p)                                              # likelihood mode
  B <- matrix(rnorm(p * p), p); V1 <- crossprod(B) + diag(p); mu1 <- rnorm(p)
  # Conditional mode under N(mu1, V1), and the likelihood Hessian there.
  P1 <- solve(V1); Si <- solve(S)
  x <- drop(solve(Si + P1, Si %*% y + P1 %*% mu1))
  lik <- EMC2:::.emc_lik_canonical(Si, x, P1, mu1)
  expect_equal(lik$lin, drop(Si %*% y))
  # A different group state: the result is that state's exact conditional.
  V2 <- diag(p) * 2; mu2 <- rnorm(p)
  eff <- structure(rep(0, p), lik = lik)
  out <- EMC2:::.emc_laplace_conditional(eff, list(mu = mu2, var = V2))
  P2 <- solve(V2)
  expect_equal(out$var, solve(Si + P2))
  expect_equal(out$mu, drop(solve(Si + P2, Si %*% y + P2 %*% mu2)))
})

test_that("an indefinite likelihood Hessian is projected, and a missing one falls back", {
  H <- diag(c(2, -1, 0.5))
  lik <- EMC2:::.emc_lik_canonical(H, c(0, 0, 0), diag(3), c(0, 0, 0))
  expect_equal(lik$prec, diag(c(2, 0, 0.5)))
  expect_null(EMC2:::.emc_laplace_conditional(c(0, 0, 0), list(mu = rep(0, 3), var = diag(3))))
  eff <- list(c(1, 2, 3))
  expect_null(attr(EMC2:::.emc_attach_lik(eff, list(lik))[[1]], "lik"))
})

test_that("the Laplace-conditional proposal leaves a conjugate posterior invariant, even when wrong", {
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
  # Exact likelihood fit, and one that is shifted and too confident: the
  # proposal changes, the stationary distribution must not.
  fits <- list(exact = list(prec = Li, lin = drop(Li %*% ybar)),
               wrong = list(prec = 3 * Li, lin = drop(3 * Li %*% (ybar + c(.4, -.3)))))
  for (lik in fits) {
    set.seed(3)
    pm <- list(list(epsilon = c(1, 1, 1), mix = c(.1, .2, .1, .6), n_particles = 5,
                    iter = 0, proposal_counts = rep(0, 4), acc_counts = rep(0, 4),
                    gd_good = FALSE))
    tune <- list(components = c(1, 1), shared_ll_idx = c(1, 1), n0 = Inf)
    eff <- structure(c(.8, -.2), lik = lik)
    a <- c(a = 0, b = 0); ll <- -0.5 * sum((a - ybar)^2 * diag(Li))
    out <- matrix(NA_real_, 20000, 2)
    for (t in seq_len(nrow(out))) {
      r <- EMC2:::new_particle(
        1, NULL, pm, eff_mu = eff, eff_var = diag(2) * .3,
        chains_mu = c(.8, -.2), chains_var = diag(2) * .2, prev_ll = ll,
        parameters = NULL, stage = "sample", type = "standard", tune = tune,
        current_alpha = a, population_mu = mu, population_var = Sig)
      a <- r$proposal; ll <- r$ll; out[t, ] <- a
    }
    for (k in 1:2) {
      x <- out[-(1:1000), k]
      expect_lt(abs(mean(x) - post_mu[k]) / bm_se(x), 4.5)
      expect_lt(abs(var(x) - post_var[k, k]) / bm_se((x - mean(x))^2), 4.5)
    }
  }
})

test_that("the subject optimiser returns the exact Gaussian likelihood in canonical form", {
  prior <- list(mu = c(-1, 0.5), var = matrix(c(2, 0.3, 0.3, 1), 2))
  observation <- c(1.2, -0.7)
  likelihood_var <- matrix(c(0.6, 0.1, 0.1, 0.4), 2)
  ll_batch <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
    d <- sweep(x, 2L, observation)
    -0.5 * rowSums((d %*% solve(likelihood_var)) * d)
  }
  result <- EMC2:::.emc_optimize_subject(c(4, 3), prior, ll_batch)
  expect_equal(result$lik$prec, solve(likelihood_var), tolerance = 1e-4)
  expect_equal(result$lik$lin, drop(solve(likelihood_var, observation)), tolerance = 1e-3)
})

test_that("production rebuild timing is read from the transition counter", {
  emc <- samples_LNR
  for (j in seq_along(emc)) { emc[[j]]$sample_proposals_next <- 200L; emc[[j]]$sample_transitions <- 150L }
  expect_false(EMC2:::.emc_sample_rebuild_due(emc))
  for (j in seq_along(emc)) emc[[j]]$sample_transitions <- 200L
  expect_true(EMC2:::.emc_sample_rebuild_due(emc))
  for (j in seq_along(emc)) emc[[j]]$sample_proposals_next <- NULL
  expect_false(EMC2:::.emc_sample_rebuild_due(emc))
})

test_that("only checked likelihood fits replace the fitted conditional proposal", {
  fit <- list(prec = diag(2), lin = c(0, 0))
  eff <- list(structure(c(1, 2), B = diag(2), given = c(0, 0), lik = fit))
  for (check in list(NULL, list(pareto_k = Inf, outside = 0),
                    list(pareto_k = NA_real_, outside = 0),
                    list(pareto_k = 0.1, outside = 0.9))) {
    bad <- fit; bad$check <- check
    out <- EMC2:::.emc_attach_lik(eff, list(bad))
    expect_null(attr(out[[1]], "lik"))
    expect_equal(EMC2:::.emc_eff_mean(out[[1]], c(3, 4)), c(4, 6))
  }
  good <- fit; good$check <- list(pareto_k = 0, outside = 0, status = "ok")
  expect_identical(attr(EMC2:::.emc_attach_lik(eff, list(good))[[1]], "lik"), good)
  expect_null(attr(EMC2:::.emc_attach_lik(eff, NULL)[[1]], "lik"))
})

test_that("likelihood checks score the actual fitted location and preserve caller RNG", {
  withr::local_seed(734)
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...) {
    -0.5 * rowSums(sweep(proposals, 2, dadm)^2)
  }, .package = "EMC2")
  fit <- list(prec = diag(2), lin = c(1, -1))
  ch <- list(type = "single", n_pars = 2L, n_subjects = 1L, par_names = c("a", "b"),
             samples = list(idx = 1L), data = list(c(1, -1)), model = NULL,
             prior = list(theta_mu_mean = c(0, 0), theta_mu_var = diag(2)),
             lik_laplace = list(fit))
  ch <- EMC2:::.emc_prepare_sampler_rng(ch)
  set.seed(sum(utf8ToInt("likelihood-check-reproducibility")))
  before <- .Random.seed
  good <- EMC2:::.emc_check_likelihood_fits(list(ch), 1)[[1]]
  expect_identical(.Random.seed, before)
  expect_equal(good$lik_laplace[[1]]$check$pareto_k, 0)
  expect_equal(good$lik_laplace[[1]]$check$outside, 0)
  ch$lik_laplace[[1]]$lin <- c(15, -15)
  bad <- EMC2:::.emc_check_likelihood_fits(list(ch), 1)[[1]]
  expect_identical(.Random.seed, before)
  expect_true(EMC2:::.emc_laplace_poor(bad$lik_laplace[[1]]$check))
  expect_identical(ch$lik_laplace[[1]]$lin, bad$lik_laplace[[1]]$lin)
})

test_that("diagnostic probes checkpoint consumed subject variates even on failure", {
  withr::local_seed(735)
  local_mocked_bindings(.emc_laplace_is_start = function(mode, ...) {
    rnorm(7L)
    if (mode[1L] > 5) stop("diagnostic failure after drawing")
    list(pareto_k = 0, outside = 0, status = "ok")
  }, .emc_burn_subject_mode = function(sampler, s, ...) {
    list(alpha = c(sampler$data[[s]][1], 0), covariance = diag(2),
         proposal_mode = c(0, 0))
  }, .package = "EMC2")
  chain <- function(lin) list(type = "single", n_pars = 2L, n_subjects = 2L,
    par_names = c("a", "b"), nuisance = FALSE,
    samples = list(idx = 1L, alpha = array(0, c(2, 2, 1))),
    data = list(c(lin[1], 0), c(10, 0)), model = NULL,
    prior = list(theta_mu_mean = c(0, 0), theta_mu_var = diag(2)),
    lik_laplace = list(list(prec = diag(2), lin = lin), NULL))
  emc <- EMC2:::.emc_prepare_chain_rng(list(chain(c(1, 0)), chain(c(20, 0))))
  caller_seed <- .Random.seed
  advanced <- function(seed) EMC2:::.emc_with_preserved_rng({
    assign(".Random.seed", seed, envir = globalenv())
    rnorm(7L)
    list(seed = .Random.seed, next_draw = runif(3L))
  })
  next_draw <- function(seed) EMC2:::.emc_with_preserved_rng({
    assign(".Random.seed", seed, envir = globalenv())
    runif(3L)
  })
  for (kind in c("likelihood", "refresh")) {
    run <- if (kind == "likelihood") EMC2:::.emc_check_likelihood_fits else
      EMC2:::.emc_laplace_refresh
    serial <- run(emc, 1L)
    if (.Platform$OS.type != "windows") expect_identical(run(emc, 2L), serial)
    expect_identical(.Random.seed, caller_seed)
    for (j in seq_along(emc)) {
      expect_identical(serial[[j]]$rng$gibbs, emc[[j]]$rng$gibbs)
      for (s in seq_len(emc[[j]]$n_subjects)) {
        if (kind == "likelihood" && s == 2L) {
          expect_identical(serial[[j]]$rng$subjects[[s]], emc[[j]]$rng$subjects[[s]])
        } else {
          expected <- advanced(emc[[j]]$rng$subjects[[s]])
          expect_identical(serial[[j]]$rng$subjects[[s]], expected$seed)
          expect_identical(next_draw(serial[[j]]$rng$subjects[[s]]), expected$next_draw)
        }
      }
    }
    if (kind == "likelihood")
      expect_identical(serial[[2]]$lik_laplace[[1]]$check$status, "error")
    else expect_identical(serial[[2]]$burn_mode_stats$laplace_status[1], "error")
  }
})
