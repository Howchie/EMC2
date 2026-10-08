test_that("the shear preserves the joint posterior with a full standard covariance", {
  skip_model_validation()
  withr::local_seed(83)
  withr::local_options(list(emc2.group_move = "shear"))
  p <- 2L
  n <- 3L
  N <- 3000L
  pn <- c("a", "b")
  P0 <- solve(matrix(c(0.5, 0.2, 0.2, 0.4), p))
  H <- matrix(c(0.5, 0.25, 0.25, 0.6), p)
  y <- matrix(c(-0.7, -0.4, 0.1, 0.2, 0.8, 0.6), p)
  v <- 10
  A <- c(0.6, 0.8)
  exact_fit <- lapply(seq_len(n), function(s) list(prec = H, lin = drop(H %*% y[, s]),
    check = list(pareto_k = 0, outside = 0, laplace_status = "ok")))
  sampler <- list(type = "standard", n_pars = p, n_subjects = n, par_names = pn,
                  nuisance = rep(FALSE, p), par_group = rep(1L, p), is_blocked = rep(TRUE, p),
                  prior = list(theta_mu_mean = rep(0, p), theta_mu_invar = P0, v = v, A = A),
                  data = lapply(seq_len(n), function(s) y[, s]), model = NULL,
                  lik_laplace = exact_fit, rng = list())
  # One round on every coordinate, along fixed directions, at a fixed step.
  plan <- list(dirs = diag(c(0.4, 0.4)), rounds = rep(1, 2L * p),
               log_step = numeric(2L * p),
               block = list(E = numeric(2L * p), n = integer(2L * p),
                            moved = integer(2L * p)))
  ll <- function(a) {
    residual <- a - y
    -colSums(residual * (H %*% residual)) / 2
  }
  # Independent rejection sampling from the generative prior, a_half
  # included. The Gaussian likelihood is bounded by one, so its value is the
  # acceptance probability.
  draws <- vector("list", N)
  k <- 0L
  while (k < N) {
    mu <- drop(t(chol(solve(P0))) %*% rnorm(p))
    a_half <- 1 / rgamma(p, shape = 1 / 2, rate = 1 / A^2)
    S <- 2 * v * diag(1 / a_half)
    V <- solve(stats::rWishart(1L, v + p - 1, solve(S))[, , 1])
    alpha <- mu + t(chol(V)) %*% matrix(rnorm(p * n), p)
    if (log(runif(1)) < sum(ll(alpha))) {
      k <- k + 1L
      draws[[k]] <- list(mu = mu, V = V, a_half = a_half, alpha = alpha, ll = ll(alpha))
    }
  }
  features <- function(mu, V, a_half, alpha) c(mu, mu[1] * mu[2], log(diag(V)), V[1, 2],
    V[1, 2] / sqrt(prod(diag(V))), log(a_half), alpha[, 1], alpha[1, 1] * alpha[2, 1],
    (alpha[1, 1] - mu[1]) * (alpha[2, 1] - mu[2]))
  before <- t(vapply(draws, function(d) features(d$mu, d$V, d$a_half, d$alpha), numeric(13)))
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...) {
    residual <- sweep(proposals, 2, dadm)
    -rowSums(residual * (residual %*% H)) / 2
  }, .package = "EMC2")
  wrong_fit <- lapply(seq_len(n), function(s) list(prec = 2 * H,
    lin = drop(2 * H %*% (y[, s] + 0.4)),
    check = list(pareto_k = 0, outside = 0, laplace_status = "ok")))
  for (fits in list(exact_fit, wrong_fit)) {
    sampler$lik_laplace <- fits
    sampler$shear <- list(plan = plan)
    after <- t(vapply(draws, function(d) {
      out <- EMC2:::.emc_shear_move(sampler,
        list(tmu = setNames(d$mu, pn), tvar = d$V, tvinv = solve(d$V), a_half = d$a_half),
        rbind(d$alpha, d$ll), "sample")
      features(out$pars$tmu, out$pars$tvar, out$pars$a_half,
               out$proposals[seq_len(p), , drop = FALSE])
    }, numeric(13)))
    difference <- after - before
    expect_gt(mean(abs(difference[, 6]) > 0), 0.01)
    expect_true(all(abs(colMeans(difference)) <= 5 * apply(difference, 2, sd) / sqrt(N) + 1e-12))
    center <- colMeans(before)
    square_difference <- sweep(after, 2, center)^2 - sweep(before, 2, center)^2
    expect_true(all(abs(colMeans(square_difference)) <
      5 * apply(square_difference, 2, sd) / sqrt(N) + 1e-12))
  }
})
