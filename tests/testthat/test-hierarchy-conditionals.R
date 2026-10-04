# Analytic checks of hierarchy full conditionals. Each multivariate normal draw
# is mocked to return its conditional mean (and zero residuals for the factor
# scores), so the coefficients a Gibbs step returns ARE the conditional means it
# requested. They are compared with the closed-form regression, written out
# independently here, for rows that mix free and fixed coefficients.

hc_mean_rmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), ...) {
  if (n > 1L) return(matrix(0, n, ncol(as.matrix(sigma))))
  matrix(as.numeric(mean), 1L)
}

hc_latent_sampler <- function(type, p, n, nf, ss = NULL, Lambda_mat = NULL) {
  pn <- letters[seq_len(p)]
  s <- list(type = type, n_pars = p, n_subjects = n, par_names = pn,
            nuisance = rep(FALSE, p), marginalised_idx = rep(FALSE, p),
            model = NULL, marginalise = NULL)
  s <- switch(type,
    factor = EMC2:::add_info_factor(s, n_factors = nf, Lambda_mat = Lambda_mat),
    infnt_factor = EMC2:::add_info_infnt_factor(s, n_factors = nf),
    SEM = EMC2:::add_info_SEM(s, sem_settings = ss))
  store <- get(paste0("sample_store_", type), asNamespace("EMC2"))
  s$samples <- store(data.frame(subjects = factor(seq_len(n))), pn, iters = 1,
                     stage = "sample", integrate = TRUE, is_nuisance = rep(FALSE, p),
                     n_factors = nf, sem_settings = ss,
                     Lambda_mat = if (type == "factor") attr(s, "Lambda_mat"))
  s$samples$idx <- 1
  start <- get(paste0("get_startpoints_", type), asNamespace("EMC2"))(s, rep(0, p), diag(p))
  s$samples <- EMC2:::fill_samples(s$samples, start, rbind(matrix(0, p, n), 0), 1, p, type)
  s
}

test_that("factor loadings regress on the response net of fixed loadings", {
  Lmat <- rbind(c(1, 0), c(Inf, 1), c(Inf, 0.5), c(Inf, Inf))
  s <- hc_latent_sampler("factor", 4L, 7L, 2L, Lambda_mat = Lmat)
  s$prior$theta_lambda_invar <- c(1, 2, 0.5, 4)   # heterogeneous by row
  set.seed(4)
  f <- matrix(rnorm(14), 7, 2) %*% matrix(c(1, 0.8, 0, 0.6), 2)  # correlated
  alpha <- t(f %*% t(rbind(c(1, 0), c(0.3, 1), c(-0.4, 0.5), c(0.7, 0.2))) +
               matrix(rnorm(28, 0, 0.3), 7))
  s$samples$lambda_untransf[, , 1] <- ifelse(Lmat == Inf, 0.2, Lmat)
  local_mocked_bindings(rmvnorm = hc_mean_rmvnorm, .package = "EMC2")
  set.seed(9)
  out <- EMC2:::gibbs_step_factor(s, alpha)
  eta <- out$eta
  yt <- sweep(t(alpha), 2, as.numeric(out$tmu))
  for (j in 2:4) {
    free <- Lmat[j, ] == Inf
    y <- yt[, j] - eta[, !free, drop = FALSE] %*% Lmat[j, !free]
    e <- out$sig_err_inv[j]
    Ef <- eta[, free, drop = FALSE]
    V <- solve(e * crossprod(Ef) + diag(s$prior$theta_lambda_invar[j], sum(free)))
    expect_equal(out$lambda_untransf[j, free], as.numeric(V %*% (e * crossprod(Ef, y))),
                 tolerance = 1e-10, label = paste("row", j))
  }
  expect_equal(out$lambda_untransf[Lmat != Inf], Lmat[Lmat != Inf])
})

test_that("infinite-factor shrinkage rates use the refreshed cumulative products", {
  s <- hc_latent_sampler("infnt_factor", 4L, 6L, 3L)
  delta0 <- c(1, 2, 3)
  s$samples$delta[, 1] <- delta0
  set.seed(2)
  alpha <- matrix(rnorm(24), 4, 6)
  draws <- c(1.5, 4, 5)
  rates <- numeric(0)
  local_mocked_bindings(rgamma = function(n, shape, rate) {
    if (n == 1L) {
      rates <<- c(rates, rate)
      return(draws[length(rates)])
    }
    stats::rgamma(n, shape, rate)
  }, .package = "EMC2")
  out <- EMC2:::gibbs_step_infnt_factor(s, alpha)
  expect_equal(out$delta, draws)
  cs <- colSums(out$psi * out$lambda^2)
  pr <- s$prior
  # delta_h's rate: prod over l <= k, l != h of delta_l, with deltas before h
  # already updated and those after h at their previous values.
  for (h in 2:3) {
    cur <- c(draws[seq_len(h - 1)], delta0[h:3])
    w <- vapply(h:3, function(k) prod(cur[setdiff(seq_len(k), h)]), 0)
    expect_equal(rates[h], pr$bd2 + 0.5 * sum(w * cs[h:3]), tolerance = 1e-10,
                 label = paste("delta", h, "rate"))
  }
})

hc_sem_settings <- function(x) {
  list(Lambda_mat = rbind(c(1, 0), c(Inf, 0.4), c(0, 1), c(Inf, Inf)),
       B_mat = matrix(c(0, Inf, 0, 0), 2), K_mat = matrix(c(Inf, 0.3, 0, Inf), 4),
       G_mat = matrix(c(Inf, 0.5), 2), factor_groups = c(1, 1), covariates = x)
}

hc_sem_state <- function(s) {
  ss <- s$sem_settings
  fill <- function(M, v) ifelse(M == Inf, v, M)
  s$samples$lambda[, , 1] <- fill(ss$Lambda_mat, 0.6)
  s$samples$K[, , 1] <- fill(ss$K_mat, -0.2)
  s$samples$G[, , 1] <- fill(ss$G_mat, 0.7)
  s$samples$B[, , 1] <- fill(ss$B_mat, 0.25)
  s$samples$delta_inv[, , 1] <- solve(matrix(c(1, 0.5, 0.5, 1.5), 2))
  s$samples$epsilon_inv[, 1] <- c(2, 3, 1.5, 2.5)
  s$samples$theta_mu[, 1] <- c(0.1, -0.1, 0.2, 0)
  s
}

test_that("SEM loading, K, G and B rows subtract fixed contributions", {
  set.seed(6)
  x <- matrix(rnorm(9), 9, 1, dimnames = list(NULL, "x"))
  s <- hc_sem_state(hc_latent_sampler("SEM", 4L, 9L, 2L, ss = hc_sem_settings(x)))
  s$prior$lambda_var <- c(0.7, 2)
  s$prior$K_var <- 1.5
  s$prior$B_var <- c(0.5, 3)
  s$prior$G_var <- 0.8
  ss <- s$sem_settings
  alpha <- matrix(rnorm(36), 4, 9)
  D0 <- s$samples$delta_inv[, , 1]
  B0 <- s$samples$B[, , 1]; G0 <- matrix(s$samples$G[, , 1], 2)
  local_mocked_bindings(rmvnorm = hc_mean_rmvnorm, .package = "EMC2")
  set.seed(1)
  out <- EMC2:::gibbs_step_SEM(s, alpha)
  eta <- out$eta
  # Loading / K rows.
  X <- cbind(x, eta)
  coef <- cbind(out$K, out$lambda)
  fixed <- cbind(ss$K_mat, ss$Lambda_mat)
  pv <- c(s$prior$K_var, s$prior$lambda_var)
  yt <- sweep(t(alpha), 2, as.numeric(out$tmu))
  for (j in c(1, 2, 4)) {
    free <- fixed[j, ] == Inf
    y <- yt[, j] - X[, !free, drop = FALSE] %*% fixed[j, !free]
    e <- out$epsilon_inv[j]
    Xf <- X[, free, drop = FALSE]
    V <- solve(e * crossprod(Xf) + diag(1 / pv[free], sum(free)))
    expect_equal(coef[j, free], as.numeric(V %*% (e * crossprod(Xf, y))),
                 tolerance = 1e-10, label = paste("loading row", j))
  }
  expect_equal(coef[fixed != Inf], fixed[fixed != Inf])
  # Structural rows, in order, with correlated disturbances (old delta_inv).
  Bc <- B0; Gc <- G0
  Z <- cbind(x, eta)
  R <- eta - eta %*% t(Bc) - x %*% t(Gc)
  sfix <- cbind(ss$G_mat, ss$B_mat)
  spv <- c(s$prior$G_var, s$prior$B_var)
  for (p in 1:2) {
    free <- sfix[p, ] == Inf
    resp <- eta[, p] - Z[, !free, drop = FALSE] %*% sfix[p, !free]
    ystar <- D0[p, p] * resp + R[, -p, drop = FALSE] %*% D0[-p, p]
    Zf <- Z[, free, drop = FALSE]
    V <- solve(D0[p, p] * crossprod(Zf) + diag(1 / spv[free], sum(free)))
    m <- as.numeric(V %*% crossprod(Zf, ystar))
    got <- c(out$G[p, ], out$B[p, ])[free]
    expect_equal(got, m, tolerance = 1e-10, label = paste("structural row", p))
    row <- c(Gc[p, ], Bc[p, ])
    row[free] <- m
    Gc[p, ] <- row[1]; Bc[p, ] <- row[-1]
    R[, p] <- eta[, p] - Z %*% row
  }
  # Subject prior = the shared moment helper.
  A <- out$lambda %*% solve(diag(2) - out$B)
  C <- out$K + A %*% out$G
  expect_equal(out$subj_mu, out$tmu + C %*% t(x), ignore_attr = TRUE)
  expect_equal(out$subj_var, A %*% solve(out$delta_inv) %*% t(A) +
                 diag(1 / out$epsilon_inv), ignore_attr = TRUE)
  expect_equal(out$tvar, out$subj_var + C %*% cov(x) %*% t(C), ignore_attr = TRUE)
})

test_that("SEM population moments keep the K x AG cross terms; implied mean adds K", {
  x <- matrix(c(0, 2), 2, 1)
  m <- EMC2:::.sem_moments(0, matrix(1), matrix(0), matrix(0.3), matrix(0.8),
                           matrix(1), 1e8, covariates = x)
  expect_equal(m$pop_mean, 1.1)
  expect_equal(as.numeric(m$pop_var - m$subj_var), (0.3 + 0.8)^2 * var(x[, 1]))
  # One parameter, one draw: dimensions survive.
  smp <- list(samples = list(theta_mu = matrix(0, 1, 1, dimnames = list("a", NULL)),
                             B = array(0, c(1, 1, 1)), G = array(0.8, c(1, 1, 1)),
                             K = array(0.3, c(1, 1, 1)), lambda = array(1, c(1, 1, 1))),
              sem_settings = list(covariates = x))
  expect_equal(EMC2:::get_mu_implied(smp, 1L), matrix(1.1, 1, 1, dimnames = list("a", NULL)))
})

test_that("SEM bridge group density is the subject-conditional density", {
  set.seed(8)
  n <- 6
  x <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "x"))
  ss <- hc_sem_settings(x)
  s <- hc_latent_sampler("SEM", 4L, n, 2L, ss = ss)
  prior <- s$prior
  mu <- c(0.1, -0.1, 0.2, 0)
  lam <- ifelse(ss$Lambda_mat == Inf, 0.6, ss$Lambda_mat)
  B <- ifelse(ss$B_mat == Inf, 0.25, ss$B_mat)
  K <- ifelse(ss$K_mat == Inf, -0.2, ss$K_mat)
  G <- ifelse(ss$G_mat == Inf, 0.7, ss$G_mat)
  eps <- c(2, 3, 1.5, 2.5)
  dinv <- solve(matrix(c(1, 0.5, 0.5, 1.5), 2))
  grp <- matrix(c(mu, lam[ss$Lambda_mat == Inf], B[ss$B_mat == Inf],
                  K[ss$K_mat == Inf], G[ss$G_mat == Inf], log(eps),
                  EMC2:::unwind_chol(dinv)), 1)
  info <- list(prior = prior, n_pars = 4L, n_factors = 2L, sem_settings = ss)
  a1 <- matrix(rnorm(4 * n), 4, n)
  a2 <- matrix(rnorm(4 * n), 4, n)
  as_list <- function(a) lapply(seq_len(n), function(k) matrix(a[, k], 1))
  b1 <- EMC2:::bridge_group_and_prior_and_jac_SEM(grp, as_list(a1), info)
  b2 <- EMC2:::bridge_group_and_prior_and_jac_SEM(grp, as_list(a2), info)
  mom <- EMC2:::.sem_moments(mu, lam, B, K, G, dinv, eps, covariates = x)
  dens <- function(a) sum(mvtnorm::dmvnorm(t(a - mom$subj_mu), sigma = mom$subj_var, log = TRUE))
  expect_equal(b1 - b2, dens(a1) - dens(a2), tolerance = 1e-8)
})

test_that("SEM rejects cyclic structural models and accepts DAGs", {
  expect_error(EMC2:::.sem_check_acyclic(matrix(c(0, Inf, Inf, 0), 2)), "cyclic")
  expect_error(EMC2:::.sem_check_acyclic(matrix(c(0, 0.2, Inf, 0), 2)), "cyclic")
  expect_error(EMC2:::.sem_check_acyclic(matrix(c(Inf, 0, 0, 0), 2)), "diagonal")
  three <- matrix(0, 3, 3); three[2, 1] <- three[3, 2] <- Inf; three[1, 3] <- 0.1
  expect_error(EMC2:::.sem_check_acyclic(three), "cyclic")
  three[1, 3] <- 0
  expect_true(EMC2:::.sem_check_acyclic(three))
  # A topological permutation is still accepted.
  perm <- c(3, 1, 2)
  expect_true(EMC2:::.sem_check_acyclic(three[perm, perm]))
})
