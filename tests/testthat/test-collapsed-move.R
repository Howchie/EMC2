local_rng_guard()
ok_check <- list(pareto_k = 0, outside = 0, laplace_status = "ok")

# A small hierarchy with an inverse-Wishart block (parameters 1-3), an
# unblocked parameter (4), subjects with and without a likelihood fit.
collapsed_fixture <- function(type) {
  set.seed(7)
  p <- 4; n <- 5
  fits <- lapply(seq_len(n - 1), function(s) {
    A <- matrix(rnorm(p * p), p)
    list(prec = crossprod(A) + diag(0.3, p), lin = rnorm(p), check = ok_check)
  })
  fits[n] <- list(NULL)
  sampler <- list(type = type, nuisance = rep(FALSE, p), n_pars = p, n_subjects = n,
                  is_blocked = c(TRUE, TRUE, TRUE, FALSE), par_group = c(1, 1, 1, 2),
                  prior = list(theta_mu_mean = rep(0, p), theta_mu_invar = diag(0.5, p),
                               v = 2, A = c(0.3, 1, 2, 0.7), shape = rep(2, p),
                               rate = rep(0.3, p)),
                  lik_laplace = fits)
  ctx <- EMC2:::.emc_collapsed_context(sampler)
  list(sampler = sampler, ctx = ctx, phi = rnorm(ctx$chart$dim, 0, 0.4),
       alpha = matrix(rnorm(p * n), p))
}

test_that("the chart is unconstrained, inverts, and has the stated Jacobian", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- collapsed_fixture(type)
    ch <- f$ctx$chart
    expect_equal(ch$dim, if (type == "standard") 4 + 6 + 1 else 8)
    st <- EMC2:::.emc_collapsed_from_phi(f$phi, ch)
    expect_equal(EMC2:::.emc_collapsed_to_phi(st$mu, st$V, ch), f$phi, tolerance = 1e-12)
    # Any phi gives a positive definite V.
    wild <- EMC2:::.emc_collapsed_from_phi(f$phi + rnorm(ch$dim, 0, 5), ch)
    expect_true(all(eigen(wild$V, only.values = TRUE)$values > 0))
    # Numerical Jacobian onto mu and the free entries of V; the stated one
    # drops the constant 2 per variance.
    free <- which(lower.tri(st$V, diag = TRUE) & st$V != 0)
    map <- function(ph) {
      s <- EMC2:::.emc_collapsed_from_phi(ph, ch)
      c(s$mu, s$V[free])
    }
    h <- 1e-6
    J <- vapply(seq_along(f$phi), function(i) {
      e <- replace(numeric(ch$dim), i, h)
      (map(f$phi + e) - map(f$phi - e)) / (2 * h)
    }, numeric(length(free) + 4L))
    expect_equal(as.numeric(determinant(J)$modulus),
                 st$log_jac + log(2) * (sum(vapply(ch$blocks, `[[`, 0L, "d")) +
                                          length(ch$gamma)), tolerance = 1e-6)
  }
})

test_that("the surrogate gradient matches finite differences", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- collapsed_fixture(type)
    sur <- function(ph) EMC2:::.emc_collapsed_surrogate(ph, f$ctx)
    g <- sur(f$phi)$grad
    h <- 1e-5
    num <- vapply(seq_along(f$phi), function(i) {
      e <- replace(numeric(length(f$phi)), i, h)
      (sur(f$phi + e)$lp - sur(f$phi - e)$lp) / (2 * h)
    }, numeric(1))
    expect_equal(g, num, tolerance = 1e-6)
  }
})

test_that("a skipped fractional collapsed move avoids rebuilding its context", {
  withr::local_options(list(emc2.group_move = "shear", emc2.collapsed_move = TRUE))
  sampler <- list(type = "standard", nuisance = rep(FALSE, 4L), gd = NULL,
                  group_designs = list(), n_subjects = 5L, n_pars = 4L,
                  collapsed = list(iterations = 0L,
                    plan = list(C = diag(11L), eps = 0.1, rounds = 0.25)))
  pars <- list()
  proposals <- matrix(0, 5L, 5L)
  result <- withr::with_seed(1L, {
    sampler$rng <- list(gibbs = get(".Random.seed", envir = globalenv()))
    EMC2:::.emc_collapsed_move(sampler, pars, proposals, "sample")
  })
  expect_identical(result$pars, pars)
  expect_identical(result$proposals, proposals)
  expect_identical(result$rounds, 0L)
  expect_equal(result$sampler$collapsed$iterations, 1L)
  expect_false(identical(result$sampler$rng$gibbs, sampler$rng$gibbs))
})

test_that("f / phat is the likelihood over its approximation at fixed residuals", {
  # f(phi) is the joint density in the (phi, r) chart: likelihood, subject
  # prior, group prior, chart Jacobian and the frame Jacobian prod_s 1/|R_s|.
  for (type in c("standard", "diagonal-gamma")) {
    f <- collapsed_fixture(type)
    ctx <- f$ctx
    p <- f$sampler$n_pars; n <- f$sampler$n_subjects
    ll <- function(A) colSums(-abs(A - 0.3)^1.5)
    frame <- function(ph) {
      s <- EMC2:::.emc_collapsed_from_phi(ph, ctx$chart)
      list(st = s, state = EMC2:::.emc_shear_state(s$mu, s$V, ctx$fits, ctx$prior, ctx$spec))
    }
    cur <- frame(f$phi)$state
    r <- vapply(seq_len(n), function(s) drop(cur$fm$R[, , s] %*% (f$alpha[, s] - cur$fm$m[, s])),
                numeric(p))
    set.seed(3)
    gap <- vapply(1:6, function(k) {
      ph <- f$phi + rnorm(ctx$chart$dim, 0, 0.3)
      fr <- frame(ph)
      A <- vapply(seq_len(n), function(s)
        fr$state$fm$m[, s] + drop(backsolve(fr$state$fm$R[, , s], r[, s])), numeric(p))
      f_val <- sum(ll(A)) + EMC2:::.emc_shear_log_hier(A, fr$state) + fr$state$lp +
        fr$st$log_jac - fr$state$logdet
      f_val - EMC2:::.emc_collapsed_surrogate(ph, ctx)$lp -
        sum(ll(A) - EMC2:::.emc_collapsed_log_lhat(A, ctx$fits))
    }, numeric(1))
    expect_lt(diff(range(gap)), 1e-9)
  }
})

test_that("the least-squares fit recovers a quadratic log likelihood", {
  set.seed(9)
  p <- 3
  H <- crossprod(matrix(rnorm(p * p), p)) + diag(0.5, p)
  lin <- rnorm(p)
  A <- matrix(rnorm(400 * p), ncol = p) %*% chol(solve(H)) + 0.4
  y <- 2 - 0.5 * rowSums((A %*% H) * A) + drop(A %*% lin)
  fit <- EMC2:::.emc_collapsed_fit_subject(A, y)
  expect_equal(fit$prec, H, tolerance = 1e-8)
  expect_equal(fit$lin, lin, tolerance = 1e-8)
  # Too few draws for the regression: no fit.
  expect_null(EMC2:::.emc_collapsed_fit_subject(A[1:10, ], y[1:10]))
  # A convex direction is projected away: the fit stays concave.
  y2 <- y + 3 * A[, 1]^2
  fit2 <- EMC2:::.emc_collapsed_fit_subject(A, y2)
  expect_true(all(eigen(fit2$prec, only.values = TRUE)$values > -1e-10))
})

test_that("empirical fit removes linear slopes in flat curvature directions", {
  set.seed(29)
  p <- 3
  Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
  H <- Q %*% diag(c(2, 1, 0)) %*% t(Q)
  lin <- drop(Q %*% c(1, -0.5, 9))
  A <- matrix(rnorm(400 * p), ncol = p)
  y <- 1 - 0.5 * rowSums((A %*% H) * A) + drop(A %*% lin)

  fit <- EMC2:::.emc_collapsed_fit_subject(A, y)
  e <- eigen(fit$prec, symmetric = TRUE)
  flat <- which(e$values <= max(e$values) * 1e-8)
  expect_length(flat, 1L)
  expect_equal(drop(crossprod(e$vectors[, flat, drop = FALSE], fit$lin)),
               0, tolerance = 1e-7)
  expect_gt(sqrt(sum(fit$lin^2)), 0.1)
})

test_that("safe base fits are cached while pooled conditioning is refreshed", {
  f <- collapsed_fixture("diagonal-gamma")
  H <- diag(c(2, 1, 0, 0))
  fit <- list(prec = H, lin = c(1, -2, 7, 8),
              cross = matrix(seq_len(8), nrow = 4), check = ok_check)
  f$sampler$collapsed <- list(fits = rep(list(fit), f$sampler$n_subjects))
  f$sampler$pooled <- list(value = c(0.5, -0.25))

  first <- EMC2:::.emc_collapsed_context(f$sampler)
  base <- first$base_fits[[1]]
  expect_true(isTRUE(base$.emc_collapsed_safe))
  expect_equal(base$lin, c(1, -2, 0, 0))
  expect_equal(base$cross, matrix(c(1, 2, 0, 0, 5, 6, 0, 0), nrow = 4))

  f$sampler$collapsed$fits <- first$base_fits
  f$sampler$pooled$value <- c(1, 2)
  second <- EMC2:::.emc_collapsed_context(f$sampler)
  expect_identical(second$base_fits[[1]], base)
  expect_true(isTRUE(second$fits[[1]]$.emc_collapsed_safe))
  expect_equal(second$fits[[1]]$lin, c(-10, -16, 0, 0))

  fallback <- f$sampler
  fallback$collapsed$fits <- NULL
  expect_false(EMC2:::.emc_collapsed_context(fallback)$cacheable)
})

test_that("safe fits project rotated cross terms into positive curvature", {
  Q <- 0.5 * rbind(c(1, 1, 1, 1), c(1, -1, 1, -1),
                   c(1, 1, -1, -1), c(1, -1, -1, 1))
  U <- Q[, 1:2, drop = FALSE]
  P <- tcrossprod(U)
  H <- U %*% diag(c(2, 1)) %*% t(U)
  cross <- matrix(seq_len(8), nrow = 4)
  fit <- EMC2:::.emc_collapsed_safe_fit(list(
    prec = H, lin = c(1, -2, 7, 8), cross = cross))

  expect_equal(fit$prec, H, tolerance = 1e-12)
  expect_equal(fit$lin, drop(P %*% c(1, -2, 7, 8)), tolerance = 1e-12)
  expect_equal(fit$cross, P %*% cross, tolerance = 1e-12)
})

# A plan as .emc_collapsed_set_plan() leaves it.
collapsed_plan <- function(p, dim, rounds, eps = 0.4, scale = 0.5) {
  zero <- list(E = matrix(0, 2 * p, 2 * p), n = 0L, acc = 0, moved = 0L)
  list(C = diag(scale, dim), eps = eps, rounds = rounds, dim = dim, last = zero,
       block = zero)
}

test_that("the collapsed move leaves the joint posterior invariant (diagonal-gamma)", {
  # The shear test's hierarchy: y_s | a_s ~ N(a_s, S), a_s ~ N(mu, V), mu ~
  # N(0, 2 I), gamma(2, 0.3) on each group precision, all diagonal, in the
  # funnel regime. Exact draws: (mu_k, log sd_k) from a fine grid of their
  # marginal posterior, then a_s.
  p <- 2; n <- 5
  S <- c(2, 1.5); m0 <- 0; P0 <- 0.5; shape <- 2; rate <- 0.3
  set.seed(21)
  y <- matrix(rnorm(p * n, 0.5, 1.5), p)
  N <- 4000
  grid_draws <- function(k) {
    mu_g <- seq(-6, 6, length.out = 301); ls_g <- seq(-5, 2.5, length.out = 301)
    lp <- outer(mu_g, ls_g, Vectorize(function(m, l) {
      tau <- exp(-2 * l)
      dnorm(m, m0, sqrt(1 / P0), log = TRUE) +
        sum(dnorm(y[k, ], m, sqrt(S[k] + exp(2 * l)), log = TRUE)) +
        dgamma(tau, shape, rate, log = TRUE) + log(2 * tau)
    }))
    w <- exp(lp - max(lp))
    cell <- sample(length(w), N, TRUE, prob = w)
    hm <- diff(mu_g[1:2]); hl <- diff(ls_g[1:2])
    cbind(mu = mu_g[row(w)[cell]] + runif(N, -hm / 2, hm / 2),
          lsd = ls_g[col(w)[cell]] + runif(N, -hl / 2, hl / 2))
  }
  g <- lapply(seq_len(p), grid_draws)
  mu <- sapply(g, function(x) x[, "mu"]); lsd <- sapply(g, function(x) x[, "lsd"])
  a <- array(NA_real_, c(N, p, n))
  for (k in seq_len(p)) {
    v <- exp(2 * lsd[, k])
    post_var <- 1 / (1 / S[k] + 1 / v)
    for (s in seq_len(n))
      a[, k, s] <- rnorm(N, post_var * (y[k, s] / S[k] + mu[, k] / v), sqrt(post_var))
  }
  before <- cbind(mu, lsd, matrix(a, N))

  Si <- diag(1 / S)
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, ...) {
      dd <- sweep(proposals, 2, dadm)
      -0.5 * rowSums((dd %*% Si) * dd)
    },
    .package = "EMC2")
  withr::local_options(list(emc2.collapsed_move = TRUE))
  pn <- c("x1", "x2")
  fits <- list(exact = lapply(seq_len(n), function(s)
                 list(prec = Si, lin = drop(Si %*% y[, s]), check = ok_check)),
               wrong = lapply(seq_len(n), function(s)
                 list(prec = 4 * Si, lin = drop(4 * Si %*% (y[, s] + 0.5)), check = ok_check)))
  prior <- list(theta_mu_mean = rep(m0, p), theta_mu_invar = diag(P0, p),
                shape = rep(shape, p), rate = rep(rate, p))
  # Exact and wrong approximations, then a broken control that drops Lhat
  # from the correction, which must fail.
  runs <- list(list(lik = fits$exact), list(lik = fits$wrong),
               list(lik = fits$wrong, broken = TRUE))
  for (run in runs) {
    if (isTRUE(run$broken)) {
      local_mocked_bindings(.emc_collapsed_log_lhat = function(alpha, fits)
        numeric(ncol(alpha)), .package = "EMC2")
    }
    sampler <- list(type = "diagonal-gamma", nuisance = rep(FALSE, p), n_subjects = n,
                    n_pars = p, par_names = pn, data = lapply(seq_len(n), function(s) y[, s]),
                    model = NULL, prior = prior, lik_laplace = run$lik, rng = list(),
                    collapsed = list(plan = collapsed_plan(p, 2 * p, 2), fits = run$lik))
    set.seed(5)
    after <- t(vapply(seq_len(N), function(i) {
      V <- diag(exp(2 * lsd[i, ]))
      al <- matrix(a[i, , ], p, dimnames = list(pn, NULL))
      ll <- vapply(seq_len(n), function(s) -0.5 * sum((al[, s] - y[, s])^2 / S), 0)
      out <- EMC2:::.emc_collapsed_move(sampler, list(tmu = setNames(mu[i, ], pn), tvar = V,
                                                      tvinv = solve(V)),
                                        rbind(al, ll), "sample")
      if (i == 1L)
        expect_true(all(vapply(out$sampler$collapsed$fits, function(f)
          is.null(f) || isTRUE(f$.emc_collapsed_safe), logical(1))))
      c(out$pars$tmu, 0.5 * log(diag(out$pars$tvar)), out$proposals[1:p, ])
    }, numeric(2 * p + p * n)))
    expect_gt(mean(after[, p + 1] != before[, p + 1]), 0.1)
    se <- apply(before, 2, sd) / sqrt(N)
    z_mean <- abs(colMeans(after) - colMeans(before)) / se
    z_var <- abs(apply(after, 2, var) / apply(before, 2, var) - 1) / sqrt(2 / N)
    if (isTRUE(run$broken)) {
      expect_gt(max(z_mean, z_var), 6)
    } else {
      expect_true(all(z_mean < 4.5))
      expect_true(all(z_var < 4.5))
    }
  }
})

test_that("the collapsed move with transport walls leaves a bounded joint posterior invariant", {
  local_test_rng(kind = "Mersenne-Twister")
  # The diagonal hierarchy above with the first parameter bounded below for
  # every subject; the move sees each subject's wall. Exact and wrong
  # approximations, then a broken control that drops the wall map's
  # Jacobian, which must fail.
  N <- 4000
  h <- diag_hierarchy_draws(N, b = -0.4)
  p <- h$p; n <- h$n; y <- h$y; Si <- h$Si
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...) h$ll(proposals, dadm),
                        .package = "EMC2")
  withr::local_options(list(emc2.collapsed_move = TRUE))
  pn <- c("x1", "x2")
  fits <- list(exact = lapply(seq_len(n), function(s)
                 list(prec = Si, lin = drop(Si %*% y[, s]), check = ok_check)),
               wrong = lapply(seq_len(n), function(s)
                 list(prec = 4 * Si, lin = drop(4 * Si %*% (y[, s] + 0.5)), check = ok_check)))
  runs <- list(list(lik = fits$exact), list(lik = fits$wrong),
               list(lik = fits$exact, broken = TRUE))
  for (run in runs) {
    if (isTRUE(run$broken)) {
      transport <- EMC2:::.emc_shear_transport
      local_mocked_bindings(.emc_shear_transport = function(...) {
        out <- transport(...)
        out$jac[] <- 0
        out
      }, .package = "EMC2")
    }
    sampler <- list(type = "diagonal-gamma", nuisance = rep(FALSE, p), n_subjects = n,
                    n_pars = p, par_names = pn, data = lapply(seq_len(n), function(s) y[, s]),
                    model = NULL, prior = h$prior, lik_laplace = run$lik, rng = list(),
                    shear = list(walls = h$walls),
                    collapsed = list(plan = collapsed_plan(p, 2 * p, 2), fits = run$lik))
    set.seed(5)
    after <- t(vapply(seq_len(N), function(i) {
      V <- diag(exp(2 * h$lsd[i, ]))
      al <- matrix(h$a[i, , ], p, dimnames = list(pn, NULL))
      ll <- vapply(seq_len(n), function(s) -0.5 * sum((al[, s] - y[, s])^2 / h$S), 0)
      out <- EMC2:::.emc_collapsed_move(sampler, list(tmu = setNames(h$mu[i, ], pn), tvar = V,
                                                      tvinv = solve(V)),
                                        rbind(al, ll), "sample")
      c(out$pars$tmu, 0.5 * log(diag(out$pars$tvar)), out$proposals[1:p, ])
    }, numeric(2 * p + p * n)))
    expect_true(all(after[, 2 * p + seq(1, p * n, by = p)] >= h$b))
    expect_gt(mean(after[, p + 1] != h$before[, p + 1]), 0.05)
    se <- apply(h$before, 2, sd) / sqrt(N)
    z_mean <- abs(colMeans(after) - colMeans(h$before)) / se
    z_var <- abs(apply(after, 2, var) / apply(h$before, 2, var) - 1) / sqrt(2 / N)
    if (isTRUE(run$broken)) {
      expect_gt(max(z_mean, z_var), 6)
    } else {
      expect_true(all(z_mean < 4.5))
      expect_true(all(z_var < 4.5))
    }
  }
})

test_that("the collapsed move keeps the joint distribution with a full covariance block", {
  local_test_rng(kind = "Mersenne-Twister")
  # Successive-conditional simulation (Geweke 2004): alternate the move with
  # exact draws of a | mu, V, y and y | a. With a valid move the chain's
  # stationary law is the prior joint, so (mu, log sd, correlation) must
  # match independent prior draws. Parameters 1-2 form an inverse-Wishart
  # block, 3 is unblocked, a_half integrated out by the move. The
  # approximations are deliberately wrong.
  p <- 3; n <- 4; v <- 2; Avec <- c(0.8, 1.2, 1); P0 <- 1
  S <- c(0.6, 0.9, 0.7); Si <- diag(1 / S)
  prior <- list(theta_mu_mean = rep(0, p), theta_mu_invar = diag(P0, p), v = v, A = Avec)
  draw_V <- function() {
    a_half <- 1 / rgamma(p, shape = 0.5, rate = 1 / Avec^2)
    V <- matrix(0, p, p)
    b <- 1:2
    V[b, b] <- solve(rWishart(1, v + 1, solve(2 * v * diag(1 / a_half[b])))[, , 1])
    V[3, 3] <- 1 / rgamma(1, v / 2, rate = v / a_half[3])
    V
  }
  feats <- function(mu, V) c(mu, 0.5 * log(diag(V)), atanh(cov2cor(V)[1, 2]))
  set.seed(31)
  M <- 6000
  ref <- t(replicate(M, feats(rnorm(p, 0, 1 / sqrt(P0)), draw_V())))
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, ...) {
      dd <- sweep(proposals, 2, dadm)
      -0.5 * rowSums((dd %*% Si) * dd)
    },
    .package = "EMC2")
  withr::local_options(list(emc2.collapsed_move = TRUE))
  pn <- paste0("x", 1:p)
  base <- list(type = "standard", nuisance = rep(FALSE, p), n_subjects = n, n_pars = p,
               par_names = pn, is_blocked = c(TRUE, TRUE, FALSE), par_group = c(1, 1, 2),
               model = NULL, prior = prior, rng = list())
  ctx_dim <- EMC2:::.emc_collapsed_chart(EMC2:::.emc_shear_v_spec(base), p)$dim
  run_chain <- function(broken) {
    if (broken) {
      local_mocked_bindings(.emc_collapsed_log_lhat = function(alpha, fits)
        numeric(ncol(alpha)), .package = "EMC2")
    }
    # The move draws from its own stream, as in the sampler; sharing the
    # test's stream would replay the draws of a and y as its momenta.
    set.seed(33)
    gibbs <- .Random.seed
    set.seed(32)
    mu <- rnorm(p, 0, 1); V <- draw_V()
    a <- mu + t(chol(V)) %*% matrix(rnorm(p * n), p)
    y <- a + sqrt(S) * matrix(rnorm(p * n), p)
    sampler <- c(base, list(collapsed = list(plan = collapsed_plan(p, ctx_dim, 1, eps = 0.35,
                                                                   scale = 0.8))))
    sampler$rng$gibbs <- gibbs
    out <- matrix(NA_real_, M, 2 * p + 1)
    for (t in seq_len(M)) {
      sampler$data <- lapply(seq_len(n), function(s) y[, s])
      sampler$lik_laplace <- lapply(seq_len(n), function(s)
        list(prec = 3 * Si, lin = drop(3 * Si %*% (y[, s] - 0.4)), check = ok_check))
      al <- matrix(a, p, dimnames = list(pn, NULL))
      ll <- vapply(seq_len(n), function(s) -0.5 * sum((al[, s] - y[, s])^2 / S), 0)
      own <- .Random.seed
      res <- EMC2:::.emc_collapsed_move(sampler, list(tmu = setNames(mu, pn), tvar = V,
                                                      tvinv = solve(V), a_half = setNames(rep(1, p), pn)),
                                        rbind(al, ll), "sample")
      assign(".Random.seed", own, envir = globalenv())
      sampler <- res$sampler
      mu <- as.numeric(res$pars$tmu); V <- unname(as.matrix(res$pars$tvar))
      Q <- solve(V) + Si
      Qi <- solve(Q)
      a <- Qi %*% (drop(solve(V, mu)) + Si %*% y) + t(chol(Qi)) %*% matrix(rnorm(p * n), p)
      y <- a + sqrt(S) * matrix(rnorm(p * n), p)
      out[t, ] <- feats(mu, V)
    }
    out
  }
  # Batch-means standard errors for the chain, iid ones for the reference.
  z_of <- function(chain) {
    nb <- 30
    bm <- apply(chain, 2, function(x) tapply(x, cut(seq_along(x), nb, labels = FALSE), mean))
    se_chain <- apply(bm, 2, sd) / sqrt(nb)
    abs(colMeans(chain) - colMeans(ref)) / sqrt(se_chain^2 + apply(ref, 2, var) / M)
  }
  good <- run_chain(FALSE)
  expect_true(all(z_of(good) < 4))
  bad <- run_chain(TRUE)
  expect_gt(max(z_of(bad)), 5)
})

test_that("the surrogate refuses states outside the cone and keeps the chain's state", {
  f <- collapsed_fixture("standard")
  bad <- replace(f$phi, f$ctx$chart$blocks[[1]]$at[1], 800)
  expect_identical(EMC2:::.emc_collapsed_surrogate(bad, f$ctx)$lp, -Inf)
  set.seed(1)
  tr <- EMC2:::collapsed_trajectories_cpp(bad, diag(f$ctx$chart$dim), 0.3, 2L, 1.5, f$ctx$cpp)
  expect_equal(as.numeric(tr$phi), bad)
  expect_identical(tr$lp, -Inf)
})

test_that("the allocation buys collapsed moves only where they pay", {
  # Two group coordinates, slow under the rest of the iteration.
  units <- list(list(kind = "fixed", rho = 0.95, cost = 0),
                list(kind = "fixed", rho = 0.9, cost = 0))
  cheap <- list(idx = 1:2, g = c(0.3, 0.4), cost = 0.05, cap = 4)
  sol <- EMC2:::.emc_alloc_solve(units, 1, joint = cheap)
  expect_gt(sol$joint, 0)
  expect_lt(sol$tau, sol$tau0)
  reached <- max(0.95 * EMC2:::.emc_alloc_factor(0.3, sol$joint),
                 0.9 * EMC2:::.emc_alloc_factor(0.4, sol$joint))
  expect_lte(reached, sol$rho + 1e-12)
  # The same move at a cost that outweighs its gain, or with no effect.
  expect_identical(EMC2:::.emc_alloc_solve(units, 1, joint = replace(cheap, "cost", list(50)))$joint, 0)
  expect_identical(EMC2:::.emc_alloc_solve(units, 1, joint = replace(cheap, "g", list(c(1, 1))))$joint, 0)
  # Without the move the solver is unchanged.
  expect_identical(EMC2:::.emc_alloc_solve(units, 1)$joint, 0)
  # A weak, cheap move: the optimum lies past the shared round cap, so the
  # cap binds until it is measured at higher rates, then the optimum does.
  weak <- list(list(kind = "fixed", rho = 0.99, cost = 0))
  wj <- list(idx = 1L, g = 0.9, cost = 0.02, cap = 4)
  expect_identical(EMC2:::.emc_alloc_solve(weak, 1, joint = wj)$joint, 4)
  expect_identical(EMC2:::.emc_alloc_solve(weak, 1, joint = replace(wj, "cap", 8))$joint, 8)
  free <- EMC2:::.emc_alloc_solve(weak, 1, joint = replace(wj, "cap", 64))$joint
  expect_gt(free, 20)
  expect_lt(free, 64)
  expect_identical(max(EMC2:::.emc_alloc_joint_grid(9.3)), 9.3)
  expect_identical(EMC2:::.emc_alloc_joint_grid(0.25), c(0.125, 0.25))
  # Caps between grid steps, at either end of the eighths.
  expect_identical(utils::tail(EMC2:::.emc_alloc_joint_grid(4.25), 2), c(4, 4.25))
  expect_identical(EMC2:::.emc_alloc_joint_grid(0.1), 0.1)
  # The cap: the shared round cap until measured, then twice the measured
  # rate where that is larger.
  ms <- list(s_rho = matrix(0.5, 1, 1), g_rho = 0.99, g_fac = NA_real_,
             g_measured = FALSE, C_rest = 1, c_round = NA_real_, coll = TRUE,
             c_fac = 0.9, c_measured = FALSE, c_coll = 0.02, c_rate = 0,
             W = list(diag(1)), s_var = matrix(1, 1, 1), c_hmc = NA_real_)
  cap_of <- function(...) EMC2:::.emc_alloc_units(modifyList(ms, list(...)), NULL,
                                                  FALSE, FALSE)$joint$cap
  expect_identical(cap_of(), EMC2:::.EMC_ALLOC_MAX_ROUNDS)
  expect_identical(cap_of(c_measured = TRUE, c_rate = 0.25), EMC2:::.EMC_ALLOC_MAX_ROUNDS)
  expect_identical(cap_of(c_measured = TRUE, c_rate = 3.9), 7.8)
  # Costs are wall time on the cores the chain holds.
  expect_equal(EMC2:::.emc_reserved_cost(0.2, 5), 1)
  expect_equal(EMC2:::.emc_reserved_cost(0.2, 0), 0.2)
})

test_that("collapsed plans pool production and retain same-stage measurements", {
  f <- collapsed_fixture("diagonal-gamma")
  p <- f$sampler$n_pars
  st <- EMC2:::.emc_collapsed_from_phi(f$phi, f$ctx$chart)
  chain <- f$sampler
  chain$samples <- list(idx = 1L, theta_mu = matrix(st$mu, p),
                        theta_var = array(st$V, c(p, p, 1)))
  chain$rng <- list()
  C <- diag(0.3, f$ctx$chart$dim)
  block <- function(n, scale) list(E = diag(scale, 2 * p), n = n, acc = n / 2, moved = n)
  rebuild <- function(chain, stage, blk) {
    chain$collapsed$plan$block <- blk
    EMC2:::.emc_collapsed_set_plan(chain, C, stage)
  }
  chain <- EMC2:::.emc_collapsed_set_plan(chain, C, "adapt")
  # The first production plan measures the adapt block alone and pools
  # nothing from it.
  chain <- rebuild(chain, "sample", block(8L, 1))
  expect_identical(chain$collapsed$plan$last$n, 8L)
  expect_identical(chain$collapsed$plan$pool$n, 0L)
  # Production blocks then add up, a thin one included, the earlier ones
  # discounted at each rebuild.
  chain <- rebuild(chain, "sample", block(20L, 2))
  chain <- rebuild(chain, "sample", block(3L, 3))
  last <- chain$collapsed$plan$last
  d <- EMC2:::.EMC_ALLOC_EVIDENCE_DECAY
  expect_equal(last$n, d * 20 + 3)
  expect_equal(last$E, diag(d * 2 + 3, 2 * p))
  expect_equal(last$acc, d * 10 + 1.5)
  # Outside production, nonempty blocks are not pooled. A thin or empty block
  # in the same stage retains the last usable measurement.
  chain <- rebuild(chain, "burn", block(6L, 1))
  expect_identical(chain$collapsed$plan$last$n, 6L)
  expect_identical(chain$collapsed$plan$pool$n, 0L)
  chain <- rebuild(chain, "burn", block(3L, 3))
  expect_identical(chain$collapsed$plan$last$n, 6L)
  chain <- rebuild(chain, "burn", block(0L, 0))
  expect_identical(chain$collapsed$plan$last$n, 6L)
  # A new stage starts with its own proposal geometry and gets a fresh probe.
  chain <- rebuild(chain, "adapt", block(0L, 0))
  expect_identical(chain$collapsed$plan$last$n, 0L)
})
