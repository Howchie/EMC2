local_rng_guard()
ok_check <- list(pareto_k = 0, outside = 0, laplace_status = "ok")

test_that("the shear plan is the centred Gibbs step's lag-1 operator", {
  set.seed(11)
  p <- 3; n <- 6
  V <- crossprod(matrix(rnorm(p * p), p)) / p + diag(0.2, p)
  fits <- lapply(seq_len(n), function(s) {
    A <- matrix(rnorm(p * p), p)
    list(prec = crossprod(A) * s / 3 + diag(0.1, p), lin = rnorm(p), check = ok_check)
  })
  P0 <- diag(0.5, p)
  plan <- EMC2:::.emc_shear_plan(fits, list(theta_mu_invar = P0), V)
  Vi <- solve(V)
  Qi <- lapply(fits, function(f) solve(f$prec + Vi))
  A <- solve(n * Vi + P0) %*% Vi %*% Reduce(`+`, Qi) %*% Vi
  expect_equal(sort(plan$lam), sort(Re(eigen(A)$values)), tolerance = 1e-10)
  # Each column moves mu by one collapsed posterior sd of its direction.
  Pc <- n * Vi + P0
  expect_equal(diag(t(plan$dirs) %*% Pc %*% plan$dirs), 1 / (1 - plan$lam), tolerance = 1e-8)
  expect_equal(plan$w, Reduce(`+`, lapply(Qi, function(q) diag(q %*% Vi))) / n, tolerance = 1e-12)
  # With a flat mean prior the mean and scale lag-1 totals agree.
  flat <- EMC2:::.emc_shear_plan(fits, list(theta_mu_invar = diag(0, p)), V)
  expect_equal(sum(flat$lam), sum(flat$w), tolerance = 1e-10)
  # No information from the data: every direction mixes as slowly as it can.
  none <- EMC2:::.emc_shear_plan(rep(list(NULL), n), list(theta_mu_invar = diag(0, p)), V)
  expect_true(all(none$lam > 1 - 1e-5))
  expect_equal(none$w, rep(1, p), tolerance = 1e-12)
})

test_that("the plan's lambda is the lag-1 autocorrelation of centred Gibbs", {
  # a_s | mu ~ N(Q^-1 (lin_s + V^-1 mu), Q^-1), mu | a ~ N(mean(a), V / n):
  # the Gaussian hierarchy the plan assumes, with V fixed.
  set.seed(3)
  p <- 2; n <- 8; N <- 40000
  V <- matrix(c(1, 0.6, 0.6, 0.8), p)
  H <- matrix(c(0.5, 0.2, 0.2, 0.3), p)
  fits <- lapply(seq_len(n), function(s) list(prec = H, lin = rnorm(p), check = ok_check))
  plan <- EMC2:::.emc_shear_plan(fits, list(theta_mu_invar = diag(0, p)), V)
  Vi <- solve(V); Qi <- solve(H + Vi); Lq <- t(chol(Qi)); Lv <- t(chol(V / n))
  lin <- sapply(fits, `[[`, "lin")
  mu <- numeric(p); draws <- matrix(0, N, p)
  for (t in seq_len(N)) {
    a <- Qi %*% (lin + drop(Vi %*% mu)) + Lq %*% matrix(rnorm(p * n), p)
    mu <- rowMeans(a) + drop(Lv %*% rnorm(p))
    draws[t, ] <- mu
  }
  z <- drop(draws %*% solve(t(plan$dirs)))
  for (i in seq_len(p))
    expect_equal(stats::acf(z[, i], lag.max = 1, plot = FALSE)$acf[2], plan$lam[i], tolerance = 0.02)
})

# A plan as .emc_shear_set_plan() leaves it, with the given expected rounds
# per coordinate.
shear_ready <- function(plan, n, rounds) {
  p <- nrow(plan$dirs)
  zero <- list(E = numeric(2 * p), n = integer(2 * p), moved = integer(2 * p))
  plan$log_step <- c(rep(0, p), rep(-0.5 * log(2 * n), p)) + log(EMC2:::.EMC_SHEAR_STEP0)
  plan$rounds <- rep_len(rounds, 2 * p)
  plan$block <- plan$last <- zero
  plan
}

# A small hierarchy with an inverse-Wishart block (parameters 1-2), an
# unblocked parameter (3), subjects with and without a likelihood fit.
shear_fixture <- function(type) {
  set.seed(7)
  p <- 3; n <- 4
  fits <- lapply(seq_len(n - 1), function(s) {
    A <- matrix(rnorm(p * p), p); H <- crossprod(A) + diag(0.3, p)
    list(prec = H, lin = rnorm(p), check = ok_check)
  })
  fits[n] <- list(NULL)
  sampler <- list(type = type, nuisance = rep(FALSE, p), n_pars = p, n_subjects = n,
                  is_blocked = c(TRUE, TRUE, FALSE), par_group = c(1, 1, 2),
                  prior = list(theta_mu_mean = rep(0, p), theta_mu_invar = diag(0.5, p),
                               v = 2, A = c(0.3, 1, 2), shape = rep(2, p), rate = rep(0.3, p)),
                  lik_laplace = fits)
  V <- if (type == "standard") {
    matrix(c(0.8, 0.3, 0, 0.3, 0.6, 0, 0, 0, 0.5), p)
  } else diag(c(0.8, 0.6, 0.5))
  list(sampler = sampler, ctx = EMC2:::.emc_shear_context(sampler), mu = rnorm(p),
       V = V, alpha = matrix(rnorm(p * n), p),
       u = c(0.3, -0.2, 0.1, 0.25, -0.15, 0.2))
}

test_that("every finite likelihood fit shapes the move, checked or not", {
  f <- shear_fixture("standard")
  expect_identical(vapply(f$ctx$fits, is.null, logical(1)), c(FALSE, FALSE, FALSE, TRUE))
  s <- f$sampler
  s$lik_laplace[[1]]$check <- NULL
  s$lik_laplace[[2]]$check <- list(pareto_k = Inf, outside = 0, laplace_status = "ok")
  s$lik_laplace[[3]]$lin[1] <- NA_real_
  expect_identical(vapply(EMC2:::.emc_shear_context(s)$fits, is.null, logical(1)),
                   c(FALSE, FALSE, TRUE, TRUE))
})

test_that("the shear map is inverted by the opposite step and has the stated Jacobian", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- shear_fixture(type)
    p <- length(f$mu); n <- ncol(f$alpha)
    fits <- f$ctx$fits; prior <- f$sampler$prior; spec <- f$ctx$spec
    state <- function(mu, V) EMC2:::.emc_shear_state(mu, V, fits, prior, spec)
    cur <- state(f$mu, f$V)
    fwd <- EMC2:::.emc_shear_propose(cur, f$u, f$alpha, fits, prior, spec)
    back <- EMC2:::.emc_shear_propose(fwd$state, -f$u, fwd$alpha, fits, prior, spec)
    expect_equal(back$alpha, f$alpha, tolerance = 1e-10)
    expect_equal(back$state$V, f$V, tolerance = 1e-12)
    expect_equal(back$state$mu, f$mu, tolerance = 1e-12)
    expect_equal(back$log_jac, -fwd$log_jac, tolerance = 1e-10)

    # Numerical Jacobian over the free coordinates: mu, the free entries of V
    # (upper triangle within blocks), and the subjects.
    free <- which(upper.tri(f$V, diag = TRUE) & f$V != 0, arr.ind = TRUE)
    pack <- function(mu, V, alpha) c(mu, V[free], alpha)
    unpack <- function(x) {
      V <- matrix(0, p, p); V[free] <- x[p + seq_len(nrow(free))]
      V[free[, 2:1, drop = FALSE]] <- x[p + seq_len(nrow(free))]
      list(mu = x[seq_len(p)], V = V,
           alpha = matrix(x[-seq_len(p + nrow(free))], p))
    }
    map <- function(x, u = f$u) {
      s <- unpack(x)
      out <- EMC2:::.emc_shear_propose(state(s$mu, s$V), u, s$alpha, fits, prior, spec)
      pack(out$state$mu, out$state$V, out$alpha)
    }
    x0 <- pack(f$mu, f$V, f$alpha)
    h <- 1e-6
    jacobian <- function(u) vapply(seq_along(x0), function(i) {
      e <- replace(numeric(length(x0)), i, h)
      (map(x0 + e, u) - map(x0 - e, u)) / (2 * h)
    }, numeric(length(x0)))
    expect_equal(as.numeric(determinant(jacobian(f$u))$modulus), fwd$log_jac, tolerance = 1e-6)

    # A mean step only moves the centres, and its Jacobian is 1.
    u_mu <- replace(f$u, p + seq_len(p), 0)
    fast <- EMC2:::.emc_shear_propose(cur, u_mu, f$alpha, fits, prior, spec)
    full <- state(f$mu + u_mu[seq_len(p)], f$V)
    expect_equal(fast$state$lp, full$lp, tolerance = 1e-12)
    expect_equal(fast$state$fm$m, full$fm$m, tolerance = 1e-12)
    expect_identical(fast$log_jac, 0)
    expect_equal(as.numeric(determinant(jacobian(u_mu))$modulus), 0, tolerance = 1e-6)

    # Maps compose: T_v T_u = T_(u+v), and so do their Jacobians, which is
    # what makes a round's candidates points of one orbit.
    v <- c(-0.1, 0.25, 0.05, -0.2, 0.1, 0.3)
    two <- EMC2:::.emc_shear_propose(fwd$state, v, fwd$alpha, fits, prior, spec)
    one <- EMC2:::.emc_shear_propose(cur, f$u + v, f$alpha, fits, prior, spec)
    expect_equal(two$alpha, one$alpha, tolerance = 1e-10)
    expect_equal(two$state$V, one$state$V, tolerance = 1e-12)
    expect_equal(fwd$log_jac + two$log_jac, one$log_jac, tolerance = 1e-10)
  }
})

test_that("mean-step candidates on the line match the full map", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- shear_fixture(type)
    p <- length(f$mu)
    fits <- f$ctx$fits; prior <- f$sampler$prior; spec <- f$ctx$spec
    cur <- EMC2:::.emc_shear_state(f$mu, f$V, fits, prior, spec)
    d <- c(0.4, -0.3, 0.2)
    t <- c(-1.3, -0.2, 0.7, 2.1)
    line <- EMC2:::.emc_shear_line(cur, d, f$alpha, t, prior)
    for (i in seq_along(t)) {
      full <- EMC2:::.emc_shear_propose(cur, c(t[i] * d, numeric(p)), f$alpha, fits, prior, spec)
      expect_equal(f$alpha + t[i] * line$Delta, full$alpha, tolerance = 1e-12)
      expect_equal(line$delta[i],
                   EMC2:::.emc_shear_log_hier(full$alpha, full$state) + full$state$lp -
                     EMC2:::.emc_shear_log_hier(f$alpha, cur) - cur$lp, tolerance = 1e-10)
    }
  }
})

# Transport walls for the fixture: subject 1 along the unblocked parameter,
# subject 3 along a mixed direction, each half-space holding at the current
# state by the given margin.
shear_walls <- function(f, margin = c(0.4, 0.05)) {
  p <- nrow(f$alpha); n <- ncol(f$alpha)
  U <- matrix(0, p, n); b <- rep(NA_real_, n); Uperp <- vector("list", n)
  has <- c(TRUE, FALSE, TRUE, FALSE)
  U[, 1] <- c(0, 0, 1); U[, 3] <- c(1, 1, 0) / sqrt(2)
  for (s in which(has)) {
    b[s] <- sum(U[, s] * f$alpha[, s]) - margin[match(s, which(has))]
    Uperp[[s]] <- qr.Q(qr(cbind(U[, s], diag(p))))[, -1L, drop = FALSE]
  }
  list(U = U, b = b, Uperp = Uperp, has = has)
}

test_that("the wall map keeps its half-space, inverts, composes and has the stated Jacobian", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- shear_fixture(type)
    p <- length(f$mu); n <- ncol(f$alpha)
    fits <- f$ctx$fits; prior <- f$sampler$prior; spec <- f$ctx$spec
    walls <- shear_walls(f)
    state <- function(mu, V) EMC2:::.emc_shear_state(mu, V, fits, prior, spec)
    propose <- function(cur, u, alpha, w = walls)
      EMC2:::.emc_shear_propose(cur, u, alpha, fits, prior, spec, walls = w)
    cur <- state(f$mu, f$V)
    plain <- propose(cur, f$u, f$alpha, NULL)
    fwd <- propose(cur, f$u, f$alpha)
    expect_equal(fwd$state$mu, plain$state$mu)
    expect_equal(fwd$state$V, plain$state$V)
    expect_equal(fwd$alpha[, !walls$has], plain$alpha[, !walls$has])
    expect_false(isTRUE(all.equal(fwd$alpha[, walls$has], plain$alpha[, walls$has])))

    # Along the orbit the walled subjects stay in their half-space where the
    # plain map leaves it.
    inside <- function(a) vapply(which(walls$has), function(s)
      sum(walls$U[, s] * a[, s]) >= walls$b[s], logical(1))
    plain_out <- FALSE
    for (k in seq(-6, 6, by = 1.5)) {
      expect_true(all(inside(propose(cur, k * f$u, f$alpha)$alpha)))
      plain_out <- plain_out || !all(inside(propose(cur, k * f$u, f$alpha, NULL)$alpha))
    }
    expect_true(plain_out)

    back <- propose(fwd$state, -f$u, fwd$alpha)
    expect_equal(back$alpha, f$alpha, tolerance = 1e-9)
    expect_equal(back$log_jac, -fwd$log_jac, tolerance = 1e-9)
    v <- c(-0.1, 0.25, 0.05, -0.2, 0.1, 0.3)
    two <- propose(fwd$state, v, fwd$alpha)
    one <- propose(cur, f$u + v, f$alpha)
    expect_equal(two$alpha, one$alpha, tolerance = 1e-9)
    expect_equal(fwd$log_jac + two$log_jac, one$log_jac, tolerance = 1e-9)

    free <- which(upper.tri(f$V, diag = TRUE) & f$V != 0, arr.ind = TRUE)
    pack <- function(mu, V, alpha) c(mu, V[free], alpha)
    unpack <- function(x) {
      V <- matrix(0, p, p); V[free] <- x[p + seq_len(nrow(free))]
      V[free[, 2:1, drop = FALSE]] <- x[p + seq_len(nrow(free))]
      list(mu = x[seq_len(p)], V = V, alpha = matrix(x[-seq_len(p + nrow(free))], p))
    }
    map <- function(x, u) {
      s <- unpack(x)
      out <- propose(state(s$mu, s$V), u, s$alpha)
      pack(out$state$mu, out$state$V, out$alpha)
    }
    x0 <- pack(f$mu, f$V, f$alpha)
    h <- 1e-6
    log_det <- function(u) {
      J <- vapply(seq_along(x0), function(i) {
        e <- replace(numeric(length(x0)), i, h)
        (map(x0 + e, u) - map(x0 - e, u)) / (2 * h)
      }, numeric(length(x0)))
      as.numeric(determinant(J)$modulus)
    }
    expect_equal(log_det(f$u), fwd$log_jac, tolerance = 1e-5)
    u_mu <- replace(f$u, p + seq_len(p), 0)
    expect_equal(log_det(u_mu), propose(cur, u_mu, f$alpha)$log_jac, tolerance = 1e-5)
    plain_mu <- propose(cur, u_mu, f$alpha, NULL)

    # A wall far below the state leaves an affine frame map: the plain
    # Jacobian, whitened residuals preserved, and for a mean step the plain
    # map itself.
    far <- walls; far$b[] <- -1e3
    far_map <- propose(cur, f$u, f$alpha, far)
    expect_equal(far_map$log_jac, plain$log_jac, tolerance = 1e-10)
    a2 <- f$alpha + 0.3
    mid <- propose(cur, f$u, (f$alpha + a2) / 2, far)$alpha
    expect_equal(mid, (far_map$alpha + propose(cur, f$u, a2, far)$alpha) / 2, tolerance = 1e-10)
    for (s in which(walls$has)) {
      f0 <- EMC2:::.emc_shear_frame(cur, s); f1 <- EMC2:::.emc_shear_frame(far_map$state, s)
      r0 <- f0$R %*% (f$alpha[, s] - f0$m)
      r1 <- f1$R %*% (far_map$alpha[, s] - f1$m)
      expect_equal(sum(r1^2), sum(r0^2), tolerance = 1e-10)
    }
    expect_equal(propose(cur, u_mu, f$alpha, far)$alpha, plain_mu$alpha, tolerance = 1e-10)
  }
})

test_that("a subject outside its wall stays outside, and the map is a bijection of the line", {
  # The reviewer's counterexample: unit precisions, frames centred at -1 and
  # +1 on the wall's normal, wall at 0, state just below it.
  f0 <- list(m = c(-1, 0), R = diag(2), logdet = 0)
  f1 <- list(m = c(1, 0), R = diag(2), logdet = 0)
  u <- c(1, 0); Uperp <- matrix(c(0, 1), 2, 1)
  map <- EMC2:::.emc_shear_wall_map
  for (a0 in list(c(-0.1, 0), c(-2, 0.3), c(0.1, 0), c(2.5, -1))) {
    fwd <- map(a0, f0, f1, u, 0, Uperp)
    expect_false(is.null(fwd))
    expect_identical(fwd$a[1] >= 0, a0[1] >= 0)
    back <- map(fwd$a, f1, f0, u, 0, Uperp)
    expect_equal(back$a, a0, tolerance = 1e-9)
    expect_equal(back$log_jac, -fwd$log_jac, tolerance = 1e-9)
    # The Jacobian of the wall coordinate's map, numerically.
    h <- 1e-6
    d <- (map(a0 + c(h, 0), f0, f1, u, 0, Uperp)$a[1] - map(a0 - c(h, 0), f0, f1, u, 0, Uperp)$a[1]) / (2 * h)
    expect_equal(log(d), fwd$log_jac, tolerance = 1e-5)
  }
  # The mirror quantile map inverts and composes like the upper one.
  q <- EMC2:::.emc_shear_qmap
  for (z in c(-3, -0.5, 0.2, 2)) {
    m1 <- q(z, 0.1, -0.4); m2 <- q(m1$z, -0.4, 0.1)
    expect_equal(m2$z, z, tolerance = 1e-9)
    expect_equal(m1$log_jac + m2$log_jac, 0, tolerance = 1e-9)
    expect_identical(m1$z >= -0.4, z >= 0.1)
  }
})

test_that("the C++ transport is the R frame map and wall map", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- shear_fixture(type)
    p <- length(f$mu); n <- ncol(f$alpha)
    fits <- f$ctx$fits; prior <- f$sampler$prior; spec <- f$ctx$spec
    walls <- shear_walls(f)
    # subject 3 starts below its wall, which takes the mirror map
    alpha <- f$alpha
    alpha[, 3] <- alpha[, 3] - 1.2 * walls$U[, 3]
    state <- function(mu, V) EMC2:::.emc_shear_state(mu, V, fits, prior, spec)
    cur <- state(f$mu, f$V)
    d <- exp(f$u[p + seq_len(p)])
    new <- state(f$mu + f$u[seq_len(p)], f$V * outer(d, d))
    tr <- EMC2:::.emc_shear_transport(alpha, cur, new, walls)
    expect_true(all(tr$ok))
    for (s in seq_len(n)) {
      f0 <- EMC2:::.emc_shear_frame(cur, s); f1 <- EMC2:::.emc_shear_frame(new, s)
      if (walls$has[s]) {
        wm <- EMC2:::.emc_shear_wall_map(alpha[, s], f0, f1, walls$U[, s], walls$b[s], walls$Uperp[[s]])
        a1 <- wm$a
        expect_equal(tr$jac[s], wm$log_jac - (f0$logdet - f1$logdet), tolerance = 1e-10)
      } else {
        a1 <- f1$m + drop(backsolve(f1$R, f0$R %*% (alpha[, s] - f0$m)))
        expect_identical(tr$jac[s], 0)
      }
      expect_equal(tr$alpha[, s], a1, tolerance = 1e-10)
      expect_equal(tr$r1sq[s], sum((f1$R %*% (a1 - f1$m))^2), tolerance = 1e-10)
    }
    expect_equal(EMC2:::.emc_shear_transport(alpha, cur, new, NULL)$alpha,
                 vapply(seq_len(n), function(s) {
                   f0 <- EMC2:::.emc_shear_frame(cur, s); f1 <- EMC2:::.emc_shear_frame(new, s)
                   f1$m + drop(backsolve(f1$R, f0$R %*% (alpha[, s] - f0$m)))
                 }, numeric(p)), tolerance = 1e-10)
    # A mean step's shifted state has the centres of the state rebuilt there.
    sh <- EMC2:::.emc_shear_shift_state(cur, f$u[seq_len(p)], prior)
    expect_equal(sh$fm$m, state(f$mu + f$u[seq_len(p)], f$V)$fm$m, tolerance = 1e-10)
  }
})

test_that("the full transport inverts and composes from below the wall too", {
  f <- shear_fixture("standard")
  p <- length(f$mu)
  fits <- f$ctx$fits; prior <- f$sampler$prior; spec <- f$ctx$spec
  walls <- shear_walls(f)
  state <- function(mu, V) EMC2:::.emc_shear_state(mu, V, fits, prior, spec)
  propose <- function(cur, u, alpha) EMC2:::.emc_shear_propose(cur, u, alpha, fits, prior, spec, walls = walls)
  cur <- state(f$mu, f$V)
  # Push every walled subject just below its wall.
  alpha <- f$alpha
  for (s in which(walls$has)) alpha[, s] <- alpha[, s] - (sum(walls$U[, s] * alpha[, s]) - walls$b[s] + 0.2) * walls$U[, s]
  below <- function(a) vapply(which(walls$has), function(s) sum(walls$U[, s] * a[, s]) < walls$b[s], logical(1))
  expect_true(all(below(alpha)))
  fwd <- propose(cur, f$u, alpha)
  expect_false(is.null(fwd))
  expect_true(all(below(fwd$alpha)))
  back <- propose(fwd$state, -f$u, fwd$alpha)
  expect_equal(back$alpha, alpha, tolerance = 1e-9)
  expect_equal(back$log_jac, -fwd$log_jac, tolerance = 1e-9)
  v <- c(-0.1, 0.25, 0.05, -0.2, 0.1, 0.3)
  two <- propose(fwd$state, v, fwd$alpha)
  one <- propose(cur, f$u + v, alpha)
  expect_equal(two$alpha, one$alpha, tolerance = 1e-9)
  expect_equal(fwd$log_jac + two$log_jac, one$log_jac, tolerance = 1e-9)
})

test_that("the C++ candidates are the R state, transport and hierarchical terms", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- shear_fixture(type)
    p <- length(f$mu)
    fits <- f$ctx$fits; prior <- f$sampler$prior; spec <- f$ctx$spec
    walls <- shear_walls(f)
    cur <- EMC2:::.emc_shear_state(f$mu, f$V, fits, prior, spec)
    cur_hier <- EMC2:::.emc_shear_log_hier(f$alpha, cur)
    t <- c(-1.3, -0.2, 0.7, 2.1)
    for (w in list(NULL, walls))
      for (v in list(f$u, replace(f$u, p + seq_len(p), 0))) {
        cand <- EMC2:::.emc_shear_candidates(cur, v, t, f$alpha, fits, prior, spec, w)
        for (i in seq_along(t)) {
          u <- t[i] * v
          d <- exp(u[p + seq_len(p)])
          mu1 <- f$mu + u[seq_len(p)]
          new <- EMC2:::.emc_shear_state(mu1, f$V * outer(d, d), fits, prior, spec)
          tr <- if (is.null(new)) NULL else EMC2:::.emc_shear_transport(f$alpha, cur, new, w)
          if (is.null(tr) || !all(tr$ok)) { expect_null(cand[[i]]); next }
          ci <- cand[[i]]
          expect_equal(ci$alpha, tr$alpha, tolerance = 1e-10)
          expect_equal(ci$state$fm, new$fm, tolerance = 1e-10)
          expect_equal(ci$state$lp, new$lp, tolerance = 1e-10)
          expect_equal(ci$hier, EMC2:::.emc_shear_log_hier(tr$alpha, new), tolerance = 1e-8)
          expect_equal(ci$log_jac, sum((spec$size + 1) * u[p + seq_len(p)]) + cur$logdet -
                         new$logdet + sum(tr$jac), tolerance = 1e-8)
        }
      }
    # Without walls a mean step's candidates lie on the line.
    d <- f$u[seq_len(p)]
    line <- EMC2:::.emc_shear_line(cur, d, f$alpha, t, prior)
    cand <- EMC2:::.emc_shear_candidates(cur, c(d, numeric(p)), t, f$alpha, fits, prior, spec)
    for (i in seq_along(t)) {
      expect_equal(cand[[i]]$alpha, f$alpha + t[i] * line$Delta, tolerance = 1e-10)
      expect_equal(cand[[i]]$hier + cand[[i]]$state$lp + cand[[i]]$log_jac - cur_hier - cur$lp,
                   line$delta[i], tolerance = 1e-8)
    }
  }
})

test_that("bounds the transform enforces leave no transport walls", {
  dat <- droplevels(forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ])
  des <- suppressMessages(design(data = dat, model = LBA, matchfun = function(d) d$S == d$lR,
                                 formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1),
                                 constants = c(sv = log(1))))
  emc <- suppressMessages(make_emc(dat, des, n_chains = 1, compress = TRUE))
  support <- EMC2:::.emc_shear_support(emc[[1]])
  expect_true(all(lengths(support) == 0L))
})

test_that("transport walls come from the model's declared bounds and the subject's closest one", {
  # The plain link leaves the declared bounds inside the transform's range.
  dat <- droplevels(forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ])
  des <- suppressMessages(design(data = dat, model = plain_link(LBA),
                                 matchfun = function(d) d$S == d$lR,
                                 formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1),
                                 constants = c(sv = log(1))))
  emc <- suppressMessages(make_emc(dat, des, n_chains = 1, compress = TRUE))
  chain <- emc[[1]]
  pn <- chain$par_names
  p <- length(pn); n <- chain$n_subjects
  support <- EMC2:::.emc_shear_support(chain)
  expect_length(support, n)
  one <- support[[1]]
  on <- vapply(one, function(h) paste(pn[h$u != 0], collapse = "+"), "")
  expect_setequal(on, c("A", "B", "B+B_Eneutral", "B+B_Eaccuracy", "t0"))
  b <- setNames(vapply(one, `[[`, numeric(1), "b"), on)
  expect_equal(b[["t0"]], log(0.05))
  expect_equal(b[["A"]], log(1e-4))
  expect_equal(b[["B+B_Eneutral"]], log(1e-4))
  expect_true(all(vapply(one, function(h) all(h$u[h$u != 0] == 1), logical(1))))

  # Window draws: subject 1 hugs the t0 floor, subject 2 the neutral B
  # threshold, subject 3 sits far from every bound.
  m <- 20
  alpha <- array(rnorm(p * n * m, 0, 0.05), c(p, n, m), dimnames = list(pn, NULL, NULL))
  alpha[c("B", "A"), , ] <- alpha[c("B", "A"), , ] + 0
  alpha["t0", 1, ] <- log(0.05) + 0.1 + alpha["t0", 1, ]
  alpha["B", 2, ] <- log(1e-4) + 0.4 + alpha["B", 2, ]
  alpha["B_Eneutral", 2, ] <- -0.3 + alpha["B_Eneutral", 2, ]
  chain$samples$alpha <- alpha
  walls <- EMC2:::.emc_shear_walls(chain, seq_len(m))
  expect_true(all(walls$has))
  expect_equal(walls$U[, 1], as.numeric(pn == "t0"))
  expect_equal(walls$b[1], log(0.05))
  expect_equal(walls$U[, 2], (pn %in% c("B", "B_Eneutral")) / sqrt(2))
  expect_equal(walls$b[2], log(1e-4) / sqrt(2))
  for (s in seq_len(n)) {
    expect_equal(crossprod(walls$Uperp[[s]]), diag(p - 1))
    expect_equal(drop(crossprod(walls$Uperp[[s]], walls$U[, s])), numeric(p - 1))
  }
  expect_true(EMC2:::.emc_shear_walls_of(c(chain, list(shear = list(walls = walls))))$has[1])
  expect_null(EMC2:::.emc_shear_walls_of(chain))
})

test_that("the wall transport leaves a bounded joint posterior invariant", {
  local_test_rng(kind = "Mersenne-Twister")
  # The diagonal hierarchy of the invariance test below with the first
  # parameter bounded below for every subject; the sampler sees each
  # subject's wall.
  N <- 4000
  h <- diag_hierarchy_draws(N, b = -0.4)
  p <- h$p; n <- h$n; y <- h$y; Si <- h$Si
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...) h$ll(proposals, dadm),
                        .package = "EMC2")
  pn <- c("x1", "x2")
  fits <- list(exact = lapply(seq_len(n), function(s)
                 list(prec = Si, lin = drop(Si %*% y[, s]), check = ok_check)),
               wrong = lapply(seq_len(n), function(s)
                 list(prec = 4 * Si, lin = drop(4 * Si %*% (y[, s] + 0.5)), check = ok_check)))
  # Walls with exact and wrong fits, the plain map on the same target, and a
  # broken control that drops the wall maps' Jacobian, which must fail.
  runs <- list(list(lik = fits$exact, walls = h$walls), list(lik = fits$wrong, walls = h$walls),
               list(lik = fits$exact, walls = NULL),
               list(lik = fits$exact, walls = h$walls, broken = TRUE))
  moved <- numeric(length(runs))
  for (r in seq_along(runs)) {
    run <- runs[[r]]
    plan <- EMC2:::.emc_shear_plan(run$lik, h$prior, diag(exp(2 * colMeans(h$lsd))))
    plan <- shear_ready(plan, n, 2)
    if (isTRUE(run$broken)) {
      candidates <- EMC2:::.emc_shear_candidates
      local_mocked_bindings(.emc_shear_candidates = function(...) {
        lapply(candidates(...), function(out) {
          if (!is.null(out)) out$log_jac <- out$log_jac - out$wall_jac
          out
        })
      }, .package = "EMC2")
    }
    sampler <- list(type = "diagonal-gamma", nuisance = rep(FALSE, p), n_subjects = n,
                    n_pars = p, par_names = pn, data = lapply(seq_len(n), function(s) y[, s]),
                    model = NULL, prior = h$prior, lik_laplace = run$lik, rng = list(),
                    shear = list(plan = plan, walls = run$walls))
    set.seed(5)
    after <- t(vapply(seq_len(N), function(i) {
      V <- diag(exp(2 * h$lsd[i, ]))
      al <- matrix(h$a[i, , ], p, dimnames = list(pn, NULL))
      ll <- vapply(seq_len(n), function(s) -0.5 * sum((al[, s] - y[, s])^2 / h$S), 0)
      out <- EMC2:::.emc_shear_move(sampler, list(tmu = setNames(h$mu[i, ], pn), tvar = V,
                                                  tvinv = solve(V)),
                                    rbind(al, ll), "sample")
      c(out$pars$tmu, 0.5 * log(diag(out$pars$tvar)), out$proposals[1:p, ])
    }, numeric(2 * p + p * n)))
    expect_true(all(after[, 2 * p + seq(1, p * n, by = p)] >= h$b))
    moved[r] <- mean(after[, 1] != h$before[, 1])
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
  expect_gt(moved[1], moved[3])
})

test_that("the shear covariance prior integrates a_half out of the Gibbs step's prior", {
  f <- shear_fixture("standard")
  v <- f$sampler$prior$v; A <- f$sampler$prior$A
  # Inverse-Wishart(v + d - 1, 2 v diag(1 / a)) on the block and a
  # gamma(v / 2, v / a) prior on the unblocked precision, a_k ~ IG(1/2, 1/A_k^2),
  # with each a_k integrated numerically.
  ig <- function(la, k) 0.5 * log(1 / A[k]^2) - 0.5 * log(pi) - 1.5 * la - exp(-la) / A[k]^2 + la
  joint <- function(V) {
    b <- 1:2; d <- 2; df <- v + d - 1; Vb <- V[b, b]; Vbi <- diag(solve(Vb))
    blk <- sum(vapply(b, function(k) log(integrate(function(la)
      exp(df / 2 * log(2 * v / exp(la)) - v * Vbi[k] / exp(la) + ig(la, k)), -30, 30,
      subdivisions = 2000L, rel.tol = 1e-10)$value), 0))
    tau <- 1 / V[3, 3]
    unb <- log(integrate(function(la) exp(dgamma(tau, v / 2, v / exp(la), log = TRUE) + ig(la, 3)),
                         -30, 30, subdivisions = 2000L, rel.tol = 1e-10)$value) + 2 * log(tau)
    -(df + d + 1) / 2 * log(det(Vb)) + blk + unb
  }
  V2 <- f$V * outer(c(1.3, 0.8, 1.1), c(1.3, 0.8, 1.1)); V2[1, 2] <- V2[2, 1] <- 0.1
  V3 <- f$V * outer(c(0.2, 3, 0.05), c(0.2, 3, 0.05))
  lp <- function(V) EMC2:::.emc_shear_log_prior_v(V, f$ctx$spec)
  expect_equal(lp(V2) - lp(f$V), joint(V2) - joint(f$V), tolerance = 1e-8)
  expect_equal(lp(V3) - lp(f$V), joint(V3) - joint(f$V), tolerance = 1e-8)
  g <- shear_fixture("diagonal-gamma")
  indep_g <- function(V) sum(dgamma(1 / diag(V), 2, 0.3, log = TRUE) + 2 * log(1 / diag(V)))
  lpg <- function(V) EMC2:::.emc_shear_log_prior_v(V, g$ctx$spec)
  expect_equal(lpg(diag(c(1, 2, 0.4))) - lpg(g$V), indep_g(diag(c(1, 2, 0.4))) - indep_g(g$V))
})

test_that("a_half is redrawn from the Gibbs step's conditional", {
  f <- shear_fixture("standard")
  set.seed(9)
  draws <- replicate(20000, EMC2:::.emc_shear_draw_a_half(f$V, f$ctx$spec))
  v <- f$sampler$prior$v; A <- f$sampler$prior$A
  prec <- c(diag(solve(f$V[1:2, 1:2])), 1 / f$V[3, 3])
  # 1 / a ~ gamma((v + d) / 2, v prec + 1 / A^2), d = 2, 2, 1.
  expect_equal(rowMeans(1 / draws), c(4, 4, 3) / 2 / (v * prec + 1 / A^2), tolerance = 0.02)
})

test_that("the shear move leaves the joint posterior invariant with the group variance free", {
  # y_s | a_s ~ N(a_s, S), a_s | mu, V ~ N(mu, V), mu ~ N(0, 2 I), and a
  # gamma(2, 0.3) prior on each group precision; S, V and the priors are
  # diagonal, so the posterior factorises over coordinates. Exact draws:
  # (mu_k, log sd_k) from a fine grid of their marginal posterior, then a_s.
  # The data are weak relative to the group spread, the funnel regime.
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
  pn <- c("x1", "x2")
  fits <- list(exact = lapply(seq_len(n), function(s)
                 list(prec = Si, lin = drop(Si %*% y[, s]), check = ok_check)),
               wrong = lapply(seq_len(n), function(s)
                 list(prec = 4 * Si, lin = drop(4 * Si %*% (y[, s] + 0.5)), check = ok_check)))
  prior <- list(theta_mu_mean = rep(m0, p), theta_mu_invar = diag(P0, p),
                shape = rep(shape, p), rate = rep(rate, p))
  # Rounds taken with the plan's lag-1 values as probabilities, then two
  # rounds on every coordinate, each with exact and with wrong fits; then a
  # broken control that drops the scale map's Jacobian, which must fail.
  runs <- list(list(rounds = NULL, lik = fits$exact), list(rounds = NULL, lik = fits$wrong),
               list(rounds = 2, lik = fits$exact), list(rounds = 2, lik = fits$wrong),
               list(rounds = 2, lik = fits$exact, broken = TRUE))
  for (run in runs) {
    lik <- run$lik
    plan <- EMC2:::.emc_shear_plan(lik, prior, diag(exp(2 * colMeans(lsd))))
    plan <- shear_ready(plan, n, if (is.null(run$rounds)) c(plan$lam, plan$w) else run$rounds)
    if (isTRUE(run$broken)) {
      candidates <- EMC2:::.emc_shear_candidates
      local_mocked_bindings(.emc_shear_candidates = function(...) {
        lapply(candidates(...), function(out) {
          if (!is.null(out)) out$log_jac <- 0
          out
        })
      }, .package = "EMC2")
    }
    sampler <- list(type = "diagonal-gamma", nuisance = rep(FALSE, p), n_subjects = n,
                    n_pars = p, par_names = pn, data = lapply(seq_len(n), function(s) y[, s]),
                    model = NULL, prior = prior, lik_laplace = lik, rng = list(),
                    shear = list(plan = plan))
    set.seed(5)
    after <- t(vapply(seq_len(N), function(i) {
      V <- diag(exp(2 * lsd[i, ]))
      al <- matrix(a[i, , ], p, dimnames = list(pn, NULL))
      ll <- vapply(seq_len(n), function(s) -0.5 * sum((al[, s] - y[, s])^2 / S), 0)
      out <- EMC2:::.emc_shear_move(sampler, list(tmu = setNames(mu[i, ], pn), tvar = V,
                                                  tvinv = solve(V)),
                                    rbind(al, ll), "sample")
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

test_that("the shear move leaves the joint posterior invariant with a full standard covariance", {
  skip_model_validation()
  withr::local_seed(83)
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
  V_bar <- Reduce(`+`, lapply(draws, `[[`, "V")) / N
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...) {
    residual <- sweep(proposals, 2, dadm)
    -rowSums(residual * (residual %*% H)) / 2
  }, .package = "EMC2")
  wrong_fit <- lapply(seq_len(n), function(s) list(prec = 2 * H,
    lin = drop(2 * H %*% (y[, s] + 0.4)),
    check = list(pareto_k = 0, outside = 0, laplace_status = "ok")))
  for (fits in list(exact_fit, wrong_fit)) {
    sampler$lik_laplace <- fits
    # Two rounds on every coordinate, from a plan built at the average V.
    plan <- EMC2:::.emc_shear_plan(fits, sampler$prior, V_bar)
    sampler$shear <- list(plan = shear_ready(plan, n, 2))
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

test_that("updates after the shear see the moved group mean as every subject's prior mean", {
  f <- shear_fixture("standard")
  p <- f$sampler$n_pars; n <- f$sampler$n_subjects
  pn <- paste0("x", seq_len(p))
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, ...) -0.5 * rowSums(sweep(proposals, 2, dadm)^2),
    .package = "EMC2")
  plan <- EMC2:::.emc_shear_plan(f$ctx$fits, f$sampler$prior, f$V)
  sampler <- c(f$sampler, list(par_names = pn, model = NULL, rng = list(),
                               data = lapply(seq_len(n), function(s) f$alpha[, s] + 0.3),
                               shear = list(plan = shear_ready(plan, n, 3))))
  al <- matrix(f$alpha, p, dimnames = list(pn, NULL))
  ll <- -0.5 * colSums(0.3^2 + 0 * al)
  pars <- list(tmu = setNames(f$mu, pn), tvar = f$V, tvinv = solve(f$V),
               a_half = setNames(rep(1, p), pn), subj_mu = matrix(f$mu, p, n))
  set.seed(4)
  out <- EMC2:::.emc_shear_move(sampler, pars, rbind(al, ll), "sample")
  expect_false(isTRUE(all.equal(out$pars$tmu, pars$tmu)))
  for (s in seq_len(n))
    expect_equal(as.numeric(EMC2:::get_group_level(out$pars, s, "standard")$mu),
                 as.numeric(out$pars$tmu))
})

test_that("worker evaluation preserves the complete shear transition and RNG", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  withr::local_options(list(emc2.worker_backend = "fork"))
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...) {
    -0.5 * rowSums(sweep(proposals, 2, dadm)^2)
  }, .package = "EMC2")
  for (type in c("standard", "diagonal-gamma")) {
    f <- shear_fixture(type)
    s <- f$sampler
    s$par_names <- paste0("x", seq_len(s$n_pars))
    s$data <- lapply(seq_len(s$n_subjects), function(j) rep(j / 3, s$n_pars))
    plan <- shear_ready(EMC2:::.emc_shear_plan(f$ctx$fits, s$prior, f$V), s$n_subjects, 1)
    s$shear <- list(plan = plan)
    pars <- list(tmu = setNames(f$mu, s$par_names), tvar = f$V,
                 tvinv = solve(f$V), a_half = c(0.7, 1.2, 0.9))
    ll <- vapply(seq_len(s$n_subjects), function(j)
      -0.5 * sum((f$alpha[, j] - s$data[[j]])^2), 0)
    props <- rbind(f$alpha, ll)
    pool <- EMC2:::.emc_wpool_start(2L, s)
    skip_if(is.null(pool), "could not start a pool here")
    tryCatch({
      for (stage in c("burn", "sample")) {
        set.seed(1729)
        serial <- EMC2:::.emc_shear_move(s, pars, props, stage)
        seed_serial <- .Random.seed
        set.seed(1729)
        parallel <- EMC2:::.emc_shear_move(s, pars, props, stage, pool,
                                           list(c(4L, 2L, 1L, 3L), integer()))
        expect_identical(.Random.seed, seed_serial)
        expect_identical(parallel[c("sampler", "pars", "proposals")],
                         serial[c("sampler", "pars", "proposals")])
        expect_true(parallel$wpool$alive)
        expect_gt(serial$sampler$shear$rounds, 0L)
      }
    }, finally = EMC2:::.emc_wpool_stop(pool))
  }
})

test_that("missing worker replies recompute candidates before the MH decision", {
  withr::local_options(list(emc2.failure_policy = "silent"))
  ctx <- list(data = list(0, 1, 2), model = NULL)
  rows <- lapply(c(2, 3, 4), function(x) matrix(x, 1))
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...)
    -(proposals[, 1] - dadm)^2, .package = "EMC2")
  pool <- list(alive = TRUE, n = 2L)
  local_mocked_bindings(
    .emc_wpool_request = function(...) TRUE,
    .emc_wpool_reply = function(pool, w, ...) {
      if (w == 1L) NULL else list(ll = -4, failures = list())
    },
    .emc_wpool_degraded = function(...) NULL,
    .package = "EMC2")
  set.seed(19); seed <- .Random.seed
  out <- EMC2:::.emc_wpool_subject_ll(pool, rows, ctx, list(c(3L, 1L), 2L), source = "shear")
  expect_identical(out$ll, rep(-4, 3))
  expect_false(out$pool$alive)
  expect_identical(.Random.seed, seed)
})

test_that("subject likelihood failures are recorded under the calling move", {
  EMC2:::.emc_reject_reset("collapsed")
  EMC2:::.emc_reject_reset("ensemble")
  # The warning is once per process; another test may have raised it.
  assign("announced", character(), envir = EMC2:::.emc_profile_state)
  withr::local_options(list(emc2.failure_policy = "silent"))
  local_mocked_bindings(calc_ll_manager = function(...) NaN, .package = "EMC2")
  ctx <- list(data = list(NULL, NULL), model = NULL)
  out <- EMC2:::.emc_wpool_subject_ll(NULL, list(matrix(0, 1), matrix(1, 1)), ctx,
                                      source = "collapsed")
  expect_identical(out$ll, c(-Inf, -Inf))
  expect_equal(EMC2:::.emc_reject_counts("collapsed")[["numerical"]], 2L)
  expect_equal(sum(EMC2:::.emc_reject_counts("ensemble")), 0L)
  withr::local_options(list(emc2.failure_policy = "strict"))
  local_mocked_bindings(calc_ll_manager = function(...) stop("broken evaluator"),
                        .package = "EMC2")
  expect_warning(expect_error(EMC2:::.emc_wpool_subject_ll(NULL, list(matrix(0, 1)),
                                          list(data = list(NULL), model = NULL),
                                          source = "shear"),
               "stopped on"), "shear update failed")
})

test_that("subject likelihood partitions cannot omit or duplicate subjects", {
  ctx <- list(data = list(NULL, NULL), model = NULL)
  rows <- list(matrix(0, 1), matrix(0, 1))
  pool <- list(alive = TRUE, n = 2L)
  for (part in list(list(1L, integer()), list(1L, 1L), list(1:2)))
    expect_error(EMC2:::.emc_wpool_subject_ll(pool, rows, ctx, part, source = "shear"),
                 "cover each subject exactly once")
})

test_that("subject likelihood workers preserve results across model families", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  for (model in c("RDM", "LBA", "DDM")) {
    fx <- audit_fixture(model, n_trials = 100L, n_particles = 3L)
    ctx <- list(data = rep(list(fx$dadm), 3), model = fx$model)
    rows <- lapply(seq_len(3), function(j) fx$prop[j, , drop = FALSE])
    serial <- EMC2:::.emc_wpool_subject_ll(NULL, rows, ctx, source = "shear")$ll
    expect_true(all(is.finite(serial)))
    pool <- EMC2:::.emc_wpool_start(2L, ctx)
    skip_if(is.null(pool), "could not start a pool here")
    tryCatch({
      for (k in seq_len(3)) {
        set.seed(21); seed <- .Random.seed
        out <- EMC2:::.emc_wpool_subject_ll(pool, rows, ctx, list(c(3L, 1L), 2L), source = "shear")
        expect_identical(out$ll, serial)
        expect_identical(.Random.seed, seed)
        expect_true(out$pool$alive)
      }
    }, finally = EMC2:::.emc_wpool_stop(pool))
  }
})

test_that("the shared subject queue gives the serial likelihoods and survives a dead worker", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  fx <- audit_fixture("LBA", n_trials = 100L, n_particles = 3L)
  ctx <- list(data = rep(list(fx$dadm), 5), model = fx$model)
  rows <- lapply(seq_len(5), function(j) fx$prop[c(j %% 3L + 1L, 1L), , drop = FALSE])
  serial <- EMC2:::.emc_wpool_subject_ll(NULL, rows, ctx, source = "shear")$ll
  pool <- EMC2:::.emc_wpool_start(3L, ctx)
  skip_if(is.null(pool), "could not start a pool here")
  on.exit(suppressWarnings(try(EMC2:::.emc_wpool_stop(pool), silent = TRUE)), add = TRUE)
  withr::local_options(list(emc2.worker_queue = TRUE))
  queue <- EMC2:::emc_wpool_subject_ll_queue_cpp
  used <- 0L
  local_mocked_bindings(emc_wpool_subject_ll_queue_cpp = function(...) {
    used <<- used + 1L
    queue(...)
  }, .package = "EMC2")
  for (k in seq_len(3)) {
    set.seed(5); seed <- .Random.seed
    out <- EMC2:::.emc_wpool_subject_ll(pool, rows, ctx, source = "shear")
    expect_identical(out$ll, serial)
    expect_identical(.Random.seed, seed)
    expect_true(out$pool$alive)
  }
  expect_identical(used, 3L)
  expect_length(list.files(pool$dir, pattern = "^shared_"), 0L)
  # A worker that dies leaves its subjects to the master.
  pid <- EMC2:::.emc_wpool_job_pid(pool$jobs[[2]])
  tools::pskill(pid, tools::SIGKILL)
  deadline <- proc.time()[["elapsed"]] + 20
  while (EMC2:::.emc_wpool_pid_alive(pid) && proc.time()[["elapsed"]] < deadline) Sys.sleep(0.02)
  out <- suppressWarnings(EMC2:::.emc_wpool_subject_ll(pool, rows, ctx, source = "shear"))
  expect_identical(out$ll, serial)
  expect_false(out$pool$alive)
})

test_that("a fit with shear and subject HMC is identical serially, statically, and queued", {
  skip_on_os("windows")
  skip_on_cran()
  force_allocation()
  dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- suppressMessages(design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1)))
  fit_on <- function(cores, queue) {
    withr::local_options(list(emc2.worker_queue = queue))
    RNGkind("L'Ecuyer-CMRG"); set.seed(23)
    emc <- suppressMessages(make_emc(dat, des, n_chains = 1, compress = TRUE))
    for (stage in c("preburn", "burn"))
      emc <- suppressMessages(run_emc(emc, stage = stage, cores_for_chains = 1,
        cores_per_chain = cores, step_size = 6, max_tries = 1, verbose = FALSE,
        stop_criteria = list(iter = 12, max_gd = Inf, min_unique = 0, min_es = 0)))
    emc[[1]]
  }
  serial <- fit_on(1, FALSE)
  # The identity checks below only mean something if both moves ran.
  expect_gt(sum(serial$shear$plan$block$moved), 0)
  expect_gt(serial$subject_hmc$totals$attempts, 0)
  expect_true(any(serial$alloc$particle$rate < 1))
  for (parallel in list(fit_on(2, FALSE), fit_on(3, TRUE))) {
    expect_identical(parallel$samples$alpha, serial$samples$alpha)
    expect_identical(parallel$samples$theta_mu, serial$samples$theta_mu)
    expect_identical(parallel$samples$theta_var, serial$samples$theta_var)
    expect_identical(parallel$shear$plan[c("rounds", "log_step", "block")],
                     serial$shear$plan[c("rounds", "log_step", "block")])
    expect_identical(parallel$rng, serial$rng)
  }
})

test_that("rebuilds carry step scales to matching directions and record the window", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  withr::local_seed(95)
  withr::local_options(list(emc2.worker_backend = "fork"))
  local_mocked_bindings(calc_ll_manager = function(proposals, ...) {
    -rowSums(proposals^2) / 2
  }, .package = "EMC2")
  subjects <- levels(forstmann$subjects)[1:3]
  dat <- droplevels(forstmann[forstmann$subjects %in% subjects, ])
  des <- suppressMessages(design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1)))
  emc <- suppressMessages(make_emc(dat, des, n_chains = 2L, type = "standard"))
  for (stage in c("preburn", "burn", "adapt"))
    emc <- suppressMessages(run_emc(emc, stage = stage, step_size = 12L, max_tries = 1L,
      verbose = FALSE, trim = FALSE, cores_for_chains = 1L, cores_per_chain = 1L,
      stop_criteria = list(iter = 24L, max_gd = Inf, selection = c("alpha", "mu"),
                           min_unique = 0, min_es = 0)))
  emc <- EMC2:::restore_duplicates(emc)
  for (ch in emc) {
    fits <- EMC2:::.emc_shear_fits(ch)
    expect_true(all(!vapply(fits, is.null, logical(1))))
    # Without fits every scale weight is exactly 1; a subject with a fit
    # carries likelihood information, so its prior share is below 1.
    expect_true(all(ch$shear$plan$w < 1))
    pl <- ch$shear$plan
    expect_true(all(is.finite(pl$rounds) & pl$rounds >= 0 &
                      pl$rounds <= EMC2:::.EMC_ALLOC_MAX_ROUNDS))
    expect_true(all(pl$log_step >= EMC2:::.EMC_SHEAR_LOG_STEP_RANGE[1] &
                      pl$log_step <= EMC2:::.EMC_SHEAR_LOG_STEP_RANGE[2]))
    # Every rebuild after the first measured a recorded window.
    h <- ch$alloc$history
    expect_gt(length(h), 2L)
    expect_false(h[[length(h)]]$stored)
    expect_gt(h[[length(h)]]$pairs[["group"]], 0)
  }
  # A rebuilt plan whose directions match keeps their scales.
  ch <- emc[[1]]
  old <- ch$shear$plan
  old$log_step[] <- log(1.7)
  ch$shear$plan <- old
  again <- EMC2:::.emc_shear_set_plan(ch, seq_len(ch$samples$idx), "burn")$shear$plan
  p <- ncol(old$dirs)
  S <- stats::cov(t(ch$samples$theta_mu[, seq_len(ch$samples$idx), drop = FALSE]))
  cw <- EMC2:::.emc_shear_carry_weights(old$dirs, again$dirs, S)
  for (k in seq_len(p)) {
    expected <- sum(cw$w[k, seq_len(p)] * (log(1.7) + log(cw$scale[k, seq_len(p)])))
    expect_equal(again$log_step[k], expected, tolerance = 1e-10)
  }
  # A log group sd always matches itself.
  expect_equal(again$log_step[p + seq_len(p)], rep(log(1.7), p))
  # Production blocks pool their rounds; the allocation reads the pool, while
  # step scales follow the last block. A log group sd's record carries over
  # unscaled.
  sd_k <- p + seq_len(p)
  rebuild <- function(ch, stage, n, moved, E) {
    pl <- ch$shear$plan
    pl$block <- list(E = rep(E, 2 * p), n = rep(n, 2 * p), moved = rep(moved, 2 * p))
    ch$shear$plan <- pl
    EMC2:::.emc_shear_set_plan(ch, seq_len(ch$samples$idx), stage)
  }
  ch$shear$plan <- again
  ch <- rebuild(ch, "sample", 12L, 0L, 1)
  expect_equal(ch$shear$plan$last$n[sd_k], rep(12L, p))
  expect_equal(ch$shear$plan$pool$n[sd_k], integer(p))
  ch <- rebuild(ch, "sample", 10L, 10L, 2)
  step <- ch$shear$plan$log_step[sd_k]
  ch <- rebuild(ch, "sample", 30L, 15L, 3)
  pl <- ch$shear$plan
  # The pool carries the earlier blocks discounted at each rebuild.
  d <- EMC2:::.EMC_ALLOC_EVIDENCE_DECAY
  expect_equal(pl$last$n[sd_k], rep(d * 10 + 30, p))
  expect_equal(pl$last$E[sd_k], rep(d * 2 + 3, p))
  expect_equal(pl$last$moved[sd_k], rep(d * 10 + 15, p))
  shift <- max(-log(2), min(log(2), 2 * (0.5 - EMC2:::.EMC_SHEAR_MOVE_TARGET)))
  expect_equal(pl$log_step[sd_k],
               pmin(EMC2:::.EMC_SHEAR_LOG_STEP_RANGE[2],
                    pmax(EMC2:::.EMC_SHEAR_LOG_STEP_RANGE[1], step + shift)))
  # Outside production the last block alone measures the rounds.
  n_burn1 <- d * (d * 10 + 30) + 7
  n_burn2 <- d * n_burn1 + 2
  ch <- rebuild(ch, "burn", 7L, 1L, 4)
  expect_equal(ch$shear$plan$last$n[sd_k], rep(n_burn1, p))
  expect_equal(ch$shear$plan$pool$n[sd_k], rep(n_burn1, p))
  ch <- rebuild(ch, "burn", 2L, 1L, 1)
  expect_equal(ch$shear$plan$last$n[sd_k], rep(n_burn2, p))
  expect_equal(ch$shear$plan$pool$n[sd_k], rep(n_burn2, p))
})

test_that("an orbit's width is measured whatever scale the grid starts at", {
  sigma <- 0.37
  local_mocked_bindings(
    .emc_shear_orbit_terms = function(ch, idx, k, grid, fits_of, n_cores = 1L)
      list(per = matrix(0, length(grid), 2), lp = -0.5 * ((grid - 0.1) / sigma)^2),
    .package = "EMC2")
  for (s0 in sigma * c(1 / 20, 1, 20)) {
    w <- EMC2:::.emc_shear_orbit_width(NULL, 1:10, 1L, NULL, s0)
    expect_equal(w, sigma, tolerance = 0.05)
  }
})

test_that("shear frames reject approximate solves instead of dropping conditional directions", {
  H <- array(0, c(2, 2, 1))
  lin <- matrix(0, 2, 1)
  mu <- c(2, -3)
  Vi <- diag(c(1, exp(96)))
  # An exact diagonal system, but below the numerical condition threshold:
  # Armadillo's approximate fallback used to lose the first centre here.
  out <- EMC2:::emc_shear_frames_cpp(H, lin, Vi, drop(Vi %*% mu))
  expect_false(out$ok)
  Vi <- matrix(c(2, 0.3, 0.3, 4), 2)
  H[, , 1] <- diag(c(3, 1))
  lin[, 1] <- c(0.5, -0.7)
  out <- EMC2:::emc_shear_frames_cpp(H, lin, Vi, drop(Vi %*% mu))
  expect_true(out$ok)
  expect_equal(drop(out$m), drop(solve(H[, , 1] + Vi, lin[, 1] + Vi %*% mu)),
               tolerance = 1e-12)
  expect_equal(crossprod(out$R[, , 1]), H[, , 1] + Vi, tolerance = 1e-12)
  lin[1, 1] <- Inf
  expect_false(EMC2:::emc_shear_frames_cpp(H, lin, Vi, drop(Vi %*% mu))$ok)
})

test_that("orbit tuning retries numerical failures on a smaller grid", {
  grids <- list()
  local_mocked_bindings(
    .emc_shear_orbit_terms = function(ch, idx, k, grid, fits_of, n_cores = 1L) {
      grids[[length(grids) + 1L]] <<- grid
      f <- if (max(abs(grid)) > 3) rep(NA_real_, length(grid)) else -0.5 * grid^2
      list(per = matrix(0, length(grid), 2), lp = f)
    }, .package = "EMC2")
  expect_equal(EMC2:::.emc_shear_orbit_width(NULL, 1:10, 1L, NULL, 4),
               1, tolerance = 0.05)
  expect_equal(max(abs(grids[[2]])), max(abs(grids[[1]])) / 4)
})

test_that("unresolved orbit widths do not replace the existing tuning", {
  grids <- list()
  local_mocked_bindings(
    .emc_shear_orbit_terms = function(ch, idx, k, grid, fits_of, n_cores = 1L) {
      grids[[length(grids) + 1L]] <<- grid
      list(per = matrix(0, length(grid), 2), lp = -1e-30 * grid^2)
    }, .package = "EMC2")
  expect_true(is.na(EMC2:::.emc_shear_orbit_width(NULL, 1:10, 1L, NULL, 1)))
  extent <- vapply(grids, function(g) max(abs(g)), numeric(1))
  expect_true(all(extent[-1] / head(extent, -1) <= 4))
})

test_that("the native shear rounds are the R rounds on the same random stream", {
  skip_on_cran()
  # The R step is the one that runs the move through the mocked binding.
  old <- options(emc2.fit_host = FALSE, emc2.native_steps = FALSE)
  on.exit(options(old), add = TRUE)
  des <- design(factors = list(subjects = 1:5, S = 1), Rlevels = 1, model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  set.seed(31)
  pars <- cbind(mu = log(0.4) + rnorm(5, 0, 0.15), sigma = log(0.05) + rnorm(5, 0, 0.3),
                tau = log(0.2) + rnorm(5, 0, 0.3))
  rownames(pars) <- 1:5
  dat <- make_data(pars, des, n_trials = 60)
  emc <- make_emc(dat, des, rt_resolution = 0.001, n_chains = 2)
  move_r <- EMC2:::.emc_shear_move
  native <- moved <- differing <- 0L
  strip <- function(x) x[setdiff(names(x), c("round_el", "build_el", "ll_el", "wpool"))]
  both <- function(sampler, pars, proposals, stage, ...) {
    seed <- get(".Random.seed", envir = globalenv())
    options(emc2.native_shear_round = FALSE)
    a <- move_r(sampler, pars, proposals, stage, ...)
    seed_a <- get(".Random.seed", envir = globalenv())
    assign(".Random.seed", seed, envir = globalenv())
    options(emc2.native_shear_round = TRUE)
    b <- move_r(sampler, pars, proposals, stage, ...)
    if (!is.null(b$rounds)) native <<- native + b$rounds
    moved <<- moved + !identical(b$pars$tvar, pars$tvar)
    if (!identical(strip(a), strip(b)) ||
        !identical(seed_a, get(".Random.seed", envir = globalenv())))
      differing <<- differing + 1L
    b
  }
  local_mocked_bindings(.emc_shear_move = both, .package = "EMC2")
  fitted <- suppressWarnings(fit(emc, verbose = FALSE, cores_for_chains = 1,
                                 stop_criteria = list(preburn = list(iter = 20),
                                                      burn = list(iter = 60),
                                                      adapt = list(iter = 60),
                                                      sample = list(iter = 100))))
  options(emc2.native_shear_round = NULL)
  # The comparison covered rounds that moved the state.
  expect_gt(native, 20L)
  expect_gt(moved, 10L)
  expect_identical(differing, 0L)
})

test_that("rebuilt mean directions inherit the old directions' records by their share", {
  set.seed(3)
  p <- 3L
  D <- matrix(rnorm(p * p), p)
  sd_k <- p + seq_len(p)
  # Same directions: every record is its own.
  cw <- EMC2:::.emc_shear_carry_weights(D, D)
  expect_equal(cw$w, diag(2L * p))
  expect_equal(cw$scale[cbind(1:p, 1:p)], rep(1, p))
  # A direction rescaled by a: its own record, moves in units of the new length.
  a <- 2
  cw <- EMC2:::.emc_shear_carry_weights(D, D %*% diag(c(a, 1, 1)))
  expect_equal(cw$w, diag(2L * p))
  expect_equal(cw$scale[1, 1], 1 / a)
  # Two directions of equal length rotated into each other: each new one
  # is half of each old one, the third keeps its own and the log sds theirs.
  D[, 2] <- D[, 2] * sqrt(sum(D[, 1]^2) / sum(D[, 2]^2))
  th <- pi / 4
  Q <- diag(p)
  Q[1:2, 1:2] <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  cw <- EMC2:::.emc_shear_carry_weights(D, D %*% Q)
  expect_equal(cw$w[1:2, 1:2], matrix(0.5, 2, 2))
  expect_equal(cw$w[3, ], replace(numeric(2L * p), 3, 1))
  expect_equal(cw$w[sd_k, sd_k], diag(p))
  # The rotated pair keeps the pair's evidence through a rebuild's carry.
  n <- c(300, 100, 50, rep(10, p))
  expect_equal(drop(cw$w %*% n)[1:2], c(200, 200))
  expect_equal(sum(drop(cw$w %*% n)[1:2]), sum(n[1:2]))
  # An old basis that cannot be inverted passes nothing to the mean directions.
  cw <- EMC2:::.emc_shear_carry_weights(matrix(0, p, p), D)
  expect_equal(rowSums(cw$w)[1:p], rep(0, p))
  expect_equal(cw$w[sd_k, sd_k], diag(p))
})

test_that("a carried record keeps its per-round factor in any new basis", {
  set.seed(5)
  p <- 4L
  A <- matrix(rnorm(p * p), p)
  S <- crossprod(A) + diag(p)
  var_in <- function(D) { Di <- solve(D); rowSums((Di %*% S) * Di) }
  # Every old coordinate moves by the same share of its variance per round.
  g0 <- 0.5
  D <- matrix(rnorm(p * p), p)
  E <- c(2 * (1 - g0) * var_in(D) * 100, numeric(p))
  n <- c(rep(100, p), numeric(p))
  g_new <- function(Dn) {
    cw <- EMC2:::.emc_shear_carry_weights(D, Dn, S)
    1 - drop((cw$w * cw$scale^2) %*% E)[1:p] / drop(cw$w %*% n)[1:p] / (2 * var_in(Dn))
  }
  # A fresh basis of unrelated lengths, and the same with one coordinate
  # pinned to a direction along which mu covaries with it.
  fresh <- matrix(rnorm(p * p), p) %*% diag(c(0.1, 1, 5, 20))
  expect_equal(g_new(fresh), rep(g0, p))
  row <- solve(D)[2, ]
  pinned <- EMC2:::.emc_shear_pin_dir(fresh, 2L, drop(S %*% row), row)
  expect_false(is.null(pinned))
  expect_equal(g_new(pinned), rep(g0, p))
  # The pinned coordinate is old coordinate 2 itself: its own record, unscaled.
  cw <- EMC2:::.emc_shear_carry_weights(D, pinned, S)
  expect_equal(cw$scale[2, 2], 1)
})

# An arbitrary curve for n subjects of p parameters: K knots of row' mu around
# `centre`, with unrelated positive-definite precisions and linear terms at
# each, so that nothing but the map's own construction makes a move valid.
arbitrary_curve <- function(p, n, row, centre, width, K = 3L, seed = 13) {
  set.seed(seed)
  cube <- array(0, c(p, p, n, K))
  lin <- array(rnorm(p * n * K), c(p, n, K))
  for (s in seq_len(n)) for (q in seq_len(K)) {
    A <- matrix(rnorm(p * p), p)
    cube[, , s, q] <- crossprod(A) * q / 2 + diag(0.2 * q, p)
  }
  list(row = row, knots = centre + width * seq(-1, 1, length.out = K), cube = cube, lin = lin)
}

test_that("a curve's fits are the knots' at the knots, between them in between, and held beyond", {
  p <- 3L; n <- 2L
  cv <- arbitrary_curve(p, n, row = c(1, -0.5, 0.25), centre = 0.3, width = 1)
  at <- function(c) EMC2:::emc_shear_curve_fits_cpp(c(c, 0, 0), cv, n)
  for (q in 1:3) {
    f <- at(cv$knots[q])
    expect_equal(f$cube, cv$cube[, , , q], tolerance = 1e-14)
    expect_equal(f$lin, cv$lin[, , q], tolerance = 1e-14)
  }
  mid <- at(cv$knots[2] + 0.25 * diff(cv$knots[2:3]))
  expect_equal(mid$cube, 0.75 * cv$cube[, , , 2] + 0.25 * cv$cube[, , , 3], tolerance = 1e-13)
  expect_equal(mid$lin, 0.75 * cv$lin[, , 2] + 0.25 * cv$lin[, , 3], tolerance = 1e-13)
  expect_equal(at(cv$knots[1] - 5)$cube, cv$cube[, , , 1])
  expect_equal(at(cv$knots[3] + 5)$lin, cv$lin[, , 3])
  # The coordinate is row' mu, whatever the rest of mu.
  expect_equal(EMC2:::emc_shear_curve_fits_cpp(c(0.3, 2, 4), cv, n), at(0.3), tolerance = 1e-13)
  expect_error(EMC2:::emc_shear_curve_fits_cpp(c(0, 0, 0), replace(cv, "knots", list(c(1, 1, 2))), n))
})

test_that("the shear map with curved frames inverts, composes and has the stated Jacobian", {
  for (type in c("standard", "diagonal-gamma")) {
    f <- shear_fixture(type)
    p <- length(f$mu); n <- ncol(f$alpha)
    prior <- f$sampler$prior; spec <- f$ctx$spec
    fits <- EMC2:::.emc_shear_pack(f$ctx$fits, p)
    row <- c(0.7, -0.4, 0.2)
    fits$curve <- arbitrary_curve(p, n, row, centre = sum(row * f$mu) + 0.05, width = 0.6)
    state <- function(mu, V) EMC2:::.emc_shear_state(mu, V, fits, prior, spec)
    cur <- state(f$mu, f$V)
    # The step moves the curve's coordinate, so the frames differ at its ends.
    expect_true(EMC2:::.emc_shear_curve_moves(fits, f$u[seq_len(p)]))
    fwd <- EMC2:::.emc_shear_propose(cur, f$u, f$alpha, fits, prior, spec)
    back <- EMC2:::.emc_shear_propose(fwd$state, -f$u, fwd$alpha, fits, prior, spec)
    expect_equal(back$alpha, f$alpha, tolerance = 1e-9)
    expect_equal(back$state$mu, f$mu, tolerance = 1e-12)
    expect_equal(back$log_jac, -fwd$log_jac, tolerance = 1e-9)
    free <- which(upper.tri(f$V, diag = TRUE) & f$V != 0, arr.ind = TRUE)
    pack <- function(mu, V, alpha) c(mu, V[free], alpha)
    unpack <- function(x) {
      V <- matrix(0, p, p); V[free] <- x[p + seq_len(nrow(free))]
      V[free[, 2:1, drop = FALSE]] <- x[p + seq_len(nrow(free))]
      list(mu = x[seq_len(p)], V = V, alpha = matrix(x[-seq_len(p + nrow(free))], p))
    }
    map <- function(x, u) {
      s <- unpack(x)
      out <- EMC2:::.emc_shear_propose(state(s$mu, s$V), u, s$alpha, fits, prior, spec)
      pack(out$state$mu, out$state$V, out$alpha)
    }
    x0 <- pack(f$mu, f$V, f$alpha)
    h <- 1e-6
    jacobian <- function(u) vapply(seq_along(x0), function(i) {
      e <- replace(numeric(length(x0)), i, h)
      (map(x0 + e, u) - map(x0 - e, u)) / (2 * h)
    }, numeric(length(x0)))
    expect_equal(as.numeric(determinant(jacobian(f$u))$modulus), fwd$log_jac, tolerance = 1e-5)
    # A mean step along the curve changes the frames: its Jacobian is not 1,
    # and it is the frames' determinant ratio.
    u_mu <- replace(f$u, p + seq_len(p), 0)
    mean_step <- EMC2:::.emc_shear_propose(cur, u_mu, f$alpha, fits, prior, spec)
    expect_gt(abs(mean_step$log_jac), 1e-3)
    expect_equal(mean_step$log_jac, cur$logdet - mean_step$state$logdet, tolerance = 1e-10)
    expect_equal(as.numeric(determinant(jacobian(u_mu))$modulus), mean_step$log_jac,
                 tolerance = 1e-5)
    full <- state(f$mu + u_mu[seq_len(p)], f$V)
    expect_equal(mean_step$state$fm$m, full$fm$m, tolerance = 1e-10)
    # A mean step that keeps the coordinate is the plain shift at the
    # current fits.
    flat <- c(0.4, 0.7, 0)
    expect_false(EMC2:::.emc_shear_curve_moves(fits, flat))
    flat_step <- EMC2:::.emc_shear_propose(cur, c(flat, 0, 0, 0), f$alpha, fits, prior, spec)
    expect_identical(flat_step$log_jac, 0)
    expect_equal(flat_step$state$fm$m, state(f$mu + flat, f$V)$fm$m, tolerance = 1e-10)
    # Maps compose across the knots.
    v <- c(-0.3, 0.25, 0.5, -0.2, 0.1, 0.3)
    two <- EMC2:::.emc_shear_propose(fwd$state, v, fwd$alpha, fits, prior, spec)
    one <- EMC2:::.emc_shear_propose(cur, f$u + v, f$alpha, fits, prior, spec)
    expect_equal(two$alpha, one$alpha, tolerance = 1e-9)
    expect_equal(fwd$log_jac + two$log_jac, one$log_jac, tolerance = 1e-9)
  }
})

test_that("the shear move with curved frames leaves the joint posterior invariant", {
  local_test_rng(kind = "Mersenne-Twister")
  N <- 4000
  h <- diag_hierarchy_draws(N)
  p <- h$p; n <- h$n; y <- h$y; Si <- h$Si
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, ...) h$ll(proposals, dadm),
                        .package = "EMC2")
  pn <- c("x1", "x2")
  lik <- lapply(seq_len(n), function(s)
    list(prec = Si, lin = drop(Si %*% y[, s]), check = ok_check))
  plan <- EMC2:::.emc_shear_plan(lik, h$prior, diag(exp(2 * colMeans(h$lsd))))
  plan <- shear_ready(plan, n, 2)
  row <- solve(plan$dirs)[1, ]
  centre <- mean(h$mu %*% row); width <- 2 * sd(h$mu %*% row)
  # The exact fits bent along the first mean direction, an unrelated curve,
  # and a broken control that drops the frames' Jacobian, which must fail.
  bent <- arbitrary_curve(p, n, row, centre, width)
  for (s in seq_len(n)) for (q in 1:3) {
    bent$cube[, , s, q] <- Si * c(0.4, 1, 2.5)[q]
    bent$lin[, s, q] <- drop(bent$cube[, , s, q] %*% (y[, s] + c(-0.6, 0, 0.8)[q]))
  }
  runs <- list(list(curve = bent), list(curve = arbitrary_curve(p, n, row, centre, width)),
               list(curve = bent, broken = TRUE))
  moved <- numeric(length(runs))
  for (r in seq_along(runs)) {
    run <- runs[[r]]
    if (isTRUE(run$broken)) {
      candidates <- EMC2:::.emc_shear_candidates
      local_mocked_bindings(.emc_shear_candidates = function(...) {
        lapply(candidates(...), function(out) {
          if (!is.null(out)) out$log_jac <- 0
          out
        })
      }, .package = "EMC2")
    }
    sampler <- list(type = "diagonal-gamma", nuisance = rep(FALSE, p), n_subjects = n,
                    n_pars = p, par_names = pn, data = lapply(seq_len(n), function(s) y[, s]),
                    model = NULL, prior = h$prior, lik_laplace = lik, rng = list(),
                    shear = list(plan = plan, curve = run$curve))
    expect_false(is.null(EMC2:::.emc_shear_context(sampler)$packed$curve))
    set.seed(5)
    after <- t(vapply(seq_len(N), function(i) {
      V <- diag(exp(2 * h$lsd[i, ]))
      al <- matrix(h$a[i, , ], p, dimnames = list(pn, NULL))
      ll <- vapply(seq_len(n), function(s) -0.5 * sum((al[, s] - y[, s])^2 / h$S), 0)
      out <- EMC2:::.emc_shear_move(sampler, list(tmu = setNames(h$mu[i, ], pn), tvar = V,
                                                  tvinv = solve(V)),
                                    rbind(al, ll), "sample")
      c(out$pars$tmu, 0.5 * log(diag(out$pars$tvar)), out$proposals[1:p, ])
    }, numeric(2 * p + p * n)))
    moved[r] <- mean(after[, 1] != h$before[, 1])
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
  expect_true(all(moved[1:2] > 0.1))
})

test_that("a window's integrated autocorrelation time pools the chains about their own means", {
  set.seed(3)
  ar1 <- function(rho, n) as.numeric(stats::filter(rnorm(n), rho, "recursive"))
  for (rho in c(0.5, 0.8, 0.9)) {
    xs <- lapply(c(0, 5, -3), function(shift) shift + ar1(rho, 20000))
    expect_equal(EMC2:::.emc_shear_window_iact(xs), (1 + rho) / (1 - rho), tolerance = 0.1)
  }
  # White noise has none to speak of, and a coordinate slower than the summed
  # lags reads as at least as slow as they resolve.
  expect_lt(EMC2:::.emc_shear_window_iact(list(rnorm(5000), rnorm(5000))), 1.3)
  expect_gt(EMC2:::.emc_shear_window_iact(list(ar1(0.999, 5000))), 200)
  expect_true(is.na(EMC2:::.emc_shear_window_iact(list(rep(1, 100), rep(2, 100)))))
  expect_true(is.na(EMC2:::.emc_shear_window_iact(list(c(1, 2)))))
})

test_that("a curve is asked for where the window shows a slow coordinate, whatever its rounds", {
  set.seed(3)
  ar1 <- function(rho, n = 400) as.numeric(stats::filter(rnorm(n), rho, "recursive"))
  cs_all <- lapply(1:3, function(j) rbind(ar1(0.97), ar1(0.5), ar1(0.9)))
  out <- EMC2:::.emc_shear_curve_slow(cs_all)
  # 1 and 3 take tens of iterations per independent draw; 2 takes about three.
  expect_identical(out$slow, c(TRUE, FALSE, TRUE))
  expect_gt(out$tau[1], out$tau[3])
  expect_lt(out$tau[2], EMC2:::.EMC_SHEAR_CURVE_TAU)
  # A window too short to measure, or one without spread, asks for nothing.
  expect_false(any(EMC2:::.emc_shear_curve_slow(lapply(cs_all, function(x) x[, 1:20]))$slow))
  expect_false(any(EMC2:::.emc_shear_curve_slow(lapply(cs_all, function(x) x[, integer(0), drop = FALSE]))$slow))
  flat <- lapply(cs_all, function(x) rbind(rep(1, ncol(x)), x[2:3, ]))
  expect_identical(EMC2:::.emc_shear_curve_slow(flat)$slow, c(FALSE, FALSE, TRUE))
})

test_that("a round at the step target refreshes the stated share of a Gaussian orbit", {
  set.seed(5)
  step <- EMC2:::.EMC_SHEAR_STEP0
  N <- EMC2:::.EMC_SHEAR_CANDIDATES
  M <- 30000L
  x <- numeric(M + 1L)
  for (i in seq_len(M)) {
    cand <- c(x[i], x[i] + rnorm(1, 0, step) + rnorm(N, 0, step))
    lw <- stats::dnorm(cand, log = TRUE)
    x[i + 1L] <- cand[sample.int(N + 1L, 1L, prob = exp(lw - max(lw)))]
  }
  expect_equal(1 - stats::cor(x[-1L], x[-(M + 1L)]), EMC2:::.EMC_SHEAR_ROUND_REFRESH,
               tolerance = 0.03)
})

test_that("a curve's worth is the fixed-frame rounds its round replaces", {
  eta <- EMC2:::.EMC_SHEAR_ROUND_REFRESH
  factor <- EMC2:::.emc_shear_round_factor
  worth <- EMC2:::.emc_shear_curve_worth
  expect_equal(factor(c(0, 0.5, 1, 1.3)), c(1, 1 - eta / 4, 1 - eta, 1 - eta))
  # Two fixed-frame rounds take off the square of one's factor.
  r_two <- sqrt((1 - factor(0.5)^2) / eta)
  expect_equal(worth(0.5, r_two), 2)
  expect_equal(worth(0.6, 0.6), 1)
  expect_gt(worth(0.4, 0.7), worth(0.4, 0.6))
  # Narrow orbits: the ratio of the shares of variance refreshed.
  expect_equal(worth(0.05, 0.1), 4, tolerance = 0.01)
  # An orbit that spans most of the margin leaves a curve nothing worth
  # having, however wide the curve's orbit is measured.
  expect_lt(worth(0.91, 0.99), 1.5)
  expect_lt(worth(0.8, 1), EMC2:::.EMC_SHEAR_CURVE_WORTH)
  expect_lt(worth(0.8, 1.2), EMC2:::.EMC_SHEAR_CURVE_WORTH)
  expect_gt(worth(0.75, 1), EMC2:::.EMC_SHEAR_CURVE_WORTH)
  expect_identical(worth(c(0.4, NA, 0, 0.4), c(0.4, 0.5, 0.5, NA)), c(1, NA, NA, NA))
  expect_equal(worth(c(0.4, 0.6), 1), c(worth(0.4, 1), worth(0.6, 1)))
})

test_that("a coordinate whose curve failed waits until the evidence has changed", {
  waits <- EMC2:::.emc_shear_curve_waits
  rec <- list(row = c(1, 0), at = 400L, len = 300L, stage = "sample", fails = 1L)
  expect_true(waits(rec, 400L, "sample"))
  expect_true(waits(rec, 699L, "sample"))
  expect_false(waits(rec, 700L, "sample"))
  # Each further failure doubles the draws it waits for.
  expect_true(waits(replace(rec, "fails", 2L), 999L, "sample"))
  expect_false(waits(replace(rec, "fails", 2L), 1000L, "sample"))
  expect_true(waits(replace(rec, "fails", 3L), 1599L, "sample"))
  # Another stage, or stored draws that no longer reach the failure, are new
  # evidence.
  expect_false(waits(rec, 450L, "adapt"))
  expect_false(waits(rec, 100L, "sample"))

  note <- EMC2:::.emc_shear_curve_note
  recs <- note(list(), c(1, 0), TRUE, 400L, 300L, "sample")
  expect_identical(recs, list(rec))
  # The same coordinate, slightly turned, counts its failures in a row.
  again <- note(recs, c(1, 0.1), TRUE, 900L, 800L, "sample")
  expect_length(again, 1L)
  expect_identical(again[[1]][c("at", "len", "fails")], list(at = 900L, len = 800L, fails = 2L))
  # Another coordinate has its own record; a curve that is kept clears its
  # coordinate's.
  both <- note(recs, c(0, 1), TRUE, 900L, 800L, "sample")
  expect_identical(vapply(both, `[[`, integer(1), "fails"), c(1L, 1L))
  expect_identical(note(both, c(-1, 0.05), FALSE, 900L, 800L, "sample"), both[2])
  # The records stay bounded.
  many <- Reduce(function(r, a) note(r, c(cos(a), sin(a)), TRUE, 1L, 1L, "sample"),
                 seq(0, 1.5, length.out = 12), list())
  expect_lte(length(many), 4L)
})

# Chains of a window of L draws whose mean coordinates (the plan's directions
# are the axes) are AR(1) series with the autocorrelations rho, for
# .emc_shear_set_curves with the orbit widths and knot fits mocked.
curve_chains <- function(rho, L = 3000L, n = 4L, stage = "sample", seed = 21) {
  set.seed(seed)
  p <- length(rho)
  lapply(1:3, function(j) {
    mu <- t(vapply(rho, function(r) as.numeric(stats::filter(rnorm(L), r, "recursive")) *
                     sqrt(1 - r^2), numeric(L)))
    list(n_pars = p, n_subjects = n, par_names = paste0("p", seq_len(p)),
         lik_laplace = rep(list(list(prec = diag(p), lin = numeric(p))), n),
         samples = list(theta_mu = mu, theta_var = array(diag(p), c(p, p, L)),
                        alpha = array(0, c(p, n, L)), stage = rep(stage, L), idx = L),
         shear = list(plan = list(dirs = diag(p), log_step = rep(0, 2L * p), stage = stage)))
  })
}

# Mocks under which a coordinate's orbit is r_fixed (r_curve under a curve) of
# its posterior sd; returns the counter of knot fits.
mock_curve_work <- function(emc, r_fixed, r_curve, env = parent.frame()) {
  sds <- apply(do.call(cbind, lapply(emc, function(ch) ch$samples$theta_mu)), 1, sd)
  count <- new.env()
  count$fits <- 0L
  count$widths <- 0L
  local_mocked_bindings(
    .emc_shear_orbit_width = function(ch, idx, k, fits_of, s0, n_cores = 1L) {
      count$widths <- count$widths + 1L
      (if (is.null(fits_of(NULL)$curve)) r_fixed else r_curve)[k] * sds[k]
    },
    .emc_shear_curve_subject = function(chain, s, mus, V, x0) {
      count$fits <- count$fits + 1L
      p <- nrow(mus)
      rep(list(list(prec = 2 * diag(p), lin = numeric(p), mode = numeric(p))), ncol(mus))
    }, .package = "EMC2", .env = env)
  count
}

test_that("a slow coordinate gets a curve when its round is worth it, and no fits when it cannot be", {
  # 1 is the slowest but its fixed-frame orbit spans the margin; 2 is slow
  # and narrow; 3 is narrow but mixes.
  emc <- curve_chains(c(0.98, 0.95, 0.3))
  count <- mock_curve_work(emc, r_fixed = c(0.9, 0.4, 0.4), r_curve = c(1, 0.7, 0.9))
  out <- EMC2:::.emc_shear_set_curves(emc, 1:3000)
  for (ch in out) {
    rec <- ch$shear$curve_rec
    expect_identical(rec$k, 2L)
    expect_true(is.list(ch$shear[["curve"]]))
    expect_identical(dim(ch$shear[["curve"]]$cube), c(3L, 3L, 4L, length(EMC2:::.EMC_SHEAR_CURVE_KNOTS)))
    expect_identical(vapply(rec$tried, `[[`, integer(1), "k"), 1:2)
    expect_identical(vapply(rec$tried, `[[`, logical(1), "fitted"), c(FALSE, TRUE))
    expect_identical(vapply(rec$tried, `[[`, logical(1), "passed"), c(FALSE, TRUE))
    expect_equal(rec$last$r_fixed, 0.4, tolerance = 1e-8)
    expect_equal(rec$last$r_curve, 0.7, tolerance = 1e-8)
    expect_equal(rec$last$worth, EMC2:::.emc_shear_curve_worth(0.4, 0.7))
    expect_gte(rec$last$tau, EMC2:::.EMC_SHEAR_CURVE_TAU)
    expect_identical(c(rec$screens, rec$builds), c(2, 1))
    # The screened coordinate is remembered, the curved one is not.
    expect_length(rec$failed, 1L)
    expect_equal(rec$failed[[1]]$row, c(1, 0, 0))
    # The curved coordinate's direction moves it by one.
    expect_equal(sum(ch$shear[["curve"]]$row * ch$shear$plan$dirs[, 2]), 1)
  }
  # Knots were fitted for the one coordinate that could pass, per subject.
  expect_identical(count$fits, 3L * 4L)
  # A curve that covers the window is kept without any work.
  count$fits <- count$widths <- 0L
  again <- EMC2:::.emc_shear_set_curves(out, 1:3000)
  expect_identical(c(count$fits, count$widths), c(0L, 0L))
  expect_identical(again, out)
})

test_that("a rebuild's curve work is bounded and failures are not repeated on the same evidence", {
  emc <- curve_chains(c(0.99, 0.985, 0.98, 0.975, 0.9, 0.9), L = 2400L)
  # Four coordinates' orbits span the margin; two are narrow but no curve
  # widens them enough.
  count <- mock_curve_work(emc, r_fixed = c(0.9, 0.9, 0.9, 0.9, 0.4, 0.4),
                           r_curve = c(1, 1, 1, 1, 0.45, 0.45))
  tried <- function(e) lapply(e, function(ch) vapply(ch$shear$curve_rec$tried, `[[`, integer(1), "k"))
  first <- EMC2:::.emc_shear_set_curves(emc, 1:1200)
  for (ks in tried(first)) expect_length(ks, EMC2:::.EMC_SHEAR_CURVE_SCREENS)
  expect_identical(count$fits, 0L)
  expect_false(any(vapply(first, function(ch) is.list(ch$shear[["curve"]]), logical(1))))
  # The plan is left as it was: a failed candidate changes no round.
  expect_identical(lapply(first, function(ch) ch$shear$plan), lapply(emc, function(ch) ch$shear$plan))
  # On the same window the failed coordinates stay out and the next ones are
  # tried, with knot fits for at most the stated number of coordinates.
  second <- EMC2:::.emc_shear_set_curves(first, 1:1200)
  for (j in 1:3) {
    expect_length(intersect(tried(second)[[j]], tried(first)[[j]]), 0L)
    expect_length(tried(second)[[j]], 2L)
    expect_length(second[[j]]$shear$curve_rec$failed, 6L)
  }
  expect_identical(count$fits, 3L * 4L * EMC2:::.EMC_SHEAR_CURVE_TRIES)
  # Every slow coordinate has failed on this window: nothing is tried.
  count$fits <- count$widths <- 0L
  third <- EMC2:::.emc_shear_set_curves(second, 1:1200)
  expect_identical(c(count$fits, count$widths), c(0L, 0L))
  expect_identical(third, second)
  # Nor while the new draws are fewer than the failures were judged on.
  expect_identical(EMC2:::.emc_shear_set_curves(second, 800:1999), second)
  expect_identical(c(count$fits, count$widths), c(0L, 0L))
  # Once they outnumber them the coordinates are tried again, and a second
  # failure waits twice as long.
  later <- EMC2:::.emc_shear_set_curves(second, 1201:2400)
  expect_gt(count$widths, 0L)
  for (ch in later) {
    again <- Filter(function(r) r$fails == 2L, ch$shear$curve_rec$failed)
    expect_gte(length(again), EMC2:::.EMC_SHEAR_CURVE_TRIES)
    expect_true(all(vapply(again, function(r) r$at == 2400L && r$len == 1200L, logical(1))))
  }
})

test_that("draws from burn-in ask for no curve, and the switch turns curves off", {
  emc <- curve_chains(c(0.98, 0.95, 0.3), stage = "burn")
  count <- mock_curve_work(emc, r_fixed = c(0.3, 0.3, 0.3), r_curve = c(0.9, 0.9, 0.9))
  expect_identical(EMC2:::.emc_shear_set_curves(emc, 1:3000), emc)
  expect_identical(count$widths, 0L)
  # Only the draws after burn-in count: too few of them measure nothing.
  mixed <- lapply(emc, function(ch) { ch$samples$stage[2971:3000] <- "adapt"; ch })
  expect_identical(EMC2:::.emc_shear_set_curves(mixed, 1:3000), mixed)
  settled <- lapply(emc, function(ch) { ch$samples$stage[1001:3000] <- "adapt"; ch })
  expect_true(all(vapply(EMC2:::.emc_shear_set_curves(settled, 1:3000),
                         function(ch) is.list(ch$shear[["curve"]]), logical(1))))
  sampled <- curve_chains(c(0.98, 0.95, 0.3))
  withr::local_options(emc2.shear_curve = FALSE)
  expect_identical(EMC2:::.emc_shear_set_curves(sampled, 1:3000), sampled)
})

test_that("a curve that no longer covers the window is refitted, and dropped once that fails twice", {
  emc <- curve_chains(c(0.98, 0.95, 0.3))
  count <- mock_curve_work(emc, r_fixed = c(0.9, 0.4, 0.4), r_curve = c(1, 0.7, 0.9))
  held <- EMC2:::.emc_shear_set_curves(emc, 1:3000)
  # The window moves off the knots along the curved coordinate.
  moved <- lapply(held, function(ch) {
    ch$samples$theta_mu[2, ] <- ch$samples$theta_mu[2, ] + 5
    ch
  })
  count$fits <- 0L
  refit <- EMC2:::.emc_shear_set_curves(moved, 1:3000)
  expect_identical(count$fits, 3L * 4L)
  for (j in 1:3) {
    expect_identical(refit[[j]]$shear$curve_rec$k, 2L)
    expect_identical(refit[[j]]$shear$curve_rec$strikes, 0L)
    expect_equal(refit[[j]]$shear[["curve"]]$knots, held[[j]]$shear[["curve"]]$knots + 5)
  }
  # A refit whose curve is no longer worth it leaves the old curve once.
  mock_curve_work(emc, r_fixed = c(0.9, 0.4, 0.4), r_curve = c(1, 0.42, 0.9))
  once <- EMC2:::.emc_shear_set_curves(moved, 1:3000)
  for (j in 1:3) {
    expect_identical(once[[j]]$shear[["curve"]], held[[j]]$shear[["curve"]])
    expect_identical(once[[j]]$shear$curve_rec[c("k", "strikes")], list(k = 2L, strikes = 1L))
    expect_identical(once[[j]]$shear$curve_rec$failed, held[[j]]$shear$curve_rec$failed)
  }
  twice <- EMC2:::.emc_shear_set_curves(once, 1:3000)
  for (ch in twice) {
    expect_null(ch$shear[["curve"]])
    expect_null(ch$shear$curve_rec$k)
    expect_equal(lapply(ch$shear$curve_rec$failed, `[[`, "row"), list(c(1, 0, 0), c(0, 1, 0)))
  }
})

test_that("a curve stays while its rounds pay in autocorrelation and cost, and is dropped once they do not", {
  emc <- curve_chains(c(0.98, 0.95, 0.3), L = 3200L)
  count <- mock_curve_work(emc, r_fixed = c(0.9, 0.4, 0.4), r_curve = c(1, 0.7, 0.9))
  held <- EMC2:::.emc_shear_set_curves(emc, 1:3000)
  for (ch in held) {
    b <- ch$shear$curve_rec$base
    expect_identical(b$at, 3000L)
    expect_gte(b$tau, EMC2:::.EMC_SHEAR_CURVE_TAU)
  }
  # The draws since the curve was built mix: it is kept, though its round
  # doubles the plan's cost.
  plan_cost <- function(e, cost, c_round) lapply(e, function(ch) {
    ch$alloc$history <- list(list(cost = cost, c_round = c_round, curved_rounds = 1))
    ch
  })
  mixing <- plan_cost(lapply(held, function(ch) {
    ch$samples$theta_mu[, 3001:3200] <- matrix(rnorm(3 * 200), 3)
    ch
  }), 2, 1)
  count$fits <- count$widths <- 0L
  expect_identical(EMC2:::.emc_shear_set_curves(mixing, 1:3200), mixing)
  expect_identical(c(count$fits, count$widths), c(0L, 0L))
  # Without a priced plan nothing is judged.
  unpriced <- lapply(mixing, function(ch) { ch$alloc$history <- NULL; ch })
  expect_identical(EMC2:::.emc_shear_set_curves(unpriced, 1:3200), unpriced)
  # The round now costs far more than the autocorrelation time fell, but too
  # few draws have run under the curve to judge it.
  dear <- plan_cost(mixing, 1e4, 1e4 - 1)
  expect_identical(EMC2:::.emc_shear_set_curves(dear, 1:3050), dear)
  # With enough of them the curve is dropped, its coordinate waits like one
  # whose curve failed, and the slow one beside it still spans its margin.
  out <- EMC2:::.emc_shear_set_curves(dear, 1:3200)
  for (ch in out) {
    expect_null(ch$shear[["curve"]])
    expect_null(ch$shear$curve_rec$k)
    expect_false(ch$shear$curve_rec$last$paid)
    rows <- lapply(ch$shear$curve_rec$failed, `[[`, "row")
    expect_true(any(vapply(rows, function(r) isTRUE(all.equal(r, c(0, 1, 0))), logical(1))))
    expect_true(EMC2:::.emc_shear_curve_waits(ch$shear$curve_rec$failed[[length(rows)]], 3200L, "sample"))
  }
  expect_identical(count$fits, 0L)
})

test_that("a curve is reused while it covers the window", {
  cv <- list(knots = 2 + 0.5 * EMC2:::.EMC_SHEAR_CURVE_KNOTS)
  set.seed(4)
  z <- rnorm(300)
  expect_true(EMC2:::.emc_shear_curve_covers(cv, 2 + 0.5 * z))
  expect_true(EMC2:::.emc_shear_curve_covers(cv, 2.5 + 0.4 * z))
  expect_false(EMC2:::.emc_shear_curve_covers(cv, 3.2 + 0.5 * z))
  expect_false(EMC2:::.emc_shear_curve_covers(cv, 2 + 1.5 * z))
  expect_false(EMC2:::.emc_shear_curve_covers(cv, 2 + 0.1 * z))
})

test_that("a coordinate that changes kernel starts its record afresh", {
  p <- 2L
  rec <- list(E = 1:4 + 0, n = 5:8 + 0, moved = 1:4 + 0, N = 1:4 + 0, iters = rep(50, 4),
              K = matrix(1, 4, 3))
  plan <- list(pool = rec, last = rec, log_step = rep(0.3, 4))
  out <- EMC2:::.emc_shear_restart_coord(plan, 2L, 0.4)
  for (nm in c("pool", "last")) {
    expect_identical(out[[nm]]$n, c(5, 0, 7, 8))
    expect_identical(out[[nm]]$E, c(1, 0, 3, 4))
    expect_identical(out[[nm]]$K[2, ], numeric(3))
    expect_identical(out[[nm]]$K[-2, ], rec$K[-2, ])
  }
  expect_equal(out$log_step, c(0.3, log(EMC2:::.EMC_SHEAR_STEP0 * 0.4), 0.3, 0.3))
  expect_identical(EMC2:::.emc_shear_restart_coord(plan, 2L)$log_step, plan$log_step)
})

test_that("a pinned direction keeps its coordinate and leaves a basis", {
  set.seed(8)
  p <- 4L
  D <- matrix(rnorm(p * p), p) + diag(2, p)
  row <- rnorm(p); dir <- rnorm(p)
  out <- EMC2:::.emc_shear_pin_dir(D, 2L, dir, row)
  # Row 2 of the inverse is the coordinate, and direction 2 moves it by one.
  expect_equal(solve(out)[2, ], row, tolerance = 1e-10)
  expect_equal(out[, 2], dir / sum(row * dir), tolerance = 1e-12)
  # The other directions lose only their part along the coordinate.
  expect_equal(drop(row %*% out[, -2]), numeric(p - 1), tolerance = 1e-12)
  expect_equal(out[, -2] + outer(out[, 2], drop(row %*% D[, -2])), D[, -2], tolerance = 1e-12)
  # Pinning a plan's own direction and coordinate changes nothing.
  expect_equal(EMC2:::.emc_shear_pin_dir(D, 3L, D[, 3], solve(D)[3, ]), D, tolerance = 1e-10)
  # A direction that does not move the coordinate cannot be pinned.
  flat <- dir - row * sum(row * dir) / sum(row^2)
  expect_null(EMC2:::.emc_shear_pin_dir(D, 2L, flat, row))
  expect_null(EMC2:::.emc_shear_pin_dir(D, 2L, dir[-1], row))
})

test_that("a curve's direction is the one the group mean covaries along", {
  set.seed(9)
  p <- 3L
  S <- crossprod(matrix(rnorm(p * p), p)) + diag(0.5, p)
  chains <- lapply(1:3, function(i) list(samples = list(
    theta_mu = t(MASS::mvrnorm(4000, rnorm(p, 0, 5), S)))))
  row <- c(1, -0.5, 0.2)
  u <- EMC2:::.emc_shear_curve_dir(chains, 1:4000, row, fallback = c(9, 9, 9))
  # Within chains, whatever their offsets: S row' scaled to move the
  # coordinate by one, along which the coordinate's line has its full width.
  expect_equal(u, drop(S %*% row) / drop(row %*% S %*% row), tolerance = 0.08)
  expect_equal(sum(row * u), 1, tolerance = 1e-12)
  expect_equal(1 / drop(u %*% solve(S) %*% u), drop(row %*% S %*% row), tolerance = 0.1)
  # No usable window: the fallback.
  one <- lapply(chains, function(ch) { ch$samples$theta_mu[] <- 1; ch })
  expect_identical(EMC2:::.emc_shear_curve_dir(one, 1:10, row, fallback = c(9, 9, 9)), c(9, 9, 9))
})

test_that("a chain without a curve tries another chain's coordinate first", {
  set.seed(10)
  p <- 4L
  D <- matrix(rnorm(p * p), p) + diag(3, p)
  chain <- function(dirs, row = NULL) list(shear = list(
    plan = list(dirs = dirs), curve = if (!is.null(row)) list(row = row)))
  # Chain 2 holds a curve on what is nearly chain 1's coordinate 3.
  near <- solve(D)[3, ] + rnorm(p, 0, 0.02)
  emc <- list(chain(D), chain(D[, c(3, 1, 2, 4)], near), chain(D))
  expect_identical(EMC2:::.emc_shear_curve_match(D, near)$k, 3L)
  expect_identical(EMC2:::.emc_shear_curve_candidates(emc, 1L, own = c(1L, 3L, 2L)),
                   c(3L, 1L, 2L))
  expect_identical(EMC2:::.emc_shear_curve_candidates(emc, 1L, own = integer(0)), 3L)
  # A coordinate that failed here and is still waited on stays out.
  expect_identical(EMC2:::.emc_shear_curve_candidates(emc, 1L, own = integer(0),
                                                      failed = list(solve(D)[3, ])), integer(0))
  # No curve anywhere: the chain's own candidates.
  none <- list(chain(D), chain(D), chain(D))
  expect_identical(EMC2:::.emc_shear_curve_candidates(none, 1L, own = c(2L, 4L)), c(2L, 4L))
  # A curve along no coordinate of this chain adds nothing.
  off <- list(chain(D), chain(D, colSums(solve(D))), chain(D))
  expect_identical(EMC2:::.emc_shear_curve_candidates(off, 1L, own = 2L), 2L)
})

test_that("what a coordinate's rounds can refresh is its orbits' mean variance over the window", {
  p <- 2L
  n <- 200L
  set.seed(4)
  chain <- function() {
    tv <- array(0, c(p, p, n))
    tv[1, 1, ] <- exp(2 * rnorm(n, 0, 0.5))
    tv[2, 2, ] <- exp(2 * rnorm(n, -1, 0.2))
    list(n_pars = p, data = list(), samples = list(idx = n, theta_mu = matrix(rnorm(p * n), p),
                                                   theta_var = tv),
         shear = list(plan = list(dirs = diag(p), log_step = rep(0, 2L * p),
                                  last = list(n = c(500, 0, 300, 10)))))
  }
  emc <- list(chain(), chain())
  seen <- new.env()
  seen$draws <- list()
  local_mocked_bindings(
    # A mean direction's orbit has sd 0.5; a log sd's is wider where it is low.
    .emc_shear_orbit_width = function(ch, idx, k, fits_of, s0, n_cores = 1L) {
      seen$draws[[length(seen$draws) + 1L]] <- c(k, ch$samples$idx)
      if (k <= p) 0.5 else exp(-0.5 * log(ch$samples$theta_var[k - p, k - p, ch$samples$idx]))
    },
    .emc_shear_context = function(ch) list(packed = NULL),
    .emc_ll_data_prepared = function(data) NULL,
    .emc_rebuild_lapply = function(X, FUN, n_cores, ...) lapply(X, FUN),
    .package = "EMC2")
  out <- EMC2:::.emc_shear_orbit_shares(emc, seq_len(n), 1L)
  m <- EMC2:::.EMC_SHEAR_SHARE_STATES
  for (j in 1:2) {
    cv <- out[[j]]$shear$plan$cond_var
    expect_identical(is.finite(cv), c(TRUE, FALSE, TRUE, FALSE))
    expect_equal(cv[1], 0.25)
    # The states are the coordinate's quantile mid-points; widths are in
    # window sds of the log sd.
    lsd <- 0.5 * log(emc[[j]]$samples$theta_var[1, 1, ])
    at <- order(lsd)[round((seq_len(m) - 0.5) / m * n)]
    expect_equal(cv[3], mean((exp(-lsd[at]) * sd(lsd))^2))
  }
  d <- do.call(rbind, seen$draws)
  expect_identical(nrow(d), 4L * m)
  expect_identical(sort(unique(d[, 1])), c(1L, 3L))
  # A state whose orbit is not resolved leaves the coordinate unmeasured.
  local_mocked_bindings(
    .emc_shear_orbit_width = function(ch, idx, k, fits_of, s0, n_cores = 1L)
      if (k == 3L && ch$samples$idx %% 2L == 0L) NA_real_ else 0.5,
    .package = "EMC2")
  out <- EMC2:::.emc_shear_orbit_shares(emc, seq_len(n), 1L)
  expect_identical(is.finite(out[[1]]$shear$plan$cond_var), c(TRUE, FALSE, FALSE, FALSE))
  # Too short a window measures nothing.
  short <- EMC2:::.emc_shear_orbit_shares(emc, seq_len(2L * m), 1L)
  expect_true(all(is.na(short[[1]]$shear$plan$cond_var)))
})
