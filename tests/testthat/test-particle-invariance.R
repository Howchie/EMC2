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
                   list("sample", c(.3, .02, .3, .38), c(1, 1, 1)))) {
    x <- run(cfg[[1]], cfg[[2]], cfg[[3]])
    expect_lt(abs(mean(x) - post_mu[1]) / bm_se(x), 4.5)
    expect_lt(abs(var(x) - post_var[1, 1]) / bm_se((x - mean(x))^2), 4.5)
  }
})


# ---- Joint group-subject move --------------------------------------------------
# .emc_group_move() translates the group mean and rescales the group SDs,
# carrying every subject with them (an interweaving step). The tests below
# check it against exact conjugate Gibbs samplers, with the real Gibbs steps.

gm_mock_ll <- function(props, data, model, marginalise = NULL) {
  r <- sweep(props, 2, data$y)
  -0.5 * rowSums((r %*% data$prec) * r)
}

gm_sampler <- function(type, p, n, par_groups = NULL, group_design = NULL) {
  pn <- letters[seq_len(p)]
  s <- list(type = type, n_pars = p, n_subjects = n, par_names = pn,
            nuisance = rep(FALSE, p), marginalised_idx = rep(FALSE, p),
            model = NULL, marginalise = NULL)
  if (type == "diagonal-gamma") return(EMC2:::add_info_diag_gamma(s, prior = NULL))
  if (is.null(group_design)) {
    return(EMC2:::add_info_standard(s, prior = NULL, par_groups = par_groups))
  }
  s$par_group <- par_groups
  s$is_blocked <- duplicated(par_groups) | duplicated(par_groups, fromLast = TRUE)
  gd <- EMC2:::add_group_design(pn, group_design, n)
  M <- sum(vapply(gd, ncol, 1L))
  s$prior <- list(theta_mu_mean = rep(0, M), theta_mu_var = diag(M), v = 2,
                  A = rep(.3, p), theta_mu_invar = diag(M))
  s$group_designs <- group_design
  s$gd <- gd
  s$XtX <- crossprod(do.call(cbind, gd))
  nc <- vapply(gd, ncol, 1L)
  s$design_row_idx <- rep(seq_along(nc), nc)
  s$design_n_cols <- nc
  s
}

# Exact Gibbs on y_s ~ N(alpha_s, obs), with the package's group-level Gibbs
# step, optionally followed by one group move per iteration.
gm_chain <- function(sampler, y, obs_prec, iters, move, V = NULL, seed = 1) {
  set.seed(seed)
  p <- nrow(y); n <- ncol(y)
  M <- length(sampler$prior$theta_mu_mean)
  sampler$data <- lapply(seq_len(n), function(s) list(y = y[, s], prec = obs_prec))
  sampler$samples <- list(theta_mu = matrix(0, M, 1),
                          theta_var = array(diag(p), c(p, p, 1)),
                          last_theta_var_inv = diag(p),
                          a_half = matrix(1, p, 1), idx = 1)
  if (!is.null(V)) {
    attr(V, "source") <- "empirical"
    sampler$group_move_cov <- V
  }
  alpha <- y
  rownames(alpha) <- sampler$par_names
  out <- matrix(NA_real_, iters, 3 * p + 1)
  for (t in seq_len(iters)) {
    pars <- EMC2:::gibbs_step(sampler, alpha, sampler$type)
    sm <- if (!is.null(pars$subj_mu)) pars$subj_mu else matrix(pars$tmu, p, n)
    Pv <- solve(pars$tvar)
    Vs <- solve(Pv + obs_prec)
    Ls <- t(chol(Vs))
    for (s in seq_len(n)) {
      alpha[, s] <- Vs %*% (Pv %*% sm[, s] + obs_prec %*% y[, s]) + Ls %*% rnorm(p)
    }
    if (move) {
      ll <- vapply(seq_len(n), function(s) {
        r <- alpha[, s] - y[, s]
        -0.5 * sum(r * (obs_prec %*% r))
      }, 0)
      sampler$rng$gibbs <- NULL
      mv <- EMC2:::.emc_group_move(sampler, pars, rbind(alpha, ll), "sample",
                                   NULL, NULL, NULL, 1L)
      sampler <- mv$sampler
      pars <- mv$pars
      alpha <- mv$proposals[seq_len(p), , drop = FALSE]
    }
    sampler$samples$theta_mu[, 1] <- pars$tmu
    sampler$samples$theta_var[, , 1] <- pars$tvar
    sampler$samples$last_theta_var_inv <- pars$tvinv
    if (!is.null(pars$a_half)) sampler$samples$a_half[, 1] <- pars$a_half
    tv <- pars$tvar
    out[t, ] <- c(pars$tmu[seq_len(p)], log(diag(tv)), alpha[, 1],
                  if (p > 1) tv[1, 2] / sqrt(tv[1, 1] * tv[2, 2]) else 0)
  }
  if (move) {
    attr(out, "move_rate") <- sampler$group_move_state$moves /
      sampler$group_move_state$attempts
  }
  out
}

gm_bm_se <- function(x, nb = 50) {
  m <- colMeans(matrix(x[seq_len(floor(length(x) / nb) * nb)], ncol = nb))
  stats::sd(m) / sqrt(nb)
}

test_that("group move leaves the joint posterior invariant for every covariance type", {
  skip_model_validation()
  local_mocked_bindings(.emc_group_move_ll_one = gm_mock_ll, .package = "EMC2")
  withr::local_options(list(emc2.group_move = TRUE, emc2.group_move_scale = TRUE))
  n <- 6
  x <- seq(-1, 1, length.out = n)
  cases <- list(
    list("standard, full IW, K = 1", "standard", 2, c(1, 1), NULL, 1L),
    list("standard, full IW, K = 4", "standard", 2, c(1, 1), NULL, 4L),
    list("standard, unblocked", "standard", 2, c(1, 2), NULL, 4L),
    list("standard, mixed blocks", "standard", 3, c(1, 1, 2), NULL, 4L),
    list("diagonal-gamma, K = 1", "diagonal-gamma", 2, NULL, NULL, 1L),
    list("diagonal-gamma, K = 4", "diagonal-gamma", 2, NULL, NULL, 4L),
    list("standard, group design", "standard", 2, c(1, 2),
         list(a = matrix(cbind(1, x), n,
                         dimnames = list(NULL, c("(Intercept)", "x")))), 4L))
  set.seed(7)
  for (cs in cases) {
    p <- cs[[3]]
    withr::local_options(list(emc2.group_move_proposals = cs[[6]]))
    smp <- gm_sampler(cs[[2]], p, n, cs[[4]], cs[[5]])
    y <- matrix(rnorm(p * n, 0.5, 1), p)
    obs_prec <- diag(1 / 1.5, p)
    ref <- gm_chain(smp, y, obs_prec, 40000L, FALSE, seed = 11)
    mv <- gm_chain(smp, y, obs_prec, 40000L, TRUE, V = diag(0.15, 2 * p), seed = 12)
    rate <- attr(mv, "move_rate")
    expect_true(rate > 0.1 && rate < 0.9, info = cs[[1]])
    for (j in seq_len(ncol(ref))) {
      if (stats::sd(ref[, j]) == 0) next
      z <- (mean(mv[, j]) - mean(ref[, j])) /
        sqrt(gm_bm_se(mv[, j])^2 + gm_bm_se(ref[, j])^2)
      expect_lt(abs(z), 5, label = sprintf("%s, column %d: |z|", cs[[1]], j))
    }
  }
})

test_that("group move spec selects intercepts, skips fixed coordinates, and weights Jacobians", {
  gd <- list(
    a = cbind(a = rep(1, 3), a_age = c(-1, 0, 1)),
    b = cbind(b_age = c(-1, 0, 1), b = rep(1, 3)))
  sampler <- list(type = "standard", n_pars = 4L, n_subjects = 3L,
                  par_names = c("a", "b", "c", "nuisance"),
                  nuisance = c(FALSE, FALSE, FALSE, TRUE),
                  marginalised_idx = c(FALSE, TRUE, FALSE, FALSE),
                  par_group = c(1, 1, 1, 2), is_blocked = c(TRUE, TRUE, TRUE, FALSE),
                  gd = c(gd, list(c = cbind(c_x = c(1, 2, 3)))))
  spec <- EMC2:::.emc_group_move_spec(sampler)
  # b is marginalised, c has no intercept column, the nuisance row is excluded.
  expect_identical(spec$trans_main, 1L)
  expect_identical(spec$mu_idx, 1L)
  expect_identical(spec$scale_main, c(1L, 3L))
  # The marginalised coordinate leaves the IW block, so a and c form a 2-block.
  expect_identical(spec$jac_weight, c(3L, 3L))
  expect_identical(spec$n_move, 3L)

  no_scale <- EMC2:::.emc_group_move_spec(sampler, scale = FALSE)
  expect_identical(no_scale$n_s, 0L)

  dg <- list(type = "diagonal-gamma", n_pars = 2L, n_subjects = 5L,
             par_names = c("a", "b"), nuisance = c(FALSE, FALSE))
  expect_identical(EMC2:::.emc_group_move_spec(dg)$jac_weight, c(2L, 2L))
  expect_null(EMC2:::.emc_group_move_spec(modifyList(dg, list(n_subjects = 1L))))
  expect_null(EMC2:::.emc_group_move_spec(modifyList(dg, list(type = "factor"))))
})

test_that("group-design move shifts intercepts, rescales residuals, and keeps slopes", {
  gd <- list(
    a = cbind(intercept = rep(1, 3), slope = c(-1, 0, 1)),
    b = cbind(slope = c(-1, 0, 1), intercept = rep(1, 3)))
  sampler <- list(type = "standard", n_pars = 2L, n_subjects = 3L,
                  par_names = c("a", "b"), nuisance = c(FALSE, FALSE),
                  par_group = c(1, 2), is_blocked = c(FALSE, FALSE), gd = gd)
  spec <- EMC2:::.emc_group_move_spec(sampler)
  tmu <- c(0.5, 2, 1, -1)
  subj_mu <- EMC2:::calculate_subject_means(gd, tmu)
  alpha <- subj_mu + matrix(c(0.3, -0.2, 0.1, 0.4, -0.5, 0.2), 2)
  state <- list(tmu = tmu, subj_mu = subj_mu, alpha = alpha,
                tvar = diag(c(0.5, 2)), tvinv = diag(c(2, 0.5)))
  theta <- c(0.7, -0.4, log(2), log(0.5))
  moved <- EMC2:::.emc_group_move_apply(theta, spec, state)
  expect_equal(moved$tmu, c(1.2, 2, 1, -1.4))
  expect_equal(moved$subj_mu, EMC2:::calculate_subject_means(gd, moved$tmu))
  expect_equal(moved$alpha - moved$subj_mu, (alpha - subj_mu) * c(2, 0.5))
  expect_equal(moved$tvar, diag(c(2, 0.5)))
  expect_equal(moved$tvinv, solve(moved$tvar))
})

gm_fake_emc <- function(chains, prior = NULL) {
  p <- nrow(chains[[1]]$mu)
  base <- list(type = "diagonal-gamma", n_pars = p, n_subjects = 10L,
               par_names = letters[seq_len(p)], nuisance = rep(FALSE, p),
               prior = if (is.null(prior)) {
                 list(theta_mu_mean = rep(0, p), theta_mu_invar = diag(1e-8, p),
                      shape = rep(2, p), rate = rep(0.3, p))
               } else prior)
  lapply(chains, function(ch) {
    tv <- array(0, c(p, p, ncol(ch$mu)))
    for (i in seq_len(ncol(ch$mu))) tv[, , i] <- diag(exp(2 * ch$logsd[, i]), p)
    modifyList(base, list(samples = list(theta_mu = ch$mu, theta_var = tv)))
  })
}

test_that("group move covariance is Cov(theta) minus the Gibbs scale, pooled within chains", {
  set.seed(3)
  N <- 4000
  p <- 2
  C_mu <- matrix(c(0.2, 0.05, 0.05, 0.15), 2)
  draw_chain <- function(offset) {
    list(mu = t(mvtnorm::rmvnorm(N, sigma = C_mu)) + offset,
         logsd = matrix(log(0.5) + rnorm(p * N, sd = 0.4), p))
  }
  # Chains centred far apart: pooling without centring would add 25 to Var(mu_a).
  emc <- gm_fake_emc(list(draw_chain(c(0, 0)), draw_chain(c(10, -5))))
  spec <- EMC2:::.emc_group_move_spec(emc[[1]])
  V <- EMC2:::.emc_group_move_covariance(emc, spec, seq_len(N))
  G <- EMC2:::.emc_group_move_gibbs_cov(spec, emc[[1]]$prior,
                                        diag(0.25 * exp(2 * 0.4^2), p))
  expect_identical(attr(V, "source"), "empirical")
  # iid draws: ESS is large, so correlations are kept (shrink = d / ESS).
  expect_lt(attr(V, "shrink"), 0.01)
  expect_equal(unname(V[1:2, 1:2]), C_mu - G[1:2, 1:2], tolerance = 0.08)
  expect_equal(unname(diag(V)[3:4]), rep(0.4^2, 2) - diag(G)[3:4], tolerance = 0.1)
})

test_that("group move covariance stays diagonal until ESS supports correlations", {
  set.seed(4)
  N <- 300
  # A slow random walk: strongly correlated coordinates, tiny ESS.
  steps <- t(mvtnorm::rmvnorm(N, sigma = matrix(c(1, 0.9, 0.9, 1), 2) * 0.01))
  mu <- t(apply(steps, 1, cumsum))
  emc <- gm_fake_emc(list(list(mu = mu, logsd = matrix(log(0.1), 2, N))))
  withr::local_options(list(emc2.group_move_scale = FALSE))
  spec <- EMC2:::.emc_group_move_spec(emc[[1]], scale = FALSE)
  V <- EMC2:::.emc_group_move_covariance(emc, spec, seq_len(N))
  expect_lt(attr(V, "ess_min"), 5 * spec$n_move)
  expect_identical(attr(V, "shrink"), 1)
  expect_equal(V[1, 2], 0)

  # Where the posterior is no wider than the Gibbs step, the floor holds V at
  # 10% of the Gibbs scale instead of collapsing it.
  flat <- gm_fake_emc(list(list(mu = matrix(rnorm(2 * N, sd = 1e-4), 2),
                                logsd = matrix(log(0.1), 2, N))))
  Vf <- EMC2:::.emc_group_move_covariance(flat, spec, seq_len(N))
  G <- EMC2:::.emc_group_move_gibbs_cov(spec, flat[[1]]$prior, diag(0.01, 2))
  expect_equal(unname(diag(Vf)), 0.1 * diag(G), tolerance = 1e-6)
})

test_that("adaptation targets are the diffusion-limit ESJD optimum", {
  # Metropolis: the Roberts-Gelman-Gilks limit, recovered from the same speed
  # function the table is built from.
  mh <- stats::optimize(function(l) -l^2 * 2 * stats::pnorm(-l / 2), c(0.5, 5))
  expect_equal(mh$minimum, 2.38, tolerance = 0.01)
  expect_equal(2 * stats::pnorm(-mh$minimum / 2), 0.234, tolerance = 0.01)
  tab <- EMC2:::.EMC_GROUP_MOVE_LIMIT
  expect_true(all(diff(tab$p_move) > 0) && all(diff(tab$l_opt) > 0))
  # The shipped table is the output of the derivation.
  skip_model_validation()
  for (K in c(2L, 4L, 8L)) {
    re <- EMC2:::.emc_group_move_limit(K, n_mc = 20000L, seed = 99L)
    row <- tab[tab$K == K, ]
    expect_equal(re[["p_move"]], row$p_move, tolerance = 0.03)
    expect_equal(re[["l_opt"]], row$l_opt, tolerance = 0.05)
  }
  expect_equal(EMC2:::.emc_group_move_target(64L, 10L),
               tab$p_move[tab$K == max(tab$K)])
})

test_that("installing a covariance keeps the adapted scale unless the shape changes a lot", {
  spec <- list(n_move = 2L)
  sampler <- list(group_move_cov = diag(2),
                  group_move_state = list(lambda = 0.4, adapt_iter = 50L, dim = 2L,
                                          K = 4L, attempts = 60L, moves = 20L))
  kept <- EMC2:::.emc_group_move_install(sampler, diag(c(1.5, 0.8)), spec, 4L)
  expect_equal(kept$group_move_state$lambda, 0.4)
  expect_identical(kept$group_move_state$adapt_iter, 50L)

  reset <- EMC2:::.emc_group_move_install(sampler, diag(c(10, 1)), spec, 4L)
  opt4 <- EMC2:::.emc_group_move_optimum(4L)
  expect_equal(reset$group_move_state$lambda, unname(opt4[["l_opt"]]) / sqrt(2))
  expect_equal(reset$group_move_state$target, unname(opt4[["p_move"]]))
  expect_identical(reset$group_move_state$adapt_iter, 0L)
  expect_identical(reset$group_move_state$attempts, 60L)

  other_k <- EMC2:::.emc_group_move_install(sampler, diag(2), spec, 1L)
  expect_equal(other_k$group_move_state$lambda, 2.38 / sqrt(2))
  expect_equal(other_k$group_move_state$target, 0.234)
})

gm_one_step_setup <- function() {
  sampler <- gm_sampler("diagonal-gamma", 2, 4)
  y <- matrix(c(0.2, 0.5, -0.1, 0.3, 0.4, -0.2, 0.1, 0.6), 2)
  sampler$data <- lapply(1:4, function(s) list(y = y[, s], prec = diag(2)))
  alpha <- y
  rownames(alpha) <- c("a", "b")
  pars <- list(tmu = c(a = 0.1, b = 0.2), tvar = diag(0.3, 2), tvinv = diag(1 / 0.3, 2),
               alpha = alpha)
  list(sampler = sampler, pars = pars, proposals = rbind(alpha, rep(0, 4)))
}

test_that("group move adapts its scale before sampling and freezes it in sample", {
  local_mocked_bindings(.emc_group_move_ll_one = gm_mock_ll, .package = "EMC2")
  withr::local_options(list(emc2.group_move_proposals = 4L))
  x <- gm_one_step_setup()
  set.seed(5)
  s <- x$sampler
  for (i in 1:30) {
    s$rng$gibbs <- NULL
    s <- EMC2:::.emc_group_move(s, x$pars, x$proposals, "burn", NULL, NULL, NULL)$sampler
  }
  expect_identical(s$group_move_state$adapt_iter, 30L)
  expect_false(isTRUE(all.equal(
    s$group_move_state$lambda,
    unname(EMC2:::.emc_group_move_optimum(4L)[["l_opt"]]) / sqrt(4))))
  lam <- s$group_move_state$lambda
  s$rng$gibbs <- NULL
  s2 <- EMC2:::.emc_group_move(s, x$pars, x$proposals, "sample", NULL, NULL, NULL)$sampler
  expect_equal(s2$group_move_state$lambda, lam)
  expect_identical(s2$group_move_state$adapt_iter, 30L)
  expect_identical(s2$group_move_state$attempts, 31L)
})

test_that("group move draws from the Gibbs stream, not the global RNG state", {
  local_mocked_bindings(.emc_group_move_ll_one = gm_mock_ll, .package = "EMC2")
  x <- gm_one_step_setup()
  set.seed(99)
  x$sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  set.seed(1)
  a <- EMC2:::.emc_group_move(x$sampler, x$pars, x$proposals, "sample", NULL, NULL, NULL)
  set.seed(2)
  b <- EMC2:::.emc_group_move(x$sampler, x$pars, x$proposals, "sample", NULL, NULL, NULL)
  expect_identical(a$proposals, b$proposals)
  expect_identical(a$sampler$rng$gibbs, b$sampler$rng$gibbs)
})

test_that("group move is skipped with one warning when the prior cannot be evaluated", {
  local_mocked_bindings(.emc_group_move_ll_one = gm_mock_ll, .package = "EMC2")
  x <- gm_one_step_setup()
  x$sampler$prior$theta_mu_invar <- diag(3)
  env <- EMC2:::.emc_group_move_env
  rm(list = ls(env), envir = env)
  withr::defer(rm(list = ls(env), envir = env))
  expect_warning(
    out <- EMC2:::.emc_group_move(x$sampler, x$pars, x$proposals, "sample",
                                  NULL, NULL, NULL),
    "group prior could not be evaluated")
  expect_identical(out$proposals, x$proposals)
  expect_no_warning(EMC2:::.emc_group_move(x$sampler, x$pars, x$proposals, "sample",
                                           NULL, NULL, NULL))
})

test_that("group move is a no-op when disabled", {
  withr::local_options(list(emc2.group_move = FALSE))
  x <- gm_one_step_setup()
  out <- EMC2:::.emc_group_move(x$sampler, x$pars, x$proposals, "burn", NULL, NULL, NULL)
  expect_identical(out$proposals, x$proposals)
  expect_null(out$sampler$group_move_state)
})

test_that("group move likelihood dispatch keeps the subject and candidate mapping", {
  context <- list(data = list(list(offset = 10), list(offset = -3)),
                  model = "mock", marginalise = NULL)
  props <- list(matrix(c(0.5, 1.5), 2), matrix(c(1.25, 2.25), 2))
  message <- list(kind = "group_move_ll", subs = c(2L, 1L), props = props)
  local_mocked_bindings(
    .emc_group_move_ll_one = function(props, data, model, marginalise) {
      props[, 1] + data$offset
    },
    .package = "EMC2")
  result <- EMC2:::.emc_wpool_compute(message, context)
  expect_equal(result$ll, cbind(c(-2.5, -1.5), c(11.25, 12.25)))

  sent <- new.env(parent = emptyenv())
  sent$messages <- vector("list", 2L)
  pool <- list(n = 2L, alive = TRUE)
  all_props <- list(matrix(c(0.5, 1), 2), matrix(c(-1, 0), 2), matrix(c(2, 3), 2))
  result <- testthat::with_mocked_bindings(
    EMC2:::.emc_wpool_group_move_ll(pool, list(c(1L, 3L), 2L), all_props,
                                    list(data = vector("list", 3L))),
    .emc_wpool_deadline = function() Inf,
    .emc_wpool_request = function(pool, w, msg, deadline) {
      sent$messages[[w]] <- msg
      TRUE
    },
    .emc_wpool_reply = function(pool, w, deadline) {
      msg <- sent$messages[[w]]
      list(ll = vapply(seq_along(msg$subs), function(i) {
        msg$props[[i]][, 1] + msg$subs[i]
      }, numeric(2)))
    },
    .package = "EMC2")
  expect_equal(result$ll, cbind(c(1.5, 2), c(1, 2), c(5, 6)))
  expect_identical(sent$messages[[1L]]$subs, c(1L, 3L))
})

test_that("fit with the group move matches a weak-data exact Gibbs reference", {
  skip_model_validation()
  skip_on_os("windows")

  data("forstmann", package = "EMC2")
  subject_names <- as.character(unique(forstmann$subjects)[1:4])
  dat <- droplevels(forstmann[forstmann$subjects %in% subject_names, ])
  y <- stats::setNames(c(-0.7, 0.1, 0.3, 1.2), subject_names)
  observation_var <- 25

  des <- design(data = dat, model = LNR,
                formula = list(m ~ 1, s ~ 1), constants = c(t0 = log(0.2)))
  fake_likelihood <- function(proposals, dadm, model, r_cores = 1L,
                              marginalise = NULL, ...) {
    subject <- as.character(unique(dadm$subjects)[1L])
    ll <- -0.5 * (proposals[, "m"] - y[[subject]])^2 / observation_var
    matrix(ll, ncol = 1L, dimnames = list(NULL, "mock_gaussian"))
  }

  set.seed(20260925L)
  fitted <- testthat::with_mocked_bindings({
    emc <- suppressMessages(make_emc(
      dat, des, type = "diagonal-gamma", compress = FALSE, n_chains = 2L))
    fit(emc, cores_for_chains = 1L,
        stop_criteria = list(
          preburn = list(iter = 20L), burn = list(iter = 50L),
          adapt = list(iter = 50L, min_unique = 0L),
          sample = list(iter = 600L)),
        verbose = FALSE, particle_factor = 10, step_size = 25L,
        max_tries = 1L, trim = FALSE)
  }, calc_ll_manager = fake_likelihood, .package = "EMC2")
  fitted <- EMC2:::restore_duplicates(fitted)

  sampler_prior <- fitted[[1L]]$prior
  prior_mean <- sampler_prior$theta_mu_mean[1L]
  prior_var <- sampler_prior$theta_mu_var[1L, 1L]
  shape <- sampler_prior$shape[1L]
  rate <- sampler_prior$rate[1L]
  set.seed(20260926L)
  n_ref <- 70000L
  burn_ref <- 10000L
  sigma2 <- 0.3
  mu <- prior_mean
  alpha <- rep(mu, length(y))
  reference <- numeric(n_ref)
  for (iter in seq_len(n_ref)) {
    alpha_var <- 1 / (1 / sigma2 + 1 / observation_var)
    alpha_mean <- alpha_var * (mu / sigma2 + y / observation_var)
    alpha <- stats::rnorm(length(y), alpha_mean, sqrt(alpha_var))
    mu_var <- 1 / (1 / prior_var + length(y) / sigma2)
    mu_mean <- mu_var * (prior_mean / prior_var + sum(alpha) / sigma2)
    mu <- stats::rnorm(1L, mu_mean, sqrt(mu_var))
    sigma2 <- 1 / stats::rgamma(
      1L, shape = shape + length(y) / 2,
      rate = rate + sum((alpha - mu)^2) / 2)
    reference[iter] <- mu
  }
  reference <- reference[-seq_len(burn_ref)]

  sample_by_chain <- lapply(fitted, function(chain) {
    idx <- which(chain$samples$stage == "sample")
    chain$samples$theta_mu["m", idx]
  })
  expect_identical(
    vapply(fitted, function(chain) sum(chain$samples$stage == "adapt"),
           integer(1)),
    rep(50L, length(fitted))
  )
  # The empirical covariance was built from posterior draws before sampling,
  # and the move was attempted every iteration and moved some of the time.
  expect_true(all(vapply(fitted, function(chain) {
    V <- chain$group_move_cov
    st <- chain$group_move_state
    !is.null(V) && identical(attr(V, "source"), "empirical") &&
      all(is.finite(V)) && st$adapt_iter > 0L &&
      st$attempts >= 720L && st$moves > 0L
  }, logical(1))))

  batch_means <- function(x, n_batch) {
    n_each <- floor(length(x) / n_batch)
    x <- x[seq_len(n_each * n_batch)]
    colMeans(matrix(x, nrow = n_each))
  }
  batch_mcse <- function(x, n_batch) {
    means <- batch_means(x, n_batch)
    stats::sd(means) / sqrt(length(means))
  }
  fit_mean <- mean(unlist(sample_by_chain))
  ref_mean <- mean(reference)
  fit_mean_se <- stats::sd(unlist(lapply(sample_by_chain, batch_means,
                                         n_batch = 20L))) /
    sqrt(20L * length(sample_by_chain))
  ref_mean_se <- batch_mcse(reference, 50L)
  expect_lt(abs(fit_mean - ref_mean) / sqrt(fit_mean_se^2 + ref_mean_se^2), 5)

  fit_second <- mean(unlist(lapply(sample_by_chain, function(x) {
    (x - ref_mean)^2
  })))
  ref_second <- mean((reference - ref_mean)^2)
  fit_second_se <- stats::sd(unlist(lapply(sample_by_chain, function(x) {
    batch_means((x - ref_mean)^2, 20L)
  }))) / sqrt(20L * length(sample_by_chain))
  ref_second_se <- batch_mcse((reference - ref_mean)^2, 50L)
  expect_lt(abs(fit_second - ref_second) /
              sqrt(fit_second_se^2 + ref_second_se^2), 5)
})
