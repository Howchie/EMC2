# Ensemble group update (R/ensemble_move.R): posterior invariance on a
# conjugate Gaussian hierarchy, using the package's own group Gibbs step.

ens_mock_ll <- function(props, data, model, marginalise = NULL) {
  r <- sweep(props, 2, data$y)
  -0.5 * rowSums((r %*% data$prec) * r)
}

ens_chain <- function(type, y, obs_prec, iters, ensemble, seed) {
  set.seed(seed)
  p <- nrow(y); n <- ncol(y); pn <- rownames(y)
  s <- list(type = type, n_pars = p, n_subjects = n, par_names = pn,
            nuisance = rep(FALSE, p), marginalised_idx = rep(FALSE, p),
            model = NULL, marginalise = NULL)
  s <- if (type == "diagonal-gamma") EMC2:::add_info_diag_gamma(s, prior = NULL) else
    EMC2:::add_info_standard(s, prior = NULL)
  s$data <- lapply(seq_len(n), function(k) list(y = y[, k], prec = obs_prec))
  s$samples <- list(theta_mu = matrix(0, p, 1), theta_var = array(diag(p), c(p, p, 1)),
                    last_theta_var_inv = diag(p), a_half = matrix(1, p, 1),
                    alpha = array(0, c(p, n, 1)), subj_ll = matrix(0, n, 1),
                    stage = "sample", idx = 1)
  # Frozen pool proposals centred on the data, independent of the group state.
  obs_var <- solve(obs_prec)
  chains_mu <- lapply(seq_len(n), function(k) y[, k])
  chains_var <- lapply(seq_len(n), function(k) obs_var)
  alpha <- y
  out <- matrix(NA_real_, iters, 2 * p + 1)
  for (t in seq_len(iters)) {
    pars <- EMC2:::gibbs_step(s, alpha, type)
    sm <- if (!is.null(pars$subj_mu)) pars$subj_mu else matrix(pars$tmu, p, n)
    Pv <- solve(pars$tvar); Vs <- solve(Pv + obs_prec); Ls <- t(chol(Vs))
    for (k in seq_len(n)) alpha[, k] <- Vs %*% (Pv %*% sm[, k] + obs_prec %*% y[, k]) +
      Ls %*% rnorm(p)
    if (ensemble) {
      ll <- vapply(seq_len(n), function(k) ens_mock_ll(t(alpha[, k]), s$data[[k]]), 0)
      s$rng$gibbs <- NULL
      mv <- EMC2:::.emc_ensemble_move(s, pars, rbind(alpha, ll), "sample", NULL, NULL,
                                      NULL, chains_mu, chains_var)
      s <- mv$sampler; pars <- mv$pars
      alpha <- mv$proposals[seq_len(p), , drop = FALSE]
    }
    s$samples$theta_mu[, 1] <- pars$tmu
    s$samples$theta_var[, , 1] <- pars$tvar
    s$samples$last_theta_var_inv <- pars$tvinv
    if (!is.null(pars$a_half)) s$samples$a_half[, 1] <- pars$a_half
    tv <- pars$tvar
    out[t, ] <- c(pars$tmu, 0.5 * log(diag(tv)), tv[1, 2] / sqrt(tv[1, 1] * tv[2, 2]))
  }
  out
}

ens_bm_se <- function(x, nb = 50) {
  m <- colMeans(matrix(x[seq_len(floor(length(x) / nb) * nb)], ncol = nb))
  stats::sd(m) / sqrt(nb)
}

test_that("ensemble group update leaves the joint posterior invariant", {
  skip_model_validation()
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, model, ...)
    ens_mock_ll(proposals, dadm, model), .package = "EMC2")
  withr::local_options(list(emc2.ensemble = TRUE, emc2.ensemble_pool = 16L,
                            emc2.ensemble_rounds = 5L))
  set.seed(3)
  y <- rbind(a = rnorm(6, 0, 0.8), b = rnorm(6, 0, 0.8))
  obs_prec <- solve(matrix(c(1.5, 1, 1, 1.5), 2))
  for (type in c("standard", "diagonal-gamma")) {
    ref <- ens_chain(type, y, obs_prec, 60000, FALSE, 1)[-(1:5000), ]
    ens <- ens_chain(type, y, obs_prec, 12000, TRUE, 2)[-(1:1000), ]
    cols <- if (type == "standard") 1:5 else 1:4
    for (k in cols) {
      z <- (mean(ens[, k]) - mean(ref[, k])) / sqrt(ens_bm_se(ens[, k])^2 + ens_bm_se(ref[, k])^2)
      expect_lt(abs(z), 4.5, label = paste(type, "column", k, "mean z"))
      zv <- (var(ens[, k]) - var(ref[, k])) /
        sqrt(ens_bm_se((ens[, k] - mean(ens[, k]))^2)^2 + ens_bm_se((ref[, k] - mean(ref[, k]))^2)^2)
      expect_lt(abs(zv), 4.5, label = paste(type, "column", k, "variance z"))
    }
  }
})

test_that("ensemble update is a no-op when disabled or unsupported", {
  withr::local_options(list(emc2.ensemble = FALSE))
  s <- list(type = "standard", n_subjects = 3L, n_pars = 2L, nuisance = c(FALSE, FALSE))
  out <- EMC2:::.emc_ensemble_move(s, list(tmu = 1), matrix(0, 3, 3), "sample",
                                   NULL, NULL, NULL, list(), list())
  expect_identical(out$pars, list(tmu = 1))
  withr::local_options(list(emc2.ensemble = TRUE))
  s$type <- "single"
  out <- EMC2:::.emc_ensemble_move(s, list(tmu = 1), matrix(0, 3, 3), "sample",
                                   NULL, NULL, NULL, list(), list())
  expect_identical(out$pars, list(tmu = 1))
})

# Latent-factor hierarchies. Subject states are drawn exactly from
#   ref : the textbook conditional a_s | h, eta_s (no ensemble), or
#   prod: the eta-marginal subject prior from get_group_level(), which is what
#         the particle step targets, optionally followed by the ensemble move
#         with either reselection weight.
# All four must target the same posterior.
ens_latent_chain <- function(type, y, obs_prec, iters, arm, seed, x = NULL) {
  set.seed(seed)
  p <- nrow(y); n <- ncol(y); pn <- rownames(y)
  nf <- if (type == "infnt_factor") 2L else 1L
  ss <- list(Lambda_mat = matrix(c(1, Inf, Inf), p, 1), B_mat = matrix(0, 1, 1),
             K_mat = matrix(0, p, 1), G_mat = matrix(Inf, 1, 1), covariates = x)
  s <- list(type = type, n_pars = p, n_subjects = n, par_names = pn,
            nuisance = rep(FALSE, p), marginalised_idx = rep(FALSE, p),
            model = NULL, marginalise = NULL)
  s <- switch(type,
    factor = EMC2:::add_info_factor(s, n_factors = nf),
    infnt_factor = EMC2:::add_info_infnt_factor(s, n_factors = nf),
    SEM = EMC2:::add_info_SEM(s, sem_settings = ss))
  store <- get(paste0("sample_store_", type), asNamespace("EMC2"))
  s$samples <- store(data.frame(subjects = factor(seq_len(n))), pn, iters = 1,
                     stage = "sample", integrate = TRUE, is_nuisance = rep(FALSE, p),
                     n_factors = nf, sem_settings = ss,
                     Lambda_mat = if (type == "factor") attr(s, "Lambda_mat"))
  s$samples$idx <- 1
  start <- get(paste0("get_startpoints_", type), asNamespace("EMC2"))(s, rep(0, p), diag(p))
  s$samples <- EMC2:::fill_samples(s$samples, start, rbind(y, 0), 1, p, type)
  s$data <- lapply(seq_len(n), function(k) list(y = y[, k], prec = obs_prec))
  draw <- function(m, V) {
    Pv <- solve(V); Vs <- solve(Pv + obs_prec); Ls <- t(chol(Vs))
    vapply(seq_len(n), function(k) as.numeric(Vs %*% (Pv %*% m[, k] + obs_prec %*% y[, k]) +
                                                Ls %*% rnorm(p)), numeric(p))
  }
  ll_of <- function(a) vapply(seq_len(n), function(k) ens_mock_ll(t(a[, k]), s$data[[k]]), 0)
  withr::local_options(list(emc2.ensemble = arm %in% c("cond", "marg"),
    emc2.ensemble_latent_weight = if (arm == "marg") "marginal" else "conditional"))
  chains_mu <- lapply(seq_len(n), function(k) y[, k])
  chains_var <- lapply(seq_len(n), function(k) solve(obs_prec))
  alpha <- y
  out <- matrix(NA_real_, iters, 2 * p + 1)
  for (t in seq_len(iters)) {
    pars <- EMC2:::gibbs_step(s, alpha, type)
    if (arm == "ref") {
      cp <- EMC2:::.emc_ensemble_subject_prior(pars, s, n, "conditional")
      alpha <- draw(cp$mean, cp$var)
    } else {
      if (arm != "prod") {
        s$rng$gibbs <- NULL
        mv <- EMC2:::.emc_ensemble_move(s, pars, rbind(alpha, ll_of(alpha)), "sample",
                                        NULL, NULL, NULL, chains_mu, chains_var)
        s <- mv$sampler; pars <- mv$pars
      }
      gl <- lapply(seq_len(n), function(k) EMC2:::get_group_level(pars, k, type))
      alpha <- draw(vapply(gl, function(g) as.numeric(g$mu), numeric(p)), gl[[1]]$var)
    }
    s$samples <- EMC2:::fill_samples(s$samples, pars, rbind(alpha, ll_of(alpha)), 1, p, type)
    tv <- pars$tvar
    out[t, ] <- c(pars$tmu, 0.5 * log(diag(tv)), tv[1, 2] / sqrt(tv[1, 1] * tv[2, 2]))
  }
  out
}

test_that("latent-factor hierarchies: particle target and ensemble weights are invariant", {
  skip_model_validation()
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, model, ...)
    ens_mock_ll(proposals, dadm, model), .package = "EMC2")
  withr::local_options(list(emc2.ensemble_pool = 16L, emc2.ensemble_rounds = 5L))
  set.seed(11)
  n <- 12
  f <- rnorm(n)
  y <- rbind(a = f, b = 0.7 * f, c = -0.6 * f) + matrix(rnorm(3 * n, 0, 0.3), 3)
  obs_prec <- diag(1 / 0.5^2, 3)
  for (type in c("factor", "infnt_factor")) {
    ref <- ens_latent_chain(type, y, obs_prec, 40000, "ref", 1)[-(1:2000), ]
    for (arm in c("prod", "cond", "marg")) {
      d <- ens_latent_chain(type, y, obs_prec, if (arm == "prod") 40000 else 10000,
                            arm, match(arm, c("prod", "cond", "marg")) + 1)[-(1:1000), ]
      for (k in seq_len(ncol(d))) {
        z <- (mean(d[, k]) - mean(ref[, k])) / sqrt(ens_bm_se(d[, k])^2 + ens_bm_se(ref[, k])^2)
        expect_lt(abs(z), 4.5, label = paste(type, arm, "column", k, "mean z"))
      }
    }
  }
})

test_that("SEM subject prior conditions on the subject's own covariates", {
  x <- matrix(c(-1, 0, 2), 3, 1, dimnames = list(NULL, "x"))
  lam <- matrix(c(1, 0.5), 2, 1)
  pars <- list(tmu = c(0.1, -0.2), lambda = lam, B = matrix(0, 1, 1),
               K = matrix(c(0.3, 0), 2, 1), G = matrix(0.8, 1, 1),
               delta_inv = matrix(2, 1, 1), epsilon_inv = c(4, 5),
               eta = matrix(0, 3, 1))
  # What gibbs_step_SEM returns for these values.
  subj_mu <- pars$tmu + (pars$K + lam %*% pars$G) %*% t(x)
  subj_var <- lam %*% (1 / pars$delta_inv) %*% t(lam) + diag(1 / pars$epsilon_inv)
  gl <- lapply(1:3, function(s) EMC2:::get_group_level_SEM(
    c(pars, list(subj_mu = subj_mu, subj_var = subj_var)), s))
  expect_equal(gl[[3]]$mu, c(0.1 + 0.3 * 2 + 0.8 * 2, -0.2 + 0.5 * 0.8 * 2))
  expect_false(isTRUE(all.equal(gl[[1]]$mu, gl[[3]]$mu)))
  expect_equal(gl[[1]]$var, subj_var)
  # And the ensemble's marginal weight uses the same prior.
  s <- list(type = "SEM", sem_settings = list(covariates = x))
  sp <- EMC2:::.emc_ensemble_subject_prior(pars, s, 3L, "marginal")
  expect_equal(sp$mean, subj_mu, ignore_attr = TRUE)
  expect_equal(sp$var, subj_var, ignore_attr = TRUE)
})
