local_rng_guard()

test_that("the composed sample iteration preserves exact bounded marginals", {
  skip_model_validation()
  N <- 2500L
  h <- diag_hierarchy_draws(N, b = -0.4, p = 2L, n = 5L)
  p <- h$p; n <- h$n; pn <- paste0("x", seq_len(p))
  fits <- lapply(seq_len(n), function(s)
    list(prec = h$Si, lin = drop(h$Si %*% h$y[, s]),
         check = list(pareto_k = 0, outside = 0, laplace_status = "ok")))
  data <- lapply(seq_len(n), function(s) data.frame(subject = s))
  attr(data, "components") <- stats::setNames(rep(1L, p), pn)
  attr(data, "shared_ll_idx") <- stats::setNames(rep(1L, p), pn)
  broken_mode <- FALSE
  hmc_round <- EMC2:::.emc_hmc_round
  local_mocked_bindings(
    calc_ll_manager = function(proposals, dadm, ...) {
      x <- as.matrix(proposals)
      s <- as.integer(dadm$subject[1L])
      ll <- h$ll(x, h$y[, s])
      ll[x[, 1L] < h$b] <- -Inf
      ll
    },
    .emc_ensemble_options = function() list(pool = 8L, rounds = 2L,
      tail_weight = 0.25, tail_df = 4, tail_excess = 3),
    .emc_hmc_round = function(plan, proposals, pars, type, streams, data, model,
      wpool = NULL, n_cores = 1L, r_cores = 1L, par_names, pick = NULL,
      ll_cache_validated = FALSE) {
      if (broken_mode) pick <- proposals[1L, ] > h$y[1L, ]
      hmc_round(plan, proposals, pars, type, streams, data, model, wpool,
                n_cores, r_cores, par_names, pick,
                ll_cache_validated = ll_cache_validated)
    },
    .package = "EMC2")

  zero_move <- list(E = numeric(2L * p), n = integer(2L * p),
                    moved = integer(2L * p), N = numeric(2L * p),
                    K = matrix(0L, 2L * p, 1L), iters = numeric(2L * p))
  collapsed_dim <- EMC2:::.emc_collapsed_chart(EMC2:::.emc_shear_v_spec(
    list(type = "diagonal-gamma", n_pars = p, n_subjects = n,
         par_group = seq_len(p), is_blocked = rep(FALSE, p), prior = h$prior)), p)$dim
  run <- function(broken = FALSE) {
    broken_mode <<- broken
    local_test_rng(if (broken) 132L else 131L)
    i <- sample.int(N, 1L)
    mu <- h$mu[i, ]; names(mu) <- pn
    V <- diag(exp(2 * h$lsd[i, ]), p)
    dimnames(V) <- list(pn, pn)
    alpha <- matrix(h$a[i, , ], p, n, dimnames = list(pn, NULL))
    ll0 <- vapply(seq_len(n), function(s)
      h$ll(matrix(alpha[, s], 1L), h$y[, s]), 0)
    ll0[alpha[1L, ] < h$b] <- -Inf
    samples <- list(idx = 1L, stage = "init",
      theta_mu = matrix(mu, p, 1L, dimnames = list(pn, NULL)),
      theta_var = array(V, c(p, p, 1L), dimnames = list(pn, pn, NULL)),
      alpha = array(alpha, c(p, n, 1L), dimnames = list(pn, NULL, NULL)),
      subj_ll = matrix(ll0, n, 1L), last_theta_var_inv = solve(V))
    pm <- lapply(seq_len(n), function(s) list(list(
      mix = c(0.05, 0.3, 0.3, 0.35), epsilon = rep(0.5, 3L),
      n_particles = 24L, iter = 100L, proposal_counts = integer(4L),
      acc_counts = integer(4L), jump_counts = numeric(4L), gd_good = FALSE,
      sel = list())))
    attr(samples, "pm_settings") <- pm
    ybar <- lapply(seq_len(n), function(s) stats::setNames(h$y[, s], pn))
    proposal_var <- diag(c(4, 3), p)
    dimnames(proposal_var) <- list(pn, pn)
    eff_var0 <- diag(h$S, p)
    dimnames(eff_var0) <- list(pn, pn)
    chains_var <- rep(list(proposal_var), n)
    hmc_subjects <- lapply(seq_len(n), function(s) list(
      eligible = TRUE, reason = "", prec = h$Si, eps = 0.7, rate = 0.5,
      h = EMC2:::.EMC_HMC_H, floor = -Inf, W = diag(p), warmup = FALSE,
      n_warm = 0L, n_tuned = 1L, calib_cpu = 0))
    sampler <- list(type = "diagonal-gamma", nuisance = rep(FALSE, p),
      n_pars = p, n_subjects = n, par_names = pn, data = data, model = NULL,
      prior = h$prior, lik_laplace = fits, samples = samples,
      eff_mu = lapply(seq_len(n), function(s) stats::setNames(numeric(p), pn)),
      eff_var = rep(list(eff_var0), n),
      chains_mu = ybar, chains_var = chains_var,
      ensemble_q = list(mu = ybar, var = chains_var), marginalise = NULL,
      shear = list(plan = list(dirs = diag(p), rounds = rep(1, 2L * p),
                               log_step = rep(log(2.5), 2L * p),
                               block = zero_move, last = zero_move,
                               stage = "sample"), walls = h$walls),
      collapsed = list(plan = list(C = diag(0.5, collapsed_dim), eps = 0.35,
        rounds = 1, dim = collapsed_dim, last = list(E = matrix(0, 2L * p, 2L * p),
          n = 0L, acc = 0, moved = 0L), block = list(E = matrix(0, 2L * p, 2L * p),
          n = 0L, acc = 0, moved = 0L)), fits = fits),
      subject_hmc = list(supported = TRUE, p = p, n = n, stage = "sample",
        subjects = hmc_subjects, evidence = vector("list", n),
        totals = list(attempts = 0L, accepted = 0L, aborts = 0L, rows = 0,
                      cpu = 0, calib_cpu = 0)),
      alloc = list(window = NULL, history = list(),
        particle = list(rate = rep(if (broken) 0 else 0.5, n)),
        last = list(rate = rep(0.5, n))))
    out <- EMC2:::run_stage(sampler, stage = "sample", iter = N,
      particles = 24L, n_cores = 1L, tune = list(), verbose = FALSE,
      verboseProgress = FALSE, r_cores = 1L)
    idx <- seq.int(1002L, N + 1L)
    t(vapply(idx, function(k) c(out$samples$theta_mu[, k],
      0.5 * log(diag(out$samples$theta_var[, , k])),
      as.numeric(out$samples$alpha[, , k])), numeric(2L * p + p * n)))
  }
  z_mean <- function(x) {
    nb <- 25L
    batches <- split(seq_len(nrow(x)), cut(seq_len(nrow(x)), nb, labels = FALSE))
    se_chain <- apply(vapply(batches, function(ii) colMeans(x[ii, , drop = FALSE]),
                             numeric(ncol(x))), 1L, stats::sd) / sqrt(nb)
    se_ref <- apply(h$before, 2L, stats::sd) / sqrt(nrow(h$before))
    abs(colMeans(x) - colMeans(h$before)) / sqrt(se_chain^2 + se_ref^2)
  }
  good <- run(FALSE)
  expect_true(all(z_mean(good) < 5))
  expect_true(all(vapply(seq_len(ncol(good)), function(k)
    suppressWarnings(stats::ks.test(good[, k], h$before[, k])$p.value) > 1e-4,
    logical(1))))
  bad <- run(TRUE)
  expect_gt(max(z_mean(bad)), 5)
})
