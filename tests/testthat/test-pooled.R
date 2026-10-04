# design(pooled = ...): coefficients shared by all subjects (R/pooled.R).

pooled_design <- function(data = NULL, pooled = "t0", ...) {
  args <- list(model = LBA, matchfun = function(d) d$S == d$lR,
               formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
               constants = c(sv = log(1)), pooled = pooled,
               report_p_vector = FALSE, ...)
  if (is.null(data)) {
    args$factors <- list(subjects = 1:6, S = c("left", "right"))
    args$Rlevels <- c("left", "right")
  } else args$data <- data
  do.call(design, args)
}

pooled_truth <- c(v = 1, v_lMTRUE = 1.5, B = log(1), A = log(.4), t0 = log(.25))

pooled_data <- function(n_trials = 80, seed = 1) {
  set.seed(seed)
  des <- pooled_design()
  re <- make_random_effects(des, pooled_truth, n_subj = 6)
  make_data(re, des, n_trials = n_trials)
}

test_that("sampled_pars leaves pooled coefficients out unless asked", {
  des <- pooled_design()
  expect_equal(names(sampled_pars(des)), c("v", "v_lMTRUE", "B", "A"))
  expect_equal(names(sampled_pars(des, include_pooled = TRUE)),
               c("v", "v_lMTRUE", "B", "A", "t0"))
  expect_error(pooled_design(pooled = "nope"), "not in sampled_pars")
  expect_error(design(model = LBA, factors = list(subjects = 1:2, S = 1:2),
                      Rlevels = 1:2, matchfun = function(d) d$S == d$lR,
                      formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                      constants = c(sv = log(1), t0 = log(.2)), pooled = "t0",
                      report_p_vector = FALSE),
               "also set as constants")
})

test_that("make_random_effects gives every subject the pooled mean", {
  set.seed(2)
  re <- make_random_effects(pooled_design(), pooled_truth, n_subj = 6)
  expect_true(all(re[, "t0"] == pooled_truth[["t0"]]))
  expect_gt(stats::sd(re[, "B"]), 0)
})

test_that("prior() gives pooled coefficients their own normal prior", {
  des <- pooled_design()
  pr <- prior(des, pmean = c(t0 = log(.3), B = 0.5), psd = c(t0 = .5))
  expect_equal(attr(pr, "pooled"),
               list(mean = c(t0 = log(.3)), var = c(t0 = .25)))
  expect_false("t0" %in% names(pr$theta_mu_mean))
  expect_equal(pr$theta_mu_mean[["B"]], 0.5)
  pr2 <- prior(des, update = pr, psd = c(t0 = 2))
  expect_equal(attr(pr2, "pooled"), list(mean = c(t0 = log(.3)), var = c(t0 = 4)))
  expect_equal(attr(prior(des), "pooled"), list(mean = c(t0 = 0), var = c(t0 = 1)))
  expect_error(prior(des, type = "single"), "hierarchical type")
})

test_that("the likelihood boundary appends pooled values only when missing", {
  dadm <- data.frame(x = 1)
  attr(dadm, "pooled") <- "t0"
  p <- matrix(1:4, 2, dimnames = list(NULL, c("a", "b")))
  expect_error(EMC2:::.emc_pooled_columns(p, dadm), "No value for pooled")
  attr(dadm, "pooled_values") <- c(t0 = -1.5)
  out <- EMC2:::.emc_pooled_columns(p, dadm, varying = c(TRUE, FALSE))
  expect_equal(out$proposals[, "t0"], c(-1.5, -1.5))
  expect_equal(out$varying, c(TRUE, FALSE, FALSE))
  supplied <- cbind(p, t0 = c(0, 1))
  expect_identical(EMC2:::.emc_pooled_columns(supplied, dadm)$proposals, supplied)
  expect_identical(EMC2:::.emc_pooled_columns(p, data.frame(x = 1))$proposals, p)
})

test_that("a pooled value scores the same as a column and as a stored value", {
  dat <- pooled_data()
  des <- pooled_design(dat)
  ml <- des$model()
  ml$pre_transform$func[["t0"]] <- "exp"
  ml$pre_transform$lower[["t0"]] <- -3
  des$model <- function() ml
  emc <- make_emc(dat, des, n_chains = 1, compress = FALSE, rt_resolution = NULL)
  d1 <- emc[[1]]$data[[1]]
  a <- c(v = 1, v_lMTRUE = 1.5, B = 0, A = log(.4))
  t0s <- log(log(c(.15, .2, .25)) + 3)  # log t0 = -3 + exp(coefficient)
  as_cols <- EMC2:::calc_ll_manager(cbind(matrix(a, 3, 4, byrow = TRUE,
                                                dimnames = list(NULL, names(a))),
                                         t0 = t0s), d1, emc[[1]]$model)
  stored <- vapply(t0s, function(x) {
    EMC2:::calc_ll_manager(t(a), EMC2:::.emc_pooled_set_one(d1, c(t0 = x)),
                           emc[[1]]$model)
  }, numeric(1))
  expect_equal(stored, as_cols, tolerance = 1e-12)
  expect_gt(stats::sd(as_cols), 1)        # the candidates are actually scored
  rows <- EMC2:::.emc_pooled_rows_ll(
    cbind(matrix(a, 3, 4, byrow = TRUE, dimnames = list(NULL, names(a))), t0 = t0s),
    EMC2:::.emc_pooled_set_one(d1, c(t0 = 99)), emc[[1]]$model)
  expect_equal(rows, as_cols, tolerance = 1e-12)
})

# Conjugate check of the pooled move: subject s contributes
# -0.5 (phi - y_s)' P (phi - y_s), so phi | alpha, y is Gaussian in closed form.
pooled_kernel_chain <- function(y, prec, prior_mean, prior_var, iters, seed) {
  set.seed(seed)
  p <- nrow(y); n <- ncol(y); pn <- rownames(y)
  local_mocked_bindings(
    .emc_pooled_rows_ll = function(rows, dadm, model) {
      r <- sweep(rows[, pn, drop = FALSE], 2, attr(dadm, "y"))
      -0.5 * rowSums((r %*% prec) * r)
    }, .package = "EMC2")
  st <- EMC2:::.emc_pooled_state(list(mean = stats::setNames(prior_mean, pn),
                                      var = stats::setNames(prior_var, pn)))
  st$cov <- solve(n * prec); st$log_eps <- log(1.5)
  data <- lapply(seq_len(n), function(k) {
    d <- data.frame(x = 0); attr(d, "y") <- y[, k]; d
  })
  s <- list(pooled = st, n_pars = 1L, n_subjects = n, par_names = "a",
            data = data, model = NULL)
  ll <- function(phi) vapply(seq_len(n), function(k) {
    r <- phi - y[, k]; -0.5 * sum(r * (prec %*% r))
  }, 0)
  props <- rbind(a = rep(0, n), ll(st$value))
  out <- matrix(NA_real_, iters, p)
  for (t in seq_len(iters)) {
    mv <- EMC2:::.emc_pooled_move(s, props, "sample", NULL, NULL, NULL, NULL,
                                  NULL, NULL, NULL)
    s <- mv$sampler; props <- mv$proposals
    out[t, ] <- s$pooled$value
  }
  list(draws = out, ll_ok = isTRUE(all.equal(props[2, ], ll(s$pooled$value))),
       sampler = s)
}

test_that("the pooled move leaves the conditional posterior invariant", {
  y <- matrix(c(0.3, -0.1, 0.5, 0.2, 0.0, 0.4, 1.0, 0.8, 1.3, 0.7, 1.1, 0.9),
              nrow = 2, byrow = TRUE, dimnames = list(c("t0", "B"), NULL))
  prec <- matrix(c(4, 1.5, 1.5, 2), 2)
  m0 <- c(1, -1); v0 <- c(0.05, 0.2)
  post_prec <- diag(1 / v0) + ncol(y) * prec
  post_var <- solve(post_prec)
  post_mean <- drop(post_var %*% (m0 / v0 + prec %*% rowSums(y)))
  res <- pooled_kernel_chain(y, prec, m0, v0, iters = 6000, seed = 11)
  d <- res$draws[-(1:200), ]
  expect_true(res$ll_ok)
  ess <- nrow(d) / 10                      # conservative for an MTM chain
  expect_lt(max(abs(colMeans(d) - post_mean) / sqrt(diag(post_var) / ess)), 4)
  expect_lt(max(abs(sqrt(diag(stats::cov(d)) / diag(post_var)) - 1)), 0.06)
  expect_lt(abs(stats::cor(d)[1, 2] - stats::cov2cor(post_var)[1, 2]), 0.06)
  expect_gt(res$sampler$pooled$moves[["sample"]] / 6000, 0.3)
})

# Shear check: subject s contributes -0.5 prec (y_s - a_s - phi)^2 (a ridge in
# a_s + phi), a_s ~ N(mu, tau2) and mu ~ N(mm, vm). Exact Gibbs draws of a and
# mu alternate with the pooled move, which carries every a_s and mu along the
# slope -1; the joint of (phi, a, mu) is Gaussian.
test_that("the sheared pooled move leaves the joint posterior invariant", {
  set.seed(5)
  n <- 5; y <- c(0.2, 0.6, -0.1, 0.4, 0.3); prec <- 50; tau2 <- 0.04
  m0 <- 0.3; v0 <- 0.5; mm <- 0; vm <- 1
  local_mocked_bindings(
    .emc_pooled_rows_ll = function(rows, dadm, model) {
      -0.5 * prec * (attr(dadm, "y") - rows[, "a"] - rows[, "t0"])^2
    }, .package = "EMC2")
  st <- EMC2:::.emc_pooled_state(list(mean = c(t0 = m0), var = c(t0 = v0)))
  st$cov <- matrix(0.05); st$log_eps <- log(1)
  st$shear <- matrix(-1)
  data <- lapply(seq_len(n), function(k) { d <- data.frame(x = 0); attr(d, "y") <- y[k]; d })
  s <- list(pooled = st, n_pars = 1L, n_subjects = n, par_names = "a", type = "standard",
            data = data, model = NULL,
            prior = list(theta_mu_mean = mm, theta_mu_var = matrix(vm),
                         theta_mu_invar = matrix(1 / vm)))
  ll <- function(a, phi) -0.5 * prec * (y - a - phi)^2
  a <- rep(0, n); mu <- 0
  iters <- 8000
  out <- matrix(NA_real_, iters, 3)
  for (t in seq_len(iters)) {
    phi <- s$pooled$value[[1]]
    pp <- prec + 1 / tau2
    a <- (prec * (y - phi) + mu / tau2) / pp + stats::rnorm(n, 0, sqrt(1 / pp))
    pm <- 1 / vm + n / tau2
    mu <- (mm / vm + sum(a) / tau2) / pm + stats::rnorm(1, 0, sqrt(1 / pm))
    pars <- list(subj_mu = matrix(mu, 1, n), tvar = matrix(tau2), tmu = mu)
    s$rng$gibbs <- NULL                    # the draws above share the global stream
    mv <- EMC2:::.emc_pooled_move(s, rbind(a, ll(a, phi)), "sample", pars, NULL,
                                  NULL, NULL, NULL, NULL, NULL)
    s <- mv$sampler; a <- mv$proposals[1, ]; mu <- mv$pars$tmu
    out[t, ] <- c(s$pooled$value[[1]], a[1], mu)
  }
  expect_equal(mv$proposals[2, ], ll(a, s$pooled$value[[1]]))
  expect_equal(mv$pars$subj_mu[1, ], rep(mu, n))
  Q <- diag(c(1 / v0 + n * prec, rep(prec + 1 / tau2, n), 1 / vm + n / tau2))
  Q[1, 2:(n + 1)] <- Q[2:(n + 1), 1] <- prec
  Q[n + 2, 2:(n + 1)] <- Q[2:(n + 1), n + 2] <- -1 / tau2
  b <- c(m0 / v0 + prec * sum(y), prec * y, mm / vm)
  V <- solve(Q); m <- drop(V %*% b)
  keep <- c(1, 2, n + 2)
  d <- out[-(1:200), ]
  sd_true <- sqrt(diag(V)[keep])
  expect_lt(max(abs(colMeans(d) - m[keep]) / (sd_true / sqrt(nrow(d) / 20))), 4)
  expect_lt(max(abs(apply(d, 2, stats::sd) / sd_true - 1)), 0.08)
  expect_gt(s$pooled$moves[["sample"]] / iters, 0.3)
})

test_that("pooled coefficients are fitted, stored and returned", {
  dat <- pooled_data(n_trials = 60)
  des <- pooled_design(dat)
  pr <- prior(des, pmean = c(t0 = log(.3)), psd = c(t0 = .5))
  emc <- make_emc(dat, des, prior_list = pr, n_chains = 1, compress = FALSE,
                  rt_resolution = NULL)
  expect_equal(emc[[1]]$par_names, c("v", "v_lMTRUE", "B", "A"))
  emc <- init_chains(emc, particles = 50, cores_for_chains = 1)
  expect_false(anyNA(emc[[1]]$samples$pooled[, 1]))
  emc <- suppressMessages(run_emc(
    emc, stage = "preburn", cores_for_chains = 1, cores_per_chain = 2,
    step_size = 15, max_tries = 1, verbose = FALSE, verboseProgress = FALSE,
    stop_criteria = list(iter = 15, max_gd = Inf, min_unique = 0, min_es = 0)))
  ch <- emc[[1]]
  j <- ch$samples$idx
  expect_equal(dim(ch$samples$pooled), c(1L, length(ch$samples$stage)))
  expect_false(anyNA(ch$samples$pooled[, seq_len(j)]))
  expect_gt(ch$pooled$tries[["preburn"]], 0)
  # Stored subject likelihoods are at the stored pooled value.
  for (s in 1:2) {
    a <- c(ch$samples$alpha[, s, j], t0 = ch$samples$pooled[, j])
    expect_equal(unname(ch$samples$subj_ll[s, j]),
                 EMC2:::calc_ll_manager(t(a), ch$data[[s]], ch$model)[[1]],
                 tolerance = 1e-10)
  }
  mu <- get_pars(emc, selection = "mu", stage = "preburn", return_mcmc = FALSE,
                 merge_chains = TRUE)
  expect_equal(rownames(mu), c("v", "v_lMTRUE", "B", "A", "t0"))
  al <- get_pars(emc, selection = "alpha", stage = "preburn", return_mcmc = FALSE,
                 merge_chains = TRUE)
  expect_equal(dimnames(al)[[1]], rownames(mu))
  expect_true(all(apply(al["t0", , ], 2, function(x) length(unique(x)) == 1L)))
  expect_equal(unname(al["t0", 1, ]), unname(mu["t0", ]))
  # The stored value is restored on resumption, even after the data cache moves.
  ch$data <- EMC2:::.emc_pooled_set_data(ch$data, c(t0 = 5))
  ch <- EMC2:::.emc_pooled_sync(ch)
  expect_equal(unname(attr(ch$data[[1]], "pooled_values")), unname(ch$samples$pooled[, j]))
  expect_message(ic <- compare(list(m = emc), stage = "preburn", print_summary = FALSE),
                 "Bayes factors are not available")
  expect_true(all(is.finite(unlist(ic[, c("DIC", "BPIC")]))))
  mapped <- get_pars(emc, selection = "alpha", stage = "preburn", map = TRUE,
                     return_mcmc = FALSE, merge_chains = TRUE)
  expect_true("t0" %in% dimnames(mapped)[[1]])
  expect_error(run_bridge_sampling(emc), "pooled")
})

test_that("pooled parameters refuse unsupported hierarchies", {
  dat <- pooled_data(n_trials = 20)
  des <- pooled_design(dat)
  expect_error(make_emc(dat, des, nuisance = "A", compress = FALSE), "nuisance")
})

test_that(".emc_pooled_mu_prec inverts matrix variance when theta_mu_invar is NULL", {
  V <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  pr <- list(theta_mu_var = V, theta_mu_invar = NULL)
  prec <- EMC2:::.emc_pooled_mu_prec(pr, 2)
  expect_equal(prec, solve(V), tolerance = 1e-12)
})

test_that("mode_start and laplace_refresh skip when pooled parameters are present", {
  dat <- pooled_data(n_trials = 20)
  des <- pooled_design(dat)
  pr <- prior(des, pmean = c(t0 = log(.3)), psd = c(t0 = .5))
  emc <- make_emc(dat, des, prior_list = pr, n_chains = 1, compress = FALSE,
                  rt_resolution = NULL)
  expect_identical(EMC2:::.emc_burn_mode_start(emc, n_cores = 1), emc)
  expect_identical(EMC2:::.emc_laplace_refresh(emc, n_cores = 1), emc)
})


test_that("pooled candidate likelihoods follow the failure policy", {
  rows <- matrix(c(0, 1, 2, 3), 4, 1, dimnames = list(NULL, "p"))
  ev <- function(f) {
    local_mocked_bindings(calc_ll_manager = function(proposals, ...) f(proposals),
                          .package = "EMC2")
    EMC2:::.emc_pooled_rows_ll(rows, NULL, NULL)
  }
  n_fail <- function(ll) length(attr(ll, "emc_failures"))
  withr::local_options(list(emc2.failure_policy = "report"))
  # Ordinary support values pass through with nothing recorded.
  ll <- ev(function(r) c(-1, -Inf, -2, -3)[seq_len(nrow(r))])
  expect_equal(as.numeric(ll), c(-1, -Inf, -2, -3)); expect_equal(n_fail(ll), 0L)
  # One bad row in a batch: retried row by row, the others survive.
  ll <- ev(function(r) if (nrow(r) > 1L || r[1, 1] == 2) stop("not positive definite") else -r[1, 1])
  expect_equal(as.numeric(ll), c(0, -1, -Inf, -3)); expect_equal(n_fail(ll), 2L)
  # NA, NaN and +Inf are recorded failures, not silent -Inf.
  ll <- ev(function(r) c(-1, NA, NaN, Inf))
  expect_equal(as.numeric(ll), c(-1, -Inf, -Inf, -Inf)); expect_equal(n_fail(ll), 1L)
  # A malformed result length is a programming failure.
  ll <- ev(function(r) if (nrow(r) > 1L) c(-1, -2) else -r[1, 1])
  expect_equal(as.numeric(ll), c(0, -1, -2, -3))
  expect_equal(attr(ll, "emc_failures")[[1]]$emc_class, "programming")
  # Counting happens in the master, under the "pooled" source.
  EMC2:::.emc_reject_reset("pooled")
  expect_warning(out <- EMC2:::.emc_pooled_record(
    ev(function(r) stop("object 'zz' not found"))), "pooled")
  expect_null(attr(out, "emc_failures"))
  expect_equal(as.numeric(out), rep(-Inf, 4))
  expect_equal(EMC2:::.emc_reject_counts("pooled")[["programming"]], 5L)
  # Strict: non-numerical failures abort, numerical ones still reject.
  withr::local_options(list(emc2.failure_policy = "strict"))
  expect_error(ev(function(r) stop("object 'zz' not found")), class = "emc_failure")
  expect_error(ev(function(r) c(-1, -2)), class = "emc_failure")
  expect_equal(as.numeric(ev(function(r) stop("not positive definite"))), rep(-Inf, 4))
  EMC2:::.emc_reject_reset("pooled")
})
