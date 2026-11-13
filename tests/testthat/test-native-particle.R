test_that("native single-subject sample blocks reproduce the R particle step", {
  skip_on_cran()
  old <- options(emc2.fit_host = FALSE)
  on.exit(options(old), add = TRUE)
  des <- design(factors = list(subjects = 1, S = 1), Rlevels = 1, model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  set.seed(11)
  dat <- make_data(c(mu = log(0.4), sigma = log(0.05), tau = log(0.2)), des,
                   n_trials = 300)
  emc <- make_emc(dat, des, type = "single", rt_resolution = 0.001)
  emc <- suppressWarnings(fit(emc, verbose = FALSE, cores_for_chains = 1,
                              stop_criteria = list(
                                preburn = list(iter = 20),
                                burn = list(iter = 60),
                                adapt = list(iter = 60),
                                sample = list(iter = 20))))
  run <- function(native) {
    options(emc2.native_particle_block = native)
    set.seed(5)
    EMC2:::run_emc(emc, stage = "sample", stop_criteria = list(iter = 150),
                   step_size = 150, verbose = FALSE, cores_for_chains = 1)
  }
  r_path <- run(FALSE)
  native <- run(TRUE)
  options(emc2.native_particle_block = NULL)
  for (j in seq_along(r_path)) {
    a <- r_path[[j]]$samples
    b <- native[[j]]$samples
    expect_identical(b$alpha, a$alpha)
    expect_identical(b$subj_ll, a$subj_ll)
    expect_identical(b$idx, a$idx)
    expect_identical(attr(b, "pm_settings"), attr(a, "pm_settings"))
  }
  # The chain moved, so the comparison covers accepted steps.
  expect_gt(length(unique(native[[1]]$samples$alpha[1, 1, ])), 50)
})

test_that("native preburn, burn and adapt blocks reproduce the R updates", {
  skip_on_cran()
  old <- options(emc2.fit_host = FALSE)
  on.exit(options(old), add = TRUE)
  des <- design(factors = list(subjects = 1, S = 1), Rlevels = 1, model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  set.seed(12)
  dat <- make_data(c(mu = log(0.4), sigma = log(0.05), tau = log(0.2)), des,
                   n_trials = 300)
  emc <- make_emc(dat, des, type = "single", rt_resolution = 0.001)
  # The first block's start is drawn outside the particle step; later blocks
  # are compared from a common state.
  options(emc2.native_particle_block = FALSE)
  emc <- suppressWarnings(EMC2:::run_emc(emc, stage = "preburn",
                                         stop_criteria = list(iter = 5),
                                         verbose = FALSE, cores_for_chains = 1))
  key <- function(e) lapply(e, function(x) list(
    alpha = x$samples$alpha, ll = x$samples$subj_ll,
    pm = attr(x$samples, "pm_settings"), rng = x$rng))
  for (st in c("preburn", "burn", "adapt")) {
    go <- function(native) {
      options(emc2.native_particle_block = native)
      set.seed(13)
      suppressWarnings(EMC2:::run_emc(emc, stage = st,
                                      stop_criteria = list(iter = 40),
                                      step_size = 40, verbose = FALSE,
                                      cores_for_chains = 1))
    }
    r_path <- go(FALSE)
    native <- go(TRUE)
    expect_identical(key(native), key(r_path), label = st)
    emc <- r_path
  }
  options(emc2.native_particle_block = NULL)
  pm <- attr(emc[[1]]$samples, "pm_settings")[[1]][[1]]
  # The comparison covered the adaptation itself, past its warm-up.
  expect_gt(pm$iter, 25)
  expect_false(is.null(pm$step_layout))
})

test_that("stored likelihood arguments reproduce the full likelihood path", {
  skip_on_cran()
  des <- design(factors = list(subjects = 1:2, S = 1), Rlevels = 1, model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  set.seed(21)
  dat <- make_data(c(mu = log(0.4), sigma = log(0.05), tau = log(0.2)), des,
                   n_trials = 80)
  emc <- make_emc(dat, des, rt_resolution = 0.001)
  data <- EMC2:::.emc_prepare_sampling_ll_data(emc[[1]]$data)
  model <- emc[[1]]$model
  P <- cbind(mu = log(0.4) + rnorm(7, 0, 0.1), sigma = log(0.05) + rnorm(7, 0, 0.1),
             tau = log(0.2) + rnorm(7, 0, 0.1))
  EMC2:::emc2_ll_cache_clear()
  full <- lapply(data, function(d) EMC2:::calc_ll_manager(P, d, model))
  expect_null(EMC2:::emc2_ll_cached(P, data[[1]], model, NULL))
  first <- lapply(data, function(d)
    EMC2:::calc_ll_manager(P, d, model, ll_cache_validated = TRUE))
  stored <- lapply(data, function(d)
    EMC2:::calc_ll_manager(P, d, model, ll_cache_validated = TRUE))
  expect_identical(first, full)
  expect_identical(stored, full)
  expect_identical(EMC2:::emc2_ll_cached(P, data[[2]], model, NULL), full[[2]])
  expect_identical(EMC2:::emc2_ll_cached_subjects(list(P, P[1:3, ]), data, model),
                   list(full[[1]], full[[2]][1:3]))
  # Repeated calls of either size reuse one parameter mapper.
  for (rows in list(1L, 1:7, 4L, 2:3, 1L)) {
    expect_identical(
      EMC2:::calc_ll_manager(P[rows, , drop = FALSE], data[[1]], model,
                             ll_cache_validated = TRUE),
      full[[1]][rows])
  }
  cand <- list(matrix(1:6 + 0.5, 3, 2), matrix(7:12 + 0.5, 3, 2))
  expect_identical(
    EMC2:::emc_candidate_rows_cpp(cand, colnames(P)),
    lapply(1:2, function(s) matrix(vapply(cand, function(a) a[, s], numeric(3)), 2,
                                   byrow = TRUE, dimnames = list(NULL, colnames(P)))))
  # Other parameter columns or another model object are not served.
  expect_null(EMC2:::emc2_ll_cached(P[, 3:1], data[[1]], model, NULL))
  other <- model
  environment(other) <- new.env(parent = environment(model))
  expect_null(EMC2:::emc2_ll_cached(P, data[[1]], other, NULL))
  expect_identical(EMC2:::emc2_ll_cached_subjects(list(P, NULL), data, model),
                   list(full[[1]], NULL))
})

test_that("native subject HMC transitions reproduce the R transitions", {
  skip_on_cran()
  des <- design(factors = list(subjects = 1:3, S = 1), Rlevels = 1, model = REXG,
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1))
  set.seed(22)
  dat <- make_data(c(mu = log(0.4), sigma = log(0.05), tau = log(0.2)), des,
                   n_trials = 120)
  emc <- make_emc(dat, des, rt_resolution = 0.001)
  data <- EMC2:::.emc_prepare_sampling_ll_data(emc[[1]]$data)
  model <- emc[[1]]$model
  pn <- emc[[1]]$par_names
  p <- length(pn)
  centre <- c(log(0.4), log(0.05), log(0.2))
  alpha <- matrix(centre + rnorm(3 * p, 0, 0.03), p, 3,
                  dimnames = list(pn, as.character(1:3)))
  ll <- vapply(1:3, function(s) EMC2:::calc_ll_manager(
    matrix(alpha[, s], 1, dimnames = list(NULL, pn)), data[[s]], model,
    ll_cache_validated = TRUE), numeric(1))
  names(ll) <- colnames(alpha)
  # A frame near the posterior scale, and a wide one whose trajectories leave
  # the bounds.
  frames <- list(list(V = diag(0.05, p), prec = diag(c(400, 300, 200)) + 20),
                 list(V = diag(25, p), prec = diag(1e-3, p)))
  ctx_for <- function(fr) list(
    data = data, model = model, r_cores = 1L, par_names = pn,
    ll_cache_validated = TRUE,
    hmc = lapply(1:3, function(s) list(
      prec = fr$prec, eps = 0.5, rate = 1, h = 1e-4,
      floor = EMC2:::.emc_hmc_ll_floor(data[[s]]),
      W = if (s < 3) diag(c(5, 4, 3)), warmup = FALSE, n_warm = 0L)))
  strip <- function(o) {
    o$tasks <- lapply(o$tasks, function(t) { t$cpu <- NULL; t })
    o
  }
  old <- options(emc2.native_particle_block = NULL)
  on.exit(options(old), add = TRUE)
  status <- character(0)
  for (fr in frames) for (eps in c(0.1, 0.9)) for (rep in 1:4) {
    ctx <- ctx_for(fr)
    set.seed(100 * rep + 1)
    seeds <- lapply(1:3, function(s) { stats::runif(1); .Random.seed })
    msg <- list(kind = "hmc", subs = 1:3, pooled = NULL, alpha = alpha, ll = ll,
                mu = matrix(centre, p, 3), var = fr$V, seeds = seeds,
                eps = rep(eps, 3))
    native <- EMC2:::.emc_native_hmc(msg, ctx)
    expect_true(all(native$done))
    options(emc2.native_particle_block = FALSE)
    r_path <- EMC2:::.emc_hmc_compute(msg, ctx)
    options(emc2.native_particle_block = TRUE)
    fast <- EMC2:::.emc_hmc_compute(msg, ctx)
    expect_identical(strip(fast), strip(r_path))
    status <- c(status, vapply(r_path$tasks, function(t)
      if (isTRUE(t$accepted)) "accepted" else t$status, ""))
  }
  # The comparison covered accepted, rejected and aborted trajectories.
  expect_true(all(c("accepted", "ok", "abort") %in% status))
})
