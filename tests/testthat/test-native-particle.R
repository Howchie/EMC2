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
