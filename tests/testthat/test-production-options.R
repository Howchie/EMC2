test_that("production checkpoints retain transition options independently of session options", {
  withr::local_options(list(emc2.group_move = "shear", emc2.ensemble_pool = 12L,
                           emc2.sample_rebuilds = numeric(0)))
  emc <- list(list(n_pars = 2L), list(n_pars = 2L))
  start <- EMC2:::.emc_record_production_options(emc)
  expect_identical(start$options$emc2.group_move, "shear")
  expect_identical(start$emc[[1]]$production_options, start$emc[[2]]$production_options)
  serialized <- unserialize(serialize(start$emc, NULL))
  withr::local_options(list(emc2.group_move = "ensemble", emc2.ensemble_pool = 64L))
  resumed <- EMC2:::.emc_record_production_options(serialized)
  expect_identical(resumed$options, start$options)
  expect_identical(getOption("emc2.group_move"), "ensemble")
  expect_false(any(c("emc2.worker_backend", "emc2.worker_queue", "emc2.worker_recycle") %in%
                     names(resumed$options)))
  withr::local_options(resumed$options)
  expect_true(EMC2:::.emc_shear_options()$enabled)
  expect_identical(EMC2:::.emc_ensemble_options()$pool, 12L)
})

test_that("frozen shear production resumes exactly under changed session options", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  withr::local_seed(94)
  withr::local_options(list(emc2.group_move = "shear",
                           emc2.worker_backend = "fork", emc2.sample_rebuilds = numeric(0),
                           emc2.independent_chains = TRUE))
  local_mocked_bindings(calc_ll_manager = function(proposals, ...) {
    -rowSums(proposals^2) / 2
  }, .package = "EMC2")
  subjects <- levels(forstmann$subjects)[1:3]
  dat <- droplevels(forstmann[forstmann$subjects %in% subjects, ])
  des <- suppressMessages(design(data = dat, model = LNR, formula = list(m ~ 1, s ~ 1, t0 ~ 1)))
  warm <- suppressMessages(make_emc(dat, des, n_chains = 2L, type = "standard"))
  run <- function(emc, stage, iter, cores = 1L) suppressMessages(run_emc(emc,
    stage = stage, step_size = 6L, max_tries = 1L, verbose = FALSE, trim = FALSE,
    cores_for_chains = 1L, cores_per_chain = cores,
    stop_criteria = list(iter = iter, max_gd = Inf, selection = c("alpha", "mu"),
                         min_unique = 0, min_es = 0)))
  for (stage in c("preburn", "burn", "adapt")) warm <- run(warm, stage, 24L)
  straight <- run(warm, "sample", 12L)
  half <- run(warm, "sample", 6L)
  checkpoint <- unserialize(serialize(EMC2:::strip_duplicates(half), NULL))
  withr::local_options(list(emc2.group_move = "ensemble",
                           emc2.sample_rebuilds = 100L))
  resumed <- EMC2:::restore_duplicates(run(checkpoint, "sample", 6L, cores = 2L))
  straight <- EMC2:::restore_duplicates(straight)
  for (j in seq_along(straight)) {
    for (field in c("alpha", "theta_mu", "theta_var", "subj_ll"))
      expect_identical(resumed[[j]]$samples[[field]], straight[[j]]$samples[[field]])
    for (field in c("rng", "shear", "production_options", "sample_proposals_next", "sample_transitions"))
      expect_identical(resumed[[j]][[field]], straight[[j]][[field]])
    expect_identical(resumed[[j]]$sample_proposals_next, Inf)
  }
  expect_identical(getOption("emc2.group_move"), "ensemble")
  expect_identical(getOption("emc2.sample_rebuilds"), 100L)
})

test_that("incompatible production configurations are not combined silently", {
  withr::local_options(list(emc2.group_move = "ensemble"))
  emc <- EMC2:::.emc_record_production_options(list(list(), list()))$emc
  emc[[2]]$production_options$emc2.group_move <- "shear"
  expect_error(EMC2:::.emc_record_production_options(emc), "different recorded production")
  emc[[2]]$production_options <- NULL
  expect_error(EMC2:::.emc_record_production_options(emc), "different recorded production")
})

test_that("continuations without recorded production options report the missing configuration", {
  expect_warning(EMC2:::.emc_record_production_options(list(
    list(sample_transitions = 5L))), "continuation uses current session settings")
  expect_warning(EMC2:::.emc_record_production_options(list(
    list(samples = list(stage = "sample")))), "continuation uses current session settings")
  expect_warning(EMC2:::.emc_record_production_options(list(
    list(samples = list(stage = "adapt")))), NA)
})
