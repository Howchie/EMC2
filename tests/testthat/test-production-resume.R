test_that("a production run split at a checkpoint matches an unsplit one", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  withr::local_seed(94)
  withr::local_options(list(emc2.worker_backend = "fork"))
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
  resumed <- EMC2:::restore_duplicates(run(checkpoint, "sample", 6L, cores = 2L))
  straight <- EMC2:::restore_duplicates(straight)
  for (j in seq_along(straight)) {
    for (field in c("alpha", "theta_mu", "theta_var", "subj_ll"))
      expect_identical(resumed[[j]]$samples[[field]], straight[[j]]$samples[[field]])
    for (field in c("rng", "shear", "sample_proposals_next", "sample_transitions"))
      expect_identical(resumed[[j]][[field]], straight[[j]][[field]])
    # Twelve production transitions stay before the first refresh.
    expect_identical(resumed[[j]]$sample_proposals_next, 100L)
  }
})
