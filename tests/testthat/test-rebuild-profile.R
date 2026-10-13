test_that("proposal profiling includes allocation and every efficient stage", {
  withr::local_options(emc2.sampler_profile = TRUE)
  called <- character()
  passthrough <- function(label) function(emc, ...) {
    called <<- c(called, label)
    emc
  }
  local_mocked_bindings(
    create_chain_proposals = passthrough("chain"),
    .emc_check_likelihood_fits = passthrough("check"),
    .emc_shear_set_plans = passthrough("shear"),
    .emc_collapsed_set_plans = passthrough("collapsed"),
    create_eff_proposals = passthrough("eff"),
    .emc_alloc_rebuild = passthrough("allocation"),
    .package = "EMC2")
  emc <- list(list(samples = list(stage = rep("burn", 201), idx = 201L),
                   nuisance = TRUE))
  add_proposals(emc, "adapt", 1L, NULL)
  expect_identical(called, c("chain", "check", "shear", "collapsed", "eff", "allocation"))
  expect_true(all(is.finite(.emc_profile_state$proposal_split)))
  expect_true(all(.emc_profile_state$proposal_split >= 0))
  called <- character()
  emc[[1L]]$sample_transitions <- 100L
  emc[[1L]]$sample_proposals_next <- 200L
  add_proposals(emc, "sample", 1L, NULL)
  expect_length(called, 0L)
  expect_identical(.emc_profile_state$proposal_split, c(chain = 0, eff = 0))
})
