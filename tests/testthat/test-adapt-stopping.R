test_that("adapt uses one min_unique default whether criteria are omitted or partial", {
  omitted <- EMC2:::get_stop_criteria("adapt", NULL, "standard")
  partial <- EMC2:::get_stop_criteria("adapt", list(iter = 100), "standard")
  expect_identical(omitted$min_unique, partial$min_unique)
})

test_that("test_adapted records why adaptation is not finished", {
  sampler <- list(n_pars = 1L, n_subjects = 2L, nuisance = FALSE, type = "standard",
                  data = structure(list(), components = 1L))
  few <- list(alpha = array(c(1, 1, 2, 2), c(1, 2, 2)))
  out <- EMC2:::test_adapted(sampler, few, min_unique = 5)
  expect_false(isTRUE(out))
  expect_identical(attr(out, "reason"), "min_unique")

  many <- list(alpha = array(rnorm(40), c(1, 2, 20)))
  out <- testthat::with_mocked_bindings(
    EMC2:::test_adapted(sampler, many, min_unique = 5),
    get_conditionals = function(...) stop("singular conditional"), .package = "EMC2")
  expect_false(isTRUE(out))
  expect_match(attr(out, "reason"), "singular conditional")
})

test_that("adaptation that cannot build proposals stops at max_tries instead of sampling", {
  emc <- structure(list(list(n_pars = 1L, n_subjects = 2L, nuisance = FALSE,
                             type = "standard",
                             data = structure(list(), components = 1L),
                             samples = list(stage = rep("adapt", 20L), idx = 20L))),
                   class = "emc")
  testthat::local_mocked_bindings(
    chain_n = function(emc) matrix(20L, 1L, 4L,
                                   dimnames = list(NULL, c("preburn", "burn", "adapt", "sample"))),
    merge_chains = function(emc) emc[[1]],
    extract_samples = function(...) list(alpha = array(rnorm(40), c(1, 2, 20))),
    get_conditionals = function(...) stop("singular conditional"),
    gd_summary = function(...) numeric(0),
    .package = "EMC2")
  crit <- EMC2:::get_stop_criteria("adapt", list(iter = 20L, min_unique = 5L), "standard")
  expect_error(
    EMC2:::check_progress(emc, "adapt", iter = 20L, stop_criteria = crit, max_tries = 1L,
                          step_size = 10L, n_cores = 1L, verbose = FALSE,
                          progress = list(iters_total = 20L, trys = 1L), n_blocks = 1L),
    "could not be built.*singular conditional")
})
