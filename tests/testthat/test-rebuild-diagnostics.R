local_rng_guard()

test_that("parallel diagnostics preserve grouped results and caller RNG", {
  skip_on_os("windows")
  set.seed(143)
  groups <- setNames(lapply(1:8, function(i) coda::mcmc.list(lapply(1:3,
    function(j) {
      m <- matrix(rnorm(6000L * 8L), 6000L, 8L)
      colnames(m) <- paste0("p", 1:8)
      coda::mcmc(m)
    }))),
    paste0("group", 1:8))
  local_mocked_bindings(get_pars = function(...) groups, .package = "EMC2")
  emc <- list(list(n_subjects = 8L))
  for (fun in list(r_hat, n_eff)) for (version in c("old", "new")) {
    expected <- get_summary_stat(emc, "alpha", fun, version = version, digits = 15,
                                 n_cores = 1L)
    seed <- .Random.seed
    actual <- get_summary_stat(emc, "alpha", fun, version = version, digits = 15,
                               n_cores = 2L)
    expect_identical(actual, expected)
    expect_identical(.Random.seed, seed)
  }
  expect_identical(
    get_summary_stat(emc, "alpha", r_hat, omit_mpsrf = FALSE, n_cores = 2L, digits = 15),
    get_summary_stat(emc, "alpha", r_hat, omit_mpsrf = FALSE, n_cores = 1L, digits = 15))

})

test_that("small summaries and custom functions keep the serial path", {
  groups <- setNames(rep(list(coda::mcmc.list(coda::mcmc(matrix(1:60, 20L)))),
                         8L), paste0("group", 1:8))
  local_mocked_bindings(
    get_pars = function(...) groups,
    .emc_rebuild_lapply = function(...) stop("unexpected fork"), .package = "EMC2")
  emc <- list(list(n_subjects = 8L))
  expect_identical(get_summary_stat(emc, "alpha", n_eff, n_cores = 9L),
                   get_summary_stat(emc, "alpha", n_eff, n_cores = 1L))
  fun <- function(x, n_cores) n_cores
  expect_identical(get_summary_stat(emc, "alpha", fun, n_cores = c(8, 9)),
                   setNames(rep(list(c(8, 9)), 8L), names(groups)))
})
