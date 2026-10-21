host_available <- function() {
  .Platform$OS.type == "unix" && EMC2:::.emc_host_wanted(2)
}

test_that("host arguments are the supplied ones, defaults left to the host", {
  f <- function(emc, a = 1, b = length(emc), ...) {
    EMC2:::.emc_host_args(f, environment(), list(...))
  }
  expect_identical(f("file.RData", a = 3, z = 4), list(emc = "file.RData", a = 3, z = 4))
})

test_that("closures ship the globals they use and the packages they need", {
  env <- globalenv()
  assign(".emc_host_test_k", 7, envir = env)
  assign(".emc_host_test_g", function(x) x + .emc_host_test_k, envir = env)
  on.exit(rm(".emc_host_test_k", ".emc_host_test_g", envir = env), add = TRUE)
  f <- function(d) .emc_host_test_g(d) + stats::sd(1:3) + median(1:3)
  environment(f) <- env
  found <- EMC2:::.emc_host_globals(list(model = list(f)))
  expect_setequal(names(found$globals), c(".emc_host_test_g", ".emc_host_test_k"))
  expect_true("stats" %in% found$packages)
})

test_that("a host call returns the value, warnings, errors and the RNG state", {
  skip_on_cran()
  skip_if_not(host_available())
  set.seed(3)
  expect_identical(EMC2:::.emc_host_call("add_defaults", list(list(a = 1), b = 2)),
                   list(a = 1, b = 2))
  set.seed(3)
  in_session <- runif(2)
  set.seed(3)
  expect_identical(EMC2:::.emc_host_call("runif", list(2)), in_session)
  expect_identical(runif(1), {set.seed(3); runif(2); runif(1)})
  expect_warning(EMC2:::.emc_host_call("warning", list("from the host")),
                 "from the host")
  expect_error(EMC2:::.emc_host_call("stop", list("host failure")), "host failure")
})

test_that("the host and everything it started die with a killed session", {
  skip_on_cran()
  skip_if_not(host_available() && dir.exists("/proc"))
  tag <- paste0("emc_host_test_", Sys.getpid(), "_", sample.int(1e6, 1))
  tagged <- function() {
    n <- 0L
    for (e in Sys.glob("/proc/[0-9]*/environ")) {
      env <- tryCatch(suppressWarnings(readBin(e, "raw", 1e6)), error = function(err) raw())
      if (!length(env)) next
      env[env == as.raw(0)] <- as.raw(10)
      if (paste0("EMC2_HOST_TEST=", tag) %in% strsplit(rawToChar(env), "\n")[[1]])
        n <- n + 1L
    }
    n
  }
  lib <- dirname(getNamespaceInfo(asNamespace("EMC2"), "path"))
  script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib)),
    "f <- get('.emc_host_call', asNamespace('EMC2'))",
    # The host forks two children that sleep, as chains would.
    "f('mclapply', list(1:2, function(i) Sys.sleep(60), mc.cores = 2))"
  ), script)
  on.exit(unlink(script), add = TRUE)
  Sys.setenv(EMC2_HOST_TEST = tag)
  on.exit(Sys.unsetenv("EMC2_HOST_TEST"), add = TRUE)
  system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
          wait = FALSE, stdout = FALSE, stderr = FALSE)
  # system2(wait = FALSE) does not return the pid; find the parent by its tag.
  deadline <- Sys.time() + 30
  while (tagged() < 4L && Sys.time() < deadline) Sys.sleep(0.2)
  expect_gte(tagged(), 4L)  # parent, host and the host's two children
  parents <- system(sprintf("pgrep -f %s", shQuote(basename(script))), intern = TRUE)
  for (p in parents) tools::pskill(as.integer(p), tools::SIGKILL)
  deadline <- Sys.time() + 10
  while (tagged() > 0L && Sys.time() < deadline) Sys.sleep(0.2)
  expect_identical(tagged(), 0L)
})
