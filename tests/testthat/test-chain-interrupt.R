local_rng_guard()

test_that("interrupting a live persistent chain reaps its workers", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(!nzchar(Sys.which("mkfifo")))
  skip_if(!dir.exists("/proc"), "needs /proc to watch child processes")

  work_dir <- tempfile("emc_chain_interrupt_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  ready_file <- file.path(work_dir, "ready.rds")

  runner <- function(sampler, ready_file) {
    pool <- EMC2:::.emc_wpool_start(
      2L, list(tag = "interrupt-test"), backend = "fork",
      retain_across_blocks = TRUE)
    if (is.null(pool)) stop("could not start retained subject workers")
    .emc_pool_state$subject_pool <- pool
    saveRDS(list(chain = Sys.getpid(), workers = vapply(
      pool$jobs, EMC2:::.emc_wpool_job_pid, integer(1))), ready_file)
    repeat Sys.sleep(0.05)
  }

  sampler <- structure(
    list(data = list(), model = identity, prior = list()), class = "pmwgs")
  pool <- EMC2:::.emc_chain_process_pool_start(list(sampler), runner)
  skip_if(is.null(pool), "could not start persistent chain process")
  stop_pool <- function() {
    if (is.null(pool)) return(invisible(NULL))
    current <- pool
    pool <<- NULL
    EMC2:::.emc_chain_process_pool_stop(current)
  }
  on.exit(stop_pool(), add = TRUE)

  parent_pid <- Sys.getpid()
  signaler <- parallel::mcparallel({
    deadline <- Sys.time() + 10
    while (!file.exists(ready_file) && Sys.time() < deadline) Sys.sleep(0.01)
    ready <- file.exists(ready_file)
    Sys.sleep(0.15)
    alive <- if (ready) {
      children <- readRDS(ready_file)
      vapply(c(children$chain, children$workers),
             EMC2:::.emc_wpool_pid_alive, logical(1))
    } else logical()
    tools::pskill(parent_pid, 2L)
    list(ready = ready, alive = alive)
  }, detached = FALSE)
  expect_false(is.null(signaler))
  on.exit(suppressWarnings(try(
    parallel::mccollect(signaler, wait = FALSE), silent = TRUE)), add = TRUE)

  caught_interrupt <- FALSE
  run_block <- function() {
    # This is the cleanup registered by run_emc() for its chain process pool.
    on.exit(stop_pool(), add = TRUE)
    EMC2:::.emc_chain_process_pool_run(
      pool, list(sampler), args = list(ready_file = ready_file), max_active = 1L)
  }
  tryCatch(run_block(), interrupt = function(e) {
    caught_interrupt <<- TRUE
  })

  signal_result <- suppressWarnings(try(
    parallel::mccollect(signaler, wait = TRUE), silent = TRUE))
  expect_true(is.list(signal_result) && length(signal_result) == 1L)
  expect_true(isTRUE(signal_result[[1L]]$ready))
  children <- readRDS(ready_file)
  expect_true(all(signal_result[[1L]]$alive))
  expect_true(caught_interrupt)
  deadline <- Sys.time() + 10
  repeat {
    alive <- vapply(c(children$chain, children$workers),
                    EMC2:::.emc_wpool_pid_alive, logical(1))
    if (!any(alive) || Sys.time() >= deadline) break
    Sys.sleep(0.02)
  }
  expect_false(any(vapply(c(children$chain, children$workers),
                          EMC2:::.emc_wpool_pid_alive, logical(1))))
})
