local_rng_guard()

test_that("an uneven block donates workers to its active chain and reclaims them", {
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")))
  work_dir <- tempfile("emc_donation_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  ready <- file.path(work_dir, "slow_ready")

  runner <- function(sampler, core_ctl, request_donation = TRUE,
                     ready_file = NULL) {
    pool <- .emc_pool_state$subject_pool
    ctx <- list(tag = "uneven-block-donation")
    if (is.null(pool)) {
      pool <- EMC2:::.emc_wpool_start(
        2L, ctx, backend = "fork", retain_across_blocks = TRUE)
      if (is.null(pool)) stop("could not start the chain's base workers")
      .emc_pool_state$subject_pool <- pool
    }
    sampler$workers_at_entry <- pool$n
    sampler$base_pids <- vapply(pool$jobs, EMC2:::.emc_wpool_job_pid,
                                integer(1))

    if (request_donation && identical(sampler$role, "fast")) {
      deadline <- Sys.time() + 5
      while (!file.exists(ready_file) && Sys.time() < deadline) Sys.sleep(0.01)
      if (!file.exists(ready_file)) stop("slow chain did not start")
      Sys.sleep(0.15)
    } else if (request_donation) {
      file.create(ready_file)
      deadline <- Sys.time() + 5
      target <- EMC2:::.emc_cores_now(core_ctl, 2L)
      while (target <= pool$n && Sys.time() < deadline) {
        Sys.sleep(0.01)
        target <- EMC2:::.emc_cores_now(core_ctl, 2L)
      }
      if (target <= pool$n) stop("no chain released capacity during the block")

      pool <- EMC2:::.emc_wpool_grow(pool, target, ctx)
      .emc_pool_state$subject_pool <- pool
      sampler$workers_after_donation <- pool$n
      sampler$donated_pids <- vapply(pool$jobs[3:4],
                                     EMC2:::.emc_wpool_job_pid, integer(1))
    }

    # run_stages() does this on exit. It first shrinks retained workers to the
    # chain's own share, then publishes completion for the current block.
    EMC2:::.emc_core_release(core_ctl)
    pool_after <- .emc_pool_state$subject_pool
    sampler$workers_after_release <- pool_after$n
    sampler$base_pids_after_release <- vapply(
      pool_after$jobs, EMC2:::.emc_wpool_job_pid, integer(1))
    sampler$donated_alive_after_release <- if (is.null(sampler$donated_pids)) {
      logical(0)
    } else vapply(sampler$donated_pids, EMC2:::.emc_wpool_pid_alive, logical(1))
    sampler
  }

  samplers <- lapply(c("fast", "slow"), function(role) structure(
    list(data = role, model = identity, prior = list(), role = role),
    class = "pmwgs"))
  pool <- EMC2:::.emc_chain_process_pool_start(samplers, runner)
  skip_if(is.null(pool), "could not start persistent chain processes")
  on.exit(EMC2:::.emc_chain_process_pool_stop(pool), add = TRUE)

  ctl <- EMC2:::.emc_core_ctl(2L, 2L, cores_for_chains = 2L)
  on.exit(unlink(ctl$dir, recursive = TRUE), add = TRUE)
  first <- EMC2:::.emc_chain_process_pool_run(
    pool, samplers, args = list(core_ctl = ctl, ready_file = ready),
    max_active = 2L)

  slow <- first[[2L]]
  expect_identical(slow$workers_at_entry, 2L)
  expect_identical(slow$workers_after_donation, 4L)
  expect_identical(slow$workers_after_release, 2L)
  expect_identical(slow$base_pids_after_release, slow$base_pids)
  expect_identical(slow$donated_alive_after_release, c(FALSE, FALSE))
  expect_identical(first[[1L]]$workers_at_entry, 2L)
  expect_identical(first[[1L]]$workers_after_release, 2L)

  # The same chain processes start the next block with only their base workers;
  # donated workers do not leak into the new block's core allocation.
  next_ctl <- EMC2:::.emc_core_ctl(2L, 2L, cores_for_chains = 2L)
  on.exit(unlink(next_ctl$dir, recursive = TRUE), add = TRUE)
  second <- EMC2:::.emc_chain_process_pool_run(
    pool, first,
    args = list(core_ctl = next_ctl, request_donation = FALSE), max_active = 2L)
  expect_identical(vapply(second, `[[`, integer(1), "workers_at_entry"),
                   c(2L, 2L))
  expect_true(EMC2:::.emc_chain_process_pool_stop(pool))
  on.exit(NULL)
})
