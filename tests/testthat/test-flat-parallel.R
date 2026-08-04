# Tests for the flattened (chain x subject) scheduler and its worker pool.
#
# The scheduler must not change what the sampler computes, only where and when
# it computes it.  The strongest available check is that a fit is identical
# whether it runs serially or over a pool, and regardless of how many workers
# the pool has -- a property the legacy nested-fork path does not have.

test_that("task chunking preserves every task and balances the deal", {
  tasks <- lapply(1:10, function(i) list(id = i))

  # One chunk when there is nothing to balance.
  expect_equal(EMC2:::.emc_chunk_tasks(tasks, n_workers = 1L), list(tasks))

  ch <- EMC2:::.emc_chunk_tasks(tasks, n_workers = 2L, chunks_per_worker = 2L)
  expect_length(ch, 4)
  ids <- vapply(unlist(ch, recursive = FALSE), `[[`, numeric(1), "id")
  expect_setequal(ids, 1:10)
  # Round-robin over a longest-first list: chunk sizes differ by at most one.
  expect_lte(diff(range(lengths(ch))), 1)

  # Never more chunks than tasks.
  expect_length(EMC2:::.emc_chunk_tasks(tasks[1:2], n_workers = 8L), 2)
})

test_that("subjects are ordered longest-processing-time first", {
  data <- list(a = data.frame(x = 1:5), b = data.frame(x = 1:50), c = data.frame(x = 1:20))
  expect_equal(EMC2:::.emc_subject_order(data), c(2L, 3L, 1L))
})

test_that("RNG streams are distinct, reproducible and leave the master ahead", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(77)
  s1 <- EMC2:::.emc_init_streams(2, 3)
  set.seed(77)
  s2 <- EMC2:::.emc_init_streams(2, 3)
  expect_identical(s1, s2)

  flat <- unlist(s1, recursive = FALSE)
  expect_length(flat, 6)
  # Every (chain, subject) gets its own stream.
  expect_length(unique(lapply(flat, paste, collapse = "-")), 6)

  # The master's own stream is moved past every stream handed out.
  set.seed(77)
  EMC2:::.emc_init_streams(2, 3)
  master <- get(".Random.seed", envir = globalenv())
  expect_false(any(vapply(flat, identical, logical(1), master)))
})

test_that("serial dispatch restores the master's RNG stream", {
  # Tasks install their own stream into .Random.seed.  In the serial path that
  # happens in the master process, so it must be put back or the Gibbs steps
  # would silently inherit a task's stream.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  clobber <- function(task) {
    assign(".Random.seed", task$seed, envir = globalenv())
    stats::runif(1)
  }
  streams <- EMC2:::.emc_init_streams(1, 2)
  assign(".Random.seed", before, envir = globalenv())
  out <- EMC2:::.emc_pool_apply(NULL, list(list(seed = streams[[1]][[1]]),
                                           list(seed = streams[[1]][[2]])), clobber)
  expect_length(out, 2)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("flat scheduling is independent of worker count and matches serial", {
  skip_on_cran()
  skip_on_os("windows")
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)

  run <- function(n_workers) {
    withr::local_options(list(emc2.flat_parallel = TRUE))
    RNGkind("L'Ecuyer-CMRG")
    set.seed(4321)
    des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
                  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                  contrasts = list(m = list(lM = ADmat)))
    emc <- suppressMessages(make_emc(dat, des, rt_resolution = 0.05, n_chains = 2))
    for (st in c("preburn", "burn")) {
      emc <- suppressMessages(run_emc(
        emc, stage = st,
        stop_criteria = list(iter = if (st == "preburn") 10 else 15,
                             max_gd = Inf, min_unique = 0, min_es = 0),
        cores_for_chains = n_workers, cores_per_chain = 1, verbose = FALSE,
        particle_factor = 20, step_size = 15, max_tries = 1))
    }
    idx <- emc[[1]]$samples$idx
    list(mu = emc[[1]]$samples$theta_mu[, idx],
         alpha = emc[[1]]$samples$alpha[, , idx],
         ll = emc[[1]]$samples$subj_ll[, idx])
  }

  serial <- run(1)
  pooled <- run(2)
  expect_equal(pooled$mu, serial$mu)
  expect_equal(pooled$alpha, serial$alpha)
  expect_equal(pooled$ll, serial$ll)
  expect_true(all(is.finite(serial$mu)))
})
