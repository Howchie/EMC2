test_that("single-subject chains spend their core budget on proposals", {
  expect_identical(
    EMC2:::.particle_core_budget(1L, n_cores = 4L, r_cores = 1L),
    list(subject = 1L, likelihood = 4L)
  )
  expect_identical(
    EMC2:::.particle_core_budget(1L, n_cores = 2L, r_cores = 4L),
    list(subject = 1L, likelihood = 4L)
  )
})

test_that("multi-subject chains retain outer parallelism", {
  expect_identical(
    EMC2:::.particle_core_budget(8L, n_cores = 4L, r_cores = 1L),
    list(subject = 4L, likelihood = 1L)
  )
})

test_that("nested likelihood workers fit inside the chain budget", {
  expect_identical(
    EMC2:::.particle_core_budget(8L, n_cores = 4L, r_cores = 4L,
                                 total_cores = 4L),
    list(subject = 1L, likelihood = 4L)
  )
  expect_identical(
    EMC2:::.particle_core_budget(8L, n_cores = 8L, r_cores = 4L,
                                 total_cores = 16L),
    list(subject = 4L, likelihood = 4L)
  )
  # An inner request larger than the chain allocation is capped rather than
  # multiplying past it.
  capped <- EMC2:::.particle_core_budget(8L, n_cores = 4L, r_cores = 8L,
                                         total_cores = 4L)
  expect_equal(capped$subject * capped$likelihood, 4L)
})

test_that("proposal likelihood splits never create empty workers", {
  idx <- EMC2:::.split_work_indices(5L, 4L)
  expect_equal(idx, c(1L, 1L, 2L, 3L, 4L))
  expect_true(all(tabulate(idx, nbins = 4L) > 0L))
})

test_that("log_likelihood_joint forwards r_cores to calc_ll_manager", {
  formals_joint <- formals(EMC2:::log_likelihood_joint)
  expect_true("r_cores" %in% names(formals_joint))
  expect_true("varying" %in% names(formals_joint))
  expect_identical(formals_joint$r_cores, 1)
})

test_that("the inner fan-out absorbs cores the outer split would truncate away", {
  pb <- EMC2:::.particle_core_budget
  # 8 cores at r_cores = 3 used to give 2 subject workers x 3 = 6 of 8.
  b <- pb(112, n_cores = 8, r_cores = 3, total_cores = 8)
  expect_equal(b$subject * b$likelihood, 8)
  expect_gte(b$likelihood, 3)   # r_cores stays a lower bound

  # The r_cores = 1 default must not gain inner workers it never asked for.
  b1 <- pb(3, n_cores = 8, r_cores = 1, total_cores = 8)
  expect_equal(b1$likelihood, 1)
  expect_equal(b1$subject, 3)

  # Widening never breaks the budget, at any combination.
  for (ns in c(1, 3, 5, 112)) for (nc in c(1, 4, 8, 16)) for (rc in c(1, 2, 3, 4, 32)) {
    b <- pb(ns, n_cores = nc, r_cores = rc, total_cores = nc)
    expect_lte(b$subject * b$likelihood, nc)
    expect_gte(b$subject, 1L)
    expect_gte(b$likelihood, 1L)
  }
})

test_that("multi-subject default r_cores stays a likelihood width of one across sizes and cores", {
  pb <- EMC2:::.particle_core_budget
  for (ns in c(2L, 8L)) {
    for (tc in c(1L, 4L, 8L, 16L, 32L)) {
      b <- pb(ns, n_cores = tc, r_cores = 1L, total_cores = tc, blas_threads = 1L)
      expect_identical(b$likelihood, 1L,
                       info = sprintf("n_subjects=%d total_cores=%d", ns, tc))
      expect_lte(b$subject * b$likelihood, tc)
      expect_gte(b$subject, 1L)
    }
  }
})

test_that("explicit r_cores stays honored across sizes and cores, budget never exceeded", {
  pb <- EMC2:::.particle_core_budget
  for (ns in c(1L, 2L, 8L)) {
    for (tc in c(1L, 4L, 8L, 16L, 32L)) {
      for (rc in c(1L, 2L)) {
        b <- pb(ns, n_cores = tc, r_cores = rc, total_cores = tc, blas_threads = 1L)
        expect_lte(b$subject * b$likelihood, tc)
        expect_gte(b$subject, 1L)
        expect_gte(b$likelihood, min(rc, tc))
      }
    }
  }
})

test_that("one subject still uses the available cores regardless of r_cores", {
  pb <- EMC2:::.particle_core_budget
  for (tc in c(1L, 4L, 8L, 16L, 32L)) {
    for (rc in c(1L, 2L)) {
      b <- pb(1L, n_cores = tc, r_cores = rc, total_cores = tc, blas_threads = 1L)
      expect_identical(b$subject, 1L)
      expect_identical(b$likelihood, min(tc, max(rc, tc)))
    }
  }
})

test_that("BLAS width cannot make the worker allocation exceed the chain budget", {
  pb <- EMC2:::.particle_core_budget
  for (ns in c(1L, 2L, 8L)) {
    for (tc in c(1L, 4L, 8L, 16L, 32L)) {
      for (bt in c(1L, 2L, 4L, 8L, 64L)) {
        b <- pb(ns, n_cores = tc, r_cores = 2L,
                total_cores = tc, blas_threads = bt,
                allow_spare = TRUE)
        # A BLAS team wider than the chain is capped to one team per chain;
        # the allocator must never create another process on top of it.
        effective_bt <- min(tc, bt)
        expect_lte(b$subject * b$likelihood * effective_bt, tc,
                   label = sprintf("subjects=%d cores=%d blas=%d", ns, tc, bt))
        expect_gte(b$subject, 1L)
        expect_gte(b$likelihood, 1L)
      }
    }
  }
})

test_that(".emc_ll_route reports serial or nested_per_call, never a silent persistent_pool on fallback", {
  route <- EMC2:::.emc_ll_route
  budget1 <- list(subject = 2L, likelihood = 1L)
  budget_wide <- list(subject = 2L, likelihood = 4L)

  # No pool at all (e.g. Windows / no FIFO support): the mcmapply fallback's
  # own core budget decides between serial and a fresh per-call fork.
  expect_identical(route(NULL, NULL, budget1, budget1), "serial")
  expect_identical(route(NULL, NULL, budget1, budget_wide), "nested_per_call")

  # A pool object exists but never came alive (single worker, or degraded to
  # the master-only fallback): this must never read as persistent_pool use,
  # only as serial or a per-call fork sized by the pool's own budget.
  dead_pool <- list(n = 1L, alive = FALSE)
  expect_identical(route(dead_pool, NULL, budget1, NULL), "serial")
  expect_identical(route(dead_pool, NULL, budget_wide, NULL), "nested_per_call")

  # A genuinely running persistent pool, subject-level or the single-subject
  # likelihood pool, is reported as such regardless of the fallback budget.
  live_pool <- list(n = 4L, alive = TRUE)
  expect_identical(route(live_pool, NULL, budget1, NULL), "persistent_pool")
  expect_identical(route(NULL, live_pool, budget1, NULL), "persistent_pool")
  expect_identical(route(dead_pool, live_pool, budget1, NULL), "persistent_pool")
  # A live handle can be stale after a request falls back in the master.  The
  # recorded arm, when available, is the diagnostic source of truth.
  expect_identical(route(NULL, live_pool, budget1, NULL,
                         list(last_route = "serial")), "serial")
  expect_identical(route(NULL, live_pool, budget1, NULL,
                         list(last_route = "persistent_pool")), "persistent_pool")
})
