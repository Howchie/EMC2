local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
# Acceptance gate 1 of the architecture-efficiency audit: "Same-build, pure
# scheduling/storage/cache changes should reproduce particle likelihood vectors
# and subject streams exactly.  Differential tests must include repeated calls,
# reordered particles, chunk boundaries, worker growth/recycling/fallback and
# restarted fits."
#
# test-worker-pool.R already covers the pool-side half of that list: a pooled
# fit repeats itself exactly, and a fit does not move when the core count does
# -- which is simultaneously the partition-invariance and the dead-pool-fallback
# test, since one core produces the synthetic dead pool that recomputes every
# subject in the master.  What is added here is the likelihood-side half, which
# nothing tested, and the restart boundary.
#
# These matter most for C11 (immutable owning native context).  The present
# mapper caches values tied to the first particle, so a context whose lifetime
# is extended across calls without explicit invalidation would give wrong
# likelihoods, and the first thing to notice would be one of these tests.

RNGkind("L'Ecuyer-CMRG")

# --- repeated and interleaved calls -----------------------------------------

test_that("repeated calls on the same fixture return the identical vector", {
  # Anything cached inside the native layer between calls shows up here first.
  for (model in c("RDM", "LBA", "DDM")) {
    fx <- audit_fixture(model, n_trials = 300L, n_particles = 6L)
    first <- audit_ll_direct(fx)
    expect_informative_ll(fx, first, model)
    for (rep in 1:3) {
      expect_bit_identical(audit_ll_direct(fx), first,
                           sprintf("%s, repeat %d", model, rep))
    }
    # And through the managed wrapper, which rebuilds the design list each call.
    expect_bit_identical(audit_ll_managed(fx), first, paste(model, "managed"))
  }
})

test_that("interleaving two fixtures does not contaminate either", {
  # Two different dadms, different designs, different widths, called
  # alternately.  A context keyed on "the last thing I saw" passes the repeated
  # call test above and fails this one.
  a <- audit_fixture("RDM", n_trials = 200L, n_particles = 4L, seed = 11L)
  b <- audit_fixture("RDM", n_trials = 300L, n_particles = 4L, width = 8L, seed = 22L)
  la <- audit_ll_direct(a)
  lb <- audit_ll_direct(b)
  expect_false(isTRUE(all.equal(length(la), 0L)))
  for (rep in 1:3) {
    expect_bit_identical(audit_ll_direct(a), la, sprintf("a after b, %d", rep))
    expect_bit_identical(audit_ll_direct(b), lb, sprintf("b after a, %d", rep))
  }
})

test_that("interleaving two different models does not contaminate either", {
  a <- audit_fixture("LBA", n_trials = 200L, n_particles = 4L)
  b <- audit_fixture("DDM", n_trials = 200L, n_particles = 4L)
  la <- audit_ll_direct(a)
  lb <- audit_ll_direct(b)
  expect_bit_identical(audit_ll_direct(a), la, "LBA after DDM")
  expect_bit_identical(audit_ll_direct(b), lb, "DDM after LBA")
  expect_bit_identical(audit_ll_direct(a), la, "LBA again")
})

# --- particle ordering and batching -----------------------------------------

test_that("reordering particles permutes the likelihood and nothing else", {
  # Every short cut in the prologue carries state forward from the template
  # particle, so the answer for particle i must not depend on which particles
  # precede it.
  for (model in c("RDM", "LBA", "DDM")) {
    fx <- audit_fixture(model, n_trials = 250L, n_particles = 8L)
    base <- audit_ll_direct(fx)
    set.seed(3)
    perm <- sample.int(nrow(fx$prop))
    shuffled <- fx
    shuffled$prop <- fx$prop[perm, , drop = FALSE]
    expect_bit_identical(audit_ll_direct(shuffled), base[perm],
                         paste(model, "permuted"))
    # Reversal too: a permutation that always moves the first particle, which
    # is the one the template lanes are primed from.
    rev_fx <- fx
    rev_fx$prop <- fx$prop[rev(seq_len(nrow(fx$prop))), , drop = FALSE]
    expect_bit_identical(audit_ll_direct(rev_fx), rev(base),
                         paste(model, "reversed"))
  }
})

test_that("chunking a particle batch does not move any value", {
  # C13 and the staged native executor both split a batch differently from the
  # way it is split today.  Whatever the chunk boundaries are, the concatenated
  # answer has to be the whole-batch answer.
  fx <- audit_fixture("RDM", n_trials = 250L, n_particles = 12L)
  base <- audit_ll_direct(fx)
  for (chunk in c(1L, 2L, 5L, 7L, 11L, 12L)) {
    pieces <- split(seq_len(nrow(fx$prop)),
                    ceiling(seq_len(nrow(fx$prop)) / chunk))
    got <- unlist(lapply(pieces, function(idx) {
      part <- fx
      part$prop <- fx$prop[idx, , drop = FALSE]
      audit_ll_direct(part)
    }), use.names = FALSE)
    expect_bit_identical(got, base, paste("chunk size", chunk))
  }
})

test_that("a single-particle call matches that particle inside a batch", {
  # The batch of one is its own boundary case: it is what the dynamic
  # likelihood queue sends, and what a fallback recomputation produces.
  for (model in c("RDM", "LBA", "DDM")) {
    fx <- audit_fixture(model, n_trials = 200L, n_particles = 5L)
    base <- audit_ll_direct(fx)
    for (i in seq_len(nrow(fx$prop))) {
      one <- fx
      one$prop <- fx$prop[i, , drop = FALSE]
      expect_bit_identical(audit_ll_direct(one), base[i],
                           sprintf("%s, particle %d alone", model, i))
    }
  }
})

test_that("the manager's inner fan-out does not move a value", {
  # calc_ll_manager() splits proposals across r_cores.  C12 adds a direct serial
  # branch to that function, so pin the split and unsplit answers together
  # before it is touched.
  skip_on_os("windows")
  fx <- audit_fixture("RDM", n_trials = 200L, n_particles = 9L)
  base <- audit_ll_managed(fx, r_cores = 1L)
  for (cores in c(2L, 3L, 4L)) {
    expect_bit_identical(audit_ll_managed(fx, r_cores = cores), base,
                         paste("r_cores =", cores))
  }
})

# --- compression ------------------------------------------------------------

test_that("compressed and uncompressed data give the same total likelihood", {
  # Compression is exact: identical trials are represented once and expanded by
  # multiplicity.  C11's note about replacing the expansion-index sum with
  # precomputed multiplicities changes the reduction order, so the equality has
  # to be pinned first.
  dat <- audit_data(300L, jitter = 0)          # no jitter: real duplicates
  packed <- audit_fixture("RDM", data = dat, n_particles = 4L)
  loose <- packed
  loose$emc <- suppressMessages(make_emc(dat, packed$des, n_chains = 1L,
                                         compress = FALSE, rt_resolution = NULL,
                                         type = "single"))
  loose$dadm <- loose$emc[[1L]]$data[[1L]]
  expect_lt(nrow(packed$dadm), nrow(loose$dadm))
  expect_equal(audit_ll_direct(packed), audit_ll_direct(loose), tolerance = 1e-10)
})

# --- subject streams --------------------------------------------------------

test_that("subject streams are independent and reproducible", {
  set.seed(99)
  a <- EMC2:::.emc_subject_streams(5L)
  set.seed(99)
  b <- EMC2:::.emc_subject_streams(5L)
  expect_identical(a, b)
  expect_length(a, 5L)
  # Distinct streams, or two subjects would draw the same proposals.
  expect_equal(length(unique(vapply(a, function(s) paste(s, collapse = ","),
                                    character(1)))), 5L)
})

test_that("block structure does not change the random streams", {
  skip_on_os("windows")
  skip_on_cran()
  dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- suppressMessages(design(data = dat, model = LNR,
                                 formula = list(m ~ 1, s ~ 1, t0 ~ 1)))
  fit <- function(step) {
    RNGkind("L'Ecuyer-CMRG"); set.seed(7)
    emc <- suppressMessages(make_emc(dat, des, n_chains = 1L, compress = TRUE))
    suppressMessages(run_emc(emc, stage = "preburn", cores_for_chains = 1,
      cores_per_chain = 2, step_size = step, max_tries = 1, verbose = FALSE,
      verboseProgress = FALSE,
      stop_criteria = list(iter = 10, max_gd = Inf, min_unique = 0, min_es = 0)))
  }
  one_block <- fit(10)
  two_blocks <- fit(5)
  expect_identical(one_block[[1]]$samples$idx, two_blocks[[1]]$samples$idx)
  expect_identical(one_block[[1]]$samples$alpha,
                   two_blocks[[1]]$samples$alpha)
  expect_identical(fit(10)[[1]]$samples$alpha, one_block[[1]]$samples$alpha)
  expect_identical(fit(5)[[1]]$samples$alpha, two_blocks[[1]]$samples$alpha)
})

# --- restarted fits ---------------------------------------------------------

test_that("a serialized and resumed fit continues bit-identically", {
  # The restart gate.  Holding the block structure fixed, a fit that is written
  # to disk, read back and continued must land exactly where the straight-
  # through run did.  C14 (appendable sample history) and C16 (worker templates
  # retained across blocks) both change what crosses this boundary.
  skip_on_os("windows")
  skip_on_cran()
  dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- suppressMessages(design(data = dat, model = LNR,
                                 formula = list(m ~ 1, s ~ 1, t0 ~ 1)))
  run <- function(emc, iter) suppressMessages(run_emc(emc, stage = "preburn",
    cores_for_chains = 1, cores_per_chain = 2, step_size = 5, max_tries = 1,
    verbose = FALSE, verboseProgress = FALSE,
    stop_criteria = list(iter = iter, max_gd = Inf, min_unique = 0, min_es = 0)))
  fresh <- function() {
    RNGkind("L'Ecuyer-CMRG"); set.seed(7)
    suppressMessages(make_emc(dat, des, n_chains = 1L, compress = TRUE))
  }
  straight <- run(fresh(), 10)
  # `stop_criteria$iter` is a target for the stage, not an increment: resuming
  # asks for 10 in total, having already stored 5.
  half <- run(fresh(), 5)
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(half, f)
  resumed <- run(readRDS(f), 10)

  expect_identical(resumed[[1]]$samples$idx, straight[[1]]$samples$idx)
  expect_identical(resumed[[1]]$samples$alpha, straight[[1]]$samples$alpha)
  expect_identical(resumed[[1]]$samples$subj_ll, straight[[1]]$samples$subj_ll)
})

# --- worker recycling -------------------------------------------------------

test_that("the recycle interval does not move a fit", {
  # Workers are re-forked periodically so they cannot drift from their shared
  # pages.  Streams live in the master, so replacing a worker must not be able
  # to move a draw -- including when it happens every single iteration.
  skip_on_os("windows")
  skip_on_cran()
  dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- suppressMessages(design(data = dat, model = LNR,
                                 formula = list(m ~ 1, s ~ 1, t0 ~ 1)))
  fit_with <- function(every) {
    op <- options(emc2.worker_recycle = every)
    on.exit(options(op), add = TRUE)
    RNGkind("L'Ecuyer-CMRG"); set.seed(5)
    emc <- suppressMessages(make_emc(dat, des, n_chains = 1L, compress = TRUE))
    suppressMessages(run_emc(emc, stage = "preburn", cores_for_chains = 1,
      cores_per_chain = 3, step_size = 8, max_tries = 1, verbose = FALSE,
      verboseProgress = FALSE,
      stop_criteria = list(iter = 8, max_gd = Inf, min_unique = 0, min_es = 0)))
  }
  every_iteration <- fit_with(1L)
  never <- fit_with(.Machine$integer.max)
  expect_identical(every_iteration[[1]]$samples$alpha, never[[1]]$samples$alpha)
  expect_identical(every_iteration[[1]]$samples$subj_ll, never[[1]]$samples$subj_ll)
})
