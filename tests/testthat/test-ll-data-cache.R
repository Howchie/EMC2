# .is_valid_ll_cache runs on every likelihood call, and only for data with a
# non-finite trial (censoring, omissions, truncation) -- the all-finite case
# returns early.  It once re-derived the whole trial partition in a per-trial
# loop containing `j0 %in% finite_set`, which is O(n_trials^2): 0.41 s per call
# on 20,000 censored trials, more than the model likelihood it was guarding.
# These tests pin both halves of the rewrite: it must still reject every
# corruption the loop rejected, and it must stay linear.

make_censored_dadm <- function(n_trials = 200, omit_every = 7L) {
  matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
  des <- design(factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = matchfun, model = LBA(),
                formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1),
                constants = c(sv = log(1)))
  p <- sampled_pars(des, doMap = FALSE)
  p["B"] <- log(1); p["t0"] <- log(0.2); p["v"] <- log(1)
  p["v_lMTRUE"] <- log(2); p["A"] <- log(0.3)
  set.seed(7)
  dat <- make_data(p, des, n_trials = n_trials)
  omit <- seq(1, nrow(dat), by = omit_every)
  dat$rt[omit] <- Inf
  dat$R[omit] <- NA
  # compress = FALSE so the row count tracks n_trials: compression collapses
  # duplicate rt/R rows, which would otherwise saturate the size sweep below.
  emc <- suppressMessages(make_emc(dat, des, type = "single", compress = FALSE))
  .cache_ll_data_attrs(emc[[1]]$data[[1]], force_rebuild = TRUE)
}

valid <- function(dadm) {
  .is_valid_ll_cache(dadm, nrow(dadm),
                     length(unique(as.integer(dadm[["lR"]]))), FALSE)
}

test_that("a freshly built ll cache validates", {
  dadm <- make_censored_dadm()
  # The early all-finite return must not be what is being exercised.
  expect_false(attr(dadm, "emc2_all_finite_trials"))
  expect_true(valid(dadm))
})

test_that("a corrupted ll cache is rejected", {
  dadm <- make_censored_dadm()

  bad <- dadm
  attr(bad, "finite_rt_mask")[which(attr(dadm, "finite_rt_mask"))[1]] <- FALSE
  expect_false(valid(bad))

  # Swap one finite trial with one non-finite trial: lengths, disjointness and
  # coverage all still hold, so only the content check can catch this.
  bad <- dadm
  fr <- attr(bad, "finite_rt_unique_trial_indices")
  ot <- attr(bad, "other_unique_trial_indices")
  attr(bad, "finite_rt_unique_trial_indices") <- c(fr[-1], ot[1])
  attr(bad, "other_unique_trial_indices") <- c(ot[-1], fr[1])
  expect_false(valid(bad))

  bad <- dadm
  attr(bad, "active_nogo_trial_mask")[1] <- !attr(bad, "active_nogo_trial_mask")[1]
  expect_false(valid(bad))

  # Data changed under a stale partition.
  bad <- dadm
  n_lR <- length(unique(as.integer(dadm[["lR"]])))
  first_finite <- attr(dadm, "finite_rt_unique_trial_indices")[1]
  bad$rt[1L + first_finite * n_lR] <- Inf
  expect_false(valid(bad))
})

test_that("the nogo branch validates and rejects", {
  dadm <- make_censored_dadm()
  levels(dadm[["lR"]]) <- c(levels(dadm[["lR"]])[1], "nogo")
  dadm <- .cache_ll_data_attrs(dadm, force_rebuild = TRUE)
  expect_true(any(attr(dadm, "active_nogo_trial_mask")))
  expect_true(valid(dadm))

  bad <- dadm
  attr(bad, "active_nogo_trial_mask")[5] <- FALSE
  expect_false(valid(bad))
})

test_that("validation is linear in trials, not quadratic", {
  small <- make_censored_dadm(n_trials = 250)
  big <- make_censored_dadm(n_trials = 4000)
  growth <- nrow(big) / nrow(small)
  expect_gt(growth, 8)
  tm <- function(dadm, reps = 20L) {
    valid(dadm)  # warm
    t1 <- proc.time()[["elapsed"]]
    for (i in seq_len(reps)) valid(dadm)
    (proc.time()[["elapsed"]] - t1) / reps
  }
  # Linear work grows like `growth`, the loop this replaced like `growth^2`.
  # The bound is generous because these are sub-millisecond timings on a shared
  # machine, but it still sits far below the quadratic cost.
  expect_lt(tm(big), 4 * growth * max(tm(small), 1e-4))
})
