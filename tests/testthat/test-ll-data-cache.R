# .is_valid_ll_cache runs on every likelihood call, and only for data with a
# non-finite trial (censoring, omissions, truncation) -- the all-finite case
# returns early.  It once re-derived the whole trial partition in a per-trial
# loop containing `j0 %in% finite_set`, which is O(n_trials^2): 0.41 s per call
# on 20,000 censored trials, more than the model likelihood it was guarding.
# These tests pin both halves of the rewrite: it must still reject every
# corruption the loop rejected, and it must stay linear.

make_censored_dadm <- function(n_trials = 200, omit_every = 7L, details = FALSE) {
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
  dadm <- EMC2:::.cache_ll_data_attrs(emc[[1]]$data[[1]], force_rebuild = TRUE)
  if (details) return(list(dadm = dadm, model = emc[[1]]$model, params = p))
  dadm
}

valid <- function(dadm) {
  EMC2:::.is_valid_ll_cache(dadm, nrow(dadm),
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
  dadm <- EMC2:::.cache_ll_data_attrs(dadm, force_rebuild = TRUE)
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

test_that("sampling setup rebuilds stale likelihood partitions and keeps list attributes", {
  dadm <- make_censored_dadm()
  first_finite <- attr(dadm, "finite_rt_unique_trial_indices")[[1L]]
  n_lR <- length(unique(as.integer(dadm[["lR"]])))
  changed <- dadm
  changed$rt[1L + first_finite * n_lR] <- Inf

  bundle <- structure(list(changed), class = "sampling_data", marker = 17L)
  prepared <- EMC2:::.emc_prepare_sampling_ll_data(bundle)

  expect_identical(attributes(prepared), attributes(bundle))
  expect_true(valid(prepared[[1L]]))
  expect_false(first_finite %in%
                 attr(prepared[[1L]], "finite_rt_unique_trial_indices"))

  finite <- dadm
  omitted <- which(!is.finite(finite$rt))
  finite$rt[omitted] <- 1
  finite$R[omitted] <- levels(finite$R)[[1L]]
  finite <- EMC2:::.cache_ll_data_attrs(finite, force_rebuild = TRUE)
  expect_true(attr(finite, "emc2_all_finite_trials"))
  finite$rt[[1L]] <- Inf
  refreshed <- EMC2:::.emc_prepare_sampling_ll_data(finite)
  expect_false(attr(refreshed, "emc2_all_finite_trials"))
  expect_true(valid(refreshed))
})

test_that("registered likelihoods reuse block validation while custom ones do not", {
  fx <- make_censored_dadm(details = TRUE)
  proposals <- matrix(fx$params, nrow = 1L,
                      dimnames = list(NULL, names(fx$params)))
  expected <- EMC2:::calc_ll_manager(proposals, fx$dadm, fx$model)
  local_mocked_bindings(.cache_ll_data_attrs = function(...) {
    stop("unexpected likelihood cache validation")
  }, .package = "EMC2")
  actual <- EMC2:::calc_ll_manager(proposals, fx$dadm, fx$model,
                                   ll_cache_validated = TRUE)
  expect_identical(actual, expected)
})

test_that("Laplace likelihood batches propagate block validation", {
  seen <- logical()
  local_mocked_bindings(calc_ll_manager = function(
      proposals, dadm, model, component = NULL, r_cores = 1L,
      marginalise = NULL, varying = NULL, ll_cache_validated = FALSE) {
    seen <<- c(seen, ll_cache_validated)
    rep(-1, nrow(proposals))
  }, .package = "EMC2")

  proposal <- matrix(0, nrow = 1L, dimnames = list(NULL, "a"))
  data <- data.frame(x = 1)
  model <- identity
  EMC2:::.emc_burn_ll_batch(proposal, data, model, "a")
  EMC2:::.emc_burn_ll_batch(
    proposal, data, model, "a", ll_cache_validated = TRUE)
  expect_identical(seen, c(FALSE, TRUE))

  sampler <- list(
    samples = list(idx = 1L,
                   alpha = array(0, dim = c(1L, 1L, 1L),
                                 dimnames = list("a", NULL, NULL))),
    n_pars = 1L, type = "single",
    prior = list(theta_mu_mean = c(a = 0),
                 theta_mu_var = matrix(1, 1L, 1L,
                                       dimnames = list("a", "a"))),
    data = list(data), model = model)
  local_mocked_bindings(.emc_optimize_subject = function(
      x0, prior, ll_batch, proposal_prior_scale = 1, blocks = NULL,
      start_covariance = FALSE) {
    ll_batch(proposal)
    list(ok = TRUE)
  }, .package = "EMC2")
  EMC2:::.emc_burn_subject_mode(
    sampler, 1L, ll_cache_validated = TRUE)
  expect_true(tail(seen, 1L))
})

test_that("retained likelihood workers inherit the block validation state", {
  seen <- NULL
  local_mocked_bindings(calc_ll_manager = function(proposals, dadm, model,
                                                   component = NULL, r_cores = 1L,
                                                   varying = NULL,
                                                   ll_cache_validated = FALSE) {
    seen <<- ll_cache_validated
    rep(-1, nrow(proposals))
  }, .package = "EMC2")
  ctx <- list(data = list(data.frame(x = 1)), model = identity,
              ll_cache_validated = TRUE)
  msg <- list(kind = "ll", proposals = matrix(0, nrow = 2L, ncol = 1L), s = 1L)

  result <- EMC2:::.emc_wpool_compute(msg, ctx)

  expect_identical(result$ll, c(-1, -1))
  expect_true(seen)
})

test_that("custom R likelihoods keep cache validation inside the manager", {
  dadm <- make_censored_dadm()
  original_cache <- EMC2:::.cache_ll_data_attrs
  cache_calls <- 0L
  local_mocked_bindings(.cache_ll_data_attrs = function(dadm, force_rebuild = FALSE) {
    cache_calls <<- cache_calls + 1L
    original_cache(dadm, force_rebuild = force_rebuild)
  }, .package = "EMC2")

  custom <- function() list(c_name = NULL, log_likelihood = function(pars, dadm, model) {
    -sum(pars^2)
  })
  out <- EMC2:::calc_ll_manager(
    matrix(c(0.2, -0.4), nrow = 1L, dimnames = list(NULL, c("a", "b"))),
    dadm, custom, ll_cache_validated = TRUE)

  expect_equal(cache_calls, 1L)
  expect_length(out, 1L)
})
