# pGuess: the uniform ("guess") outlier mixture.
#
# pContaminant is a Bernoulli OMISSION mixture -- it adds mass only at rt == +Inf
# and does nothing for observed RTs.  pGuess is the standard uniform-outlier
# mixture (Ratcliff & Tuerlinckx 2002; HDDM's w_outlier): it contributes a flat
# density directly to observed RT densities.  The two are nested:
#
#   P(omission) = pC,  P(guess) = (1 - pC) * pG,  P(process) = (1 - pC)(1 - pG)
#
# so at pG = 0 every expression reduces exactly to the pre-pGuess arithmetic.
# See src/contaminant_mixture.h and resolve_guess_window() in R/design.R.

lI_single <- function(d) factor(rep(1, nrow(d)), levels = 1)

# ---------------------------------------------------------------------------
# A test-local reference implementation of the two-case mixture.  The R
# likelihood paths deliberately still implement pContaminant alone, so they are
# no longer an oracle; this is.
# ---------------------------------------------------------------------------
ref_mix <- function(ll_proc, pC, pG, log_g, is_omission) {
  v <- ll_proc
  if (pG > 0) {
    v <- if (is.infinite(log_g) && log_g < 0) log(1 - pG) + v else
      log((1 - pG) * exp(v) + pG * exp(log_g))
  }
  if (pC > 0) {
    v <- if (is_omission) log(pC + (1 - pC) * exp(v)) else log(1 - pC) + v
  }
  v
}

# Per-trial log-likelihoods from the compiled kernel.
trial_lls <- function(dadm, model, p) {
  p_types <- names(model$p_types)
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  as.numeric(EMC2:::calc_ll_oo_pw(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = constants, designs = designs, type = model$c_name,
    bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = p_types,
    min_ll = log(1e-10), trend = model$trend))
}

total_ll <- function(dadm, model, p) {
  p_types <- names(model$p_types)
  designs <- list()
  for (nm in p_types) {
    designs[[nm]] <- attr(dadm, "designs")[[nm]][
      attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  as.numeric(EMC2:::calc_ll_oo(
    matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = constants, designs = designs, type = model$c_name,
    bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = p_types,
    min_ll = log(1e-10), trend = model$trend))
}


test_that("pGuess = 0 leaves every family's likelihood untouched", {
  matchfun <- function(d) d$S == d$lR
  cases <- list(
    LBA = list(model = LBA, formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1)),
    RDM = list(model = RDM, formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1)),
    LNR = list(model = LNR, formula = list(m ~ lM, s ~ 1, t0 ~ 1))
  )
  for (nm in names(cases)) {
    set.seed(7)
    des <- design(factors = list(subjects = 1, S = c("left", "right")),
                  Rlevels = c("left", "right"), matchfun = matchfun,
                  formula = cases[[nm]]$formula, model = cases[[nm]]$model)
    p <- sampled_pars(des); p[] <- 0.3
    dat <- make_data(p, des, n_trials = 25)
    dadm <- EMC2:::design_model(dat, des, compress = TRUE, rt_resolution = NULL)
    model <- attr(dadm, "model")()

    pv <- c(p, pContaminant = -Inf, pGuess = -Inf)
    pv <- pv[names(pv) %in% names(model$p_types) | names(pv) %in% names(p)]
    base <- total_ll(dadm, model, p)
    expect_true(is.finite(base), info = nm)

    # An explicitly zero pGuess must be indistinguishable from an absent one.
    des0 <- design(factors = list(subjects = 1, S = c("left", "right")),
                   Rlevels = c("left", "right"), matchfun = matchfun,
                   formula = c(cases[[nm]]$formula, list(pGuess ~ 1)),
                   model = cases[[nm]]$model)
    dadm0 <- EMC2:::design_model(dat, des0, compress = TRUE, rt_resolution = NULL)
    p0 <- c(p, pGuess = qnorm(0))
    expect_equal(total_ll(dadm0, attr(dadm0, "model")(), p0), base,
                 tolerance = 0, info = nm)
  }
})


test_that("the compiled mixture matches the reference on every trial kind", {
  # Finite RTs plus a +Inf omission, a -Inf left-censored trial and an NA.
  des <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1, pContaminant ~ 1, pGuess ~ 1),
    functions = list(lI = lI_single), model = REXG)
  dat <- data.frame(
    subjects = factor(1), S = factor(1), R = factor(1, levels = 1),
    rt = c(0.45, 0.55, 0.70, 0.95, Inf, -Inf, NA),
    LT = 0.20, LC = 0.35, UC = 2, UT = 3)
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()

  gw <- attr(dadm, "guess_window")
  n_resp <- attr(dadm, "guess_n_resp")
  expect_equal(gw, c(0.35, 2))          # max(LT, LC), min(UC, UT)
  expect_identical(n_resp, 1L)
  log_g <- -log(n_resp * diff(gw))

  base_p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18), t0 = log(0.05))
  pC <- 0.08; pG <- 0.15

  # Process-only per-trial log-likelihoods.
  ll_proc <- trial_lls(dadm, model,
                       c(base_p, pContaminant = qnorm(0), pGuess = qnorm(0)))
  got <- trial_lls(dadm, model,
                   c(base_p, pContaminant = qnorm(pC), pGuess = qnorm(pG)))

  rt <- dat$rt
  want <- vapply(seq_along(rt), function(i) {
    lg <- if (is.finite(rt[i]) && rt[i] > 0) log_g else -Inf
    ref_mix(ll_proc[i], pC, pG, lg, identical(rt[i], Inf))
  }, numeric(1))
  expect_equal(got, want, tolerance = 1e-10)

  # The point of the mixture: every observed RT now has a likelihood FLOOR at
  # the guess mass, however badly the process fits it.
  finite <- is.finite(rt) & rt > 0
  expect_true(all(got[finite] >= log((1 - pC) * pG) + log_g - 1e-12))
})


test_that("the mixture is proper: it integrates to 1 over the window", {
  des <- design(
    factors = list(subjects = 1, S = 1), Rlevels = c("a", "b"),
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, pGuess ~ 1),
    model = RDM, constants = c(s = log(1)))
  pG <- 0.2
  p <- c(v = log(2), B = log(1), A = log(0.3), t0 = log(0.2), pGuess = qnorm(pG))

  # Pin the window so the quadrature grid covers it exactly; otherwise the
  # missing slice of uniform mass shows up as a spurious propriety failure.
  des$TC <- list(guess_window = c(0, 12))
  # A small trapezoid grid is enough for the package regression contract.  The
  # dense grid used while developing the mixture belongs in the validation
  # scripts, not in every package test run.
  grid <- seq(0, 12, length.out = 241)
  dat <- do.call(rbind, lapply(c("a", "b"), function(r) data.frame(
    subjects = factor(1), S = factor(1),
    R = factor(r, levels = c("a", "b")), rt = grid)))
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()
  expect_equal(attr(dadm, "guess_window"), c(0, 12))
  expect_identical(attr(dadm, "guess_n_resp"), 2L)

  dens <- exp(trial_lls(dadm, model, p))
  # Trapezoid over the grid, summed across both responses.
  mass <- sum(vapply(split(dens, rep(1:2, each = length(grid))), function(d)
    sum(diff(grid) * (utils::head(d, -1) + utils::tail(d, -1)) / 2), numeric(1)))
  expect_equal(mass, 1, tolerance = 5e-3)
})


test_that("a two-response DDM on [0, 5] reproduces HDDM's w_outlier = 0.1", {
  des <- design(factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"),
                formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1, pGuess ~ 1),
                model = DDM, constants = c(s = log(1)))
  # A guess window of exactly [0, 5] over two responses is HDDM's setup.
  dat <- data.frame(subjects = factor(1), S = factor("left", levels = c("left", "right")),
                    R = factor(c("left", "right"), levels = c("left", "right")),
                    rt = c(0.4, 0.4), LT = 0, UT = Inf, LC = 0, UC = Inf)
  des$TC <- list(guess_window = c(0, 5))
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  expect_equal(attr(dadm, "guess_window"), c(0, 5))
  expect_identical(attr(dadm, "guess_n_resp"), 2L)
  model <- attr(dadm, "model")()

  # The guess component's per-response density is 1 / (n_resp * W) = 0.1,
  # exactly HDDM's w_outlier.  Recover it from a half-and-half mixture:
  #   d_mix = (1 - pG) * d_proc + pG * 0.1
  base <- c(v = 1, a = log(1), t0 = log(0.2), Z = qnorm(0.5))
  d_proc <- exp(trial_lls(dadm, model, c(base, pGuess = qnorm(0))))
  pG <- 0.5
  d_mix <- exp(trial_lls(dadm, model, c(base, pGuess = qnorm(pG))))
  implied <- (d_mix - (1 - pG) * d_proc) / pG
  expect_equal(implied, c(0.1, 0.1), tolerance = 1e-6)

  # w_outlier is accepted as the HDDM spelling and converts to the same window.
  des2 <- des; des2$TC <- list(w_outlier = 0.1)
  dadm2 <- EMC2:::design_model(dat, des2, compress = FALSE, rt_resolution = NULL)
  expect_equal(attr(dadm2, "guess_window"), c(0, 5))
})


test_that("a guess window narrower than the data is refused, not silently widened", {
  # HDDM's w_outlier = 0.1 means a 5 s window on two responses.  A 9 s RT sits
  # outside it, where the uniform guess density is zero -- so that trial, the
  # slow outlier pGuess exists to catch, would get no guess component while its
  # neighbours did.  This used to be applied silently: every trial, in or out of
  # the window, received density 1 / (n_resp * (UG - LG)).
  des <- design(factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"),
                formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1, pGuess ~ 1),
                model = DDM, constants = c(s = log(1)))
  dat <- data.frame(subjects = factor(1),
                    S = factor("left", levels = c("left", "right")),
                    R = factor("left", levels = c("left", "right")),
                    rt = c(0.4, 9.0), LT = 0, UT = Inf, LC = 0, UC = Inf)

  des$TC <- list(w_outlier = 0.1)
  expect_error(EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL),
               "does not cover the observed RTs")
  # The message names the value that would work.
  expect_error(EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL),
               "w_outlier <=")

  # The explicit spelling is held to the same standard ...
  des$TC <- list(guess_window = c(0, 5))
  expect_error(EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL),
               "does not cover the observed RTs")

  # ... and a window that does cover the data is accepted unchanged.
  des$TC <- list(guess_window = c(0, 10))
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  expect_equal(attr(dadm, "guess_window"), c(0, 10))

  # The derived window cannot fail the check: it is built from the same bounds
  # the RTs were already validated against, and widens past the slowest RT when
  # there is no finite upper edge.
  des$TC <- NULL
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  expect_gte(attr(dadm, "guess_window")[2], max(dat$rt))
})


test_that("expand and compressed branches agree with pGuess active", {
  # Mirrors test-pcontaminant-compressed.R: on all-unique trials the two
  # summation branches must produce the same total.
  des <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1,
    formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, t0 ~ 1, pContaminant ~ 1, pGuess ~ 1),
    functions = list(lI = lI_single), model = REXG)
  dat <- data.frame(
    subjects = factor(1), S = factor(1), R = factor(1, levels = 1),
    rt = c(0.45, 0.55, 0.70, 0.95, -Inf, Inf),
    LT = 0.20, LC = 0.35, UC = Inf, UT = Inf)
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()
  p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18), t0 = log(0.05),
         pContaminant = qnorm(0.05), pGuess = qnorm(0.2))

  expect_true(length(attr(dadm, "expand")) > 0)
  ll_expand <- total_ll(dadm, model, p)
  stripped <- dadm
  attr(stripped, "expand") <- integer(0)
  expect_true(is.finite(ll_expand))
  expect_equal(total_ll(stripped, model, p), ll_expand, tolerance = 1e-10)
})


test_that("resolve_guess_window follows its documented resolution order", {
  d <- data.frame(rt = c(0.5, 1.2, 7.3), R = factor(c("a", "b", "a")),
                  LT = 0, UT = Inf, LC = 0, UC = Inf)
  # No finite upper edge -> next whole second above max(rt), floored at 5.
  expect_equal(EMC2:::resolve_guess_window(d)$window, c(0, 8))
  d2 <- d; d2$rt <- c(0.5, 1.2, 2.3)
  expect_equal(EMC2:::resolve_guess_window(d2)$window, c(0, 5))
  # Bounds win over the fallback.
  d3 <- d; d3$LT <- 0.2; d3$LC <- 0.4; d3$UC <- 3; d3$UT <- 2.5
  expect_equal(EMC2:::resolve_guess_window(d3)$window, c(0.4, 2.5))
  # TC wins over everything.  The window still has to cover the RTs -- d3's
  # fastest is 0.5, so a lower edge of 1 is now refused rather than silently
  # leaving that trial without a guess component.
  expect_equal(EMC2:::resolve_guess_window(d3, list(guess_window = c(0.3, 9)))$window,
               c(0.3, 9))
  expect_error(EMC2:::resolve_guess_window(d3, list(guess_window = c(1, 9))),
               "does not cover the observed RTs")
  expect_error(EMC2:::resolve_guess_window(d3, list(guess_window = c(0, 5))),
               "does not cover the observed RTs")
  # nogo and time are not overt responses and do not count.
  d4 <- d; d4$R <- factor(c("a", "b", "a"), levels = c("a", "b", "nogo", "time"))
  expect_identical(EMC2:::resolve_guess_window(d4)$n_resp, 2L)
})


test_that("make_missing draws guesses inside the window, after truncation", {
  set.seed(3)
  n <- 500
  d <- data.frame(subjects = factor(1), R = factor(sample(c("a", "b"), n, TRUE)),
                  rt = runif(n, 0.3, 1.2), LT = 0, UT = Inf, LC = 0, UC = Inf)
  out <- make_missing(d, pGuess = 0.25, guess_window = c(2, 4),
                      rt_resolution = NULL)
  guessed <- out$rt > 1.5
  expect_lt(abs(mean(guessed) - 0.25), 0.08)
  expect_true(all(out$rt[guessed] >= 2 & out$rt[guessed] <= 4))
  # Nested with the omission: P(guess) = (1 - pC) * pG, and no trial is both.
  out2 <- make_missing(d, pContaminant = 0.2, pGuess = 0.25,
                       guess_window = c(2, 4), rt_resolution = NULL)
  expect_equal(mean(is.infinite(out2$rt)), 0.2, tolerance = 0.1)
  expect_lt(abs(mean(is.finite(out2$rt) & out2$rt > 1.5) - 0.8 * 0.25), 0.1)
})


test_that("pGuess is refused when the model has no compiled likelihood", {
  # A free (or non-zero) pGuess on an R-only likelihood would be sampled and
  # silently ignored; design() turns that into a loud error.
  fake <- function() {
    m <- LNR()
    m$c_name <- NULL
    m
  }
  expect_error(
    design(factors = list(subjects = 1, S = c("left", "right")),
           Rlevels = c("left", "right"),
           matchfun = function(d) d$S == d$lR,
           formula = list(m ~ lM, s ~ 1, t0 ~ 1, pGuess ~ 1), model = fake),
    "compiled likelihood")
  # A pGuess pinned to its inert default is fine.
  expect_no_error(suppressMessages(utils::capture.output(
    design(factors = list(subjects = 1, S = c("left", "right")),
           Rlevels = c("left", "right"),
           matchfun = function(d) d$S == d$lR,
           formula = list(m ~ lM, s ~ 1, t0 ~ 1), model = fake))))
})


test_that("guess_exclude validation rejects bad inputs", {
  d <- data.frame(subjects = factor(1), rt = c(0.5, 1.2),
                  R = factor(c("a", "b"), levels = c("a", "b", "c")),
                  LT = 0, UT = Inf, LC = 0, UC = Inf)

  # Non-character
  expect_error(EMC2:::resolve_guess_window(d, list(guess_exclude = 123)),
               "character vector")
  expect_error(make_missing(d, pGuess = 0.2, guess_exclude = 123),
               "character vector")

  # NAs or empty strings
  expect_error(EMC2:::resolve_guess_window(d, list(guess_exclude = c("a", NA))),
               "character vector")
  expect_error(EMC2:::resolve_guess_window(d, list(guess_exclude = c("a", ""))),
               "character vector")

  # Duplicates
  expect_error(EMC2:::resolve_guess_window(d, list(guess_exclude = c("a", "a"))),
               "duplicate")

  # Unknown response levels (nogo and time are allowed)
  expect_error(EMC2:::resolve_guess_window(d, list(guess_exclude = "unknown")),
               "not present in R")
  expect_no_error(EMC2:::resolve_guess_window(d, list(guess_exclude = "nogo")))

  # Excluding all overt response levels
  expect_error(EMC2:::resolve_guess_window(d, list(guess_exclude = c("a", "b", "c"))),
               "at least one overt response level eligible")
})


test_that("guess_exclude removes excluded levels during simulation", {
  local_test_rng(42)
  n <- 500
  d <- data.frame(subjects = factor(1),
                  R = factor(sample(c("a", "b", "c"), n, replace = TRUE),
                             levels = c("a", "b", "c")),
                  rt = runif(n, 0.4, 0.8), LT = 0, UT = Inf, LC = 0, UC = Inf)

  # pGuess = 1 with "b" excluded: all trials are guesses, but "b" is never produced
  out <- make_missing(d, pGuess = 1.0, guess_exclude = "b",
                      guess_window = c(2, 4), rt_resolution = NULL)
  expect_true(all(out$rt >= 2 & out$rt <= 4))
  expect_false(any(out$R == "b"))
  expect_true(all(out$R %in% c("a", "c")))
  # Distribution across remaining levels should be roughly equal (uniform)
  prop_a <- mean(out$R == "a")
  expect_gt(prop_a, 0.4)
  expect_lt(prop_a, 0.6)

  # make_data also inherits and respects TC$guess_exclude
  matchfun <- function(d) d$S == d$lR
  des <- design(factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = matchfun,
                formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, pGuess ~ 1),
                model = RDM, constants = c(s = log(1)),
                TC = list(guess_exclude = "right", guess_window = c(2, 4)))
  p <- c(v = log(2), B = log(1), A = log(0.3), t0 = log(0.2), pGuess = qnorm(0.999))
  sim_dat <- make_data(p, des, n_trials = 200)
  slow_guesses <- sim_dat[sim_dat$rt >= 2 & sim_dat$rt <= 4, ]
  # About 160 of the 400 trials, sd about 10.
  expect_gt(nrow(slow_guesses), 120)
  expect_true(all(slow_guesses$R == "left"))
})


test_that("guess_exclude modifies guess density and sets floor only on eligible responses in DDM", {
  des <- design(factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"),
                formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1, pGuess ~ 1),
                model = DDM, constants = c(s = log(1)),
                TC = list(guess_window = c(0, 5), guess_exclude = "right"))
  dat <- data.frame(subjects = factor(1), S = factor("left", levels = c("left", "right")),
                    R = factor(c("left", "right"), levels = c("left", "right")),
                    rt = c(0.4, 0.4), LT = 0, UT = Inf, LC = 0, UC = Inf)
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()

  # With 'right' excluded, n_resp is 1, so guess density for 'left' is 1 / (1 * 5) = 0.2
  expect_identical(attr(dadm, "guess_n_resp"), 1L)
  expect_identical(attr(dadm, "guess_response_levels"), "left")
  expect_identical(attr(dadm, "guess_response_codes"), 1L)

  base <- c(v = 1, a = log(1), t0 = log(0.2), Z = qnorm(0.5))
  pG <- 0.4
  ll_proc <- trial_lls(dadm, model, c(base, pGuess = qnorm(0)))
  ll_mix  <- trial_lls(dadm, model, c(base, pGuess = qnorm(pG)))

  # Eligible trial ('left'): mixture with g = 0.2
  expect_equal(exp(ll_mix[1]),
               (1 - pG) * exp(ll_proc[1]) + pG * 0.2,
               tolerance = 1e-6)

  # Ineligible trial ('right'): no guess component, downweighted by (1 - pG)
  expect_equal(exp(ll_mix[2]),
               (1 - pG) * exp(ll_proc[2]),
               tolerance = 1e-6)
})


test_that("guess_exclude works in race models across fast and slow paths", {
  matchfun <- function(d) d$S == d$lR
  des <- design(factors = list(subjects = 1, S = c("r1", "r2", "r3")),
                Rlevels = c("r1", "r2", "r3"), matchfun = matchfun,
                formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, pGuess ~ 1),
                model = RDM, constants = c(s = log(1)),
                TC = list(guess_window = c(0, 10), guess_exclude = "r2"))
  dat <- data.frame(
    subjects = factor(1),
    S = factor("r1", levels = c("r1", "r2", "r3")),
    R = factor(c("r1", "r2", "r3"), levels = c("r1", "r2", "r3")),
    rt = c(0.5, 0.5, 0.5), LT = 0, UT = Inf, LC = 0, UC = Inf)

  dadm_fast <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm_fast, "model")()
  expect_identical(attr(dadm_fast, "guess_n_resp"), 2L)
  expect_equal(sort(attr(dadm_fast, "guess_response_levels")), c("r1", "r3"))

  pG <- 0.3
  p_base <- c(v = log(2), B = log(1), A = log(0.3), t0 = log(0.2))
  p_0 <- c(p_base, pGuess = qnorm(0))
  p_mix <- c(p_base, pGuess = qnorm(pG))

  ll_proc <- trial_lls(dadm_fast, model, p_0)
  ll_mix  <- trial_lls(dadm_fast, model, p_mix)

  # Eligible trials ('r1', 'r3'): g = 1 / (2 * 10) = 0.05
  expect_equal(exp(ll_mix[1]), (1 - pG) * exp(ll_proc[1]) + pG * 0.05, tolerance = 1e-6)
  expect_equal(exp(ll_mix[3]), (1 - pG) * exp(ll_proc[3]) + pG * 0.05, tolerance = 1e-6)

  # Ineligible trial ('r2'): g = 0
  expect_equal(exp(ll_mix[2]), (1 - pG) * exp(ll_proc[2]), tolerance = 1e-6)

  # Also test slow path (with truncation LT > 0)
  dat$LT <- 0.1
  dadm_slow <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  ll_slow_proc <- trial_lls(dadm_slow, model, p_0)
  ll_slow_mix  <- trial_lls(dadm_slow, model, p_mix)
  expect_equal(exp(ll_slow_mix[1]), (1 - pG) * exp(ll_slow_proc[1]) + pG * 0.05, tolerance = 1e-6)
  expect_equal(exp(ll_slow_mix[2]), (1 - pG) * exp(ll_slow_proc[2]), tolerance = 1e-6)
  expect_equal(exp(ll_slow_mix[3]), (1 - pG) * exp(ll_slow_proc[3]) + pG * 0.05, tolerance = 1e-6)
})


test_that("guess_exclude preserves propriety: total mass across responses integrates to 1", {
  matchfun <- function(d) d$S == d$lR
  des <- design(factors = list(subjects = 1, S = c("a", "b")),
                Rlevels = c("a", "b"), matchfun = matchfun,
                formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, pGuess ~ 1),
                model = RDM, constants = c(s = log(1)),
                TC = list(guess_window = c(0, 10), guess_exclude = "b"))
  pG <- 0.35
  p <- c(v = log(2), B = log(1), A = log(0.3), t0 = log(0.2), pGuess = qnorm(pG))

  grid <- seq(0, 10, length.out = 201)
  dat <- do.call(rbind, lapply(c("a", "b"), function(r) data.frame(
    subjects = factor(1), S = factor("a", levels = c("a", "b")),
    R = factor(r, levels = c("a", "b")), rt = grid,
    LT = 0, UT = Inf, LC = 0, UC = Inf)))
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()

  dens <- exp(trial_lls(dadm, model, p))
  mass <- sum(vapply(split(dens, rep(1:2, each = length(grid))), function(d)
    sum(diff(grid) * (utils::head(d, -1) + utils::tail(d, -1)) / 2), numeric(1)))
  expect_equal(mass, 1, tolerance = 5e-3)
})


test_that("guess_exclude works in LogicalRules", {
  matchfun <- function(d) d$lR %in% c("A", "B", "n_A", "n_B")
  template <- data.frame(
    subjects = factor(1), S = factor("AB", levels = "AB"),
    LogicalRule = factor("OR", levels = "OR"),
    R = factor("yes", levels = c("yes", "no")),
    rt = 0.5, LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  des <- design(
    data = template,
    Rlevels = c("yes", "no"),
    fixed_accumulator_roles = factor(c("A", "B", "n_A", "n_B"),
                                     levels = c("A", "B", "n_A", "n_B")),
    matchfun = matchfun,
    model = LogicalRulesRDMSWTN,
    formula = list(v ~ 0 + GoA + GoB + NegA + NegB,
                   B ~ 0 + GoA + GoB + NegA + NegB,
                   t0 ~ 1, A ~ 1, s ~ 1, pGuess ~ 1),
    constants = c(s = 0, mG = 0, mK = 0, sv = 0),
    functions = list(
      GoA = function(d) ifelse(d$lR == "A", 1, 0),
      GoB = function(d) ifelse(d$lR == "B", 1, 0),
      NegA = function(d) ifelse(d$lR == "n_A", 1, 0),
      NegB = function(d) ifelse(d$lR == "n_B", 1, 0)
    ),
    TC = list(guess_window = c(0, 10), guess_exclude = "no")
  )
  dat <- data.frame(
    subjects = factor(c(1, 1)), S = factor(c("AB", "AB"), levels = "AB"),
    LogicalRule = factor(c("OR", "OR"), levels = "OR"),
    R = factor(c("yes", "no"), levels = c("yes", "no")),
    rt = c(0.5, 0.5), LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()
  expect_identical(attr(dadm, "guess_n_resp"), 1L)
  expect_identical(attr(dadm, "guess_response_levels"), "yes")

  p <- sampled_pars(des); p[] <- 0
  vals <- c(v_GoA = 1.6, v_GoB = 1.1, v_NegA = 0.95, v_NegB = 1.25,
            B_GoA = log(0.75), B_GoB = log(0.95), B_NegA = log(0.85), B_NegB = log(0.70),
            t0 = log(0.2), A = log(0.3))
  for (nm in names(vals)) if (nm %in% names(p)) p[[nm]] <- vals[[nm]]
  pG <- 0.25
  ll_proc <- trial_lls(dadm, model, c(p, pGuess = qnorm(0)))
  ll_mix  <- trial_lls(dadm, model, c(p, pGuess = qnorm(pG)))

  # Eligible trial ('yes'): guess density = 1 / (1 * 10) = 0.1
  expect_equal(exp(ll_mix[1]), (1 - pG) * exp(ll_proc[1]) + pG * 0.1, tolerance = 1e-6)
  # Ineligible trial ('no'): no guess density
  expect_equal(exp(ll_mix[2]), (1 - pG) * exp(ll_proc[2]), tolerance = 1e-6)
})


test_that("guess_exclude works in BAwLcorr", {
  matchfun <- function(d) d$S == d$lR
  dat <- data.frame(
    subjects = factor(rep(1, 2)),
    S = factor(c("left", "left"), levels = c("left", "right")),
    R = factor(c("left", "right"), levels = c("left", "right")),
    rt = c(0.5, 0.5), LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  m_corr <- function() BAwL(correlated = TRUE)
  des <- design(
    data = dat, Rlevels = c("left", "right"), matchfun = matchfun,
    formula = list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1, rho ~ 1, pGuess ~ 1),
    constants = c(sv = log(1), k = log(0), mG = 0, mK = 0), model = m_corr,
    TC = list(guess_window = c(0, 10), guess_exclude = "right")
  )
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()
  expect_identical(attr(dadm, "guess_n_resp"), 1L)
  expect_identical(attr(dadm, "guess_response_levels"), "left")

  p_base <- sampled_pars(des); p_base[] <- 0
  p_base["v_lRleft"] <- log(1.5); p_base["v_lRright"] <- log(1.0)
  p_base["B"] <- log(1); p_base["A"] <- log(0.2); p_base["t0"] <- log(0.15)
  p_base["rho"] <- qnorm((0.5 + 1) / 2)

  pG <- 0.3
  p_0 <- c(p_base, pGuess = qnorm(0))
  p_mix <- c(p_base, pGuess = qnorm(pG))

  ll_proc <- trial_lls(dadm, model, p_0)
  ll_mix  <- trial_lls(dadm, model, p_mix)

  # Eligible ('left'): g = 1 / (1 * 10) = 0.1
  expect_equal(exp(ll_mix[1]), (1 - pG) * exp(ll_proc[1]) + pG * 0.1, tolerance = 1e-6)
  # Ineligible ('right'): g = 0
  expect_equal(exp(ll_mix[2]), (1 - pG) * exp(ll_proc[2]), tolerance = 1e-6)
})


test_that("guess_exclude is respected in the BAwLcorr quadrature fallback", {
  matchfun <- function(d) d$S == d$lR
  response_levels <- c("a", "b", "pm")
  dat <- data.frame(
    subjects = factor(rep(1, 3)),
    S = factor(rep("a", 3), levels = response_levels),
    R = factor(response_levels, levels = response_levels),
    rt = c(0.5, 0.6, 0.7), LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  m_corr <- function() BAwL(correlated = TRUE)
  des <- design(
    data = dat, Rlevels = response_levels, matchfun = matchfun,
    formula = list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1,
                   rho ~ 1, pGuess ~ 1),
    constants = c(sv = log(1), k = 0, mG = 0, mK = 0), model = m_corr,
    TC = list(guess_window = c(0, 10), guess_exclude = "pm")
  )
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()
  p_base <- sampled_pars(des); p_base[] <- 0
  p_base[grep("^v_", names(p_base))] <- log(1.2)
  p_base["B"] <- log(1); p_base["A"] <- log(0.2); p_base["t0"] <- log(0.15)
  p_base["rho"] <- qnorm(0.9)

  # Three loaded racers use the shared-factor quadrature fallback.
  old <- Sys.getenv("EMC2_BAWLCORR_COUNTERS", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("EMC2_BAWLCORR_COUNTERS")
    else Sys.setenv(EMC2_BAWLCORR_COUNTERS = old)
  }, add = TRUE)
  Sys.setenv(EMC2_BAWLCORR_COUNTERS = "1")
  EMC2:::bawl_corr_counters_reset()
  total_ll(dadm, model, c(p_base, pGuess = qnorm(0)))
  expect_gt(EMC2:::bawl_corr_counter_values()[["gh_no_clock_trials"]], 0)

  pG <- 0.3
  ll_proc <- trial_lls(dadm, model, c(p_base, pGuess = qnorm(0)))
  ll_mix <- trial_lls(dadm, model, c(p_base, pGuess = qnorm(pG)))
  # a and b are eligible (g = 1 / (2 * 10)); pm has no guess component.
  expect_equal(exp(ll_mix[1]), (1 - pG) * exp(ll_proc[1]) + pG * 0.05,
               tolerance = 1e-6)
  expect_equal(exp(ll_mix[2]), (1 - pG) * exp(ll_proc[2]) + pG * 0.05,
               tolerance = 1e-6)
  expect_equal(exp(ll_mix[3]), (1 - pG) * exp(ll_proc[3]),
               tolerance = 1e-6)
})


test_that("guess_exclude works in RDMSWTNcorr", {
  dat <- data.frame(
    subjects = factor(rep(1, 2)),
    R = factor(c("a", "b"), levels = c("a", "b")),
    rt = c(0.7, 0.8),
    LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  des <- design(
    data = dat, Rlevels = c("a", "b"),
    formula = list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1,
                   mG ~ 1, mK ~ 1, rho ~ 1, pGuess ~ 1),
    constants = c(s = 0, mG = 0, mK = 0),
    model = RDMSWTNcorr,
    TC = list(guess_window = c(0, 10), guess_exclude = "b")
  )
  dadm <- EMC2:::design_model(dat, des, compress = FALSE, rt_resolution = NULL)
  model <- attr(dadm, "model")()
  expect_identical(attr(dadm, "guess_n_resp"), 1L)
  expect_identical(attr(dadm, "guess_response_levels"), "a")

  p <- sampled_pars(des); p[] <- 0
  p["v_lRa"] <- log(1.2); p["v_lRb"] <- log(0.8)
  p["B"] <- log(1); p["A"] <- log(0.2); p["t0"] <- log(0.1); p["sv"] <- log(0.4)
  p["rho"] <- qnorm((0.5 + 1) / 2)

  pG <- 0.25
  ll_proc <- trial_lls(dadm, model, c(p, pGuess = qnorm(0)))
  ll_mix  <- trial_lls(dadm, model, c(p, pGuess = qnorm(pG)))

  # Eligible trial ('a'): g = 1 / (1 * 10) = 0.1
  expect_equal(exp(ll_mix[1]), (1 - pG) * exp(ll_proc[1]) + pG * 0.1, tolerance = 1e-6)
  # Ineligible trial ('b'): g = 0
  expect_equal(exp(ll_mix[2]), (1 - pG) * exp(ll_proc[2]), tolerance = 1e-6)
})


test_that("guess_exclude applies 1/W guess density to missing responses (R = NA)", {
  # 1. Race model (RDM)
  des_race <- design(
    factors = list(subjects = 1, S = c("s1", "s2")),
    Rlevels = c("r1", "r2", "r3"),
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, pGuess ~ 1),
    model = RDM,
    constants = c(s = log(1)),
    TC = list(guess_window = c(0, 10), guess_exclude = "r2")
  )
  dat_race <- data.frame(
    subjects = factor(rep(1, 3)),
    S = factor(rep("s1", 3), levels = c("s1", "s2")),
    R = factor(c("r1", "r2", NA), levels = c("r1", "r2", "r3")),
    rt = c(0.5, 0.5, 0.5),
    LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  dadm_race <- EMC2:::design_model(dat_race, des_race, compress = FALSE, rt_resolution = NULL)
  model_race <- attr(dadm_race, "model")()

  pG <- 0.3
  p_base <- c(v = log(2), B = log(1), A = log(0.3), t0 = log(0.2))
  ll_proc_race <- trial_lls(dadm_race, model_race, c(p_base, pGuess = qnorm(0)))
  ll_mix_race  <- trial_lls(dadm_race, model_race, c(p_base, pGuess = qnorm(pG)))

  # r1 is eligible (g = 1 / (2 * 10) = 0.05)
  expect_equal(exp(ll_mix_race[1]), (1 - pG) * exp(ll_proc_race[1]) + pG * 0.05, tolerance = 1e-6)
  # r2 is excluded (g = 0)
  expect_equal(exp(ll_mix_race[2]), (1 - pG) * exp(ll_proc_race[2]), tolerance = 1e-6)
  # R = NA takes marginal 1/W = 1/10 = 0.10
  expect_equal(exp(ll_mix_race[3]), (1 - pG) * exp(ll_proc_race[3]) + pG * 0.10, tolerance = 1e-6)

  # 2. DDM
  des_ddm <- design(
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1, pGuess ~ 1),
    model = DDM,
    constants = c(s = log(1)),
    TC = list(guess_window = c(0, 5), guess_exclude = "right")
  )
  dat_ddm <- data.frame(
    subjects = factor(rep(1, 3)),
    S = factor(rep("left", 3), levels = c("left", "right")),
    R = factor(c("left", "right", NA), levels = c("left", "right")),
    rt = c(0.4, 0.4, 0.4),
    LT = 0, UT = Inf, LC = 0, UC = Inf
  )
  dadm_ddm <- EMC2:::design_model(dat_ddm, des_ddm, compress = FALSE, rt_resolution = NULL)
  model_ddm <- attr(dadm_ddm, "model")()

  p_ddm_base <- c(v = 1, a = log(1), t0 = log(0.2), Z = qnorm(0.5))
  pG_ddm <- 0.2
  ll_proc_ddm <- trial_lls(dadm_ddm, model_ddm, c(p_ddm_base, pGuess = qnorm(0)))
  ll_mix_ddm  <- trial_lls(dadm_ddm, model_ddm, c(p_ddm_base, pGuess = qnorm(pG_ddm)))

  # left is eligible (g = 1 / (1 * 5) = 0.2)
  expect_equal(exp(ll_mix_ddm[1]), (1 - pG_ddm) * exp(ll_proc_ddm[1]) + pG_ddm * 0.2, tolerance = 1e-6)
  # right is excluded (g = 0)
  expect_equal(exp(ll_mix_ddm[2]), (1 - pG_ddm) * exp(ll_proc_ddm[2]), tolerance = 1e-6)
  # R = NA takes marginal 1/W = 1/5 = 0.2
  expect_equal(exp(ll_mix_ddm[3]), (1 - pG_ddm) * exp(ll_proc_ddm[3]) + pG_ddm * 0.2, tolerance = 1e-6)
})
