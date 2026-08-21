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

  base_p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18), t0 = log(0))
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
  grid <- seq(0, 12, length.out = 8000)
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
  p <- c(mu = log(0.35), sigma = log(0.12), tau = log(0.18), t0 = log(0),
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
  n <- 4000
  d <- data.frame(subjects = factor(1), R = factor(sample(c("a", "b"), n, TRUE)),
                  rt = runif(n, 0.3, 1.2), LT = 0, UT = Inf, LC = 0, UC = Inf)
  out <- make_missing(d, pGuess = 0.25, guess_window = c(2, 4),
                      rt_resolution = NULL)
  guessed <- out$rt > 1.5
  expect_equal(mean(guessed), 0.25, tolerance = 0.08)
  expect_true(all(out$rt[guessed] >= 2 & out$rt[guessed] <= 4))
  # Nested with the omission: P(guess) = (1 - pC) * pG, and no trial is both.
  out2 <- make_missing(d, pContaminant = 0.2, pGuess = 0.25,
                       guess_window = c(2, 4), rt_resolution = NULL)
  expect_equal(mean(is.infinite(out2$rt)), 0.2, tolerance = 0.1)
  expect_equal(mean(is.finite(out2$rt) & out2$rt > 1.5), 0.8 * 0.25,
               tolerance = 0.1)
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
