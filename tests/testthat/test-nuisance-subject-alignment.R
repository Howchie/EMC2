# Per-subject nuisance parameters (pContaminant, pGuess) must reach the trials
# of the subject they belong to.
#
# Regression: make_data() hands make_missing() a ONE-PER-TRIAL vector carrying
# the pars row names ("1", "2", ...).  get_missing() classified a named vector
# as a subject lookup whenever every subject level appeared among its names,
# which is automatically true for numerically labelled subjects.  Every subject
# was then given the value of the data row whose name matched its label -- all
# of them rows belonging to the first subject -- so simulated data ignored the
# fitted per-subject omission/guess rates entirely while the likelihood used
# them.  Posterior predictives silently disagreed with the fit.

test_that("get_missing keeps a per-trial vector per-trial", {
  data <- data.frame(subjects = factor(rep(1:3, each = 4)), rt = 1)
  # As built by make_data(): one value per trial, named by data row.
  supplied <- setNames(rep(c(0.1, 0.5, 0.9), each = 4), as.character(1:12))

  out <- EMC2:::get_missing(supplied, data, "pContaminant", 0, "numeric")

  expect_false(attr(out, "subjectwise"))
  expect_equal(as.numeric(out), as.numeric(supplied))
  expect_equal(as.numeric(tapply(out, data$subjects, mean)), c(0.1, 0.5, 0.9))
})

test_that("get_missing still honours a genuine subject lookup", {
  data <- data.frame(subjects = factor(rep(1:3, each = 4)), rt = 1)
  supplied <- c("1" = 0.1, "2" = 0.5, "3" = 0.9)

  out <- EMC2:::get_missing(supplied, data, "pContaminant", 0, "numeric")

  expect_true(attr(out, "subjectwise"))
  expect_equal(as.numeric(tapply(out, data$subjects, mean)), c(0.1, 0.5, 0.9))

  # Named by exactly the subject levels, so still a lookup even where the
  # per-trial and per-subject lengths coincide (one trial per subject).
  data1 <- data.frame(subjects = factor(1:3), rt = 1)
  out1 <- EMC2:::get_missing(supplied[c(3, 1, 2)], data1, "pContaminant", 0, "numeric")
  expect_true(attr(out1, "subjectwise"))
  expect_equal(as.numeric(out1), c(0.1, 0.5, 0.9))
})

test_that("get_missing recycles a scalar and reads a data column", {
  data <- data.frame(subjects = factor(rep(1:3, each = 2)), rt = 1, LT = 0.2)

  out <- EMC2:::get_missing(0.3, data, "LT", 0, "numeric")
  expect_false(attr(out, "subjectwise"))
  expect_equal(as.numeric(out), rep(0.3, 6))

  out <- EMC2:::get_missing(NULL, data, "LT", 0, "numeric")
  expect_equal(as.numeric(out), rep(0.2, 6))
})

test_that("make_missing applies pContaminant and pGuess subject by subject", {
  set.seed(11)
  n <- 4000
  data <- data.frame(
    subjects = factor(rep(1:3, each = n)),
    R = factor("Yes", levels = "Yes"),
    rt = runif(3 * n, 0.3, 0.8)
  )
  pC <- rep(c(0.05, 0.4, 0.2), each = n)
  pG <- rep(c(0.3, 0.05, 0.1), each = n)
  names(pC) <- names(pG) <- as.character(seq_len(3 * n))

  out <- EMC2:::make_missing(data, UC = 3, pContaminant = pC, pGuess = pG,
                             guess_window = c(0.1, 3))

  omit <- tapply(is.na(out$R), out$subjects, mean)
  expect_equal(as.numeric(omit), c(0.05, 0.4, 0.2), tolerance = 0.03)

  # A guess lands outside the process window [0.3, 0.8]; P(guess) = (1-pC)*pG.
  guess <- tapply(out$rt < 0.3 | (out$rt > 0.8 & is.finite(out$rt)),
                  out$subjects, function(x) mean(x, na.rm = TRUE))
  expect_true(all(guess > 0))
  expect_equal(order(as.numeric(guess)), c(2, 3, 1))
})

test_that("simulated omission rates track a per-subject pContaminant", {
  skip_on_cran()
  set.seed(12)
  subjects <- 1:3
  des <- design(
    factors = list(subjects = subjects, S = c("left", "right")),
    Rlevels = c("left", "right"),
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, pContaminant ~ 1),
    constants = c(s = log(1)),
    model = RDM
  )
  pC <- c(0.02, 0.45, 0.25)
  pars <- matrix(
    rep(c(v = log(2), B = log(1), A = log(0.5), t0 = log(0.2)), each = 3),
    nrow = 3, dimnames = list(as.character(subjects),
                              c("v", "B", "A", "t0"))
  )
  pars <- cbind(pars, pContaminant = qnorm(pC))

  dat <- make_data(pars, des, n_trials = 1500, TC = list(UC = 3))

  omit <- tapply(is.na(dat$R), dat$subjects, mean)
  expect_equal(as.numeric(omit), pC, tolerance = 0.04)
})
