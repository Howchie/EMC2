test_that("make_missing truncates before contamination", {
  set.seed(1)

  dat <- data.frame(
    subjects = factor(c("s1", "s1", "s1")),
    rt = c(0.4, 0.8, 1.2),
    R = factor(rep("left", 3), levels = "left")
  )

  out <- make_missing(dat, UT = 1, pContaminant = 1, verbose = FALSE, rt_resolution = NULL)

  expect_equal(nrow(out), 2)
  expect_true(all(is.infinite(out$rt)))
  expect_true(all(is.na(out$R)))
})

test_that("make_data filters intrinsic omissions as the design declares", {
  lbaio_design <- function(TC = NULL) design(
    factors = list(subjects = "s1", S = "x"),
    Rlevels = c("yes", "no"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    model = LBA(posdrift = FALSE),
    formula = list(v ~ 1),
    constants = c(sv = log(0.01), B = log(1), A = log(0.5), t0 = log(0.2)),
    TC = TC
  )

  # Default: the never-finish outcome is retained, as the likelihood assumes.
  set.seed(1)
  retained <- make_data(c(v = -10), lbaio_design(), n_trials = 20,
                        TC = list(UT = 4), rt_resolution = NULL)
  expect_equal(nrow(retained), 20L)
  expect_true(all(is.infinite(retained$rt)))
  expect_true(all(is.na(retained$R)))

  # design(TC = list(filter_defective = TRUE)) governs fit and simulation.
  filtering <- lbaio_design(list(UT = 4, filter_defective = TRUE))
  set.seed(1)
  expect_equal(nrow(make_data(c(v = -10), filtering, n_trials = 20,
                              rt_resolution = NULL)), 0L)
  set.seed(1)
  expect_equal(nrow(make_data(c(v = -10), lbaio_design(), n_trials = 20,
                              TC = list(UT = 4, filter_defective = TRUE),
                              rt_resolution = NULL)), 0L)
  # An explicit argument still wins.
  set.seed(1)
  expect_equal(nrow(make_data(c(v = -10), filtering, n_trials = 20,
                              rt_resolution = NULL, filter_defective = FALSE)), 20L)
})

test_that("make_missing can opt into defective omission truncation", {
  dat <- data.frame(
    subjects = factor(rep("s1", 2)),
    R = factor(c(NA, "yes"), levels = "yes"),
    rt = c(Inf, 0.5)
  )

  retained <- make_missing(dat, UT = 4, rt_resolution = NULL)
  filtered <- make_missing(dat, UT = 4, rt_resolution = NULL,
                           filter_defective = TRUE)
  expect_equal(nrow(retained), 2L)
  expect_equal(nrow(filtered), 1L)
  expect_equal(filtered$rt, 0.5)
})

test_that("censoring and truncation bounds are read from columns, not attributes", {
  dat <- data.frame(
    subjects = factor(rep("s1", 2)),
    R = factor(rep("left", 2), levels = "left"),
    lR = factor(rep("left", 2), levels = "left"),
    rt = c(0.4, 0.4)
  )

  # These are legacy representations and must not affect any downstream
  # bound resolution.  The ordinary attributes on data frames remain valid;
  # only the four censoring/truncation attributes are intentionally ignored.
  attr(dat, "LT") <- c(0.1, 0.1)
  attr(dat, "LC") <- c(0.2, 0.2)
  attr(dat, "UT") <- c(0.8, 0.8)
  attr(dat, "UC") <- c(0.7, 0.7)

  window <- EMC2:::resolve_guess_window(dat, TC = list())$window
  expect_equal(window, c(0, 5))

  dm <- matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "v"))
  attr(dm, "expand") <- c(1L, 1L)
  compressed <- EMC2:::compress_dadm(dat, list(dm), NULL, NULL)
  expect_equal(nrow(compressed), 1L)
})
