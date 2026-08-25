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
