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
