test_that("plot_cdf smooth is opt-in and uses the requested RT grid", {
  observed <- data.frame(
    rt = c(0.2, 0.3, 0.4, 0.5),
    subjects = factor("s1"),
    R = factor(c("a", "a", "b", "b")),
    UC = Inf
  )

  raw <- plot_cdf(observed, to_plot = "data", return_cdf = TRUE)
  smooth <- plot_cdf(observed, to_plot = "data", smooth = TRUE,
                     rt_resolution = 0.1, layout = c(1, 1),
                     legendpos = c(NA, NA), return_cdf = TRUE)

  expect_length(raw$data[["All Data"]][["a"]][, "x"], 100L)
  expect_equal(smooth$data[["All Data"]][["a"]][, "x"],
               seq(0, 0.5, by = 0.1))
  expect_true(all(vapply(smooth$data[["All Data"]], function(cdf) {
    all(diff(cdf[, "y"]) >= -sqrt(.Machine$double.eps)) &&
      all(cdf[, "y"] >= 0)
  }, logical(1))))
  expect_false(identical(raw$data[["All Data"]][["a"]][, "x"],
                         smooth$data[["All Data"]][["a"]][, "x"]))
})


test_that("smoothed CDF grid includes finite upper censoring bounds", {
  observed <- data.frame(
    rt = c(0.2, Inf, Inf),
    subjects = factor("s1"),
    R = factor(c("a", "b", "c")),
    UC = c(Inf, 1.0, 2.0),
    UT = c(Inf, 1.5, Inf)
  )
  smooth <- plot_cdf(observed, to_plot = "data", smooth = TRUE,
                     rt_resolution = 0.25, return_cdf = TRUE)

  expect_equal(max(smooth$data[["All Data"]][["a"]][, "x"]), 2)
})


test_that("smoothed CDF plotting handles posterior draws and percentile markers", {
  observed <- data.frame(
    rt = c(0.2, 0.3, 0.4, 0.5),
    subjects = factor("s1"),
    R = factor(c("a", "a", "b", "b")),
    UC = Inf
  )
  predictive <- do.call(rbind, lapply(seq_len(2), function(postn) {
    out <- observed
    out$rt <- out$rt + postn / 100
    out$postn <- postn
    out
  }))

  png(tempfile(fileext = ".png"))
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot_cdf(
    observed, post_predict = predictive, smooth = TRUE,
    rt_resolution = 0.1, layout = c(1, 1), legendpos = c(NA, NA)
  ))
})
