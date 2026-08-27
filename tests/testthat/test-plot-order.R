test_that("data factor levels determine plot panel order", {
  observed <- data.frame(
    rt = rep(c(0.2, 0.3, 0.4, 0.5), 3),
    subjects = factor(rep("s1", 12)),
    R = factor(rep(c("error", "correct"), 6),
               levels = c("error", "correct")),
    L = factor(rep(c("third", "first", "second"), each = 4),
               levels = c("second", "third", "first"))
  )
  predictive <- observed
  predictive$L <- factor(as.character(predictive$L))
  predictive$rt <- predictive$rt + 0.01
  expected <- c("L=second", "L=third", "L=first")

  with_png <- function(expr) {
    png(tempfile(fileext = ".png"))
    on.exit(dev.off(), add = TRUE)
    force(expr)
  }

  density_result <- with_png(plot_density(
    observed, post_predict = predictive, factors = "L", defective_factor = "R",
    layout = c(1, 3), legendpos = c(NA, NA)
  ))
  expect_equal(levels(density_result$datasets$data$group_key), expected)
  expect_equal(levels(density_result$datasets$posterior$group_key),
               c("L=first", "L=second", "L=third"))

  stat_result <- with_png(plot_stat(
    observed, post_predict = predictive, factors = "L",
    stat_fun = function(data) mean(data$rt), layout = c(1, 3),
    legendpos = c(NA, NA)
  ))
  expect_equal(as.character(stat_result$data$L),
               c("second", "third", "first"))

  cdf_result <- plot_cdf(
    observed, post_predict = predictive, factors = "L", defective_factor = "R",
    layout = c(1, 3), legendpos = c(NA, NA), return_cdf = TRUE
  )
  expect_equal(attr(cdf_result, "unique_group_keys"), expected)
})
