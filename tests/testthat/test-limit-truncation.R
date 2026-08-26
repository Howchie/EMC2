# Behavioral contract tests for xlim handling in plot_density/prep_data_plot.
#
# Pins the agreed contract:
#   1. Row preservation: predictives with genuine mass beyond the observed data
#      are never subset to the axis window by prep_data_plot.
#   2. Tail mass: compute_def_dens evaluates on a shared grid spanning the
#      display window extended by 5% of its span per side, and returns
#      list(x, y) per defective level, so curves exist past the observed data.
#   3. Defective proportions: level densities are scaled by proportions
#      computed on the untruncated frame (nrow(level)/nrow(group) * p_finite).
#   4. User-supplied xlim wins over the quantile-derived default window.
#
# Note: parts of this contract are implemented in parallel in R/plot_data.R;
# until that change lands, failures of exactly those expectations are expected
# and must not be "fixed" by weakening the assertions here.

set.seed(123)

# Deterministic observed data capped below 1.0 s (columns required by
# check_data_plot: rt, subjects, plus defective_factor R and panel factor L).
make_observed <- function(n_per_level = 100L) {
  u <- (seq_len(n_per_level) - 0.5) / n_per_level
  data.frame(
    rt       = c(0.30 + 0.60 * u, 0.35 + 0.55 * u),
    subjects = "s1",
    R        = factor(rep(c("1", "2"), each = n_per_level)),
    L        = factor("A")
  )
}

# Posterior predictive: same columns plus postn. Bulk below ~0.95 s plus a
# genuine slow tail reaching ~2.87 s in both levels. The tail occupies ~3% of
# rows so that the quantile-derived default window stays below 2.5 s.
make_pred <- function(n_postn = 4L) {
  ub <- (seq_len(120L) - 0.5) / 120L
  ut <- (seq_len(4L) - 0.5) / 4L
  pool <- c(0.30 + 0.65 * ub, 1.10 + 1.768 * ut)  # bulk + tail to ~2.868 s
  do.call(rbind, lapply(c("1", "2"), function(lev) {
    do.call(rbind, lapply(seq_len(n_postn), function(p) {
      data.frame(
        rt       = pool,
        subjects = "s1",
        R        = factor(rep(lev, length(pool))),
        L        = factor("A"),
        postn    = p
      )
    }))
  }))
}

# Run expr with a png device open; device is closed even on expectation failure.
with_png_device <- function(expr) {
  png(tempfile(fileext = ".png"))
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

test_that("prep_data_plot preserves every posterior predictive row", {
  obs  <- make_observed()
  pred <- make_pred()
  expect_lt(max(obs$rt), 1.0)
  expect_gte(max(pred$rt), 2.4)

  check <- EMC2:::prep_data_plot(obs,
    post_predict = pred, prior_predict = NULL,
    to_plot = c("data", "posterior"), limits = c("data", "posterior"),
    factors = "L", defective_factor = "R", subject = NULL,
    n_cores = 1, n_post = 10, functions = NULL)

  expect_equal(nrow(check$datasets$posterior), nrow(pred))
  expect_equal(nrow(check$datasets$data), nrow(obs))
})

test_that("plot_density keeps predictive rows and evaluates density past observed range", {
  obs  <- make_observed()
  pred <- make_pred()

  res <- with_png_device(
    plot_density(obs, post_predict = pred, defective_factor = "R",
                 factors = "L", xlim = c(0, 2.5))
  )

  expect_true(is.list(res) && all(c("datasets", "sources", "xlim") %in% names(res)))
  # (a) returned invisible structure carries the untruncated posterior source
  expect_equal(nrow(res$datasets$posterior), nrow(pred))
  # user xlim forwarded through dots becomes the resolved window
  expect_equal(res$xlim, c(0, 2.5))

  # (b) density evaluation grid spans the display window extended by 5% of its
  # span per side (contract: NOT the full pooled RT range), i.e. past 2.4 s.
  win <- c(0, 2.5)
  pad <- 0.05 * diff(win)
  dd <- EMC2:::compute_def_dens(res$datasets$posterior, "R", list(),
                                from = win[1] - pad, to = win[2] + pad)
  expect_true(is.list(dd[["1"]]) && identical(names(dd[["1"]]), c("x", "y")))
  expect_length(dd[["1"]]$x, 512)
  expect_gte(max(dd[["1"]]$x), 2.4)
})

test_that("defective densities scale by untruncated-frame proportions", {
  win <- c(0, 2.5)
  pad <- 0.05 * diff(win)
  from <- win[1] - pad
  to   <- win[2] + pad

  # (i) Exact cross-level ratio: both levels share an identical multiset of
  # rts (bulk + slow tail) replicated 3 vs 2 times. With a fixed bandwidth the
  # unscaled density estimates are identical, so max(y) ratios equal the
  # untruncated count proportions exactly.
  ub <- (seq_len(120L) - 0.5) / 120L
  ut <- (seq_len(6L) - 0.5) / 6L
  pool <- c(0.30 + 0.65 * ub, 1.10 + 1.768 * ut)
  pred_ratio <- data.frame(
    rt       = c(rep(pool, 3), rep(pool, 2)),
    subjects = "s1",
    R        = factor(rep(c("1", "2"),
                          times = c(3 * length(pool), 2 * length(pool)))),
    L        = factor("A")
  )
  dd <- EMC2:::compute_def_dens(pred_ratio, "R", list(bw = 0.08),
                                from = from, to = to)
  expect_true(is.list(dd[["1"]]) && identical(names(dd[["1"]]), c("x", "y")))
  n1 <- sum(pred_ratio$R == "1")
  n2 <- sum(pred_ratio$R == "2")
  expect_equal(max(dd[["1"]]$y) / max(dd[["2"]]$y), n1 / n2,
               tolerance = 1e-6)

  # (ii) Per-level recomputation on an asymmetric frame (slow tail almost
  # exclusively in level 2): max(y) equals the level's share of the full
  # untruncated frame times the unscaled density maximum.
  pred_asym <- data.frame(
    rt       = c(pool[seq_len(120L)], pool),
    subjects = "s1",
    R        = factor(rep(c("1", "2"),
                          times = c(length(pool) - 6L, length(pool)))),
    L        = factor("A")
  )
  dda <- EMC2:::compute_def_dens(pred_asym, "R", list(bw = 0.08),
                                 from = from, to = to)
  for (lev in c("1", "2")) {
    rtv <- pred_asym$rt[pred_asym$R == lev]
    ref <- stats::density(rtv[is.finite(rtv)], bw = 0.08,
                          from = from, to = to, n = 512)
    p_exp <- mean(pred_asym$R == lev)  # share on the FULL frame; p_finite = 1
    expect_equal(max(dda[[lev]]$y), p_exp * max(ref$y), tolerance = 1e-6)
  }
})

test_that("user xlim overrides the quantile-derived window in prep_data_plot", {
  obs  <- make_observed()
  pred <- make_pred()
  args <- list(obs, post_predict = pred, prior_predict = NULL,
               to_plot = c("data", "posterior"),
               limits = c("data", "posterior"), factors = "L",
               defective_factor = "R", subject = NULL, n_cores = 1,
               n_post = 10, functions = NULL)

  default_win <- do.call(EMC2:::prep_data_plot, args)
  # precondition: the quantile-derived window is genuinely narrower than 2.5
  expect_lt(max(default_win$xlim), 2.5)

  user_win <- do.call(EMC2:::prep_data_plot, c(args, list(user_xlim = c(0, 2.5))))
  expect_equal(user_win$xlim, c(0, 2.5))

  rev_win <- do.call(EMC2:::prep_data_plot, c(args, list(user_xlim = c(2.5, 0))))
  expect_equal(rev_win$xlim, c(0, 2.5))  # supplied xlim is sorted
})
