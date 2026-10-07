# Floored response times: design_model() stores a point inside each
# rt_resolution bin and the race likelihood averages the density over the bin
# (FlooredRtBins in src/likelihood_shared.h).

bin_design <- function(model, formula) {
  design(factors = list(subjects = 1, S = c("left", "right")),
         Rlevels = c("left", "right"),
         matchfun = function(d) d$S == d$lR,
         formula = formula, model = model, verbose = FALSE)
}

bin_dadm <- function(dat, des, rt_resolution, compress = FALSE) {
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   rt_resolution = rt_resolution, compress = compress))
  emc[[1]]$data[[1]]
}

bin_ll <- function(p, dadm, des) {
  as.numeric(EMC2:::calc_ll_manager(matrix(p, 1, dimnames = list(NULL, names(p))),
                                    dadm, des$model))
}

gauss_legendre <- function(n) {
  b <- seq_len(n - 1) / sqrt(4 * seq_len(n - 1)^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(1:(n - 1), 2:n)] <- b
  J[cbind(2:n, 1:(n - 1))] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

# log of the bin-averaged race density of each trial, integrated in
# z = sqrt(t - t0) over the part of the bin above t0.
bin_reference <- function(dat, des, p, h, n_nodes = 24) {
  g <- gauss_legendre(n_nodes)
  t0 <- exp(p[["t0"]])
  vapply(seq_len(nrow(dat)), function(i) {
    r <- h * floor(dat$rt[i] / h + 1e-9)
    lo <- max(r, t0); hi <- r + h
    zl <- sqrt(lo - t0); zh <- sqrt(hi - t0)
    z <- (zh - zl) / 2 * g$x + (zh + zl) / 2
    w <- (zh - zl) / 2 * g$w
    d <- bin_dadm(dat[i, ], des, NULL)
    f <- vapply(t0 + z^2, function(t) { d$rt[] <- t; exp(bin_ll(p, d, des)) }, 0)
    log(sum(w * 2 * z * f) / h)
  }, 0)
}

rdm_case <- list(model = RDM, formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
                 p = c(v = log(2), v_lMTRUE = log(1.6), B = 0, A = log(0.3),
                       t0 = log(0.25), s = 0))
lba_case <- list(model = LBA, formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                 p = c(v = 1, v_lMTRUE = 0.8, sv = 0, B = log(0.5), A = log(0.4),
                       t0 = log(0.25)))
lnr_case <- list(model = LNR, formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                 p = c(m = -0.8, m_lMTRUE = -0.4, s = log(0.4), t0 = log(0.2)))

test_that("design_model stores bin midpoints and records the bin width", {
  des <- bin_design(rdm_case$model, rdm_case$formula)
  set.seed(1)
  dat <- make_data(rdm_case$p, des, n_trials = 20, rt_resolution = NULL)
  h <- 0.05
  dadm <- bin_dadm(dat, des, h)
  expect_equal(attr(dadm, "rt_resolution"), h)
  q <- dadm$rt / h
  expect_true(all(abs(q - floor(q) - 0.5) < 1e-9))
  # Re-binning a binned rt is a no-op, so simulated data floored by make_data
  # land on the same midpoints.
  expect_equal(EMC2:::.rt_bin_midpoint(dadm$rt, h), dadm$rt)
  expect_null(attr(bin_dadm(dat, des, NULL), "rt_resolution"))
})

test_that("bin midpoints stay inside the censoring and truncation bounds", {
  x <- c(0.40, 1.47, 1.49)
  mid <- EMC2:::.rt_bin_midpoint(x, 0.05, lower = 0.42, upper = 1.48)
  expect_equal(mid, c((0.42 + 0.45) / 2, (1.45 + 1.48) / 2, (1.45 + 1.48) / 2))
  expect_true(all(mid >= 0.42 & mid < 1.48))
})

test_that("the binned race likelihood matches the bin-averaged density", {
  h <- 1 / 60
  for (case in list(rdm_case, lba_case, lnr_case)) {
    des <- bin_design(case$model, case$formula)
    set.seed(2)
    dat <- make_data(case$p, des, n_trials = 6, rt_resolution = NULL)
    ref <- bin_reference(dat, des, case$p, h)
    got <- vapply(seq_len(nrow(dat)), function(i) bin_ll(case$p, bin_dadm(dat[i, ], des, h), des), 0)
    edge <- vapply(seq_len(nrow(dat)), function(i)
      bin_ll(case$p, bin_dadm(transform(dat[i, ], rt = h * floor(rt / h + 1e-9)), des, NULL), des), 0)
    expect_lt(max(abs(got - ref)), 0.02)
    expect_lt(sum(abs(got - ref)), sum(abs(edge - ref)) / 5)
  }
})

test_that("t0 inside the earliest bin softens the likelihood instead of walling it", {
  h <- 0.05
  des <- bin_design(lba_case$model, lba_case$formula)
  set.seed(3)
  dat <- make_data(lba_case$p, des, n_trials = 10, rt_resolution = NULL)
  dadm <- bin_dadm(dat, des, h)
  lo <- h * floor(min(dat$rt) / h + 1e-9)
  ll_at <- function(t0) { p <- lba_case$p; p[["t0"]] <- log(t0); bin_ll(p, dadm, des) }
  # Continuous where t0 crosses the bin's lower edge.
  expect_equal(ll_at(lo - 1e-10), ll_at(lo + 1e-10), tolerance = 1e-6)
  # Above the stored midpoint the earliest trial still has probability, and
  # the likelihood falls as less of its bin remains feasible.
  above <- vapply(lo + c(0.6, 0.75, 0.9) * h, ll_at, 0)
  expect_true(all(is.finite(above)))
  expect_true(all(diff(above) < 0))
})

test_that("compression and the censored-data path apply the same bin rule", {
  h <- 1 / 60
  des <- bin_design(rdm_case$model, rdm_case$formula)
  set.seed(4)
  dat <- make_data(rdm_case$p, des, n_trials = 40, rt_resolution = NULL)
  p <- rdm_case$p
  p[["t0"]] <- log(h * floor(min(dat$rt) / h + 1e-9) + 0.4 * h)
  expect_equal(bin_ll(p, bin_dadm(dat, des, h, compress = TRUE), des),
               bin_ll(p, bin_dadm(dat, des, h), des), tolerance = 1e-10)
  # A censored trial sends the data down the general race path; its own
  # contribution does not depend on binning.
  cens <- dat[1, ]
  cens$rt <- Inf
  both <- rbind(dat, cens)
  both$trials <- seq_len(nrow(both))
  expect_equal(bin_ll(p, bin_dadm(both, des, h), des) - bin_ll(p, bin_dadm(dat, des, h), des),
               bin_ll(p, bin_dadm(both, des, NULL), des) - bin_ll(p, bin_dadm(dat, des, NULL), des),
               tolerance = 1e-8)
})

test_that("every subject's data keeps the bin width", {
  des <- design(factors = list(subjects = 1:2, S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
                formula = rdm_case$formula, model = RDM, verbose = FALSE)
  set.seed(5)
  dat <- make_data(rdm_case$p, des, n_trials = 10)
  emc <- suppressMessages(make_emc(dat, des, n_chains = 1, rt_resolution = 1 / 60))
  for (d in emc[[1]]$data) expect_equal(attr(d, "rt_resolution"), 1 / 60)
})
