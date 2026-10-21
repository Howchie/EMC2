local_rng_guard()  # see helper-rng.R: keep this file's RNG changes inside it
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

test_that("compare", {
  expect_snapshot(
    compare(list(samples_LNR), cores_for_props = 1)
  )
})

test_that("compare can opt into LOO without WAIC or Bayes factors", {
  out <- compare(
    list(samples_LNR),
    WAIC = FALSE,
    LOO = TRUE,
    BayesFactor = FALSE,
    print_summary = FALSE
  )

  expect_true(all(c("LOO", "wLOO", "DIC", "wDIC", "BPIC", "wBPIC") %in% names(out)))
  expect_false(any(c("WAIC", "wWAIC", "MD", "wMD") %in% names(out)))

  out_lower <- compare(
    list(samples_LNR),
    WAIC = FALSE,
    loo = TRUE,
    BayesFactor = FALSE,
    print_summary = FALSE
  )
  expect_equal(out_lower$LOO, out$LOO)
})

test_that("compare supports subject-level pointwise aggregation for LOO", {
  out_trial <- compare(
    list(samples_LNR),
    WAIC = FALSE,
    LOO = TRUE,
    pointwise = "trial",
    BayesFactor = FALSE,
    print_summary = FALSE
  )
  out_subject <- compare(
    list(samples_LNR),
    WAIC = FALSE,
    LOO = TRUE,
    pointwise = "subject",
    BayesFactor = FALSE,
    print_summary = FALSE
  )

  expect_true(all(c("LOO", "wLOO") %in% names(out_subject)))
  expect_false(isTRUE(all.equal(out_trial$LOO, out_subject$LOO)))
})

test_that("balance_by gives every cell equal weight in WAIC and LOO", {
  dat <- get_data(samples_LNR)
  ll <- EMC2:::.ll_matrix_pooled(samples_LNR)
  pw_waic <- suppressWarnings(loo::waic(ll))$pointwise[, "waic"]
  pw_loo <- suppressWarnings(loo::loo(ll))$pointwise[, "looic"]
  n_c <- table(dat$R)
  w <- nrow(dat) / (length(n_c) * as.numeric(n_c[as.character(dat$R)]))
  expect_equal(sum(w), nrow(dat))
  expect_equal(waic_pooled(samples_LNR, balance_by = "R"), sum(w * pw_waic))
  expect_equal(loo_pooled(samples_LNR, balance_by = "R"), sum(w * pw_loo))

  out <- compare(list(samples_LNR), WAIC = TRUE, LOO = TRUE, balance_by = "R",
                 BayesFactor = FALSE, print_summary = FALSE)
  expect_equal(out$WAIC, sum(w * pw_waic))
  expect_equal(out$LOO, sum(w * pw_loo))
  expect_identical(attr(out, "balance_by"), "R")

  sub <- names(samples_LNR[[1]]$data)[1]
  keep <- dat$subjects == sub
  n_s <- table(dat$R[keep])
  w_s <- sum(keep) / (length(n_s) * as.numeric(n_s[as.character(dat$R[keep])]))
  expect_equal(waic_subject(samples_LNR, subject = sub, balance_by = "R"),
               sum(w_s * pw_waic[keep]))

  expect_error(waic_pooled(samples_LNR, balance_by = "nope"), "not in the data")
  expect_error(loo_pooled(samples_LNR, balance_by = "R", pointwise = "subject"),
               "pointwise")

  plain <- compare(list(samples_LNR), WAIC = TRUE, BayesFactor = FALSE,
                   print_summary = FALSE)
  expect_error(compare(list(out, plain), print_summary = FALSE), "balance_by")
})

test_that("compare combines saved outputs and recomputes all metric weights", {
  saved <- data.frame(
    MD = c(10, 12), wMD = c(.99, .01),
    WAIC = c(9, 13), wWAIC = c(.8, .2),
    LOO = c(8, 11), wLOO = c(.01, .99),
    DIC = c(7, 9), wDIC = c(.5, .5),
    BPIC = c(6, 10), wBPIC = c(.2, .8),
    EffectiveN = c(2, 3), label = c("first", "second"),
    row.names = c("model_a", "model_b")
  )
  out <- compare(list(saved[1, , drop = FALSE], saved[2, , drop = FALSE]),
                 print_summary = FALSE)

  expect_identical(rownames(out), rownames(saved))
  expect_identical(out[, c("MD", "LOO", "DIC", "BPIC", "EffectiveN", "label")],
                   saved[, c("MD", "LOO", "DIC", "BPIC", "EffectiveN", "label")])
  for (metric in c("MD", "WAIC", "LOO", "DIC", "BPIC")) {
    expected <- exp(-(saved[[metric]] - min(saved[[metric]])) / 2)
    expected <- expected / sum(expected)
    expect_equal(out[[paste0("w", metric)]], expected)
  }
  expect_null(attr(out, "pw_ll"))
})

test_that("compare combines LOO-only saved outputs without marginal deviance", {
  a <- data.frame(LOO = 4, wLOO = 1, DIC = 5, wDIC = 1,
                  BPIC = 6, wBPIC = 1, row.names = "a")
  b <- data.frame(LOO = 6, wLOO = 1, DIC = 7, wDIC = 1,
                  BPIC = 8, wBPIC = 1, row.names = "b")
  out <- compare(list(a, b), print_summary = FALSE)
  expect_false(any(c("WAIC", "MD", "wWAIC", "wMD") %in% names(out)))
  expect_equal(sum(out$wLOO), 1)
  expect_equal(sum(out$wDIC), 1)
  expect_equal(sum(out$wBPIC), 1)
})

test_that("compare rejects malformed or mixed saved outputs", {
  good <- data.frame(DIC = 1, wDIC = 1, BPIC = 2, wBPIC = 1)
  loo <- data.frame(LOO = 1, wLOO = 1, DIC = 1, wDIC = 1, BPIC = 2, wBPIC = 1)
  bad_value <- good
  bad_value$DIC <- NA_real_

  expect_error(compare(list(good, loo), print_summary = FALSE), "same information-criterion")
  expect_error(compare(list(good, bad_value), print_summary = FALSE), "finite numeric")
  expect_error(compare(list(list(not_a_fit = TRUE), good), print_summary = FALSE),
               "Cannot mix")
})

test_that("fit-based compare does not retain pointwise likelihood matrices", {
  out <- compare(list(samples_LNR), WAIC = TRUE, LOO = FALSE,
                 BayesFactor = FALSE, print_summary = FALSE)
  expect_null(attr(out, "pw_ll"))
})


test_that("savage-dickey", {
  expect_snapshot(
    round(hypothesis(samples_LNR, parameter = "m", do_plot = F, H0 = -1), 2))
  expect_snapshot(
    round(hypothesis(samples_LNR, fun = function(d) d["m"] - d["m_lMd"],
                  H0 = -0.5, do_plot = F), 2))
})
