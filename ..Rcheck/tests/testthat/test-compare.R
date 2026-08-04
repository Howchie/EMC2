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

test_that("savage-dickey", {
  expect_snapshot(
    round(hypothesis(samples_LNR, parameter = "m", do_plot = F, H0 = -1), 2))
  expect_snapshot(
    round(hypothesis(samples_LNR, fun = function(d) d["m"] - d["m_lMd"],
                  H0 = -0.5, do_plot = F), 2))
})
