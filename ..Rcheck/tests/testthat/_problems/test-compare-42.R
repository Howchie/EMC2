# Extracted from test-compare.R:42

# prequel ----------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# test -------------------------------------------------------------------------
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
