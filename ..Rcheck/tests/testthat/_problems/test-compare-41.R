# Extracted from test-compare.R:41

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
