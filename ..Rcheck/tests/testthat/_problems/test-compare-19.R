# Extracted from test-compare.R:19

# prequel ----------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# test -------------------------------------------------------------------------
out <- compare(
    list(samples_LNR),
    WAIC = FALSE,
    LOO = TRUE,
    BayesFactor = FALSE,
    print_summary = FALSE
  )
expect_true(all(c("LOO", "wLOO", "DIC", "wDIC", "BPIC", "wBPIC") %in% names(out)))
