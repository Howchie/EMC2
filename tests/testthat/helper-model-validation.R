# Development-only model validation is run explicitly with
# `EMC2_TEST_LEVEL=full`.  The normal package pass keeps the small contract
# suite in test-model-family-contracts.R.
skip_model_validation <- function() {
  if (!identical(Sys.getenv("EMC2_TEST_LEVEL"), "full")) {
    testthat::skip("development model validation (set EMC2_TEST_LEVEL=full)")
  }
}
