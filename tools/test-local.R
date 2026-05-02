#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
filter <- if (length(args) > 0) args[[1]] else NULL

if (is.null(filter) || !nzchar(filter)) {
  testthat::test_local(stop_on_failure = TRUE)
} else {
  testthat::test_local(filter = filter, stop_on_failure = TRUE)
}
