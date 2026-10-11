# Reproduce the allocation before the incremental collapsed-gain check.
# All kernel definitions, movement measurements and other allocation rules stay
# those of the loaded package; only that additional decision is removed.
previous_budget_solver <- function() {
  env <- new.env(parent = asNamespace("EMC2"))
  solve <- EMC2:::.emc_alloc_solve
  statements <- as.list(body(solve))
  gate <- which(vapply(statements, function(x)
    is.call(x) && identical(x[[1L]], as.name("if")) &&
      grepl("without_joint$objective", paste(deparse(x[[2L]]), collapse = " "), fixed = TRUE),
    logical(1)))
  stopifnot(length(gate) == 1L)
  body(solve)[[gate]] <- quote(NULL)
  environment(solve) <- env
  env$.emc_alloc_solve <- solve
  units <- EMC2:::.emc_alloc_solve_units
  statements <- as.list(body(units))
  gate <- which(vapply(statements, function(x)
    is.call(x) && identical(x[[1L]], as.name("if")) &&
      grepl("best$joint > 0", paste(deparse(x[[2L]]), collapse = " "), fixed = TRUE),
    logical(1)))
  stopifnot(length(gate) == 1L)
  body(units)[[gate]] <- quote(NULL)
  environment(units) <- env
  units
}
