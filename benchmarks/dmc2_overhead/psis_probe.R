# Invoke in a fresh process per observation under the benchmark lock.
args <- commandArgs(TRUE)
if (length(args) != 1L || !args[[1]] %in% c('cold', 'warm'))
  stop('usage: psis_probe.R cold|warm')
set.seed(101)
before <- .Random.seed
weights <- seq(-2, 2, length.out = 100)
if (args[[1]] == 'warm') invisible(suppressWarnings(loo::psis(weights, r_eff = 1)))
parent_unchanged <- identical(before, .Random.seed)
timing <- system.time(result <- parallel::mclapply(1:3, function(i) {
  start <- proc.time()[['elapsed']]
  value <- suppressWarnings(loo::psis(weights, r_eff = 1))
  c(elapsed = proc.time()[['elapsed']] - start,
    k = as.numeric(value$diagnostics$pareto_k))
}, mc.cores = 3L, mc.set.seed = FALSE))
stopifnot(parent_unchanged, identical(before, .Random.seed),
          all(vapply(result, is.numeric, logical(1))))
cat(args[[1]], 'wall', unname(timing[['elapsed']]), 'children',
    paste(vapply(result, `[[`, 0, 'elapsed'), collapse = ','), 'k',
    paste(vapply(result, `[[`, 0, 'k'), collapse = ','), '\n')
