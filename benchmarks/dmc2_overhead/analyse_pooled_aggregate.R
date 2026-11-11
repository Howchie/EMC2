# Combine independent seeded fits as chains for an aggregate posterior check.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_pooled_aggregate.R results-dir')
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
seeds <- as.integer(strsplit(Sys.getenv(
  "POOLED_SEEDS", "101 202 303 404 505 606 707"), "\\s+")[[1L]])
probs <- c(.05, .5, .95)
chain_matrix <- function(ch) {
  a <- aperm(ch$samples$alpha, c(3L, 1L, 2L))
  dim(a) <- c(dim(a)[1L], dim(a)[2L] * dim(a)[3L])
  colnames(a) <- unlist(lapply(dimnames(ch$samples$alpha)[[2L]], function(s)
    paste(dimnames(ch$samples$alpha)[[1L]], s, sep = '__')), use.names = FALSE)
  pooled <- t(ch$samples$pooled)
  colnames(pooled) <- rownames(ch$samples$pooled)
  cbind(a, pooled)
}
load_arm <- function(arm) {
  out <- list(); times <- numeric()
  for (seed in seeds) {
    f <- if (seed == 101L) paste0('pooled_', arm, '.rds') else
      paste0('pooled_', seed, '_', arm, '.rds')
    x <- readRDS(file.path(root, f))
    if (is.null(out$hash)) out$hash <- x$data_md5
    stopifnot(identical(out$hash, x$data_md5))
    out$chains <- c(out$chains, lapply(x$fit, chain_matrix))
    out$times <- c(out$times, x$timing[['elapsed']])
  }
  out
}
base <- load_arm('base'); candidate <- load_arm('new')
stopifnot(identical(base$hash, candidate$hash))
make_array <- function(chains) {
  vars <- colnames(chains[[1L]])
  draws <- array(NA_real_, c(nrow(chains[[1L]]), length(chains), length(vars)),
                 dimnames = list(NULL, NULL, vars))
  for (j in seq_along(chains)) draws[, j, ] <- chains[[j]][, vars, drop = FALSE]
  as_draws_array(draws)
}
summary <- function(draws) summarise_draws(draws, 'mean', 'mcse_mean', 'ess_bulk', 'ess_tail', 'rhat',
  ~quantile2(.x, probs = probs),
  ~setNames(vapply(probs, function(p) mcse_quantile(.x, probs = p), numeric(1)),
            c('mcse_q05', 'mcse_q50', 'mcse_q95')))
da <- make_array(base$chains); db <- make_array(candidate$chains)
sa <- summary(da); sb <- summary(db)
diff <- cbind(abs(sa$q5 - sb$q5) / sqrt(sa$mcse_q05^2 + sb$mcse_q05^2),
              abs(sa$q50 - sb$q50) / sqrt(sa$mcse_q50^2 + sb$mcse_q50^2),
              abs(sa$q95 - sb$q95) / sqrt(sa$mcse_q95^2 + sb$mcse_q95^2))
colnames(diff) <- c('q05', 'q50', 'q95')
cat('common data hash', base$hash, 'seeds', paste(seeds, collapse = ','),
    'chains per arm', length(base$chains), 'draws per chain', dim(da)[1L], '\n')
cat('base/candidate total seconds', sum(base$times), sum(candidate$times),
    'paired mean seconds', mean(base$times), mean(candidate$times),
    'mean paired reduction percent', 100 * mean((base$times - candidate$times) / base$times), '\n')
for (kind in c('bulk', 'tail')) {
  col <- paste0('ess_', kind)
  cat('minimum', kind, 'ESS per draw base/candidate', min(sa[[col]]) / prod(dim(da)[1:2]),
      min(sb[[col]]) / prod(dim(db)[1:2]), 'ESS per second',
      min(sa[[col]]) / sum(base$times), min(sb[[col]]) / sum(candidate$times), '\n')
}
cat('pooled-chain quantile checks beyond four combined MCSE', sum(diff > 4), '/', length(diff),
    'largest standardized difference', max(diff), '\n')
cat('largest per-variable quantile differences:\n')
print(data.frame(variable = sa$variable, max_mcse = apply(diff, 1L, max)), row.names = FALSE)
cat('\nBaseline aggregate diagnostics:\n'); print(sa, width = Inf)
cat('\nCandidate aggregate diagnostics:\n'); print(sb, width = Inf)
