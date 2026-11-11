# Compare state fields for fixture pairs; timing evidence is reported separately.
args <- commandArgs(TRUE)
if (length(args) < 2) stop('usage: analyse_pairs.R results-dir pair-prefix...')
for (prefix in args[-1]) {
  path <- file.path(args[[1]], prefix)
  a <- readRDS(paste0(path, '_base.rds')); b <- readRDS(paste0(path, '_new.rds'))
  cat('\n', prefix, 'elapsed baseline', a$timing[['elapsed']], 'candidate', b$timing[['elapsed']], '\n')
  if (!is.null(a$data_md5) || !is.null(b$data_md5)) cat('data checksum equal', identical(a$data_md5, b$data_md5), '\n')
  if (!is.null(a$fixture_md5) || !is.null(b$fixture_md5)) cat('fixture checksum equal', identical(a$fixture_md5, b$fixture_md5), '\n')
  cat('caller RNG equal', identical(a$seed, b$seed), '\n')
  for (j in seq_along(a$fit)) {
    x <- a$fit[[j]]; y <- b$fit[[j]]
    fields <- setdiff(union(names(x), names(y)), c('model', 'prior', 'data'))
    different <- fields[!vapply(fields, function(f) identical(x[[f]], y[[f]]), logical(1))]
    cat('chain', j, 'samples equal', identical(x$samples, y$samples),
        'RNG equal', identical(x$rng, y$rng), 'differing state fields', paste(different, collapse = ','), '\n')
    for (f in setdiff(different, c('samples', 'rng'))) {
      delta <- all.equal(x[[f]], y[[f]])
      if (!isTRUE(delta)) cat(' ', f, ':', paste(delta, collapse = '; '), '\n')
    }
  }
}
