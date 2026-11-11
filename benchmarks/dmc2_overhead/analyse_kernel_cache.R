args <- commandArgs(TRUE)
if (length(args) < 1L) stop('usage: analyse_kernel_cache.R results-dir [prefix]')
root <- args[[1L]]
prefix <- if (length(args) > 1L) args[[2L]] else 'rexg_cache'
a <- readRDS(file.path(root, paste0(prefix, '_base.rds')))
b <- readRDS(file.path(root, paste0(prefix, '_candidate.rds')))
stopifnot(identical(a$data_md5, b$data_md5), length(a$result) == length(b$result))
d <- abs(a$result - b$result)
cat('data hash equal:', identical(a$data_md5, b$data_md5), a$data_md5, '\n')
cat('likelihood vector bit-identical:', identical(a$result, b$result), '\n')
cat('maximum absolute difference:', max(d), '\n')
cat('median baseline/candidate seconds per 256-particle call:',
    median(a$timing), median(b$timing), '\n')
cat('median speedup:', median(a$timing) / median(b$timing), '\n')
