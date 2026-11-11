args <- commandArgs(TRUE)
if (length(args) != 1L) stop('usage: analyse_profiles.R results-dir')
for (case in c('EXG', 'LNR')) {
  files <- file.path(args[[1]], paste0(case, '_', c('base', 'new'), '_profile.rds'))
  a <- readRDS(files[[1]]); b <- readRDS(files[[2]])
  cols <- c('stage', 'block', 'block_total', 'block_prepare', 'block_proposals',
            'block_sample', 'block_concat', 'block_check_progress')
  stopifnot(identical(a$blocks[, c('stage', 'block')], b$blocks[, c('stage', 'block')]))
  x <- a$blocks[, cols]; y <- b$blocks[, cols[-c(1,2)]]
  names(x)[-(1:2)] <- paste0(names(x)[-(1:2)], '_base')
  names(y) <- paste0(names(y), '_new')
  cat('\n', case, '\n'); print(cbind(x, y), row.names = FALSE)
  for (arm in c('base', 'new')) {
    z <- if (arm == 'base') a else b
    it <- z$iterations
    cat(arm, 'rows', nrow(it), 'wire bytes', if (all(is.na(it$wire_bytes))) NA else sum(it$wire_bytes, na.rm = TRUE),
        'peak sampled private memory', if (all(is.na(it$private_memory))) NA else max(it$private_memory, na.rm = TRUE), '\n')
  }
}
