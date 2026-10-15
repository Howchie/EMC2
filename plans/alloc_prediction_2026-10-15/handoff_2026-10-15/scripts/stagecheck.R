suppressPackageStartupMessages(library(posterior))
tau_of <- function(G, idx) { a <- simplify2array(lapply(G, function(g) g[idx, , drop = FALSE]))
  apply(a, 2, function(x) nrow(x) * ncol(x) / ess_basic(x)) }
for (f in commandArgs(TRUE)) {
  G <- readRDS(f)$draws
  run <- tau_of(G, 101:400); slow <- order(-run)[1:8]
  cat(basename(f), " realized slowest 8:", round(run[slow]), "\n")
  for (b in 1:3) { st <- tau_of(G, 1:(100 * b))
    cat(sprintf("  stage to %d: slow8 pred %s | median ratio slow8 %.2f all %.2f spearman %.2f\n", 100 * b,
      paste(round(st[slow]), collapse = " "), median(st[slow] / run[slow]), median(st / run), cor(st, run, method = "spearman"))) }
}
