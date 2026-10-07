# Shear move vs ensemble update, from results_shear/<case or ''><move>/.
#   Rscript benchmarks/production-adaptation/compare_shear.R
root <- "benchmarks/production-adaptation/results_shear"
dirs <- list.dirs(root, recursive = FALSE)
dirs <- dirs[grepl("(shear|ensemble)$", dirs)]
rows <- list()
for (d in dirs) for (f in list.files(d, "^prod_.*\\.rds$", full.names = TRUE)) {
  r <- readRDS(f); s <- r$summary; pv <- r$per_var
  w <- r$warmup
  emc_f <- sub("/prod_", "/emc_", f)
  acc <- NA_real_
  if (file.exists(emc_f)) {
    e <- readRDS(emc_f)
    sh <- lapply(e, `[[`, "shear")
    if (!is.null(sh[[1]])) acc <- mean(vapply(sh, function(x) x$accepted / max(1, x$proposed), 0))
  }
  g <- function(what, col) s[[col]][s$what == what]
  rows[[length(rows) + 1L]] <- data.frame(
    case = r$case, seed = r$seed, move = sub(".*_", "", r$tag),
    burn = w$burn_iter, burn_ok = w$burn_converged, sample_conv = r$stop$converged,
    transitions = r$stop$transitions, kept = r$stop$retained,
    rhat_a = r$stop$rhat_alpha_max, rhat_mu = r$stop$rhat_mu_max,
    a_min = g("alpha", "ess_bulk_min"), a_med = g("alpha", "ess_bulk_med"),
    mu_min = g("mu", "ess_bulk_min"), mu_med = g("mu", "ess_bulk_med"),
    lsd_min = g("logsd", "ess_bulk_min"),
    warm_cpu = w$burn_cpu + w$adapt_cpu, prod_cpu = r$prod_cpu,
    mu_min_per_kcpu = 1000 * g("mu", "ess_bulk_min") / r$prod_cpu,
    shear_acc = acc)
}
tab <- do.call(rbind, rows)
print(tab[order(tab$case, tab$seed, tab$move), ], digits = 3, row.names = FALSE)
