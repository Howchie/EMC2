# Summary of continue_arm.R results: one row per run, then per label the
# median and range over seeds of ESS per kCPU-second on the slowest fixed
# direction, and each label's ratio to a reference label.
#   Rscript benchmarks/subject-hmc/summarise_cont.R <dir> [reference_label]
args <- commandArgs(TRUE)
dir <- args[1]; ref <- if (length(args) > 1) args[2] else NULL
files <- list.files(dir, "\\.rds$", full.names = TRUE)
rows <- lapply(files, function(f) {
  r <- readRDS(f)
  k <- which.min(r$dirs$ess_bulk)
  data.frame(label = r$label, seed = r$seed, iters = r$iters, scored = r$scored,
             cpu = round(r$cpu), s_iter_chain = round(r$cpu_per_iter_chain, 2),
             particles = paste(unique(r$n_particles), collapse = "-"),
             slow_dir = k, ess_slow = round(r$dirs$ess_bulk[k]),
             ess_slow_h1 = round(r$half_dirs[[1]][k]), ess_slow_h2 = round(r$half_dirs[[2]][k]),
             slow_per_kcpu = round(1000 * r$dirs$ess_bulk[k] / r$cpu, 2),
             dir_med_per_kcpu = round(1000 * median(r$dirs$ess_bulk) / r$cpu, 2),
             mu_min_per_kcpu = round(1000 * min(r$mu$ess_bulk) / r$cpu, 2),
             lsd_min_per_kcpu = round(1000 * min(r$lsd$ess_bulk) / r$cpu, 2),
             alpha_min_per_kcpu = round(1000 * min(r$alpha$ess_bulk) / r$cpu, 2),
             rhat_max = round(max(r$mu$rhat, r$lsd$rhat, r$alpha$rhat), 3),
             hmc = round(mean(r$moves["hmc", ])), shear = round(mean(r$moves["shear", ])),
             collapsed = round(mean(r$moves["collapsed", ])))
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$label, tab$seed), ]
print(tab, row.names = FALSE)
metrics <- c("slow_per_kcpu", "dir_med_per_kcpu", "mu_min_per_kcpu", "lsd_min_per_kcpu", "alpha_min_per_kcpu")
agg <- do.call(rbind, lapply(split(tab, tab$label), function(t) {
  data.frame(label = t$label[1], seeds = nrow(t),
             do.call(cbind, lapply(metrics, function(m) {
               d <- data.frame(median(t[[m]]), min(t[[m]]), max(t[[m]]))
               names(d) <- paste0(m, c("_med", "_min", "_max")); d })))
}))
cat("\nPer label (median, min, max over seeds):\n")
print(agg[, c("label", "seeds", grep("slow_per_kcpu|mu_min", names(agg), value = TRUE))], row.names = FALSE)
if (!is.null(ref) && ref %in% agg$label) {
  r0 <- agg[agg$label == ref, ]
  cat(sprintf("\nRatio to %s (median over seeds; 'separated' = the seed ranges do not overlap):\n", ref))
  for (l in setdiff(agg$label, ref)) {
    a <- agg[agg$label == l, ]
    for (m in metrics) {
      sep <- a[[paste0(m, "_min")]] > r0[[paste0(m, "_max")]] || a[[paste0(m, "_max")]] < r0[[paste0(m, "_min")]]
      cat(sprintf("  %-20s %-20s %5.2fx %s\n", l, m, a[[paste0(m, "_med")]] / r0[[paste0(m, "_med")]],
                  if (sep) "separated" else "overlapping"))
    }
  }
}
