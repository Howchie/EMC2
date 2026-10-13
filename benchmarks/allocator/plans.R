# Plan summary per rebuild and chain for a calibrate.R run, with the
# coordinate that binds the minimax and how it realised.
#   EMC_LIB=<lib> Rscript benchmarks/allocator/plans.R <run_dir>
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
run_dir <- commandArgs(TRUE)[1]
tab <- readRDS(file.path(run_dir, "calib.rds"))
ks <- sort(unique(tab$k))
for (k in ks) {
  cap <- readRDS(file.path(run_dir, sprintf("alloc_%02d.rds", k)))
  for (j in seq_along(cap$records)) {
    r <- cap$records[[j]]
    t <- tab[tab$k == k & tab$chain == j, ]
    b <- which.max(t$pred)
    cat(sprintf("k%d ch%d B=%d tau %.2f tau0 %.2f regime %.1f rounds(mu) %.2f rounds(lsd) %.2f coll %.2f hmc>0 %d/%d mean r %.2f mean q %.2f | bind %s %s pred %.3f pt %.3f real %.3f | real max %.3f (%s)\n",
      k, j, r$B, r$tau, r$tau0, r$regime, sum(r$rounds[seq_len(length(r$rounds) / 2)]),
      sum(r$rounds[-seq_len(length(r$rounds) / 2)]), r$collapsed_rounds,
      sum(r$rate > 0), length(r$rate), mean(r$rate), mean(r$particle_rate),
      t$kind[b], paste(t$coord[b], t$subject[b]), t$pred[b], t$pred_pt[b], t$real[b],
      max(t$real), t$kind[which.max(t$real)]))
  }
}
