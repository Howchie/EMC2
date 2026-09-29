suppressMessages(library(EMC2))
load("/data/work/EMC2_dev_oo/plans/convergence/samples_lba_rates_d45_snapshot_0218.RData")
## The checkpoint keeps only adapt + sample once sampling starts; use all stored
## iterations (mode-start ll is no longer available).
st <- emc[[1]]$samples$stage; burn <- seq_along(st); last <- tail(burn, 200); bi <- 1
cat("stored stages:", paste(names(table(st)), table(st)), "\n")
subs <- c("52176", "52091", "53422", "51874", "52534")
for (sj in subs) {
  cat("\n== subject", sj, "\n")
  pooled_sd <- apply(do.call(cbind, lapply(emc, function(ch) ch$samples$alpha[, sj, last])), 1, sd)
  for (k in seq_along(emc)) {
    a <- emc[[k]]$samples$alpha[, sj, ]
    ll <- emc[[k]]$samples$subj_ll[sj, ]
    d <- diff(t(a[, burn]))
    acc <- mean(rowSums(abs(d)) > 0)
    msj <- mean(rowSums(sweep(d, 2, pooled_sd, "/")^2)) / nrow(a)   # per-coordinate standardized sq jump
    worst <- names(which.max(abs(rowMeans(a[, last]) - rowMeans(do.call(cbind, lapply(emc, function(ch) ch$samples$alpha[, sj, last]))))))
    cat(sprintf("chain %d: accept %.2f | std sq jump/coord %.3f | ll first %.1f, last-500 mean %.1f (max %.1f) | %s start %.3f -> mean %.3f\n",
        k, acc, msj, ll[bi], mean(ll[last]), max(ll[burn]), worst, a[worst, bi], mean(a[worst, last])))
  }
}
