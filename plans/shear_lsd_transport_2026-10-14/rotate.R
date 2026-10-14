suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(commandArgs(TRUE)[1])), quiet = TRUE)
for (j in 1:3) {
  ch0 <- emc[[j]]; cv <- ch0$shear$curve; k <- ch0$n_pars + cv$l
  idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch0$samples$idx), "sample")
  nodes <- which(vapply(ch0$shear$atlas$nodes, Negate(is.null), TRUE)); kept <- which(!vapply(cv$B, is.null, TRUE))
  cat(sprintf("chain %d: subjects with nodes %d, kept curves %d\n", j, length(nodes), length(kept)))
  for (d in c(1001, 1000, 980, 950, 900, 800, 700)) {
    ch <- ch0; ch$samples$idx <- d
    ctx <- ns$.emc_shear_context(ch)
    sp <- ns$.emc_shear_orbit_spread(ch, idx, k, list(function(st) ns$.emc_shear_fits_at(ctx, st$mu, st$V), function(st) ctx$packed), 1L)
    o <- order(-sp[, 2])[1:6]
    cat(sprintf("  draw %4d: total fixed %7.1f curved %7.1f | share of fixed spread in subjects with nodes %.0f%% | top: %s\n", d, sum(sp[, 2]), sum(sp[, 1]),
      100 * sum(sp[nodes, 2]) / sum(sp[, 2]), paste(sprintf("%d%s:%.0f", o, ifelse(o %in% kept, "*", ""), sp[o, 2]), collapse = " ")))
  }
}
