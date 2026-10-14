suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(commandArgs(TRUE)[1])), quiet = TRUE)
grid <- c(-2, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2)
who <- list(c(39, 21, 75, 58, 2), c(75, 60, 58, 2), c(53, 35, 75, 2))
for (j in seq_along(emc)) {
  ch <- emc[[j]]; cv <- ch$shear$curve; k <- ch$n_pars + cv$l
  idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch$samples$idx), "sample")
  ctx <- ns$.emc_shear_context(ch)
  sets <- list(fixed = function(st) ctx$packed, curved = function(st) ns$.emc_shear_fits_at(ctx, st$mu, st$V))
  for (nm in names(sets)) {
    tm <- ns$.emc_shear_orbit_terms(ch, idx, k, grid, sets[[nm]], 1L)
    cat(sprintf("chain %d %s frames, grid %s\n", j, nm, paste(grid, collapse = " ")))
    for (s in who[[j]]) cat(sprintf("   subj %2d: %s\n", s, paste(sprintf("%7.1f", tm$per[, s] - tm$per[which(grid == 0), s]), collapse = " ")))
    if (nm == "fixed") cat(sprintf("   lp (total): %s\n", paste(sprintf("%7.1f", tm$lp - tm$lp[which(grid == 0)]), collapse = " ")))
  }
  w <- ch$shear$walls; if (!is.null(w)) cat("   walls on subjects:", paste(which(w$has)[which(w$has) %in% who[[j]]], collapse = ","), "\n")
}
