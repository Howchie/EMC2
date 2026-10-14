# Ceiling for Laplace-following frames: every subject refitted at every grid state.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(commandArgs(TRUE)[1])), quiet = TRUE)
draws <- as.integer(strsplit(commandArgs(TRUE)[2], ",")[[1]]); j <- as.integer(commandArgs(TRUE)[3])
grid <- ns$.EMC_SHEAR_PILOT_GRID
ch0 <- emc[[j]]; p <- ch0$n_pars; k <- p + ch0$shear$curve$l; n <- ch0$n_subjects
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch0$samples$idx), "sample")
for (d in draws) {
  ch <- ch0; ch$samples$idx <- d
  ctx <- ns$.emc_shear_context(ch)
  mu <- as.numeric(ch$samples$theta_mu[, d]); V <- ch$samples$theta_var[, , d]
  path <- ns$.emc_shear_coord_path(ch, idx, k, mu, V, grid)
  t0 <- proc.time()[3]
  oracle <- lapply(path, function(st) {
    fits <- parallel::mclapply(seq_len(n), function(s) ns$.emc_shear_node_fit(ch, s, st$mu, st$V, as.numeric(ch$samples$alpha[, s, d])), mc.cores = 5)
    nf <- sum(vapply(fits, is.null, TRUE))
    fits <- lapply(seq_len(n), function(s) if (is.null(fits[[s]])) ctx$fits[[s]] else fits[[s]])
    list(packed = ns$.emc_shear_pack(fits, p), fail = nf, key = c(st$mu, diag(as.matrix(st$V))))
  })
  okey <- function(st) { g <- which.min(vapply(oracle, function(o) sum(abs(o$key - c(st$mu, diag(as.matrix(st$V))))), 1)); oracle[[g]]$packed }
  sp <- ns$.emc_shear_orbit_spread(ch, idx, k, list(function(st) ctx$packed, function(st) ns$.emc_shear_fits_at(ctx, st$mu, st$V), okey), 5L)
  o <- order(-sp[, 1])[1:4]
  cat(sprintf("chain %d draw %d (%.0f s, failed fits %s): total spread fixed %.1f | current curve %.1f | oracle %.1f ; top fixed subjects %s -> oracle %s\n",
    j, d, proc.time()[3] - t0, paste(vapply(oracle, `[[`, 1, "fail"), collapse = "/"), sum(sp[, 1]), sum(sp[, 2]), sum(sp[, 3]),
    paste(sprintf("%d:%.0f", o, sp[o, 1]), collapse = " "), paste(sprintf("%.1f", sp[o, 3]), collapse = " ")))
}
