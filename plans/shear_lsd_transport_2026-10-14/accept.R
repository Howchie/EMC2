# Acceptance of a jump along the lsd coordinate to each grid offset: exp of the orbit log-density change.
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
ns <- asNamespace("EMC2")
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(commandArgs(TRUE)[1])), quiet = TRUE)
j <- as.integer(commandArgs(TRUE)[3]); draws <- as.integer(strsplit(commandArgs(TRUE)[2], ",")[[1]])
grid <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2); NC <- 5L
ch0 <- emc[[j]]; p <- ch0$n_pars; l <- ch0$shear$curve$l; k <- p + l; n <- ch0$n_subjects
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch0$samples$idx), "sample")
ctx0 <- ns$.emc_shear_context(ch0)
cmode <- function(lin, P, Vi, mu) solve(P + Vi, lin + Vi %*% mu)
refit_all <- function(ch, mu, V) { Vi <- solve(V)
  fits <- parallel::mclapply(seq_len(n), function(s) ns$.emc_shear_node_fit(ch, s, mu, V, drop(cmode(ctx0$packed$lin[, s], ctx0$packed$cube[, , s], Vi, mu))), mc.cores = NC)
  ns$.emc_shear_pack(lapply(seq_len(n), function(s) if (is.null(fits[[s]])) ctx0$fits[[s]] else fits[[s]]), p) }
cat(sprintf("chain %d, coordinate log sd(%s); step in window sds; window sd of that log sd = %.3f\n", j, rownames(ch0$samples$theta_mu)[l], ch0$shear$curve$lsd_scale))
acc <- list(fixed = NULL, refit = NULL)
for (d in draws) {
  ch <- ch0; ch$samples$idx <- d; ctx <- ns$.emc_shear_context(ch)
  memo <- new.env(); rf <- function(st) { key <- paste(signif(diag(as.matrix(st$V)), 10), collapse = ","); if (is.null(memo[[key]])) memo[[key]] <- refit_all(ch, st$mu, st$V); memo[[key]] }
  for (nm in names(acc)) {
    tm <- ns$.emc_shear_orbit_terms(ch, idx, k, grid, if (nm == "fixed") function(st) ctx$packed else rf, NC)
    lpo <- rowSums(tm$per) + tm$lp; lpo <- lpo - lpo[grid == 0]
    acc[[nm]] <- rbind(acc[[nm]], pmin(1, exp(lpo)))
  }
}
for (nm in names(acc)) { a <- acc[[nm]]; colnames(a) <- grid; rownames(a) <- draws; cat("\n", nm, "frames: acceptance of a jump to each offset (rows = draws)\n"); print(round(a, 2)); cat(" mean:", round(colMeans(a), 2), "\n") }
