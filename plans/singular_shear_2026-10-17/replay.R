# Read-only replay of the singular-solve investigation. Run with either the
# original library or an isolated patched installation:
# R_LIBS=/tmp/emc2-shear-check Rscript plans/singular_shear_2026-10-17/replay.R checkpoint.rds
library(EMC2)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
emc <- readRDS(args[1])
cat("Library:", find.package("EMC2"), "\n")
emc <- EMC2:::restore_custom_kernel_pointers(emc, quiet = TRUE)
ch <- emc[[1]]
p <- nrow(ch$samples$theta_mu)
i <- ch$samples$idx
mu <- ch$samples$theta_mu[, i]
V <- ch$samples$theta_var[, , i]
packed <- EMC2:::.emc_shear_pack(ch$lik_laplace, p)
frames <- function(V) {
  Vi <- chol2inv(chol(V))
  EMC2:::emc_shear_frames_cpp(packed$cube, packed$lin, Vi, drop(Vi %*% mu))
}
cat("Ordinary frames ok:", frames(V)$ok, "\n")
scale <- rep(1, p)
scale[match("sv_LoadLow", rownames(V))] <- exp(-48)
extreme <- frames(V * outer(scale, scale))
cat("Extreme frames ok:", extreme$ok, "\n")
if (isTRUE(extreme$ok)) {
  # Compare with triangular substitution, which does not drop directions.
  s <- 1L
  Vi <- chol2inv(chol(V * outer(scale, scale)))
  R <- chol(packed$cube[, , s] + Vi)
  ref <- backsolve(R, forwardsolve(t(R), packed$lin[, s] + drop(Vi %*% mu)))
  cat("Extreme centre max absolute error:", max(abs(ref - extreme$m[, s])), "\n")
}
trace(".emc_shear_orbit_terms", where = asNamespace("EMC2"), print = FALSE,
      tracer = quote(cat("  probe extent:", max(abs(grid)), "window SDs\n")))
idx <- seq_len(i)
for (name in c("v.s_low", "v.g_lure", "sv_LoadLow")) {
  j <- match(name, ch$par_names)
  unit <- sd(EMC2:::.emc_shear_window_lsd(ch, idx)[j, ])
  s0 <- exp(ch$shear$plan$log_step[p + j]) / EMC2:::.EMC_SHEAR_STEP0 / unit
  cat("Coordinate:", name, "initial width:", s0, "\n")
  timing <- system.time(w <- EMC2:::.emc_shear_orbit_width(
    ch, idx, p + j, function(st) packed, s0))
  cat("Measured width:", w, "log-SD units:", w * unit,
      "elapsed:", timing[["elapsed"]], "\n")
}
cat("Stress search: v.s_low, initial width 100 window SDs\n")
j <- match("v.s_low", ch$par_names)
stress <- tryCatch(EMC2:::.emc_shear_orbit_width(
  ch, idx, p + j, function(st) packed, 100), error = function(e) conditionMessage(e))
cat("Stress search width:", stress, "\n")
untrace(".emc_shear_orbit_terms", where = asNamespace("EMC2"))
