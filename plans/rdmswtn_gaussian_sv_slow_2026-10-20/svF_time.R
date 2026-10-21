suppressMessages(library(EMC2))
e <- new.env(); load("Fits/Match_Control_Exp1_rdm7r_blk_gL_noA_noS_sepA_sepL_sepV_sv_t0L_mt_gaussian.RData", e)
ch <- e$emc[[1]]; a <- ch$samples$alpha; ni <- dim(a)[3]; model <- ch$model
invisible(EMC2:::calc_ll_manager(t(a[, 1, ni - 0:1]), ch$data[[1]], model))
tm <- function(j, svf, s = NULL) { P <- t(a[, j, ni - 0:9]); P[, "sv_lMFALSE"] <- log(svf); if (!is.null(s)) P[, "s"] <- log(s)
  system.time(EMC2:::calc_ll_manager(P, ch$data[[j]], model))[["elapsed"]] / 10 * 1000 }
cat("ms per subject-likelihood call (10 draws averaged), by sv_lMFALSE\n")
for (j in c(1, 30, 60, 90)) cat(sprintf("subj %-6s at draws: %5.1f | svF=1e-2: %5.1f  1e-3: %5.1f  1e-4: %5.1f  1e-5: %5.1f  1e-6: %5.1f | s=.1,svF=1e-4: %5.1f\n",
  ch$subjects[j], { P <- t(a[, j, ni - 0:9]); system.time(EMC2:::calc_ll_manager(P, ch$data[[j]], model))[["elapsed"]] / 10 * 1000 },
  tm(j, 1e-2), tm(j, 1e-3), tm(j, 1e-4), tm(j, 1e-5), tm(j, 1e-6), tm(j, 1e-4, 0.1)))
