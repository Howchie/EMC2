suppressMessages(library(EMC2)); SP <- Sys.getenv("SP"); s <- as.integer(Sys.getenv("S", "22"))
e <- new.env(); load("Samples/samples_control_Exp1_rdm7r_blk_gL_noA_noS_sepA_sepL_sepV_sv_t0L_mt_gaussian.RData", e)
emc <- e$emc; for (j in seq_along(emc)) emc[[j]] <- EMC2:::.emc_pooled_sync(emc[[j]])
emc <- EMC2:::.emc_prepare_chain_rng(emc); smp <- emc[[1]]
orig <- EMC2:::.emc_burn_ll_batch; log <- new.env(); log$n <- 0; log$slow <- list(); t0 <- proc.time()[[3]]
wrap <- function(x, ...) { st <- proc.time()[[3]]; out <- orig(x, ...); el <- proc.time()[[3]] - st
  log$n <- log$n + 1; X <- if (is.matrix(x)) x else matrix(x, 1)
  if (el > 0.5) { cat(sprintf("call %d: %d rows, %.1f s (t=%.0f)\n", log$n, nrow(X), el, proc.time()[[3]] - t0))
    log$slow[[length(log$slow) + 1]] <- list(x = X, el = el, ll = out); saveRDS(log$slow, file.path(SP, sprintf("slow_s%d.rds", s))) }
  if (proc.time()[[3]] - t0 > 150) stop("trace budget")
  out }
assignInNamespace(".emc_burn_ll_batch", wrap, "EMC2")
assign(".Random.seed", smp$rng$subjects[[s]], envir = globalenv())
r <- tryCatch(EMC2:::.emc_burn_subject_mode(smp, s, proposal_prior_scale = EMC2:::.emc_laplace_prior_scale, start_covariance = TRUE), error = function(e) conditionMessage(e))
cat("calls", log$n, "\n")
