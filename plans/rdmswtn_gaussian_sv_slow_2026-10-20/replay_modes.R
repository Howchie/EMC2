suppressMessages(library(EMC2)); SP <- Sys.getenv("SP")
e <- new.env(); load("Samples/samples_control_Exp1_rdm7r_blk_gL_noA_noS_sepA_sepL_sepV_sv_t0L_mt_gaussian.RData", e)
emc <- e$emc; print(table(emc[[1]]$samples$stage))
for (j in seq_along(emc)) emc[[j]] <- EMC2:::.emc_pooled_sync(emc[[j]])
emc <- EMC2:::.emc_prepare_chain_rng(emc)
orig <- EMC2:::.emc_burn_ll_batch
invisible(orig(t(emc[[1]]$samples$alpha[, 1, 2:3]), emc[[1]]$data[[1]], emc[[1]]$model, rownames(emc[[1]]$samples$alpha)))  # warm the model before forking
jobs <- expand.grid(chain = seq_along(emc), subject = seq_len(emc[[1]]$n_subjects))
res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
  smp <- emc[[jobs$chain[k]]]; s <- jobs$subject[k]; tag <- sprintf("c%d_s%03d", jobs$chain[k], s)
  t0 <- proc.time()[[3]]
  wrap <- function(x, ...) { st <- proc.time()[[3]]; out <- orig(x, ...); el <- proc.time()[[3]] - st
    if (el > 1) saveRDS(list(x = x, el = el, ll = out), file.path(SP, "slow2", sprintf("%s_%.0f.rds", tag, st - t0)))
    if (proc.time()[[3]] - t0 > 120) stop("budget")
    out }
  assignInNamespace(".emc_burn_ll_batch", wrap, "EMC2")
  assign(".Random.seed", smp$rng$subjects[[s]], envir = globalenv())
  r <- tryCatch(EMC2:::.emc_burn_subject_mode(smp, s, proposal_prior_scale = EMC2:::.emc_laplace_prior_scale, start_covariance = TRUE), error = function(e) conditionMessage(e))
  c(k = k, t = proc.time()[[3]] - t0)
}, mc.cores = 4, mc.preschedule = FALSE)
tt <- do.call(rbind, res); o <- order(-tt[, "t"]); jb <- jobs[tt[o, "k"], ]
print(cbind(jb, t = round(tt[o, "t"], 1))[1:12, ])
