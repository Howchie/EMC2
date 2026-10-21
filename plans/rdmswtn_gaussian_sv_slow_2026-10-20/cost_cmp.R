suppressMessages(library(EMC2))
f <- c(gaussian = "Fits/Match_Control_Exp1_rdm7r_blk_gL_noA_noS_sepA_sepL_sepV_sv_t0L_mt_gaussian.RData",
       gamma = "Samples/samples_control_Exp1_rdm7r_blk_gL_noA_noS_sepA_sepL_sepV_sv_t0L_mt_gamma.RData")
for (nm in names(f)) { e <- new.env(); load(f[nm], e); emc <- e$emc
  cat("=====", nm, "\n"); print(table(emc[[1]]$samples$stage))
  for (c in seq_along(emc)) { x <- emc[[c]]; h <- x$subject_hmc; a <- h$alloc
    tot <- unlist(h$totals); w <- x$alloc$window
    cat(sprintf("chain %d | hmc attempts %d acc %.2f aborts %d rows %.0f cpu %.0f s | C_rest %.1f c_round %.2f hmc_fixed %.2f tau %.0f cost %.1f | window cpu_rest %.0f, sum cpu_hmc %.0f, n %s\n",
      c, tot["attempts"], tot["accepted"]/tot["attempts"], tot["aborts"], tot["rows"], tot["cpu"], a$C_rest, a$c_round, a$hmc_fixed, a$tau, a$cost,
      w$cpu_rest, sum(w$cpu_hmc), w$n))
    cp <- a$c_part; ch <- a$c_hmc
    cat(sprintf("   per-subject c_part median %.3f max %.3f (subj %s) | c_hmc median %.3f max %.3f (subj %s)\n",
      median(cp), max(cp), x$subjects[which.max(cp)], median(ch), max(ch), x$subjects[which.max(ch)]))
    if (!is.null(x$collapsed)) cat("   collapsed:", paste(names(unlist(x$collapsed))[1:min(6,length(unlist(x$collapsed)))], collapse=" "), "\n")
  }
  s <- emc[[1]]$samples; i <- s$stage == "sample"
  if (any(i)) for (p in c("s","sv_lMFALSE","sv_lMTRUE")) cat(sprintf("  %s exp() quantiles sample: %s\n", p, paste(signif(exp(quantile(s$alpha[p,,i], c(0,.01,.5,.99,1))),2), collapse=" ")))
}
