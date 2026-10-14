# Compare bridge displacement with the existing eight-candidate round.
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2"); args <- commandArgs(TRUE)
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(readRDS(args[1])), quiet = TRUE)
ch <- emc[[1]]; p <- ch$n_pars; l <- match("B_Espeed", ch$par_names)
idx <- ns$.emc_chain_proposal_window(emc, seq_len(ch$samples$idx), "sample")
wsd <- sd(0.5 * log(ch$samples$theta_var[l, l, idx]))
draws <- c(1001, 1000, 980, 950, 900, 800, 700, 600)
RNGkind("L'Ecuyer-CMRG"); set.seed(sum(utf8ToInt("eight-candidate-displacement")))
results <- parallel::mclapply(draws, function(d) {
  c <- ch; c$samples$idx <- d; ctx <- ns$.emc_shear_context(c)
  aa <- matrix(c$samples$alpha[, , d], p)
  cur <- ns$.emc_shear_state(c$samples$theta_mu[, d], c$samples$theta_var[, , d], ctx$packed, c$prior, ctx$spec)
  x <- list(cur = cur, cur_hier = ns$.emc_shear_log_hier(aa, cur), alpha = aa,
            ll = c$samples$subj_ll[, d], rounds = 0L, build_el = 0, ll_el = 0,
            wpool = NULL, wpool_part = NULL)
  direction <- replace(numeric(2 * p), p + l, 1)
  out <- NULL
  for (scale in c(0.05, 0.1, 0.2, 0.4, 0.8, 1.6)) for (rep in 1:12) {
    tt <- proc.time(); o <- ns$.emc_shear_round(x, direction, scale * wsd, c, ctx)
    out <- rbind(out, data.frame(draw = d, scale = scale, E = o$E / wsd^2,
                 cpu = sum((proc.time() - tt)[1:2]), build_el = o$build_el, ll_el = o$ll_el))
  }
  out
}, mc.cores = 3L, mc.set.seed = TRUE)
results <- do.call(rbind, results)
saveRDS(list(results = results, wsd = wsd, idx = idx), args[2])
print(aggregate(cbind(E, cpu, build_el, ll_el) ~ scale, results, mean))
