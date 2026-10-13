# Group movement per ensemble sweep as a function of the pool size M and the
# rounds R, on saved sample-stage states, without fits. For each state, M and
# replicate: fresh pools from the chain's frozen proposals, their likelihoods,
# then single rounds chained up to RMAX (identical to one call with R rounds),
# recording after each round the standardised squared jump of the group means
# and log sds, the share of subjects whose selection changed, the pools'
# weight ESS, and the CPU of the likelihood rows and of the rounds.
#   EMC_LIB=<lib> EMC_IN=<fit.rds> STATES=8 REPS=4 MS=8,16,32,64,128 RMAX=20 CORES=4 Rscript sweep_mr.R out.rds
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
cat("EMC2 from", find.package("EMC2"), "\n")
ns <- asNamespace("EMC2")
emc <- readRDS(Sys.getenv("EMC_IN"))
emc <- ns$restore_custom_kernel_pointers(ns$restore_duplicates(emc), quiet = TRUE)
s <- emc[[1]]
type <- s$type; p <- s$n_pars; n <- s$n_subjects
n_states <- as.integer(Sys.getenv("STATES", "8")); reps <- as.integer(Sys.getenv("REPS", "4"))
Ms <- as.integer(strsplit(Sys.getenv("MS", "8,16,32,64,128"), ",")[[1]]); Rmax <- as.integer(Sys.getenv("RMAX", "20"))
cores <- as.integer(Sys.getenv("CORES", "4"))
cat(sprintf("type %s p %d subjects %d; states %d reps %d M %s rounds up to %d\n", type, p, n, n_states, reps, paste(Ms, collapse = ","), Rmax))
# Posterior scale of the group means and log sds over every chain's sample draws.
draws <- function(f) do.call(cbind, lapply(emc, function(ch) {
  w <- which(ch$samples$stage == "sample" & seq_along(ch$samples$stage) <= ch$samples$idx)
  vapply(w, function(t) f(ch$samples, t), numeric(p))
}))
sd_mu <- apply(draws(function(sm, t) sm$theta_mu[, t]), 1, sd)
sd_lsd <- apply(draws(function(sm, t) 0.5 * log(diag(as.matrix(sm$theta_var[, , t])))), 1, sd)
sm <- s$samples
w_sample <- which(sm$stage == "sample" & seq_along(sm$stage) <= sm$idx)
states <- w_sample[round(seq(1, length(w_sample), length.out = n_states + 2))[-c(1, n_states + 2)]]
pars_at <- function(t) {
  tmu <- sm$theta_mu[, t]; tvar <- as.matrix(sm$theta_var[, , t])
  if (type == "standard") list(tmu = tmu, tvar = tvar, tvinv = solve(tvar), a_half = sm$a_half[, t],
                               subj_mu = matrix(tmu, p, n), alpha = sm$alpha[, , t]) else
    list(tmu = tmu, tvar = tvar, tvinv = diag(1 / diag(tvar), p), alpha = sm$alpha[, , t])
}
ens_q <- s$ensemble_q
stopifnot(!is.null(ens_q))
cfg <- ns$.emc_ensemble_options()
cpu_now <- function() { x <- proc.time(); sum(x[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE) }
ll_rows <- function(rows) {
  out <- parallel::mclapply(seq_len(n), function(k)
    ns$.emc_ensemble_ll_candidate(rows[[k]], s$data[[k]], s$model), mc.cores = cores)
  matrix(unlist(out), nrow = nrow(rows[[1]]))
}
rows_out <- list()
for (j in seq_along(states)) {
  t <- states[j]
  pars0 <- pars_at(t)
  alpha0 <- matrix(sm$alpha[, , t], p, n, dimnames = list(s$par_names, NULL)); ll0 <- sm$subj_ll[, t]
  lsd0 <- 0.5 * log(diag(pars0$tvar))
  for (M in Ms) for (rep in seq_len(reps)) {
    set.seed(1000 * j + 10 * M + rep)
    pools <- withr::with_options(list(emc2.ensemble_pool = M),
                                 ns$.emc_ensemble_draw_pools(s, "sample", ens_q$mu, ens_q$var))
    c0 <- cpu_now(); e0 <- proc.time()[["elapsed"]]
    ll_new <- ll_rows(pools$rows)
    cpu_ll <- cpu_now() - c0; el_ll <- proc.time()[["elapsed"]] - e0
    k0 <- vapply(seq_len(n), function(k) sample.int(M, 1L), 1L)
    st <- ns$emc_ensemble_stack_cpp(pools$rows, pools$lq, ll_new, rbind(alpha0, ll0), k0,
                                    pools$L, pools$mu, pools$logdet, pools$c, cfg$tail_weight, cfg$tail_df)
    stack <- st$pool_stack; colnames(stack) <- s$par_names
    cfg1 <- cfg; cfg1$pool <- M; cfg1$rounds <- 1L
    pars <- pars0; sel <- k0; cpu_rounds <- 0
    for (r in seq_len(Rmax)) {
      c0 <- Sys.time()
      res <- ns$.emc_ensemble_rounds_cpp(s, pars, stack, st$base_w, sel, M, cfg1)
      cpu_rounds <- cpu_rounds + as.numeric(Sys.time() - c0, units = "secs")
      prev <- sel; sel <- res$selected; pars <- res$pars
      rows_out[[length(rows_out) + 1L]] <- data.frame(
        state = j, M = M, rep = rep, round = r,
        jump_mu = mean(((pars$tmu - pars0$tmu) / sd_mu)^2),
        jump_lsd = mean(((0.5 * log(diag(pars$tvar)) - lsd0) / sd_lsd)^2),
        jump_mu_max = max(((pars$tmu - pars0$tmu) / sd_mu)^2),
        moved = mean(sel != k0), changed = mean(sel != prev),
        wess = mean(res$wess) / M, wess_min = min(res$wess) / M,
        eff = mean(res$eff), eff_med = median(res$eff), eff_q10 = unname(quantile(res$eff, 0.1)),
        cpu_ll = cpu_ll, el_ll = el_ll, cpu_rounds = cpu_rounds)
    }
  }
  cat(sprintf("state %d/%d done\n", j, length(states)))
}
out <- do.call(rbind, rows_out)
res <- list(table = out, type = type, p = p, n = n, sd_mu = sd_mu, sd_lsd = sd_lsd,
            states = states, lib = find.package("EMC2"), fixture = Sys.getenv("EMC_IN"))
saveRDS(res, commandArgs(TRUE)[1])
agg <- aggregate(cbind(jump_mu, jump_lsd, moved, changed, wess, eff, eff_med, cpu_ll, cpu_rounds) ~ M + round, out, mean)
print(agg[agg$round %in% c(1, 2, 3, 5, 10, 20), ], digits = 3, row.names = FALSE)
