# Isolate the exact subject-likelihood cost of one rejected collapsed round.
# This deliberately uses the serial fallback (no worker pool), repeats the
# same proposal from the same checkpoint, and alternates scales to reduce
# order/cache bias. It measures round cost, not whole-sampler ESS/sec.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
options(emc2.failure_policy = "strict")

prepared <- readRDS("/tmp/collapsed_transport_followup/safe_fit_ess_start.rds")
chain <- prepared$start[[2L]]
ctx <- EMC2:::.emc_collapsed_context(chain)
idx <- chain$samples$idx
p <- chain$n_pars
cur <- EMC2:::.emc_shear_state(chain$samples$theta_mu[, idx],
  as.matrix(chain$samples$theta_var[, , idx]), ctx$fits, chain$prior, ctx$spec)
phi0 <- EMC2:::.emc_collapsed_to_phi(cur$mu, cur$V, ctx$chart)
x0 <- list(cur = cur, alpha = matrix(chain$samples$alpha[, , idx], p),
  ll = chain$samples$subj_ll[, idx], rounds = 0L, wpool = NULL,
  wpool_part = NULL)

make_plan <- function(scale, t0_mean_scale, t0_hyper_scale) {
  plan <- chain$collapsed$plan
  plan$scale <- scale
  plan$t0_ratio <- t0_mean_scale / scale
  plan$t0_mean_ratio <- plan$t0_ratio
  plan$t0_hyper_ratio <- t0_hyper_scale / scale
  plan$mean_idx <- integer()
  plan$mean_ratio <- 1
  plan$length <- 3
  C <- EMC2:::.emc_collapsed_preconditioner(plan$C, plan$scale,
    plan$t0_ratio, plan$t0_idx, plan$t0_hyper_ratio, plan$t0_mean_idx)
  plan$eps <- EMC2:::.emc_with_preserved_rng(
    EMC2:::.emc_collapsed_tune_eps(phi0, C, plan$eps, ctx, plan$length))
  plan
}

plans <- list(
  unit = make_plan(1, 1, 1),
  fixed = make_plan(0.0441941738241592, 0.0441941738241592,
                    0.0110485434560398))

ns <- asNamespace("EMC2")
binding <- ".emc_wpool_subject_ll"
original <- get(binding, envir = ns, inherits = FALSE)
active_arm <- NULL
measure_likelihood <- FALSE
records <- list()
unlockBinding(binding, ns)
assign(binding, function(pool, rows, ctx, part = NULL, source) {
  cpu0 <- proc.time()
  wall0 <- proc.time()[["elapsed"]]
  out <- original(pool, rows, ctx, part, source)
  elapsed <- proc.time()[["elapsed"]] - wall0
  cpu <- proc.time() - cpu0
  if (measure_likelihood && identical(source, "collapsed")) {
    records[[length(records) + 1L]] <<- data.frame(
      arm = active_arm, subjects = length(rows),
      user_seconds = unname(cpu[["user.self"]]),
      system_seconds = unname(cpu[["sys.self"]]),
      elapsed_seconds = elapsed)
  }
  out
}, envir = ns)
lockBinding(binding, ns)

seed <- 20261010L
run_one <- function(arm, measure = TRUE) {
  active_arm <<- arm
  measure_likelihood <<- measure
  set.seed(seed)
  x <- x0
  x$moved <- FALSE
  started <- proc.time()
  out <- EMC2:::.emc_collapsed_round(x, plans[[arm]], chain, ctx)
  elapsed <- unname((proc.time() - started)[["elapsed"]])
  data.frame(arm = arm, candidate_acceptance_probability = out$acc,
    moved = isTRUE(out$moved), round_elapsed_seconds = elapsed)
}

# Warm each path once, then alternate measurement order on the same proposal.
invisible(run_one("unit", measure = FALSE))
invisible(run_one("fixed", measure = FALSE))
rounds <- 12L
out <- vector("list", 2L * rounds)
pos <- 0L
for (k in seq_len(rounds)) {
  arms <- if (k %% 2L) c("unit", "fixed") else c("fixed", "unit")
  for (arm in arms) {
    pos <- pos + 1L
    out[[pos]] <- run_one(arm)
  }
}
round_summary <- do.call(rbind, out)
eval_summary <- do.call(rbind, records)
summary <- do.call(rbind, lapply(split(round_summary, round_summary$arm),
  function(d) {
    e <- eval_summary[eval_summary$arm == d$arm[1L], , drop = FALSE]
    data.frame(arm = d$arm[1L], rounds = nrow(d), subjects_per_round = chain$n_subjects,
      exact_subject_rows = sum(e$subjects),
      same_rng_candidate_repeated = TRUE,
      candidate_exact_acceptance_probability = mean(d$candidate_acceptance_probability),
      observed_accepted_rounds = sum(d$moved), mean_exact_eval_wall_seconds = mean(e$elapsed_seconds),
      total_exact_eval_wall_seconds = sum(e$elapsed_seconds),
      mean_exact_eval_cpu_seconds = mean(e$user_seconds + e$system_seconds),
      total_exact_eval_cpu_seconds = sum(e$user_seconds + e$system_seconds),
      mean_round_wall_seconds = mean(d$round_elapsed_seconds),
      total_round_wall_seconds = sum(d$round_elapsed_seconds),
      step_size = plans[[d$arm[1L]]]$eps)
  }))
summary <- summary[match(c("unit", "fixed"), summary$arm), ]
write.csv(summary,
  "plans/collapsed_transport_adaptation_2026-10-09/zero_acceptance_round_cost.csv",
  row.names = FALSE)
print(summary, row.names = FALSE)
