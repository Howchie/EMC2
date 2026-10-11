# Compare subject-HMC trajectory lengths on frozen LBA subject conditionals.
# Each arm starts from the same retained state and subject RNG stream. The
# group state is fixed, so this isolates the subject conditional kernel.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)

snapshot <- Sys.getenv("HMC_PATH_SNAPSHOT",
  "/tmp/lba_sampling_audit/calibrated_snapshot.rds")
out_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
raw <- readRDS(snapshot)
emc <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(raw), quiet = TRUE)

transition_length <- function(a0, ll0, fr, ll_fn, eps_base, h, floor,
                              trajectory_length) {
  p <- fr$p
  U <- stats::runif(1L, 0.9, 1.1)
  v0 <- stats::rnorm(p)
  u_acc <- stats::runif(1L)
  eps <- eps_base * U
  L <- max(1L, as.integer(ceiling(trajectory_length / eps)))
  out <- list(alpha = a0, ll = ll0, accepted = FALSE, acc_prob = 0,
    disp = numeric(p), U = U, eps = eps, L = L, steps = 0L,
    rows = 0L, status = "abort", energy = NA_real_)
  cur <- EMC2:::.emc_hmc_eval(drop(fr$R %*% (a0 - fr$mu)), fr, ll_fn,
    h, floor, ll_centre = ll0)
  out$rows <- cur$rows
  if (!cur$ok) return(out)
  v <- v0 + 0.5 * eps * cur$g
  z <- cur$z
  for (l in seq_len(L)) {
    z <- z + eps * v
    new <- EMC2:::.emc_hmc_eval(z, fr, ll_fn, h, floor)
    out$rows <- out$rows + new$rows
    out$steps <- l
    if (!new$ok) return(out)
    v <- v + (if (l < L) eps else 0.5 * eps) * new$g
  }
  energy <- (-new$lp + 0.5 * sum(v^2)) -
    (-cur$lp + 0.5 * sum(v0^2))
  out$energy <- energy
  out$status <- "ok"
  out$acc_prob <- if (is.finite(energy)) min(1, exp(-energy)) else 0
  out$disp <- new$a - a0
  if (u_acc < out$acc_prob) {
    out$alpha <- new$a
    out$ll <- new$ll
    out$accepted <- TRUE
  }
  out
}

warmup_moves <- as.integer(Sys.getenv("HMC_PATH_WARMUP", "0"))
path_lengths <- if (warmup_moves > 0L) {
  c(short = pi / 4, default = pi / 2, long = 3 * pi / 4)
} else {
  c(very_short = pi / 8, short = pi / 4,
    medium = 3 * pi / 8, default = pi / 2)
}
subjects <- as.integer(Sys.getenv("HMC_PATH_SUBJECTS", "14"))
n_transitions <- as.integer(Sys.getenv("HMC_PATH_DRAWS", "400"))
chains <- as.integer(strsplit(Sys.getenv("HMC_PATH_CHAINS", "2,3"), ",",
                              fixed = TRUE)[[1L]])
rows <- list()
draws_out <- list()

for (j in chains) for (s in subjects) for (arm in names(path_lengths)) {
  ch <- emc[[j]]
  i <- ch$samples$idx
  prior <- EMC2:::.emc_burn_subject_prior(ch, s)
  frame <- EMC2:::.emc_hmc_frame(ch$lik_laplace[[s]]$prec,
                                  prior$mu, prior$var)
  alpha <- as.numeric(ch$samples$alpha[, s, i])
  ll <- as.numeric(ch$samples$subj_ll[s, i])
  data <- ch$data[[s]]
  par_names <- ch$par_names
  ll_fn <- function(X) {
    colnames(X) <- par_names
    as.numeric(calc_ll_manager(X, dadm = data, model = ch$model,
                               r_cores = 1L))
  }
  eps <- ch$subject_hmc$subjects[[s]]$eps
  eps_start <- eps
  seed <- ch$rng$subjects[[s]]
  assign(".Random.seed", seed, envir = globalenv())
  trajectory_length <- path_lengths[[arm]]
  alpha_draws <- matrix(NA_real_, nrow = n_transitions,
                        ncol = length(alpha),
                        dimnames = list(NULL, par_names))
  accepted <- logical(n_transitions)
  acc_prob <- numeric(n_transitions)
  abort_steps <- integer(n_transitions)
  n_rows <- numeric(n_transitions)
  n_steps <- integer(n_transitions)
  warmup_accept <- logical(warmup_moves)
  warmup_abort <- logical(warmup_moves)
  started <- proc.time()[["elapsed"]]
  warm_state <- list(eps = eps, n_warm = 0L, warmup = TRUE)
  if (warmup_moves > 0L) for (k in seq_len(warmup_moves)) {
    tr <- if (identical(arm, "default")) {
      EMC2:::.emc_hmc_transition(alpha, ll, frame, ll_fn, warm_state$eps,
        1e-4, EMC2:::.emc_hmc_ll_floor(data))
    } else {
      transition_length(alpha, ll, frame, ll_fn, warm_state$eps, 1e-4,
        EMC2:::.emc_hmc_ll_floor(data), trajectory_length)
    }
    warmup_accept[k] <- tr$accepted
    warmup_abort[k] <- tr$status != "ok"
    if (tr$accepted) {
      alpha <- tr$alpha
      ll <- tr$ll
    }
    warm_state <- EMC2:::.emc_hmc_warmup_step(warm_state, tr)
  }
  eps <- warm_state$eps
  for (k in seq_len(n_transitions)) {
    tr <- if (identical(arm, "default")) {
      EMC2:::.emc_hmc_transition(alpha, ll, frame, ll_fn, eps, 1e-4,
        EMC2:::.emc_hmc_ll_floor(data))
    } else {
      transition_length(alpha, ll, frame, ll_fn, eps, 1e-4,
        EMC2:::.emc_hmc_ll_floor(data), trajectory_length)
    }
    if (tr$accepted) {
      alpha <- tr$alpha
      ll <- tr$ll
    }
    alpha_draws[k, ] <- alpha
    accepted[k] <- tr$accepted
    acc_prob[k] <- tr$acc_prob
    abort_steps[k] <- if (tr$status == "ok") NA_integer_ else tr$steps
    n_rows[k] <- tr$rows
    n_steps[k] <- tr$L
  }
  elapsed <- proc.time()[["elapsed"]] - started
  arr <- array(alpha_draws, dim = c(n_transitions, 1L, length(par_names)),
               dimnames = list(NULL, NULL, par_names))
  ess_summary <- posterior::summarise_draws(
    posterior::as_draws_array(arr), posterior::ess_bulk)
  ess <- stats::setNames(ess_summary[[2L]], ess_summary$variable)
  t0_names <- startsWith(names(ess), "t0")
  rows[[length(rows) + 1L]] <- data.frame(
    chain = j, subject = s, arm = arm,
    requested_length = trajectory_length,
    eps_start = eps_start, eps = eps, warmup_moves = warmup_moves,
    warmup_acceptance = if (warmup_moves) mean(warmup_accept) else NA_real_,
    warmup_abort_rate = if (warmup_moves) mean(warmup_abort) else NA_real_,
    transitions = n_transitions,
    seconds = elapsed,
    rows = sum(n_rows), mean_rows = mean(n_rows),
    mean_leapfrog_steps = mean(n_steps),
    accepted = mean(accepted), mean_acceptance_probability = mean(acc_prob),
    abort_rate = mean(!is.na(abort_steps)),
    abort_step_median = if (all(is.na(abort_steps))) NA_real_ else
      median(abort_steps, na.rm = TRUE),
    median_parameter_ess_per_second = median(ess / elapsed),
    median_t0_ess_per_second = if (any(t0_names))
      median(ess[t0_names] / elapsed) else NA_real_,
    sum_parameter_ess_per_second = sum(ess) / elapsed,
    stringsAsFactors = FALSE)
  draws_out[[paste(j, s, arm, sep = "_")]] <-
    list(draws = alpha_draws, ess = ess)
}

summary <- do.call(rbind, rows)
suffix <- if (warmup_moves > 0L) "warmup" else "refined"
suffix <- paste0(suffix, "_s", paste(subjects, collapse = "_"))
write.csv(summary, file.path(out_dir,
  paste0("hmc_path_length_", suffix, "_summary.csv")),
          row.names = FALSE)
saveRDS(list(summary = summary, draws = draws_out),
        file.path(out_dir, paste0("hmc_path_length_", suffix, "_probe.rds")))
print(summary, row.names = FALSE)
