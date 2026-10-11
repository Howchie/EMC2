pkgload::load_all(".", quiet = TRUE, compile = FALSE)

artifact_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
snapshot <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")),
  quiet = TRUE)
n_draws <- 400L
warmup_moves <- as.integer(Sys.getenv(
  "HMC_T0_WARMUP_MOVES", as.character(EMC2:::.EMC_HMC_WARMUP_MOVES)))
stopifnot(length(warmup_moves) == 1L, is.finite(warmup_moves), warmup_moves >= 1L)
subject_id <- "52091"
coupled_names <- c("t0_LoadLow", "t0_LoadHigh", "t0_itemTypeTarget",
  "B", "B_lR_target", "B_LoadLow", "B.w_prac", "sv_lMTRUE", "A")

conditional_prior <- function(mu, V, block, alpha) {
  rest <- setdiff(seq_along(mu), block)
  if (!length(rest)) return(list(mu = mu[block], var = V[block, block, drop = FALSE]))
  Vrr <- V[rest, rest, drop = FALSE]
  Vbr <- V[block, rest, drop = FALSE]
  Vrb <- V[rest, block, drop = FALSE]
  list(
    mu = drop(mu[block] + Vbr %*% solve(Vrr, alpha[rest] - mu[rest])),
    var = V[block, block, drop = FALSE] - Vbr %*% solve(Vrr, Vrb)
  )
}

run_arm <- function(ch, subject, blocks, seed, arm, chain_id) {
  p <- ch$n_pars
  pn <- ch$par_names
  idx <- ch$samples$idx
  hmc <- ch$subject_hmc$subjects[[subject]]
  prior <- EMC2:::.emc_burn_subject_prior(ch, subject)
  H <- hmc$prec
  V <- as.matrix(prior$var)
  mu <- as.numeric(prior$mu)
  alpha0 <- as.numeric(ch$samples$alpha[, subject, idx])
  ll0 <- as.numeric(ch$samples$subj_ll[subject, idx])
  floor <- EMC2:::.emc_hmc_ll_floor(ch$data[[subject]])
  h <- hmc$h

  make_frame <- function(block, alpha) {
    cp <- conditional_prior(mu, V, block, alpha)
    EMC2:::.emc_hmc_frame(H[block, block, drop = FALSE], cp$mu, cp$var)
  }
  make_ll <- function(block, alpha) function(Xb) {
    X <- matrix(rep(alpha, times = nrow(Xb)), nrow = nrow(Xb), byrow = TRUE)
    X[, block] <- Xb
    colnames(X) <- pn
    as.numeric(calc_ll_manager(X, dadm = ch$data[[subject]], model = ch$model,
                               r_cores = 1L))
  }

  set.seed(seed)
  started <- proc.time()[["elapsed"]]
  alpha <- alpha0
  ll <- ll0
  eps <- numeric(length(blocks))
  warmup_rows <- 0L
  for (b in seq_along(blocks)) {
    block <- blocks[[b]]
    ps <- list(eps = EMC2:::.emc_hmc_clip_eps(length(block)^(-1 / 4)),
               n_warm = 0L, warmup = TRUE)
    for (step in seq_len(warmup_moves)) {
      fr <- make_frame(block, alpha)
      ll_fn <- make_ll(block, alpha)
      tr <- EMC2:::.emc_hmc_transition(alpha[block], ll, fr, ll_fn,
                                       ps$eps, h, floor)
      warmup_rows <- warmup_rows + tr$rows
      if (isTRUE(tr$accepted)) {
        alpha[block] <- tr$alpha
        ll <- tr$ll
      }
      ps <- EMC2:::.emc_hmc_warmup_step(ps, list(acc_prob = tr$acc_prob))
      # The production warmup normally ends after 20 transitions.  This
      # benchmark can extend the same adaptation rule to test whether that
      # short window leaves support-sensitive subjects badly tuned.
      ps$warmup <- step < warmup_moves
      ps$rate <- if (ps$warmup) 1 else 0
    }
    eps[b] <- ps$eps
  }

  draws <- matrix(NA_real_, n_draws, p,
    dimnames = list(NULL, pn))
  accepted <- aborted <- rows <- numeric(length(blocks))
  for (iter in seq_len(n_draws)) {
    for (b in seq_along(blocks)) {
      block <- blocks[[b]]
      fr <- make_frame(block, alpha)
      ll_fn <- make_ll(block, alpha)
      tr <- EMC2:::.emc_hmc_transition(alpha[block], ll, fr, ll_fn,
                                       eps[b], h, floor)
      accepted[b] <- accepted[b] + as.integer(isTRUE(tr$accepted))
      aborted[b] <- aborted[b] + as.integer(!identical(tr$status, "ok"))
      rows[b] <- rows[b] + tr$rows
      if (isTRUE(tr$accepted)) {
        alpha[block] <- tr$alpha
        ll <- tr$ll
      }
    }
    draws[iter, ] <- alpha
  }
  elapsed <- proc.time()[["elapsed"]] - started
  arr <- posterior::as_draws_array(array(draws,
    dim = c(n_draws, 1L, p), dimnames = list(NULL, NULL, pn)))
  ess_tbl <- posterior::summarise_draws(arr, posterior::ess_bulk)
  ec <- grep("ess_bulk", names(ess_tbl), value = TRUE)[1L]
  ess <- stats::setNames(ess_tbl[[ec]], ess_tbl$variable)
  data.frame(arm = arm, chain = chain_id,
    subject = subject_id, draws = n_draws, block_sizes =
      paste(lengths(blocks), collapse = "+"), eps = paste(signif(eps, 4), collapse = "+"),
    block_acceptance = paste(round(accepted / n_draws, 3), collapse = "+"),
    block_aborts = paste(aborted, collapse = "+"), mean_rows_per_iter =
    (sum(rows) + warmup_rows) / n_draws, tuning_rows = warmup_rows,
    warmup_moves = warmup_moves,
    seconds = elapsed, sum_ess_per_second = sum(ess) / elapsed,
    median_ess_per_second = median(ess) / elapsed,
    min_ess_per_second = min(ess) / elapsed,
    median_t0_ess_per_second = median(ess[startsWith(names(ess), "t0_")]) / elapsed,
    max_t0_ess_per_second = max(ess[startsWith(names(ess), "t0_")]) / elapsed,
    setNames(as.list(ess), paste0("ess_", names(ess))),
    check.names = FALSE)
}

all_rows <- list()
for (j in c(2L, 3L)) {
  ch <- snapshot[[j]]
  subject <- match(subject_id, dimnames(ch$samples$alpha)[[2L]])
  stopifnot(!is.na(subject))
  p <- ch$n_pars
  coupled <- match(coupled_names, ch$par_names)
  stopifnot(!anyNA(coupled))
  other <- setdiff(seq_len(p), coupled)
  arms <- list(full = list(seq_len(p)), t0_coupled = list(coupled, other))
  for (arm in names(arms)) {
    cat("chain", j, "arm", arm, "\n")
    all_rows[[paste(j, arm, sep = "_")]] <- run_arm(
      ch, subject, arms[[arm]], seed = 20261010L + j, arm = arm,
      chain_id = j)
    flush.console()
  }
}
summary <- do.call(rbind, all_rows)
suffix <- if (warmup_moves == EMC2:::.EMC_HMC_WARMUP_MOVES) "" else
  paste0("_warmup", warmup_moves)
write.csv(summary, file.path(artifact_dir,
  paste0("hmc_t0_coupled_block_summary", suffix, ".csv")), row.names = FALSE)
saveRDS(summary, file.path(artifact_dir,
  paste0("hmc_t0_coupled_block_summary", suffix, ".rds")))
print(summary[, c("arm", "chain", "block_sizes", "eps", "block_acceptance",
                  "block_aborts", "mean_rows_per_iter", "seconds",
                  "sum_ess_per_second", "median_ess_per_second",
                  "min_ess_per_second", "median_t0_ess_per_second")],
      row.names = FALSE)
