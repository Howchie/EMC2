# Check whether a small exact pilot ranks subject-HMC path lengths by the
# slowest standardized coordinate movement per unit CPU.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
snapshot <- Sys.getenv("HMC_PATH_SNAPSHOT",
  "/tmp/lba_sampling_audit/calibrated_snapshot.rds")
out_dir <- "plans/collapsed_transport_adaptation_2026-10-09"
raw <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(snapshot)), quiet = TRUE)
paths <- c(short = pi / 4, default = pi / 2, long = 3 * pi / 4)
chains <- as.integer(strsplit(Sys.getenv("HMC_PATH_CHAINS", "2,3"), ",",
                              fixed = TRUE)[[1L]])
subjects <- as.integer(strsplit(Sys.getenv("HMC_PATH_SUBJECTS", "14,48"), ",",
                                fixed = TRUE)[[1L]])
n_pilots <- as.integer(Sys.getenv("HMC_PATH_PILOTS", "7"))

transition_length <- function(a0, ll0, fr, ll_fn, eps_base, h, floor,
                              trajectory_length) {
  p <- fr$p
  U <- stats::runif(1L, 0.9, 1.1)
  v0 <- stats::rnorm(p)
  u_acc <- stats::runif(1L)
  eps <- eps_base * U
  L <- max(1L, as.integer(ceiling(trajectory_length / eps)))
  out <- list(alpha = a0, ll = ll0, accepted = FALSE, acc_prob = 0,
              disp = numeric(p), rows = 0L, status = "abort")
  cur <- EMC2:::.emc_hmc_eval(drop(fr$R %*% (a0 - fr$mu)), fr,
    ll_fn, h, floor, ll_centre = ll0)
  out$rows <- cur$rows
  if (!cur$ok) return(out)
  v <- v0 + 0.5 * eps * cur$g
  z <- cur$z
  for (l in seq_len(L)) {
    z <- z + eps * v
    new <- EMC2:::.emc_hmc_eval(z, fr, ll_fn, h, floor)
    out$rows <- out$rows + new$rows
    if (!new$ok) return(out)
    v <- v + (if (l < L) eps else 0.5 * eps) * new$g
  }
  energy <- (-new$lp + 0.5 * sum(v^2)) -
    (-cur$lp + 0.5 * sum(v0^2))
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

rows <- list()
for (j in chains) for (s in subjects) {
  ch <- raw[[j]]
  i <- ch$samples$idx
  prior <- EMC2:::.emc_burn_subject_prior(ch, s)
  frame <- EMC2:::.emc_hmc_frame(ch$lik_laplace[[s]]$prec,
                                  prior$mu, prior$var)
  alpha0 <- as.numeric(ch$samples$alpha[, s, i])
  ll0 <- as.numeric(ch$samples$subj_ll[s, i])
  data <- ch$data[[s]]
  par_names <- ch$par_names
  ll_fn <- function(X) {
    colnames(X) <- par_names
    as.numeric(calc_ll_manager(X, dadm = data, model = ch$model,
                               r_cores = 1L))
  }
  eps <- ch$subject_hmc$subjects[[s]]$eps
  posterior_idx <- seq.int(max(1L, i - 199L), i)
  variances <- apply(ch$samples$alpha[, s, posterior_idx, drop = FALSE], 1L,
                     stats::var)
  variances[!is.finite(variances) | variances <= .Machine$double.eps] <- 1
  t0 <- startsWith(par_names, "t0")
  for (arm in names(paths)) {
    assign(".Random.seed", ch$rng$subjects[[s]], envir = globalenv())
    alpha <- alpha0
    ll <- ll0
    Mdiag <- numeric(length(alpha))
    cpu <- rows_used <- 0
    accepted <- aborted <- 0L
    for (k in seq_len(n_pilots)) {
      started <- proc.time()
      tr <- if (identical(arm, "default")) {
        EMC2:::.emc_hmc_transition(alpha, ll, frame, ll_fn, eps, 1e-4,
          EMC2:::.emc_hmc_ll_floor(data))
      } else {
        transition_length(alpha, ll, frame, ll_fn, eps, 1e-4,
          EMC2:::.emc_hmc_ll_floor(data), paths[[arm]])
      }
      spent <- proc.time() - started
      cpu <- cpu + unname(spent[["user.self"]] + spent[["sys.self"]])
      rows_used <- rows_used + tr$rows
      Mdiag <- Mdiag + tr$acc_prob * tr$disp^2
      accepted <- accepted + as.integer(tr$accepted)
      aborted <- aborted + as.integer(tr$status != "ok")
      if (tr$accepted) {
        alpha <- tr$alpha
        ll <- tr$ll
      }
    }
    coord <- Mdiag / (n_pilots * variances)
    rows[[length(rows) + 1L]] <- data.frame(
      chain = j, subject = s, arm = arm, path_length = paths[[arm]],
      eps = eps, pilots = n_pilots, cpu_seconds = cpu, rows = rows_used,
      acceptance = accepted / n_pilots, abort_rate = aborted / n_pilots,
      min_coordinate_score_per_cpu = min(coord) / max(cpu, .Machine$double.eps),
      min_t0_score_per_cpu = if (any(t0))
        min(coord[t0]) / max(cpu, .Machine$double.eps) else NA_real_,
      min_coordinate_score_per_row = min(coord) / max(rows_used, 1),
      min_t0_score_per_row = if (any(t0))
        min(coord[t0]) / max(rows_used, 1) else NA_real_,
      stringsAsFactors = FALSE)
  }
}
summary <- do.call(rbind, rows)
write.csv(summary, file.path(out_dir, "hmc_path_length_pilot_scores.csv"),
          row.names = FALSE)
print(summary, row.names = FALSE)
