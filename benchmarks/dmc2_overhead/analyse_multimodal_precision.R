args <- commandArgs(TRUE)
if (length(args) != 1L)
  stop("usage: analyse_multimodal_precision.R results-dir")
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
seeds <- as.integer(strsplit(Sys.getenv("SEEDS", "101 202 303"), "\\s+")[[1L]])
arms <- c("base", "new")
fits <- setNames(vector("list", length(arms)), arms)
seed_fits <- setNames(vector("list", length(arms)), arms)
times <- setNames(vector("list", length(arms)), arms)
seed_times <- setNames(vector("list", length(arms)), arms)
mode_draws <- setNames(vector("list", length(arms)), arms)
theta_draws <- setNames(vector("list", length(arms)), arms)
switch_counts <- setNames(vector("list", length(arms)), arms)

for (arm in arms) {
  for (seed in seeds) {
    path <- file.path(root, sprintf("multimodal_%d_%s.rds", seed, arm))
    x <- readRDS(path)
    stopifnot(identical(x$sampler_seed, seed),
              identical(as.numeric(x$cache_max_reuses),
                        if (arm == "base") 0 else Inf),
              all(vapply(x$fit, function(ch)
                dim(ch$samples$alpha)[3L] == x$sample_draws, logical(1))))
    fits[[arm]] <- c(fits[[arm]], x$fit)
    seed_fits[[arm]][[as.character(seed)]] <- x$fit
    times[[arm]] <- c(times[[arm]], x$timing[["elapsed"]])
    seed_times[[arm]][[as.character(seed)]] <- x$timing[["elapsed"]]
  }
  chains <- fits[[arm]]
  mode_by_chain <- lapply(chains, function(ch)
    as.numeric(ch$samples$alpha[1L, 1L, ] > 0))
  theta_by_chain <- lapply(chains, function(ch)
    as.numeric(ch$samples$alpha[1L, 1L, ]))
  n <- min(lengths(mode_by_chain))
  mode_arr <- array(NA_real_, c(n, length(chains), 1L),
    dimnames = list(NULL, NULL, "right_mode"))
  theta_arr <- array(NA_real_, c(n, length(chains), 1L),
    dimnames = list(NULL, NULL, "theta"))
  for (j in seq_along(chains)) {
    mode_arr[, j, 1L] <- mode_by_chain[[j]][seq_len(n)]
    theta_arr[, j, 1L] <- theta_by_chain[[j]][seq_len(n)]
  }
  mode_draws[[arm]] <- as_draws_array(mode_arr)
  theta_draws[[arm]] <- as_draws_array(theta_arr)
  switch_counts[[arm]] <- vapply(mode_by_chain, function(z) sum(diff(z) != 0), integer(1))
}

mode_stats <- lapply(mode_draws, function(x)
  summarise_draws(x, "mean", "mcse_mean", "ess_bulk", "rhat"))
theta_stats <- lapply(theta_draws, function(x)
  summarise_draws(x, "ess_bulk", "ess_tail",
    ~quantile2(.x, probs = c(.05, .5, .95)),
    ~setNames(vapply(c(.05, .5, .95), function(p)
      mcse_quantile(.x, probs = p), numeric(1)),
      c("mcse_q05", "mcse_q50", "mcse_q95"))))

for (arm in arms) {
  mode <- mode_stats[[arm]]
  theta <- theta_stats[[arm]]
  cat(arm, "chains", length(fits[[arm]]), "draws/chain",
      dim(mode_draws[[arm]])[1L], "total elapsed", sum(times[[arm]]),
      "mean elapsed", mean(times[[arm]]), "\n")
  cat("  right-mode mass", mode$mean, "MCSE", mode$mcse_mean,
      "distance from exact 0.5 in MCSE", abs(mode$mean - 0.5) / mode$mcse_mean,
      "ESS", mode$ess_bulk, "Rhat", mode$rhat, "\n")
  cat("  theta 5%/95% quantiles", theta$q5, theta$q95,
      "bulk/tail ESS", theta$ess_bulk, theta$ess_tail, "\n")
  cat("  mode switches per chain", paste(switch_counts[[arm]], collapse = ","), "\n")
}

mode_ess_ratios <- theta_bulk_ratios <- theta_tail_ratios <- numeric(length(seeds))
for (seed_i in seq_along(seeds)) {
  seed <- seeds[seed_i]
  get_seed_mode <- function(fit) {
    chains <- lapply(fit, function(ch) as.numeric(ch$samples$alpha[1L, 1L, ] > 0))
    n <- min(lengths(chains))
    arr <- array(NA_real_, c(n, length(chains), 1L),
      dimnames = list(NULL, NULL, "right_mode"))
    for (j in seq_along(chains)) arr[, j, 1L] <- chains[[j]][seq_len(n)]
    summarise_draws(as_draws_array(arr), "mean", "mcse_mean", "ess_bulk")
  }
  get_seed_theta <- function(fit) {
    chains <- lapply(fit, function(ch) as.numeric(ch$samples$alpha[1L, 1L, ]))
    n <- min(lengths(chains))
    arr <- array(NA_real_, c(n, length(chains), 1L),
      dimnames = list(NULL, NULL, "theta"))
    for (j in seq_along(chains)) arr[, j, 1L] <- chains[[j]][seq_len(n)]
    summarise_draws(as_draws_array(arr), "ess_bulk", "ess_tail")
  }
  a <- get_seed_mode(seed_fits$base[[as.character(seed)]])
  b <- get_seed_mode(seed_fits$new[[as.character(seed)]])
  z <- (b$mean - a$mean) / sqrt(a$mcse_mean^2 + b$mcse_mean^2)
  theta_a <- get_seed_theta(seed_fits$base[[as.character(seed)]])
  theta_b <- get_seed_theta(seed_fits$new[[as.character(seed)]])
  mode_ess_ratios[seed_i] <- b$ess_bulk / a$ess_bulk
  theta_bulk_ratios[seed_i] <- theta_b$ess_bulk / theta_a$ess_bulk
  theta_tail_ratios[seed_i] <- theta_b$ess_tail / theta_a$ess_tail
  cat("seed", seed, "right-mode mass base/cache", a$mean, b$mean,
      "paired difference in combined MCSE", z,
      "mode ESS ratio", b$ess_bulk / a$ess_bulk,
      "theta bulk/tail ESS ratios", theta_b$ess_bulk / theta_a$ess_bulk,
      theta_b$ess_tail / theta_a$ess_tail,
      "elapsed base/cache", seed_times$base[[as.character(seed)]],
      seed_times$new[[as.character(seed)]], "\n")
}
ratio_summary <- function(x) {
  z <- log(x)
  c(geomean = exp(mean(z)),
    lower = exp(mean(z) - qt(.975, length(z) - 1L) * sd(z) / sqrt(length(z))),
    upper = exp(mean(z) + qt(.975, length(z) - 1L) * sd(z) / sqrt(length(z))))
}
cat("across-seed geometric mean and 95% t intervals, mode-indicator ESS ratio:",
    paste(signif(ratio_summary(mode_ess_ratios), 4), collapse = ","), "\n")
cat("across-seed geometric mean and 95% t intervals, theta bulk ESS ratio:",
    paste(signif(ratio_summary(theta_bulk_ratios), 4), collapse = ","),
    "tail:", paste(signif(ratio_summary(theta_tail_ratios), 4), collapse = ","), "\n")

base_mode <- mode_stats$base
new_mode <- mode_stats$new
mode_z <- abs(base_mode$mean - new_mode$mean) /
  sqrt(base_mode$mcse_mean^2 + new_mode$mcse_mean^2)
base_theta <- theta_stats$base
new_theta <- theta_stats$new
q_names <- c("q5", "q95")
mcse_names <- c("mcse_q05", "mcse_q95")
quantile_z <- vapply(seq_along(q_names), function(i)
  (new_theta[[q_names[i]]] - base_theta[[q_names[i]]]) /
    sqrt(new_theta[[mcse_names[i]]]^2 + base_theta[[mcse_names[i]]]^2), numeric(1))
cat("cache minus fresh right-mode mass, combined MCSE units", mode_z, "\n")
cat("theta quantile z scores (5%, 95%)",
    paste(signif(quantile_z, 4), collapse = ","), "\n")
cat("candidate/baseline theta bulk and tail ESS-per-draw ratios",
    new_theta$ess_bulk / base_theta$ess_bulk,
    new_theta$ess_tail / base_theta$ess_tail, "\n")
speedups <- times$base / times$new
log_speedups <- log(speedups)
speed_ci <- exp(mean(log_speedups) + c(-1, 1) *
  qt(.975, length(speedups) - 1L) * sd(log_speedups) / sqrt(length(speedups)))
cat("paired speedups", paste(signif(speedups, 4), collapse = ","),
    "geometric mean", exp(mean(log_speedups)),
    "95% t interval", paste(signif(speed_ci, 4), collapse = " to "), "\n")
