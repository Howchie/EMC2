args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: analyse_lnr_warm_start.R results-dir")
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
probs <- c(.05, .5, .95)
seeds <- as.integer(strsplit(Sys.getenv("SEEDS", "101 202 303"), "\\s+")[[1L]])

draws_by_chain <- function(fit) {
  lapply(fit, function(ch) {
    keep <- which(ch$samples$stage == "sample")
    alpha <- ch$samples$alpha[, 1L, keep, drop = FALSE]
    out <- t(alpha[, 1L, ])
    colnames(out) <- dimnames(ch$samples$alpha)[[1L]]
    out
  })
}

summarise_fit <- function(fit) {
  chains <- draws_by_chain(fit)
  n <- min(vapply(chains, nrow, integer(1)))
  vars <- colnames(chains[[1L]])
  draws <- array(NA_real_, c(n, length(chains), length(vars)),
                 dimnames = list(NULL, NULL, vars))
  for (j in seq_along(chains)) draws[, j, ] <- chains[[j]][seq_len(n), vars]
  stats <- summarise_draws(as_draws_array(draws), "ess_bulk", "ess_tail",
    ~quantile2(.x, probs = probs),
    ~setNames(vapply(probs, function(p) mcse_quantile(.x, probs = p), numeric(1)),
              c("mcse_q05", "mcse_q50", "mcse_q95")))
  list(chains = chains, n = n, stats = stats)
}

base_times <- new_times <- time_ratios <- bulk_ratios <- tail_ratios <- numeric()
for (seed in seeds) {
  prefix <- file.path(root, paste0("LNR_", seed))
  base <- readRDS(paste0(prefix, "_base_precision.rds"))
  new <- readRDS(paste0(prefix, "_new_precision.rds"))
  stopifnot(identical(base$data_md5, new$data_md5),
            identical(base$initial_seed, new$initial_seed))
  sb <- summarise_fit(base$fit)
  sn <- summarise_fit(new$fit)
  diff <- do.call(rbind, lapply(seq_along(probs), function(i) {
    q <- c("q5", "q50", "q95")[[i]]
    mcse <- c("mcse_q05", "mcse_q50", "mcse_q95")[[i]]
    data.frame(variable = sb$stats$variable,
      quantile = c("q5", "q50", "q95")[[i]],
      standardized_difference = abs(sb$stats[[q]] - sn$stats[[q]]) /
        sqrt(sb$stats[[mcse]]^2 + sn$stats[[mcse]]^2))
  }))
  fields <- c("alpha", "stage", "subj_ll", "idx")
  samples_equal <- vapply(seq_along(base$fit), function(j)
    all(vapply(fields, function(f)
      identical(base$fit[[j]]$samples[[f]], new$fit[[j]]$samples[[f]]), logical(1))),
    logical(1))
  rng_equal <- vapply(seq_along(base$fit), function(j)
    identical(base$fit[[j]]$rng, new$fit[[j]]$rng), logical(1))
  bulk_ratios <- c(bulk_ratios,
    (min(sn$stats$ess_bulk) / (sn$n * length(sn$chains))) /
      (min(sb$stats$ess_bulk) / (sb$n * length(sb$chains))))
  tail_ratios <- c(tail_ratios,
    (min(sn$stats$ess_tail) / (sn$n * length(sn$chains))) /
      (min(sb$stats$ess_tail) / (sb$n * length(sb$chains))))
  base_times <- c(base_times, base$timing[["elapsed"]])
  new_times <- c(new_times, new$timing[["elapsed"]])
  time_ratios <- c(time_ratios, base$timing[["elapsed"]] /
                                  new$timing[["elapsed"]])
  cat("\nLNR seed", seed,
      "elapsed base/candidate", base_times[[length(base_times)]],
      new_times[[length(new_times)]],
      "draws/chain base/candidate", sb$n, sn$n,
      "sample equality", samples_equal,
      "chain RNG equality", rng_equal,
      "caller RNG equality", identical(base$final_seed, new$final_seed),
      "precision target base/candidate",
      min(sb$stats$ess_bulk) > 800, min(sn$stats$ess_bulk) > 800,
      "quantile failures", sum(diff$standardized_difference > 4),
      "largest quantile difference", max(diff$standardized_difference), "\n")
  cat("  min bulk ESS/draw base/candidate",
      min(sb$stats$ess_bulk) / (sb$n * length(sb$chains)),
      min(sn$stats$ess_bulk) / (sn$n * length(sn$chains)),
      "tail",
      min(sb$stats$ess_tail) / (sb$n * length(sb$chains)),
      min(sn$stats$ess_tail) / (sn$n * length(sn$chains)), "\n")
  cat("  min bulk ESS/s base/candidate",
      min(sb$stats$ess_bulk) / base$timing[["elapsed"]],
      min(sn$stats$ess_bulk) / new$timing[["elapsed"]],
      "tail",
      min(sb$stats$ess_tail) / base$timing[["elapsed"]],
      min(sn$stats$ess_tail) / new$timing[["elapsed"]], "\n")
}
cat("\nMean elapsed base/candidate", mean(base_times), mean(new_times),
    "speedup", mean(base_times) / mean(new_times), "\n")
time_ci <- exp(mean(log(time_ratios)) + c(-1, 1) *
  qt(.975, df = length(time_ratios) - 1L) * sd(log(time_ratios)) /
    sqrt(length(time_ratios)))
cat("paired time speedup geometric mean", exp(mean(log(time_ratios))),
    "95% t interval", time_ci[[1L]], time_ci[[2L]],
    "seed ratios", time_ratios, "\n")
for (kind in c("bulk", "tail")) {
  x <- if (kind == "bulk") bulk_ratios else tail_ratios
  ci <- exp(mean(log(x)) + c(-1, 1) *
    qt(.975, df = length(x) - 1L) * sd(log(x)) / sqrt(length(x)))
  cat(kind, "ESS minimum ratio geometric mean", exp(mean(log(x))),
      "95% t interval", ci[[1L]], ci[[2L]], "seed ratios", x, "\n")
}
