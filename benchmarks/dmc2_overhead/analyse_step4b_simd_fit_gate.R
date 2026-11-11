args <- commandArgs(TRUE)
if (length(args) != 1L)
  stop("usage: analyse_step4b_simd_fit_gate.R results-dir")
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
sample_fields <- c("alpha", "stage", "subj_ll", "idx")

as_draw_array <- function(fit) {
  arr <- simplify2array(lapply(fit, function(ch) ch$samples$alpha[, 1, ]))
  arr <- aperm(arr, c(2, 3, 1))
  dimnames(arr) <- list(NULL, NULL, fit[[1L]]$par_names)
  arr
}

as_one_chain_array <- function(chain) {
  x <- aperm(chain$samples$alpha, c(3, 2, 1))
  dimnames(x) <- list(NULL, NULL, chain$par_names)
  x
}

draw_summary <- function(arr) {
  summarise_draws(as_draws_array(arr), "ess_bulk", "ess_tail",
    ~quantile2(.x, probs = c(.05, .5, .95)),
    ~setNames(vapply(c(.05, .5, .95), function(prob)
      mcse_quantile(.x, probs = prob), numeric(1)),
      c("mcse_q05", "mcse_q50", "mcse_q95")))
}

quantile_z <- function(base, candidate) {
  mcse_names <- c("mcse_q05", "mcse_q50", "mcse_q95")
  q_names <- c("q05", "q50", "q95")
  unlist(lapply(seq_along(q_names), function(i) {
    (candidate[[q_names[i]]] - base[[q_names[i]]]) /
      sqrt(candidate[[mcse_names[i]]]^2 + base[[mcse_names[i]]]^2)
  }))
}

seeds <- as.integer(strsplit(Sys.getenv("SEEDS", "101 202 303"), "\\s+")[[1L]])
summaries <- list()
for (seed in seeds) {
  prefix <- file.path(root, paste("EXG", seed, sep = "_"))
  base <- readRDS(paste0(prefix, "_base.rds"))
  new <- readRDS(paste0(prefix, "_new.rds"))
  stopifnot(identical(base$data_md5, new$data_md5),
            identical(base$initial_seed, new$initial_seed))
  draws_equal <- vapply(seq_along(base$fit), function(j)
    all(vapply(sample_fields, function(field)
      identical(base$fit[[j]]$samples[[field]], new$fit[[j]]$samples[[field]]),
      logical(1))), logical(1))
  rng_equal <- vapply(seq_along(base$fit), function(j)
    identical(base$fit[[j]]$rng, new$fit[[j]]$rng), logical(1))
  baseline <- draw_summary(as_draw_array(base$fit))
  candidate <- draw_summary(as_draw_array(new$fit))
  z <- quantile_z(baseline, candidate)
  chain_z <- vapply(seq_along(base$fit), function(j) {
    max(abs(quantile_z(draw_summary(as_one_chain_array(base$fit[[j]])),
                       draw_summary(as_one_chain_array(new$fit[[j]])))))
  }, numeric(1))
  chain_ess_ratio <- function(field) vapply(seq_along(base$fit), function(j) {
    old <- draw_summary(as_one_chain_array(base$fit[[j]]))[[field]]
    candidate <- draw_summary(as_one_chain_array(new$fit[[j]]))[[field]]
    min(candidate) / min(old)
  }, numeric(1))
  chain_bulk_ratio <- chain_ess_ratio("ess_bulk")
  chain_tail_ratio <- chain_ess_ratio("ess_tail")
  base_bulk <- min(baseline$ess_bulk) / length(base$fit[[1L]]$samples$stage)
  new_bulk <- min(candidate$ess_bulk) / length(new$fit[[1L]]$samples$stage)
  base_tail <- min(baseline$ess_tail) / length(base$fit[[1L]]$samples$stage)
  new_tail <- min(candidate$ess_tail) / length(new$fit[[1L]]$samples$stage)
  speedup <- base$timing[["elapsed"]] / new$timing[["elapsed"]]
  summaries[[as.character(seed)]] <- c(speedup = speedup,
    bulk_ratio = new_bulk / base_bulk, tail_ratio = new_tail / base_tail,
    min_chain_bulk_ratio = min(chain_bulk_ratio),
    min_chain_tail_ratio = min(chain_tail_ratio))
  cat("EXG", seed, "elapsed baseline/candidate",
      base$timing[["elapsed"]], new$timing[["elapsed"]],
      "speedup", speedup, "draws equal", paste(draws_equal, collapse = ","),
      "chain RNG equal", paste(rng_equal, collapse = ","),
      "caller RNG equal", identical(base$final_seed, new$final_seed),
      "max |quantile z| combined/per chain", max(abs(z)),
      paste(signif(chain_z, 4), collapse = ","),
      "minimum bulk ESS/draw ratio", new_bulk / base_bulk,
      "minimum tail ESS/draw ratio", new_tail / base_tail,
      "per-chain minimum bulk/tail ratios",
      paste(signif(chain_bulk_ratio, 4), collapse = ","), "/",
      paste(signif(chain_tail_ratio, 4), collapse = ","), "\n")
}
summaries <- do.call(rbind, summaries)
cat("Geometric mean speedup:", exp(mean(log(summaries[, "speedup"]))), "\n")
log_speed <- log(summaries[, "speedup"])
speed_ci <- exp(mean(log_speed) + c(-1, 1) *
  qt(.975, df = length(log_speed) - 1) * sd(log_speed) / sqrt(length(log_speed)))
cat("Geometric mean speedup 95% t interval:",
    paste(signif(speed_ci, 4), collapse = " to "), "\n")
cat("Paired speedup ratios:", paste(signif(summaries[, "speedup"], 4),
    collapse = ", "), "\n")
cat("Minimum ESS/draw ratios, bulk:",
    paste(signif(summaries[, "bulk_ratio"], 4), collapse = ", "),
    "tail:", paste(signif(summaries[, "tail_ratio"], 4), collapse = ", "), "\n")
cat("Per-chain minimum ESS ratios, bulk:",
    paste(signif(summaries[, "min_chain_bulk_ratio"], 4), collapse = ", "),
    "tail:", paste(signif(summaries[, "min_chain_tail_ratio"], 4), collapse = ", "), "\n")
for (field in c("bulk_ratio", "tail_ratio", "min_chain_bulk_ratio",
                "min_chain_tail_ratio")) {
  x <- summaries[, field]
  ci <- mean(x) + c(-1, 1) * qt(.975, df = length(x) - 1) * sd(x) / sqrt(length(x))
  cat(field, "mean and 95% t interval:", mean(x),
      paste(signif(ci, 4), collapse = " to "), "\n")
}
