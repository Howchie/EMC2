# Posterior and efficiency gates for the Step 2 conditional-fit cache.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: analyse_likelihood_cache.R results-dir")
suppressPackageStartupMessages(library(posterior))
root <- args[[1L]]
probs <- c(.05, .5, .95)
default_seeds <- Sys.getenv("SEEDS", "101 202 303")
cases <- strsplit(Sys.getenv("CASES", "EXG LNR"), "\\s+")[[1L]]
candidate_tag <- Sys.getenv("CANDIDATE_TAG", "default")
case_seeds <- lapply(cases, function(case)
  as.integer(strsplit(Sys.getenv(paste0(case, "_SEEDS"), default_seeds), "\\s+")[[1L]]))
names(case_seeds) <- cases

chain_draws <- function(fit) {
  lapply(fit, function(ch) {
    alpha <- ch$samples$alpha
    draws <- t(alpha[, 1L, , drop = FALSE][, 1L, ])
    colnames(draws) <- dimnames(alpha)[[1L]]
    draws
  })
}

summarise_chains <- function(chains) {
  vars <- colnames(chains[[1L]])
  stopifnot(all(vapply(chains, function(x) identical(colnames(x), vars), logical(1))))
  n <- min(vapply(chains, nrow, integer(1)))
  draws <- array(NA_real_, c(n, length(chains), length(vars)),
                 dimnames = list(NULL, NULL, vars))
  for (j in seq_along(chains)) draws[, j, ] <- chains[[j]][seq_len(n), vars, drop = FALSE]
  stats <- summarise_draws(as_draws_array(draws), "mean", "ess_bulk", "ess_tail", "rhat",
    ~quantile2(.x, probs = probs),
    ~setNames(vapply(probs, function(p) mcse_quantile(.x, probs = p), numeric(1)),
              c("mcse_q05", "mcse_q50", "mcse_q95")))
  list(draws = draws, stats = stats, n = n)
}

gate <- function(a, b) {
  stopifnot(identical(a$stats$variable, b$stats$variable))
  z <- do.call(rbind, lapply(seq_along(probs), function(i) {
    q <- c("q5", "q50", "q95")[[i]]
    mcse <- c("mcse_q05", "mcse_q50", "mcse_q95")[[i]]
    data.frame(variable = a$stats$variable, quantile = q,
      standardized_difference = abs(a$stats[[q]] - b$stats[[q]]) /
        sqrt(a$stats[[mcse]]^2 + b$stats[[mcse]]^2))
  }))
  z
}

ratios <- setNames(lapply(cases, function(x)
  list(bulk = numeric(), tail = numeric())), cases)
all_base_chains <- all_new_chains <- setNames(vector("list", length(cases)), cases)
all_base_times <- all_new_times <- setNames(vector("list", length(cases)), cases)
for (case in cases) for (seed in case_seeds[[case]]) {
  prefix <- file.path(root, paste(case, seed, sep = "_"))
  base <- readRDS(paste0(prefix, "_base_default.rds"))
  new <- readRDS(paste0(prefix, "_new_", candidate_tag, ".rds"))
  stopifnot(identical(base$data_md5, new$data_md5),
            identical(base$initial_seed, new$initial_seed))
  base_chains <- chain_draws(base$fit)
  new_chains <- chain_draws(new$fit)
  sb <- summarise_chains(base_chains)
  sn <- summarise_chains(new_chains)
  diffs <- gate(sb, sn)
  ratios[[case]]$bulk <- c(ratios[[case]]$bulk,
                            min(sn$stats$ess_bulk) / min(sb$stats$ess_bulk))
  ratios[[case]]$tail <- c(ratios[[case]]$tail,
                            min(sn$stats$ess_tail) / min(sb$stats$ess_tail))
  all_base_chains[[case]] <- c(all_base_chains[[case]], base_chains)
  all_new_chains[[case]] <- c(all_new_chains[[case]], new_chains)
  all_base_times[[case]] <- c(all_base_times[[case]], base$timing[["elapsed"]])
  all_new_times[[case]] <- c(all_new_times[[case]], new$timing[["elapsed"]])
  cat("\n", case, seed,
      "draws/chain base,candidate", sb$n, sn$n,
      "fit elapsed base,candidate", base$timing[["elapsed"]], new$timing[["elapsed"]],
      "caller RNG identical", identical(base$final_seed, new$final_seed),
      "quantile failures", sum(diffs$standardized_difference > 4L), "/", nrow(diffs),
      "largest combined-MCSE difference", max(diffs$standardized_difference), "\n")
  cat("minimum ESS/draw base,candidate bulk",
      min(sb$stats$ess_bulk) / (sb$n * dim(sb$draws)[2L]),
      min(sn$stats$ess_bulk) / (sn$n * dim(sn$draws)[2L]), "tail",
      min(sb$stats$ess_tail) / (sb$n * dim(sb$draws)[2L]),
      min(sn$stats$ess_tail) / (sn$n * dim(sn$draws)[2L]), "\n")
  cat("quantile standardized differences:\n"); print(diffs, row.names = FALSE)
}

for (case in cases) {
  cat("\n", case, "mean elapsed base/candidate",
      mean(all_base_times[[case]]), mean(all_new_times[[case]]),
      "mean candidate/base ratio",
      mean(all_new_times[[case]] / all_base_times[[case]]), "\n")
  for (kind in names(ratios[[case]])) {
    log_ratio <- log(ratios[[case]][[kind]])
    ci <- exp(mean(log_ratio) + c(-1, 1) *
      qt(.975, df = length(log_ratio) - 1L) * sd(log_ratio) / sqrt(length(log_ratio)))
    cat(case, kind, "minimum ESS geometric-mean ratio", exp(mean(log_ratio)),
        "95% t interval", ci[[1L]], ci[[2L]],
        "seed ratios", ratios[[case]][[kind]], "\n")
  }
  base_all <- summarise_chains(all_base_chains[[case]])
  new_all <- summarise_chains(all_new_chains[[case]])
  aggregate_diffs <- gate(base_all, new_all)
  cat("\n", case, "combined", length(all_base_chains[[case]]),
      "chains; common draws per chain", base_all$n, new_all$n, "\n",
      "Aggregate quantile failures:", sum(aggregate_diffs$standardized_difference > 4L),
      "/", nrow(aggregate_diffs), "largest standardized difference:",
      max(aggregate_diffs$standardized_difference), "\n")
  cat("Aggregate minimum ESS/draw base,candidate bulk",
      min(base_all$stats$ess_bulk) / (base_all$n * length(all_base_chains[[case]])),
      min(new_all$stats$ess_bulk) / (new_all$n * length(all_new_chains[[case]])), "tail",
      min(base_all$stats$ess_tail) / (base_all$n * length(all_base_chains[[case]])),
      min(new_all$stats$ess_tail) / (new_all$n * length(all_new_chains[[case]])), "\n")
  cat("Aggregate per-parameter ESS/draw and Rhat:\n")
  print(data.frame(variable = base_all$stats$variable,
    bulk_base = base_all$stats$ess_bulk / (base_all$n * length(all_base_chains[[case]])),
    bulk_candidate = new_all$stats$ess_bulk / (new_all$n * length(all_new_chains[[case]])),
    tail_base = base_all$stats$ess_tail / (base_all$n * length(all_base_chains[[case]])),
    tail_candidate = new_all$stats$ess_tail / (new_all$n * length(all_new_chains[[case]])),
    rhat_base = base_all$stats$rhat, rhat_candidate = new_all$stats$rhat))
  cat("Aggregate quantile standardized differences:\n"); print(aggregate_diffs, row.names = FALSE)
}
