args <- commandArgs(TRUE)
if (length(args) != 1L)
  stop("usage: analyse_step6c_particle_weights_gate.R results-dir")
root <- args[[1L]]
case_seeds <- list(
  EXG = as.integer(strsplit(Sys.getenv("EXG_SEEDS", "101 202 303"), "\\s+")[[1L]]),
  LNR = as.integer(strsplit(Sys.getenv("LNR_SEEDS", "101"), "\\s+")[[1L]]))

ratios <- list()
sample_values <- function(samples) {
  # pm_settings records the active option, so that metadata necessarily differs
  # between the two arms even when all simulated values match.
  attr(samples, "pm_settings") <- NULL
  samples
}
for (case in names(case_seeds)) {
  ratios[[case]] <- numeric()
  for (seed in case_seeds[[case]]) {
    prefix <- file.path(root, paste(case, seed, sep = "_"))
    off <- readRDS(paste0(prefix, "_off.rds"))
    on <- readRDS(paste0(prefix, "_on.rds"))
    stopifnot(identical(off$data_md5, on$data_md5),
              identical(off$initial_seed, on$initial_seed),
              identical(off$final_seed, on$final_seed),
              identical(off$options[["emc2.particle_weights_cpp"]], FALSE),
              identical(on$options[["emc2.particle_weights_cpp"]], TRUE))
    same_samples <- vapply(seq_along(off$fit), function(j)
      identical(sample_values(off$fit[[j]]$samples),
                sample_values(on$fit[[j]]$samples)), logical(1))
    same_rng <- vapply(seq_along(off$fit), function(j)
      identical(off$fit[[j]]$rng, on$fit[[j]]$rng), logical(1))
    stopifnot(all(same_samples), all(same_rng))
    off_s <- unname(off$timing[["elapsed"]])
    on_s <- unname(on$timing[["elapsed"]])
    ratios[[case]] <- c(ratios[[case]], off_s / on_s)
    cat(case, "seed", seed, "data", off$data_md5,
        "elapsed R/C++", off_s, on_s, "speedup", off_s / on_s,
        "draws identical", paste(same_samples, collapse = ","),
        "chain RNG identical", paste(same_rng, collapse = ","),
        "caller RNG identical", identical(off$final_seed, on$final_seed), "\n")
  }
  cat(case, "geometric mean speedup", exp(mean(log(ratios[[case]]))),
      "paired ratios", paste(signif(ratios[[case]], 4), collapse = ", "), "\n")
}
