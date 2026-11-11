args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: analyse_step6b_context_fit_gate.R results-dir")
root <- args[[1L]]
seed <- as.integer(Sys.getenv("SEED", "101"))

for (case in c("EXG", "LNR")) {
  prefix <- file.path(root, paste(case, seed, sep = "_"))
  off <- readRDS(paste0(prefix, "_off.rds"))
  on <- readRDS(paste0(prefix, "_on.rds"))
  stopifnot(identical(off$data_md5, on$data_md5),
            identical(off$initial_seed, on$initial_seed),
            identical(off$final_seed, on$final_seed))
  stopifnot(identical(off$options[["emc2.native_pt_context"]], FALSE),
            identical(on$options[["emc2.native_pt_context"]], TRUE))
  same_samples <- vapply(seq_along(off$fit), function(j)
    identical(off$fit[[j]]$samples, on$fit[[j]]$samples), logical(1))
  same_rng <- vapply(seq_along(off$fit), function(j)
    identical(off$fit[[j]]$rng, on$fit[[j]]$rng), logical(1))
  state_fields <- setdiff(union(names(off$fit[[1L]]), names(on$fit[[1L]])),
                          c("data", "prior", "model"))
  same_state <- vapply(state_fields, function(field)
    all(vapply(seq_along(off$fit), function(j)
      identical(off$fit[[j]][[field]], on$fit[[j]][[field]]), logical(1))),
    logical(1))
  stopifnot(all(same_samples), all(same_rng))
  off_s <- unname(off$timing[["elapsed"]])
  on_s <- unname(on$timing[["elapsed"]])
  cat(case, "seed", seed, "data", off$data_md5,
      "elapsed context-off/on", off_s, on_s,
      "speedup", off_s / on_s,
      "draws identical per chain", paste(same_samples, collapse = ","),
      "chain RNG identical", paste(same_rng, collapse = ","),
      "caller RNG identical", identical(off$final_seed, on$final_seed),
      "differing continuation fields", paste(state_fields[!same_state], collapse = ","),
      "\n")
}
