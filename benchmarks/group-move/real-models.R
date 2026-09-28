# Compare the original sampler, retained group-move methods, in short real
# RDM and LBA fits. This measures integration and sampling efficiency, not
# posterior accuracy.
suppressPackageStartupMessages(pkgload::load_all(quiet = TRUE))
if (!requireNamespace("posterior", quietly = TRUE))
  stop("Install posterior for benchmark diagnostics")

data("forstmann", package = "EMC2")
subjects <- unique(forstmann$subjects)[1:2]
dat <- droplevels(do.call(rbind, lapply(subjects, function(id) {
  head(forstmann[forstmann$subjects == id, ], 50L)
})))
matchfun <- function(d) d$S == d$lR
AD <- matrix(c(-0.5, 0.5), ncol = 1L, dimnames = list(NULL, "d"))

cases <- list(
  RDM = list(model = RDM, formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
             constants = c(s = 0)),
  LBA = list(model = LBA, formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1),
             constants = c(sv = 0))
)
methods <- c("original", "legacy", "am")
warmup <- 250L
settle <- 100L
preburn <- 50L
burn <- 50L
draws_n <- 1000L
seed <- 260926L
rows <- list()

for (model_name in names(cases)) for (method in methods) {
  cfg <- cases[[model_name]]
  des <- design(data = dat, model = cfg$model, matchfun = matchfun,
    formula = cfg$formula, contrasts = list(v = list(lM = AD)),
    constants = cfg$constants)
  proposal_count <- if (method == "legacy") 4L else 1L
  configured_method <- if (method == "original") "am" else method
  options(emc2.group_move_method = configured_method,
          emc2.group_move = method != "original",
          emc2.group_move_proposals = proposal_count, emc2.group_move_scale = TRUE,
          emc2.group_move_warmup = warmup,
          emc2.group_move_settle = settle)
  set.seed(seed)
  emc <- suppressMessages(make_emc(dat, des, type = "standard", n_chains = 2L,
                                   compress = TRUE, rt_resolution = 0.05))
  started <- proc.time()[["elapsed"]]
  emc <- fit(emc, cores_for_chains = 1L, cores_per_chain = 1L,
    particle_factor = 10L, step_size = 50L, max_tries = 20L,
    verbose = FALSE, trim = FALSE,
    stop_criteria = list(preburn = list(iter = preburn),
      burn = list(iter = burn),
      adapt = list(iter = warmup + settle, min_unique = 0L),
      sample = list(iter = draws_n)))
  elapsed <- proc.time()[["elapsed"]] - started

  n_mu <- nrow(emc[[1]]$samples$theta_mu)
  n_var <- dim(emc[[1]]$samples$theta_var)[1L]
  z <- array(NA_real_, c(draws_n, length(emc), n_mu + n_var))
  for (ch in seq_along(emc)) {
    s <- emc[[ch]]$samples
    idx <- which(s$stage == "sample")
    if (length(idx) != draws_n)
      stop(model_name, "/", method, " retained ", length(idx),
           " production draws; expected ", draws_n)
    z[, ch, seq_len(n_mu)] <- t(s$theta_mu[, idx, drop = FALSE])
    z[, ch, n_mu + seq_len(n_var)] <- t(vapply(idx, function(i) {
      log(sqrt(diag(s$theta_var[, , i])))
    }, numeric(n_var)))
  }
  rhat <- posterior::rhat(z)
  ess <- posterior::ess_bulk(z)
  gm <- group_move_diagnostics(emc)
  rows[[length(rows) + 1L]] <- data.frame(
    model = model_name, method = method, seed = seed,
    max_rhat = max(rhat, na.rm = TRUE), min_bulk_ess = min(ess, na.rm = TRUE),
    ess_per_second = min(ess, na.rm = TRUE) / elapsed,
    elapsed_seconds = elapsed, group_dimension = gm$dimension[1L],
    mean_adaptation = mean(gm$n_adapt, na.rm = TRUE),
    group_move_likelihood_evaluations =
      if (method == "original") 0L else
        (preburn + burn + warmup + settle + draws_n) * length(emc) *
          length(subjects) * proposal_count,
    mean_sample_acceptance = mean(gm$sample_acceptance, na.rm = TRUE),
    stringsAsFactors = FALSE)
  cat(model_name, method, "max Rhat", round(max(rhat, na.rm = TRUE), 3),
      "min bulk ESS", round(min(ess, na.rm = TRUE)),
      "ESS/s", round(min(ess, na.rm = TRUE) / elapsed, 2),
      "seconds", round(elapsed, 1), "\n")
}

result <- do.call(rbind, rows)
outdir <- "benchmarks/group-move/results/real-models"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(outdir, "metrics.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.txt"))
print(result, row.names = FALSE)
