pkgload::load_all(".", quiet = TRUE, compile = FALSE)

input <- "/tmp/collapsed_transport_followup/fresh_retuned_calibrated.rds"
fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(input)), quiet = TRUE)
chain <- fit[[1L]]
old <- chain$collapsed$plan
old$calibrate <- TRUE
old$calibration_screen <- FALSE
old$independent_scales <- FALSE
old$scale_updates <- 0L
old$t0_updates <- 0L
old$t0_hyper_updates <- 0L
old$length_updates <- 0L
chain$collapsed$plan <- old

rank_and_null_slope <- function(fit) {
  if (is.null(fit)) return(c(rank = NA_real_, null_slope = NA_real_))
  H <- (fit$prec + t(fit$prec)) / 2
  e <- eigen(H, symmetric = TRUE)
  tol <- max(abs(e$values), .Machine$double.xmin) * 1e-10
  keep <- which(e$values > tol)
  flat <- setdiff(seq_along(e$values), keep)
  null_slope <- if (length(flat))
    sqrt(sum(crossprod(e$vectors[, flat, drop = FALSE], fit$lin)^2)) else 0
  c(rank = length(keep), null_slope = null_slope)
}

raw <- lapply(chain$collapsed$fits, rank_and_null_slope)
ctx <- EMC2:::.emc_collapsed_context(chain)
safe <- lapply(ctx$fits, rank_and_null_slope)
raw <- do.call(rbind, raw)
safe <- do.call(rbind, safe)

retained_before <- chain$samples
started <- proc.time()[["elapsed"]]
chain <- EMC2:::.emc_collapsed_calibrate(chain, n_cores = 4L)
elapsed <- proc.time()[["elapsed"]] - started
stopifnot(identical(chain$samples, retained_before))
plan <- chain$collapsed$plan

fit_summary <- data.frame(
  stage = rep(c("stored", "sanitized"), each = nrow(raw)),
  subject = rep(seq_len(nrow(raw)), 2L),
  rank = c(raw[, "rank"], safe[, "rank"]),
  null_slope = c(raw[, "null_slope"], safe[, "null_slope"])
)
write.csv(fit_summary,
  "plans/collapsed_transport_adaptation_2026-10-09/safe_fit_rank_audit.csv",
  row.names = FALSE)

summary <- data.frame(
  start_scale = old$scale,
  start_t0_mean_scale = old$scale * old$t0_ratio,
  start_t0_hyper_scale = old$scale * old$t0_hyper_ratio,
  selected_non_t0_scale = plan$scale,
  selected_t0_mean_scale = plan$scale * plan$t0_ratio,
  selected_t0_hyper_scale = plan$scale * plan$t0_hyper_ratio,
  path_length = plan$length,
  eps = plan$eps,
  validation_acceptance = plan$calibration$acceptance,
  marginal_cost = plan$calibration$cost,
  setup_cost = plan$calibration$fixed_cost,
  seconds = elapsed,
  pilots = plan$calibration$pilots
)
write.csv(summary,
  "plans/collapsed_transport_adaptation_2026-10-09/safe_fit_calibration_summary.csv",
  row.names = FALSE)
saveRDS(chain,
  "/tmp/collapsed_transport_followup/safe_fit_independent_calibrated_chain1.rds")
print(summary, row.names = FALSE)
cat("Stored fit ranks:", paste(sort(unique(raw[, "rank"])), collapse = ","),
    "subjects with nonzero null slope:", sum(raw[, "null_slope"] > 1e-8),
    "of", nrow(raw), "\n")
cat("Sanitized fit ranks:", paste(sort(unique(safe[, "rank"])), collapse = ","),
    "subjects with nonzero null slope:", sum(safe[, "null_slope"] > 1e-8),
    "of", nrow(safe), "\n")
