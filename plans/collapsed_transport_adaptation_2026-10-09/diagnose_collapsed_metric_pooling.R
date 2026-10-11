pkgload::load_all(".", quiet = TRUE, compile = FALSE)
snapshot <- Sys.getenv("COLL_SCREEN_SNAPSHOT",
  "/tmp/lba_sampling_audit/calibrated_snapshot.rds")
fit <- EMC2:::restore_custom_kernel_pointers(
  EMC2:::restore_duplicates(readRDS(snapshot)), quiet = TRUE)
idx <- EMC2:::.emc_chain_proposal_window(fit, stage = "sample")
chart <- EMC2:::.emc_collapsed_context(fit[[1L]])$chart
phi <- lapply(fit, function(ch) t(vapply(idx, function(i)
  EMC2:::.emc_collapsed_to_phi(ch$samples$theta_mu[, i],
    as.matrix(ch$samples$theta_var[, , i]), chart), numeric(chart$dim))))
all_phi <- do.call(rbind, phi)
S_all <- stats::cov(all_phi)
N <- nrow(all_phi)
p <- chart$dim
within_df <- sum(vapply(phi, nrow, integer(1))) - length(phi)
S_within <- Reduce(`+`, lapply(phi, function(a) (nrow(a) - 1) * stats::cov(a))) /
  within_df
grand <- colMeans(all_phi)
S_between <- Reduce(`+`, lapply(phi, function(a)
  nrow(a) * tcrossprod(colMeans(a) - grand))) / (N - 1)
metric_summary <- function(S, n) {
  lambda <- min(1, p / n)
  M <- (1 - lambda) * S + lambda *
    diag(pmax(diag(S), 1e-8), p)
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  c(lambda = lambda, trace = sum(diag(M)), logdet = sum(log(ev)),
    condition = max(ev) / min(ev), min_eigenvalue = min(ev),
    max_eigenvalue = max(ev))
}
metrics <- rbind(
  pooled = metric_summary(S_all, N),
  chain1 = metric_summary(stats::cov(phi[[1L]]), nrow(phi[[1L]])),
  chain2 = metric_summary(stats::cov(phi[[2L]]), nrow(phi[[2L]])),
  chain3 = metric_summary(stats::cov(phi[[3L]]), nrow(phi[[3L]])))
between_trace_share <- sum(diag(S_between)) / sum(diag(S_all))
between_vs_within <- sum(diag(S_between)) / sum(diag(S_within))
out <- data.frame(metric = rownames(metrics), metrics,
  between_trace_share = between_trace_share,
  between_vs_within_trace = between_vs_within,
  window_draws = nrow(phi[[1L]]), chart_dim = p)
write.csv(out,
  "plans/collapsed_transport_adaptation_2026-10-09/collapsed_metric_pooling.csv",
  row.names = FALSE)
print(out, row.names = FALSE)
