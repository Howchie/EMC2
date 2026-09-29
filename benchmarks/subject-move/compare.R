# Rscript compare.R <candidate-root> <reference-root> [reference-fallback-root]
args <- commandArgs(TRUE)
summarize <- function(r) {
  group <- Filter(Negate(is.null), r$draws)
  group <- lapply(group, function(x) x[, , apply(x, 3, sd) > 0, drop = FALSE])
  rh <- unlist(lapply(group, function(x) apply(x, 3, posterior::rhat)))
  ess <- unlist(lapply(group, function(x) apply(x, 3, posterior::ess_bulk)))
  c(warmup_seconds = sum(r$timing[c("preburn", "burn", "adapt")]),
    sample_seconds = r$timing[["sample"]], max_group_rhat = max(rh),
    max_subject_rhat = if (is.null(r$alpha_rhat)) NA_real_ else max(r$alpha_rhat),
    min_subject_ess = min(r$alpha_ess), median_subject_ess = median(r$alpha_ess),
    min_group_ess = min(ess),
    min_subject_ess_sec = min(r$alpha_ess) / r$timing[["sample"]],
    median_subject_ess_sec = median(r$alpha_ess) / r$timing[["sample"]],
    min_group_ess_sec = min(ess) / r$timing[["sample"]])
}
paths <- list.files(args[1], pattern = "-seed[0-9]+[.]rds$", full.names = TRUE)
rows <- lapply(paths, function(path) {
  r <- readRDS(path); candidate <- summarize(r)
  refs <- file.path(args[-1], basename(path)); refs <- refs[file.exists(refs)]
  if (!length(refs)) return(NULL)
  reference <- summarize(readRDS(refs[1]))
  rate <- c("min_subject_ess_sec", "median_subject_ess_sec", "min_group_ess_sec")
  out <- data.frame(case = r$case, arm = r$arm, seed = r$seed, t(candidate))
  for (nm in rate) out[[paste0(nm, "_ratio")]] <- candidate[[nm]] / reference[[nm]]
  out$reference <- refs[1]
  out
})
out <- do.call(rbind, rows)
print(out, row.names = FALSE, digits = 4)
if (!is.null(out)) write.csv(out, file.path(args[1], "comparison.csv"), row.names = FALSE)
