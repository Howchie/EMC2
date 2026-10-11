# Summary table of fit_arm.R results.
#   Rscript benchmarks/subject-hmc/summarise.R <dir> [out.csv]
args <- commandArgs(TRUE)
dir <- args[1]
`%||%` <- function(a, b) if (is.null(a)) b else a
files <- list.files(dir, "^[A-Za-z0-9]+_(base|shear|hmc|fullc?)_[0-9]+\\.rds$", full.names = TRUE)
rows <- lapply(files, function(f) {
  r <- readRDS(f)
  log <- readLines(sub("\\.rds$", ".log", f))
  u <- grep("^usage_usec", log, value = TRUE)
  cpu <- if (length(u)) as.numeric(sub("usage_usec ", "", u[1])) / 1e6 else NA_real_
  prod_cpu <- if (!is.null(r$stage_cpu$sample)) r$stage_cpu$sample[["cpu"]] else NA_real_
  cls <- names(r$diag)
  min_bulk <- min(vapply(r$diag, `[[`, 0, "min_ess_bulk"))
  min_tail <- min(vapply(r$diag, `[[`, 0, "min_ess_tail"))
  h <- Filter(Negate(is.null), r$hmc)
  tot <- function(k) if (length(h)) sum(vapply(h, function(x) as.numeric(x$totals[[k]] %||% 0), 0)) else NA_real_
  collapsed <- Filter(Negate(is.null), r$moves)
  collapsed <- Filter(function(x) !is.null(x$collapsed), collapsed)
  coll_tot <- function(k) if (length(collapsed))
    sum(vapply(collapsed, function(x) as.numeric(x$collapsed[[k]] %||% 0), 0)) else 0
  coll_attempts <- coll_tot("attempts")
  coll_accepted <- coll_tot("accepted_rounds")
  coll_expected <- coll_tot("acceptance_probability_sum")
  data.frame(
    case = r$case, arm = r$arm, seed = r$seed, stop_ok = r$stop_ok, try = r$sample_try,
    burn = r$stage_iters[["burn"]], adapt = r$stage_iters[["adapt"]], sample = r$stage_iters[["sample"]],
    kept = r$n_draws, cpu = round(cpu), prod_cpu = round(prod_cpu),
    rhat_max = round(max(vapply(r$diag, `[[`, 0, "max_rhat")), 3),
    n_rhat_1.01 = sum(vapply(r$diag, `[[`, 0L, "n_over_1.01")),
    rhat_mu = round(r$diag$mu$max_rhat, 3), rhat_alpha = round(r$diag$alpha$max_rhat, 3),
    ess_bulk_min = round(min_bulk), ess_tail_min = round(min_tail),
    ess_mu = round(r$diag$mu$min_ess_bulk), ess_alpha = round(r$diag$alpha$min_ess_bulk),
    ess_per_kcpu = round(1000 * min(min_bulk, min_tail) / cpu, 2),
    ess_per_kcpu_prod = round(1000 * min(min_bulk, min_tail) / prod_cpu, 2),
    hmc_attempts = tot("attempts"), hmc_accept = round(tot("accepted") / max(1, tot("attempts")), 2),
    collapsed_attempts = coll_attempts, collapsed_accepted_rounds = coll_accepted,
    collapsed_accept = if (coll_attempts > 0) round(coll_accepted / coll_attempts, 3) else NA_real_,
    collapsed_expected_accept = if (coll_attempts > 0) round(coll_expected / coll_attempts, 3) else NA_real_,
    hmc_row_share = round(100 * (tot("rows") + tot("calib_rows")) /
                            (tot("rows") + tot("calib_rows") + tot("base_rows")), 1),
    calib_share = round(100 * tot("calib_rows") / (tot("rows") + tot("calib_rows") + tot("base_rows")), 1),
    max_z_mu = if (is.null(r$recovery)) NA else round(max(abs(r$recovery$z_mu)), 2),
    max_z_lsd = if (is.null(r$recovery)) NA else round(max(abs(r$recovery$z_lsd)), 2))
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$case, tab$seed, tab$arm), ]
print(tab, row.names = FALSE)
if (length(args) > 1) utils::write.csv(tab, args[2], row.names = FALSE)
