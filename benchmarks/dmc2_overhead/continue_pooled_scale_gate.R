# Continue a fixed-work pooled run until the DMC2 hardfit R-hat/ESS subset
# passes, checking every 100 sample draws. Likelihood-gap and subject-spread
# rules are not included; see RESULTS.md.
args <- commandArgs(TRUE)
if (length(args) != 3L)
  stop('usage: continue_pooled_scale_gate.R input.rds output.rds max-sample-draws')
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv('EMC_LIB')))
suppressPackageStartupMessages(library(posterior))
input <- args[[1L]]
output <- args[[2L]]
max_sample_draws <- as.integer(args[[3L]])
if (!is.finite(max_sample_draws) || max_sample_draws < 1000L ||
    max_sample_draws %% 100L != 0L || max_sample_draws > 4000L)
  stop('max-sample-draws must be a multiple of 100 from 1000 through 4000')

initial <- readRDS(input)
emc <- initial$fit
if (!is.null(initial$seed)) assign('.Random.seed', initial$seed, envir = .GlobalEnv)
options(emc2.fit_host = FALSE, emc2.sampler_profile = TRUE)

sample_indices <- function(chain) {
  s <- chain$samples
  which(s$stage[seq_len(s$idx)] == 'sample')
}

sample_draw_count <- function(x) min(vapply(x, function(ch)
  length(sample_indices(ch)), integer(1)))

prefix_fit <- function(x, n) {
  out <- lapply(x, function(ch) {
    s <- ch$samples
    sample_idx <- sample_indices(ch)
    if (length(sample_idx) < n) stop('prefix exceeds available sample draws')
    end <- sample_idx[[n]]
    for (nm in names(s)) {
      if (nm == 'idx') next
      z <- s[[nm]]
      d <- dim(z)
      if (!is.null(d) && length(d) && tail(d, 1L) == s$idx) {
        subs <- c(rep(list(TRUE), length(d) - 1L), list(seq_len(end)),
                  list(drop = FALSE))
        s[[nm]] <- do.call(`[`, c(list(z), subs))
      } else if (is.null(d) && length(z) == s$idx) {
        s[[nm]] <- z[seq_len(end)]
      }
    }
    s$idx <- end
    ch$samples <- s
    ch
  })
  class(out) <- class(x)
  out
}

group_draws <- function(x) {
  n <- sample_draw_count(x)
  chains <- lapply(x, function(ch) {
    idx <- sample_indices(ch)
    mu <- t(ch$samples$theta_mu[, idx, drop = FALSE])
    colnames(mu) <- paste('mu', rownames(ch$samples$theta_mu), sep = '__')
    v <- ch$samples$theta_var[, , idx, drop = FALSE]
    vn <- dimnames(ch$samples$theta_var)[[1L]]
    if (is.null(vn)) vn <- rownames(ch$samples$alpha)
    sigma2 <- t(vapply(seq_len(n), function(i) diag(v[, , i]), numeric(length(vn))))
    colnames(sigma2) <- paste('sigma2', vn, sep = '__')
    log_sd <- 0.5 * log(sigma2)
    colnames(log_sd) <- paste('log_sd', vn, sep = '__')
    cbind(mu, log_sd)
  })
  vars <- colnames(chains[[1L]])
  draws <- array(NA_real_, c(n, length(chains), length(vars)),
                 dimnames = list(NULL, NULL, vars))
  for (j in seq_along(chains)) draws[, j, ] <- chains[[j]][, vars, drop = FALSE]
  as_draws_array(draws)
}

hardfit_subset <- function(x) {
  n <- sample_draw_count(x)
  get_rhat <- function(selection) unlist(EMC2::gd_summary(
    x, selection = selection, stage = 'sample', stat = NULL, version = 'old'))
  alpha <- get_rhat('alpha')
  mu <- get_rhat('mu')
  sigma2 <- get_rhat('sigma2')
  ess <- summarise_draws(group_draws(x), 'ess_bulk')
  min_group_log_sd_ess <- min(ess$ess_bulk[grepl('^(mu|log_sd)__', ess$variable)])
  alpha_q99 <- unname(quantile(alpha, .99))
  mu_max <- max(mu)
  sigma2_max <- max(sigma2)
  pass <- n >= 1000L && is.finite(alpha_q99) && is.finite(mu_max) &&
    is.finite(sigma2_max) && is.finite(min_group_log_sd_ess) &&
    alpha_q99 < 1.1 && mu_max < 1.1 && sigma2_max < 1.1 &&
    min_group_log_sd_ess >= 100
  list(sample_draws = n, alpha_rhat_q99 = alpha_q99, mu_rhat_max = mu_max,
       sigma2_rhat_max = sigma2_max, group_log_sd_ess_bulk_min = min_group_log_sd_ess,
       pass = pass)
}

profile_time_at <- function(x, n) {
  b <- x$blocks
  if (is.null(b) || !all(c('stage', 'block_iterations', 'block_total') %in% names(b)))
    stop('input lacks per-block timing records')
  rows <- which(b$stage == 'sample')
  cumulative <- cumsum(b$block_iterations[rows])
  j <- which(cumulative >= n)[1L]
  if (is.na(j)) stop('input profile does not cover the requested sample prefix')
  sum(b$block_total[seq_len(rows[[j]])], na.rm = TRUE)
}

checks <- list()
total_elapsed <- as.numeric(initial$timing[['elapsed']])
target <- NULL
diagnostic_elapsed <- 0
initial_n <- sample_draw_count(emc)
for (n in seq.int(100L, min(initial_n, max_sample_draws), by = 100L)) {
  diagnostic_started <- proc.time()[['elapsed']]
  m <- hardfit_subset(prefix_fit(emc, n))
  diagnostic_elapsed <- diagnostic_elapsed +
    proc.time()[['elapsed']] - diagnostic_started
  m$elapsed_seconds <- profile_time_at(initial, n) + diagnostic_elapsed
  checks[[length(checks) + 1L]] <- m
  if (n >= 1000L && isTRUE(m$pass)) {
    target <- m
    break
  }
}
initial_diagnostic_elapsed <- diagnostic_elapsed

extension_elapsed <- 0
while (is.null(target) && sample_draw_count(emc) < max_sample_draws) {
  remaining <- max_sample_draws - sample_draw_count(emc)
  block_iter <- min(100L, remaining)
  started <- proc.time()[['elapsed']]
  emc <- fit(emc, stage = 'sample', stop_criteria = list(iter = block_iter),
             step_size = block_iter, verbose = FALSE, cores_per_chain = 2,
             cores_for_chains = 2, max_tries = 1, thin = FALSE,
             particle_factor = 50, n_blocks = 1, r_cores = 1)
  block_elapsed <- proc.time()[['elapsed']] - started
  diagnostic_started <- proc.time()[['elapsed']]
  m <- hardfit_subset(emc)
  diagnostic_elapsed <- diagnostic_elapsed +
    proc.time()[['elapsed']] - diagnostic_started
  extension_elapsed <- extension_elapsed + block_elapsed +
    (proc.time()[['elapsed']] - diagnostic_started)
  total_elapsed <- as.numeric(initial$timing[['elapsed']]) +
    initial_diagnostic_elapsed + extension_elapsed
  m$elapsed_seconds <- total_elapsed
  checks[[length(checks) + 1L]] <- m
  if (isTRUE(m$pass)) target <- m
  saveRDS(list(input = input, fit = emc, checks = checks,
               extension_elapsed = extension_elapsed), paste0(output, '.checkpoint.rds'))
}

result <- list(input = input, data_md5 = initial$data_md5,
               target = target, checks = checks,
               reached = !is.null(target), max_sample_draws = max_sample_draws,
               total_elapsed_seconds = if (is.null(target)) total_elapsed else target$elapsed_seconds,
               initial_elapsed_seconds = as.numeric(initial$timing[['elapsed']]),
               extension_elapsed_seconds = extension_elapsed,
               fit = emc)
saveRDS(result, output)
unlink(paste0(output, '.checkpoint.rds'))
cat('sample draws', sample_draw_count(emc), 'reached target', result$reached,
    'elapsed seconds', result$total_elapsed_seconds, '\n')
if (!is.null(target)) {
  print(target)
}
