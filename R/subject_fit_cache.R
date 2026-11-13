# Durable, chain-local cache for conditional subject fits and their
# importance diagnostics. The configuration fingerprint is built at run
# setup, not in the likelihood evaluation path.
.emc_likelihood_config_fingerprint <- function(sampler) {
  strip_pooled <- function(x) {
    attr(x, "pooled_values") <- NULL
    if (is.list(x) && !is.data.frame(x))
      for (i in seq_along(x)) x[[i]] <- strip_pooled(x[[i]])
    x
  }
  data <- strip_pooled(sampler$data)
  path <- tempfile("emc-likelihood-config-")
  con <- file(path, open = "wb")
  on.exit({
    if (!is.null(con)) close(con)
    unlink(path)
  }, add = TRUE)
  ok <- tryCatch({
    serialize(list(data = data, model = sampler$model,
                   par_names = sampler$par_names, n_pars = sampler$n_pars,
                   type = sampler$type), con, version = 3L)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  close(con)
  con <- NULL
  unname(tools::md5sum(path))
}

.emc_prepare_subject_fit_cache <- function(emc, validate = FALSE) {
  for (j in seq_along(emc)) {
    chain <- emc[[j]]
    old <- chain$likelihood_cache_version
    config <- if (isTRUE(validate) || !is.list(old) || is.null(old$config))
      .emc_likelihood_config_fingerprint(chain) else old$config
    if (is.null(config)) {
      chain$likelihood_cache_version <- NULL
      chain$subject_fit_cache <- NULL
      chain$laplace_check_cache <- NULL
      emc[[j]] <- chain
      next
    }
    version <- list(schema = 1L, config = config)
    config_changed <- is.list(old) && !is.null(old$config) &&
      !identical(old$config, config)
    if (!identical(old, version)) {
      chain$subject_fit_cache <- vector("list", chain$n_subjects)
      chain$laplace_check_cache <- vector("list", chain$n_subjects)
      chain$likelihood_cache_version <- version
    } else {
      if (!is.list(chain$subject_fit_cache) ||
          length(chain$subject_fit_cache) != chain$n_subjects)
        chain$subject_fit_cache <- vector("list", chain$n_subjects)
      if (!is.list(chain$laplace_check_cache) ||
          length(chain$laplace_check_cache) != chain$n_subjects)
        chain$laplace_check_cache <- vector("list", chain$n_subjects)
    }
    if (config_changed) {
      chain$lik_laplace <- vector("list", chain$n_subjects)
      chain$burn_mode_mu <- NULL
      chain$burn_mode_var <- NULL
      chain$burn_mode_stats <- NULL
    }
    emc[[j]] <- chain
  }
  emc
}

# REXG is the only registered native model enabled by default. LNR failed the
# tail-ESS gate. The 18-seed multimodal custom fixture passed, but custom R
# likelihoods can capture mutable state that the cache key cannot describe, so
# keep them opt-in through the explicit override while more cache contracts are
# evaluated.
.emc_subject_fit_cache_max_reuses <- function(sampler) {
  override <- getOption("emc2.subject_fit_cache_max_reuses", NULL)
  if (is.numeric(override) && length(override) == 1L &&
      !is.na(override) && override >= 0) return(override)
  model_name <- tryCatch({
    if (is.function(sampler$model)) sampler$model()$c_name else NULL
  }, error = function(e) NULL)
  if (identical(model_name, "REXG")) Inf else 0L
}

# Refit the exact cached conditioning state, using its previous mode only as
# the optimizer start. The Laplace diagnostic still runs with fresh draws.
.emc_subject_fit_cache_warm_start <- function(sampler) {
  isTRUE(getOption("emc2.subject_fit_cache_warm_start", FALSE))
}

.emc_subject_fit_key <- function(chain, subject, prior, proposal_prior_scale) {
  version <- chain$likelihood_cache_version
  if (!is.list(version) || is.null(version$config)) return(NULL)
  list(config = version$config, subject = as.integer(subject),
       prior_mu = as.numeric(prior$mu),
       prior_var = unname(as.matrix(prior$var)),
       pooled = chain$pooled$value,
       proposal_prior_scale = as.numeric(proposal_prior_scale))
}

.emc_subject_fit_refresh_due <- function(chain, subject, key, max_reuses = Inf) {
  if (!is.numeric(max_reuses) || length(max_reuses) != 1L ||
      is.na(max_reuses) || max_reuses < 0 || is.infinite(max_reuses) ||
      is.null(key) || length(chain$subject_fit_cache) < subject) return(FALSE)
  entry <- chain$subject_fit_cache[[subject]]
  if (!is.list(entry) || !identical(entry$key, key) ||
      !.emc_usable_mode(entry$result)) return(FALSE)
  reuse_count <- entry$reuse_count
  if (!is.numeric(reuse_count) || length(reuse_count) != 1L ||
      is.na(reuse_count) || reuse_count < 0) reuse_count <- 0L
  reuse_count >= max_reuses
}

.emc_subject_fit_lookup <- function(chain, subject, key, force_refresh = FALSE,
                                    warm_start_on_refresh = FALSE) {
  if (is.null(key) || length(chain$subject_fit_cache) < subject)
    return(list(hit = FALSE, result = NULL, warm_start = NULL,
                reuse_count = 0L))
  entry <- chain$subject_fit_cache[[subject]]
  if (!is.list(entry) || !is.list(entry$key) || !is.list(entry$result))
    return(list(hit = FALSE, result = NULL, warm_start = NULL,
                reuse_count = 0L))
  reuse_count <- entry$reuse_count
  if (!is.numeric(reuse_count) || length(reuse_count) != 1L ||
      is.na(reuse_count) || reuse_count < 0) reuse_count <- 0L
  exact <- identical(entry$key, key) && .emc_usable_mode(entry$result)
  if (exact && !isTRUE(force_refresh)) {
    result <- entry$result
    # Gain measures the improvement from the original optimizer start, so it
    # is not a current-run statistic when the fit itself is reused.
    result$gain <- NA_real_
    return(list(hit = TRUE, result = result,
                warm_start = entry$result$alpha,
                reuse_count = reuse_count))
  }
  warm_start <- NULL
  if (isTRUE(force_refresh) && isTRUE(warm_start_on_refresh) && exact) {
    warm_start <- entry$result$alpha
  } else if (!isTRUE(force_refresh) &&
             identical(entry$key$config, key$config) &&
             .emc_usable_mode(entry$result) &&
             length(entry$result$alpha) == length(key$prior_mu)) {
    warm_start <- entry$result$alpha
  }
  list(hit = FALSE, result = NULL, warm_start = warm_start,
       reuse_count = if (exact) reuse_count else 0L)
}

.emc_subject_fit_store <- function(chain, subject, key, result,
                                   reuse_count = 0L) {
  if (is.null(key) || !.emc_usable_mode(result) ||
      length(chain$subject_fit_cache) < subject) return(chain)
  chain$subject_fit_cache[[subject]] <- list(
    key = key, result = result, reuse_count = as.integer(reuse_count))
  chain
}

.emc_laplace_check_key <- function(chain, subject, mode, covariance, prior) {
  version <- chain$likelihood_cache_version
  if (!is.list(version) || is.null(version$config)) return(NULL)
  .emc_laplace_check_key_values(
    version$config, subject, mode, covariance, prior,
    if (is.null(chain$pooled)) NULL else chain$pooled$value)
}

.emc_laplace_check_key_values <- function(config, subject, mode, covariance,
                                           prior, pooled = NULL) {
  list(config = config, subject = as.integer(subject),
       mode = as.numeric(mode), covariance = unname(as.matrix(covariance)),
       prior_mu = as.numeric(prior$mu),
       prior_var = unname(as.matrix(prior$var)),
       pooled = pooled)
}

.emc_laplace_check_lookup <- function(chain, subject, key) {
  if (is.null(key) || length(chain$laplace_check_cache) < subject)
    return(NULL)
  entries <- chain$laplace_check_cache[[subject]]
  if (!is.list(entries) || !length(entries)) return(NULL)
  for (i in rev(seq_along(entries)))
    if (identical(entries[[i]]$key, key)) return(entries[[i]]$check)
  NULL
}

.emc_laplace_check_store <- function(chain, subject, key, check,
                                     max_entries = 8L, draws = NULL) {
  if (is.null(key) || is.null(check) || identical(check$status, "error") ||
      length(chain$laplace_check_cache) < subject) return(chain)
  entries <- chain$laplace_check_cache[[subject]]
  if (!is.list(entries)) entries <- list()
  same <- Filter(function(entry) identical(entry$key, key), entries)
  entries <- Filter(function(entry) !identical(entry$key, key), entries)
  # Only the latest entry keeps importance draws for pooling; an exact repeat
  # keeps the draws it already had.
  if (is.null(draws) && length(same)) draws <- same[[length(same)]]$draws
  entries <- lapply(entries, function(entry) { entry$draws <- NULL; entry })
  entry <- list(key = key, check = check)
  if (!is.null(draws)) entry$draws <- draws
  entries[[length(entries) + 1L]] <- entry
  if (length(entries) > max_entries)
    entries <- entries[seq.int(length(entries) - max_entries + 1L,
                               length(entries))]
  chain$laplace_check_cache[[subject]] <- entries
  chain
}
