# Package-level environment for storing suppressed loo warnings
.loo_warn_store <- new.env(parent = emptyenv())
.loo_warn_store$msgs <- character(0)

# Helper: run thunk, capture any loo/pareto-k warnings, stash for later
.loo_call <- function(thunk) {
  withCallingHandlers(thunk(), warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("pareto|p_waic|waic|loo", msg, ignore.case = TRUE)) {
      .loo_warn_store$msgs <- c(.loo_warn_store$msgs, msg)
      invokeRestart("muffleWarning")
    }
  })
}

#' Retrieve suppressed WAIC/loo warnings
#'
#' WAIC computation suppresses Pareto-k and p_waic diagnostic warnings
#' to keep output clean. Call this function after \code{compare()} or
#' \code{compare_subject()} to inspect any suppressed messages.
#'
#' @return Character vector of warning messages (invisibly), printed to console.
#' @export
waic_warnings <- function() {
  msgs <- .loo_warn_store$msgs
  if (length(msgs) == 0) {
    message("No WAIC warnings stored.")
  } else {
    message("Stored WAIC warnings:")
    for (m in msgs) message("  ", m)
  }
  invisible(msgs)
}

# Compute the pointwise log-likelihood matrix for one subject.
.ll_matrix_subject <- function(emc, stage = "sample", filter = 0, subject) {
  # Try to get cached first
  pw_ll <- get_pars(emc, selection = "pw_ll", stage = stage, filter = filter,
                    merge_chains = TRUE, return_mcmc = FALSE, subject = subject)
  if (!is.null(pw_ll[[1]])) {
    return(t(pw_ll[[1]]))
  }
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  proposals <- do.call(rbind, alpha[[subject]])   # [n_iter x n_pars]
  calc_ll_pw(proposals, emc[[1]]$data[[subject]], emc[[1]]$model)
}

# Concatenate per-trial pointwise log-likelihood matrices across subjects.
# Returns [n_iter x total_trials].
.ll_matrix_pooled <- function(emc, stage = "sample", filter = 0, cores = 1) {
  # Try to get cached first
  pw_ll <- get_pars(emc, selection = "pw_ll", stage = stage, filter = filter,
                    merge_chains = TRUE, return_mcmc = FALSE)
  if (!is.null(pw_ll)) {
    if (is.list(pw_ll)) pw_ll <- pw_ll[[1]]
    return(t(pw_ll))   # [n_iter x total_trials]
  }
  has_any <- any(sapply(emc, function(x) is.list(x) && !is.null(x$samples$pw_ll)))
  if (has_any)
    message("pw_ll cache incomplete (not all chains have it) — recomputing. ",
            "Check that add_pw_ll() or save_pw_ll=TRUE completed without errors.")
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  data  <- emc[[1]]$data
  model <- emc[[1]]$model
  ll_list <- auto_mclapply(names(alpha), function(sub) {
    proposals <- do.call(rbind, alpha[[sub]])     # [n_iter x n_pars]
    calc_ll_pw(proposals, data[[sub]], model)     # [n_iter x n_trials_sub]
  }, mc.cores = cores)
  do.call(cbind, ll_list)                         # [n_iter x total_trials]
}

# Marginal subject log-likelihood matrix for hierarchical leave-one-subject-out LOO.
#
# The conditional approach (summing per-trial LLs within a subject) gives
# importance weights with Pareto-k >> 1 because removing hundreds of trials at
# once is too influential for PSIS. The correct quantity for hierarchical LOO is
# the MARGINAL subject likelihood p(y_j | psi) = E_{alpha_j ~ p(alpha|psi)}[p(y_j | alpha_j)],
# which integrates out subject parameters and varies only with the group-level
# parameters (theta_mu, theta_var). Its lower iteration-to-iteration variance
# gives well-conditioned PSIS importance weights.
#
# For each posterior draw of (theta_mu, theta_var), K proposals are drawn from
# MVN(theta_mu, theta_var) and the marginal is estimated via log-mean-exp.
# Falls back to conditional subject sums for single-level (no group) models.
#
# Returns [n_iter x n_subjects].
.marg_ll_matrix <- function(emc, stage = "sample", filter = 0, K = 200) {
  theta_mu  <- get_pars(emc, selection = "mu",    stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)

  if (is.null(theta_mu) || is.null(theta_var)) {
    # Single-level model: marginalisation isn't possible; sum conditional LLs.
    ll_mat    <- .ll_matrix_pooled(emc, stage = stage, filter = filter)
    subjects  <- names(emc[[1]]$data)
    trial_counts <- sapply(emc[[1]]$data[subjects], nrow)
    cum_trials   <- c(0, cumsum(trial_counts))
    out <- matrix(nrow = nrow(ll_mat), ncol = length(subjects))
    colnames(out) <- subjects
    for (i in seq_along(subjects))
      out[, i] <- rowSums(ll_mat[, (cum_trials[i]+1):cum_trials[i+1], drop = FALSE])
    return(out)
  }

  # theta_mu:  [n_pars x n_iter]
  # theta_var: [n_pars x n_pars x n_iter]
  n_iter   <- ncol(theta_mu)
  subjects <- names(emc[[1]]$data)
  n_subj   <- length(subjects)
  data_list <- emc[[1]]$data
  model     <- emc[[1]]$model
  log_K     <- log(K)

  out <- matrix(NA_real_, nrow = n_iter, ncol = n_subj)
  colnames(out) <- subjects

  for (iter in seq_len(n_iter)) {
    props <- tryCatch(
      MASS::mvrnorm(K, theta_mu[, iter], theta_var[,, iter]),
      error = function(e) NULL
    )
    if (is.null(props)) next   # non-PD Sigma; leave row NA

    for (s in seq_along(subjects)) {
      ll_mat <- calc_ll_pw(props, data_list[[subjects[s]]], model)  # [K x n_trials]
      subj_lls <- rowSums(ll_mat)                                    # [K]
      out[iter, s] <- matrixStats::logSumExp(subj_lls) - log_K
    }
  }
  out
}

# Per-subject WAIC: [n_iter x n_trials] log-likelihood matrix -> scalar WAIC
waic_from_ll <- function(ll_mat) {
  .loo_call(function() loo::waic(ll_mat)$estimates["waic", "Estimate"])
}

# Pooled or subject-level PSIS-LOO from a pointwise log-likelihood matrix.
loo_from_ll <- function(ll_mat, cores = 1) {
  .loo_call(function() loo::loo(ll_mat, cores = cores)$estimates["looic", "Estimate"])
}

# Per-subject WAIC: [n_iter x n_trials] log-likelihood matrix -> scalar WAIC
waic_subject <- function(emc, stage = "sample", filter = 0, subject) {
  ll_mat <- .ll_matrix_subject(emc, stage = stage, filter = filter, subject = subject)
  waic_from_ll(ll_mat)
}

# Pooled WAIC: per-trial LLs for "trial"; marginal subject LLs for "subject".
waic_pooled <- function(emc, stage = "sample", filter = 0,
                        pointwise = c("trial", "subject"), K = 200, cores = 1) {
  pointwise <- match.arg(pointwise)
  ll_all <- if (pointwise == "trial")
    .ll_matrix_pooled(emc, stage = stage, filter = filter, cores = cores)
  else
    .marg_ll_matrix(emc, stage = stage, filter = filter, K = K)
  waic_from_ll(ll_all)
}

# Per-subject PSIS-LOO: [n_iter x n_trials] log-likelihood matrix -> scalar LOOIC
loo_subject <- function(emc, stage = "sample", filter = 0, subject) {
  ll_mat <- .ll_matrix_subject(emc, stage = stage, filter = filter, subject = subject)
  loo_from_ll(ll_mat)
}

# Pooled PSIS-LOO: per-trial LLs for "trial"; marginal subject LLs for "subject".
loo_pooled <- function(emc, stage = "sample", filter = 0,
                       pointwise = c("trial", "subject"), K = 200, cores = 1) {
  pointwise <- match.arg(pointwise)
  ll_all <- if (pointwise == "trial")
    .ll_matrix_pooled(emc, stage = stage, filter = filter, cores = cores)
  else
    .marg_ll_matrix(emc, stage = stage, filter = filter, K = K)
  loo_from_ll(ll_all, cores = cores)
}
