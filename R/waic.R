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
    message("pw_ll cache incomplete (not all chains have it) - recomputing. ",
            "Check that add_pw_ll() or save_pw_ll=TRUE completed without errors.")
  # Use the raw array form here.  The mcmc-list form is not defined for
  # single-level fits, while the array is uniform across single and
  # hierarchical objects: [parameter, subject, iteration].
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE,
                    return_mcmc = FALSE)
  data  <- emc[[1]]$data
  model <- emc[[1]]$model
  if (is.array(alpha) && length(dim(alpha)) >= 3L) {
    subjects <- dimnames(alpha)[[2]]
    if (is.null(subjects)) subjects <- names(data)
    alpha_by_subject <- setNames(lapply(seq_along(subjects), function(j) {
      a <- alpha[, j, , drop = TRUE]
      if (!is.matrix(a)) a <- matrix(a, nrow = dim(alpha)[1L])
      a
    }), subjects)
  } else if (is.list(alpha)) {
    alpha_by_subject <- alpha
    subjects <- names(alpha_by_subject)
    if (is.null(subjects)) subjects <- names(data)
    names(alpha_by_subject) <- subjects
  } else {
    stop("Could not extract subject-level alpha samples for pointwise likelihood")
  }

  # `mc.preschedule = FALSE` gives every subject a fresh child. A child
  # computes one subject/model matrix, returns it, and exits; only this master
  # process retains the matrices for comparison. This is the same lifecycle
  # as the sampler's short-lived fallback and avoids a worker accumulating a
  # subject queue across the model.
  ll_list <- auto_mclapply(subjects, function(sub) {
    proposals <- if (is.matrix(alpha_by_subject[[sub]]))
      t(alpha_by_subject[[sub]]) else do.call(rbind, alpha_by_subject[[sub]])
    calc_ll_pw(proposals, data[[sub]], model)
  }, mc.cores = cores, mc.preschedule = FALSE)
  names(ll_list) <- subjects
  if (any(!vapply(ll_list, is.matrix, logical(1)))) {
    bad <- subjects[!vapply(ll_list, is.matrix, logical(1))]
    stop("Pointwise likelihood failed for subject(s): ",
         paste(bad, collapse = ", "))
  }
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
.marg_ll_matrix <- function(emc, stage = "sample", filter = 0, K = 200,
                            cores = 1) {
  theta_mu  <- get_pars(emc, selection = "mu",    stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)

  if (is.null(theta_mu) || is.null(theta_var)) {
    # Single-level model: marginalisation isn't possible; sum conditional LLs.
    ll_mat    <- .ll_matrix_pooled(emc, stage = stage, filter = filter,
                                   cores = cores)
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

  # Draw the proposal sets once so serial and parallel execution use the same
  # Monte Carlo draws.  The likelihood work is independent across subjects.
  proposals <- lapply(seq_len(n_iter), function(iter) {
    tryCatch(
      MASS::mvrnorm(K, theta_mu[, iter], theta_var[,, iter]),
      error = function(e) NULL
    )
  })

  subject_ll <- function(s) {
    vapply(seq_len(n_iter), function(iter) {
      props <- proposals[[iter]]
      if (is.null(props)) return(NA_real_)
      ll_mat <- calc_ll_pw(props, data_list[[subjects[s]]], model)
      subj_lls <- rowSums(ll_mat)
      matrixStats::logSumExp(subj_lls) - log_K
    }, numeric(1))
  }

  n_cores <- max(1L, min(as.integer(cores), n_subj))
  # Each subject is also an independent one-shot job for the marginal path.
  ll_list <- auto_mclapply(seq_along(subjects), subject_ll,
                           mc.cores = n_cores, mc.preschedule = FALSE)
  out <- do.call(cbind, ll_list)
  colnames(out) <- subjects
  out
}

# Per-subject WAIC: [n_iter x n_trials] log-likelihood matrix -> scalar WAIC
waic_from_ll <- function(ll_mat) {
  .loo_call(function() loo::waic(ll_mat)$estimates["waic", "Estimate"])
}

# Pooled or subject-level PSIS-LOO from a pointwise log-likelihood matrix.
loo_from_ll <- function(ll_mat, cores = 1) {
  .loo_call(function() {
    fit <- if (cores > 1L)
      .emc_wpool_run_loo(ll_mat, cores = cores)
    else
      NULL
    # If a clean wrapper cannot be created, retain correctness and avoid
    # forking from compare()'s dirty frame by running PSIS serially.
    if (is.null(fit)) fit <- loo::loo(ll_mat, cores = 1L)
    fit$estimates["looic", "Estimate"]
  })
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
    .marg_ll_matrix(emc, stage = stage, filter = filter, K = K, cores = cores)
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
    .marg_ll_matrix(emc, stage = stage, filter = filter, K = K, cores = cores)
  loo_from_ll(ll_all, cores = cores)
}
