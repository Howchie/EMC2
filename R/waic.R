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
  emc <- restore_custom_kernel_pointers(emc, quiet = TRUE)
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
  emc <- restore_custom_kernel_pointers(emc, quiet = TRUE)
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
  emc <- restore_custom_kernel_pointers(emc, quiet = TRUE)
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

# Per-trial weights that give every cell of the crossing of `balance_by` the
# same total weight. Weights are aligned with the columns of the pointwise
# log-likelihood matrix and average 1, so a weighted criterion equals the
# ordinary one when the cells are balanced.
.balance_weights <- function(emc, balance_by, subject = NULL) {
  dat <- get_data(emc)
  if (!is.data.frame(dat))
    stop("`balance_by` is not supported for joint models.")
  if (!is.character(balance_by) || length(balance_by) == 0L || anyNA(balance_by))
    stop("`balance_by` must be a character vector of data column names.")
  missing <- setdiff(balance_by, names(dat))
  if (length(missing) > 0L)
    stop("`balance_by` column(s) not in the data: ", paste(missing, collapse = ", "))
  if (!is.null(subject)) {
    if (is.numeric(subject)) subject <- names(emc[[1]]$data)[subject]
    dat <- dat[as.character(dat$subjects) == as.character(subject), , drop = FALSE]
  }
  cell <- interaction(dat[balance_by], drop = TRUE, sep = ":")
  n_c <- tabulate(cell)
  nrow(dat) / (nlevels(cell) * n_c[as.integer(cell)])
}

# Weighted total of a pointwise criterion; plain total without weights.
.weighted_ic <- function(pointwise, weights) {
  if (is.null(weights)) return(sum(pointwise))
  if (length(weights) != length(pointwise))
    stop("`balance_by` weights do not match the pointwise log-likelihood (",
         length(weights), " trials vs ", length(pointwise), " columns).")
  sum(weights * pointwise)
}

# Per-subject WAIC: [n_iter x n_trials] log-likelihood matrix -> scalar WAIC
waic_from_ll <- function(ll_mat, weights = NULL) {
  .loo_call(function() {
    fit <- loo::waic(ll_mat)
    if (is.null(weights)) fit$estimates["waic", "Estimate"]
    else .weighted_ic(fit$pointwise[, "waic"], weights)
  })
}

# Pooled or subject-level PSIS-LOO from a pointwise log-likelihood matrix.
loo_from_ll <- function(ll_mat, cores = 1, weights = NULL) {
  .loo_call(function() {
    fit <- if (cores > 1L)
      .emc_wpool_run_loo(ll_mat, cores = cores)
    else
      NULL
    # If a clean wrapper cannot be created, retain correctness and avoid
    # forking from compare()'s dirty frame by running PSIS serially.
    if (is.null(fit)) fit <- loo::loo(ll_mat, cores = 1L)
    if (is.null(weights)) fit$estimates["looic", "Estimate"]
    else .weighted_ic(fit$pointwise[, "looic"], weights)
  })
}

#' WAIC and PSIS-LOO for a single fitted model
#'
#' Compute WAIC or PSIS-LOO for one `emc` fit, either for a single subject or
#' pooled over all subjects. These are the per-model building blocks behind the
#' `WAIC` and `LOO` columns of \code{\link{compare}()}; use them when you only
#' need one model's criterion. Both are on the deviance scale (lower is better).
#'
#' Pointwise log-likelihoods are read from the fit if they were stored during
#' sampling, and recomputed from the posterior otherwise. Suppressed Pareto-k
#' and p_waic diagnostics can be retrieved with \code{\link{waic_warnings}()}.
#'
#' @param emc An emc object.
#' @param stage Character. The sampling stage to use. Defaults to `"sample"`.
#' @param filter Integer or numeric vector. Iterations to remove from the start
#' of `stage`, or a vector of iterations to keep, as in \code{\link{get_pars}()}.
#' @param subject Subject name or index to compute the criterion for.
#' @param pointwise Character string, `"trial"` (default) or `"subject"`.
#' `"trial"` treats each trial as one pointwise unit. `"subject"` gives the
#' hierarchical leave-one-subject-out criterion: each subject's log-likelihood
#' is marginalised over the group distribution using `K` importance draws per
#' posterior iteration. Single-level fits fall back to summed conditional
#' subject log-likelihoods.
#' @param K Integer (default 200). Importance draws from the group distribution
#' per posterior iteration when `pointwise = "subject"`.
#' @param cores Integer (default 1). Cores for the pointwise log-likelihood
#' computation and, for LOO, for PSIS.
#' @param balance_by `NULL` (default) or a character vector of data column
#' names. If supplied, every cell of the crossing of these columns receives
#' equal total weight: trial \eqn{i} in cell \eqn{c} is weighted
#' \eqn{N / (C n_c)}, where \eqn{N} is the number of trials, \eqn{C} the number
#' of observed cells and \eqn{n_c} the trials in cell \eqn{c}. The criterion
#' then estimates predictive accuracy for a design in which all cells are
#' equally frequent, on the same scale as the unweighted criterion (the two
#' coincide when cells are balanced). Requires `pointwise = "trial"`.
#'
#' @return A single number: the WAIC estimate (`waic_*`) or the LOOIC estimate
#' (`loo_*`).
#' @seealso \code{\link{compare}()}, \code{\link{waic_warnings}()}
#' @examples
#' \dontrun{
#' waic_pooled(samples_LNR)
#' loo_subject(samples_LNR, subject = 1)
#' loo_pooled(samples_LNR, pointwise = "subject", cores = 4)
#' loo_pooled(samples_LNR, balance_by = "S")
#' }
#' @name waic_loo
NULL

#' @rdname waic_loo
#' @export
waic_subject <- function(emc, stage = "sample", filter = 0, subject,
                         balance_by = NULL) {
  ll_mat <- .ll_matrix_subject(emc, stage = stage, filter = filter, subject = subject)
  weights <- if (is.null(balance_by)) NULL else
    .balance_weights(emc, balance_by, subject = subject)
  waic_from_ll(ll_mat, weights = weights)
}

#' @rdname waic_loo
#' @export
waic_pooled <- function(emc, stage = "sample", filter = 0,
                        pointwise = c("trial", "subject"), K = 200, cores = 1,
                        balance_by = NULL) {
  pointwise <- match.arg(pointwise)
  if (!is.null(balance_by) && pointwise == "subject")
    stop("`balance_by` requires `pointwise = \"trial\"`.")
  ll_all <- if (pointwise == "trial")
    .ll_matrix_pooled(emc, stage = stage, filter = filter, cores = cores)
  else
    .marg_ll_matrix(emc, stage = stage, filter = filter, K = K, cores = cores)
  weights <- if (is.null(balance_by)) NULL else .balance_weights(emc, balance_by)
  waic_from_ll(ll_all, weights = weights)
}

#' @rdname waic_loo
#' @export
loo_subject <- function(emc, stage = "sample", filter = 0, subject,
                        balance_by = NULL) {
  ll_mat <- .ll_matrix_subject(emc, stage = stage, filter = filter, subject = subject)
  weights <- if (is.null(balance_by)) NULL else
    .balance_weights(emc, balance_by, subject = subject)
  loo_from_ll(ll_mat, weights = weights)
}

#' @rdname waic_loo
#' @export
loo_pooled <- function(emc, stage = "sample", filter = 0,
                       pointwise = c("trial", "subject"), K = 200, cores = 1,
                       balance_by = NULL) {
  pointwise <- match.arg(pointwise)
  if (!is.null(balance_by) && pointwise == "subject")
    stop("`balance_by` requires `pointwise = \"trial\"`.")
  ll_all <- if (pointwise == "trial")
    .ll_matrix_pooled(emc, stage = stage, filter = filter, cores = cores)
  else
    .marg_ll_matrix(emc, stage = stage, filter = filter, K = K, cores = cores)
  weights <- if (is.null(balance_by)) NULL else .balance_weights(emc, balance_by)
  loo_from_ll(ll_all, cores = cores, weights = weights)
}
