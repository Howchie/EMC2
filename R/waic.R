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
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  proposals <- do.call(rbind, alpha[[subject]])   # [n_iter x n_pars]
  calc_ll_manager_pw(proposals, emc[[1]]$data[[subject]], emc[[1]]$model)
}

# Concatenate pointwise log-likelihood matrices across subjects.
.ll_matrix_pooled <- function(emc, stage = "sample", filter = 0,
                              pointwise = c("trial", "subject")) {
  pointwise <- match.arg(pointwise)
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  data  <- emc[[1]]$data
  model <- emc[[1]]$model
  ll_list <- lapply(names(alpha), function(sub) {
    proposals <- do.call(rbind, alpha[[sub]])     # [n_iter x n_pars]
    calc_ll_manager_pw(proposals, data[[sub]], model)  # [n_iter x n_trials_sub]
  })
  if (pointwise == "trial") {
    return(do.call(cbind, ll_list))              # [n_iter x total_trials]
  }
  ll_by_subject <- lapply(ll_list, rowSums)
  out <- do.call(cbind, ll_by_subject)           # [n_iter x n_subjects]
  colnames(out) <- names(alpha)
  out
}

# Per-subject WAIC: [n_iter x n_trials] log-likelihood matrix -> scalar WAIC
waic_from_ll <- function(ll_mat) {
  .loo_call(function() loo::waic(ll_mat)$estimates["waic", "Estimate"])
}

# Pooled or subject-level PSIS-LOO from a pointwise log-likelihood matrix.
loo_from_ll <- function(ll_mat) {
  .loo_call(function() loo::loo(ll_mat)$estimates["looic", "Estimate"])
}

# Per-subject WAIC: [n_iter x n_trials] log-likelihood matrix -> scalar WAIC
waic_subject <- function(emc, stage = "sample", filter = 0, subject) {
  ll_mat <- .ll_matrix_subject(emc, stage = stage, filter = filter, subject = subject)
  waic_from_ll(ll_mat)
}

# Pooled WAIC: concatenate per-trial LLs across all subjects, then waic()
waic_pooled <- function(emc, stage = "sample", filter = 0,
                        pointwise = c("trial", "subject")) {
  ll_all <- .ll_matrix_pooled(emc, stage = stage, filter = filter, pointwise = pointwise)
  waic_from_ll(ll_all)
}

# Per-subject PSIS-LOO: [n_iter x n_trials] log-likelihood matrix -> scalar LOOIC
loo_subject <- function(emc, stage = "sample", filter = 0, subject) {
  ll_mat <- .ll_matrix_subject(emc, stage = stage, filter = filter, subject = subject)
  loo_from_ll(ll_mat)
}

# Pooled PSIS-LOO: concatenate per-trial LLs across all subjects, then loo()
loo_pooled <- function(emc, stage = "sample", filter = 0,
                       pointwise = c("trial", "subject")) {
  ll_all <- .ll_matrix_pooled(emc, stage = stage, filter = filter, pointwise = pointwise)
  loo_from_ll(ll_all)
}
