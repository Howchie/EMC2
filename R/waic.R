# Package-level environment for storing suppressed loo warnings
.loo_warn_store <- new.env(parent = emptyenv())
.loo_warn_store$msgs <- character(0)

# Helper: run thunk, capture any loo/pareto-k warnings, stash for later
.waic_call <- function(thunk) {
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

# Per-subject WAIC: [n_iter x n_trials] log-likelihood matrix -> scalar WAIC
waic_subject <- function(emc, stage = "sample", filter = 0, subject, r_cores = 1) {
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  proposals <- do.call(rbind, alpha[[subject]])   # [n_iter x n_pars]
  ll_mat <- calc_ll_manager_pw(proposals, emc[[1]]$data[[subject]], emc[[1]]$model, r_cores = r_cores)
  .waic_call(function() loo::waic(ll_mat)$estimates["waic", "Estimate"])
}

# Pooled WAIC: concatenate per-trial LLs across all subjects, then waic()
waic_pooled <- function(emc, stage = "sample", filter = 0, r_cores = 1) {
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  data  <- emc[[1]]$data
  model <- emc[[1]]$model
  ll_list <- auto_mclapply(names(alpha), function(sub) {
    proposals <- do.call(rbind, alpha[[sub]])     # [n_iter x n_pars]
    calc_ll_manager_pw(proposals, data[[sub]], model, r_cores = 1)  # [n_iter x n_trials_sub]
  }, mc.cores = r_cores)
  ll_all <- do.call(cbind, ll_list)              # [n_iter x total_trials]
  .waic_call(function() loo::waic(ll_all)$estimates["waic", "Estimate"])
}
