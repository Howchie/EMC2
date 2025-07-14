calc_ll_R <- function(p_vector, model, dadm){
  if(!is.null(model$transform)){
    pars <- get_pars_matrix(p_vector, dadm, model)
  } else{
    pars <- p_vector
  }
  ll <- model$log_likelihood(pars, dadm, model)
  return(ll)
}


log_likelihood_race <- function(pars,dadm,model,min_ll=log(1e-10))
  # Race model summed log likelihood
{
  if (any(names(dadm)=="RACE")){# Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA
  }
  
  if (is.null(attr(pars,"ok"))){
    ok <- !logical(dim(pars)[1])
  } else ok <- attr(pars,"ok")

  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(model$dfun(rt=dadm$rt[dadm$winner],
                                     pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <- log(1-model$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[!ok] <- min_ll
  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- rowSums(cbind(ll, lds[!dadm$winner]), na.rm=TRUE)
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,function(x){sum(x,na.rm = TRUE)})
    }
    ll[is.na(ll)] <- min_ll
    return(sum(pmax(min_ll,ll[attr(dadm,"expand")]),na.rm = TRUE))
  } else return(sum(pmax(min_ll,lds[attr(dadm,"expand")]),na.rm = TRUE))
}


# log_likelihood_race_cens_trunc <- function(pars, dadm, model, min_ll = log(1e-10)) {
#   # Hardcoded censoring and truncation values for now
#   LT <- dadm$LT;UT <- dadm$UT;LC <- dadm$LC;UC <- dadm$UC
#   cn=colnames(pars)
#   if (any(names(dadm) == "RACE")) {
#     pars[as.numeric(dadm$lR) > as.numeric(as.character(dadm$RACE)), ] <- NA
#   }
#   
#   ok_params <- if (!is.null(attr(pars, "ok"))) attr(pars, "ok") else rep(TRUE, dim(pars)[1])
#   n_acc <- length(levels(dadm$R))
#   n_unique_trials <- dim(dadm)[1]/n_acc
#   
#   ll_unique <- numeric(n_unique_trials)
#   
#   f_race_integrand <- function(t, p_trial_this_winner_first, model) {
#     p_winner <- p_trial_this_winner_first[1, , drop = FALSE]
#     
#     pdf_winner <- model$dfun(rt = t, pars = p_winner)
#     pdf_winner[is.na(pdf_winner) | !is.finite(pdf_winner) | pdf_winner < 0] <- 0
#     
#     if (nrow(p_trial_this_winner_first) > 1) {
#       survivor_losers <- 1
#       for (i in 2:nrow(p_trial_this_winner_first)) {
#         p_loser_i <- p_trial_this_winner_first[i, , drop = FALSE]
#         s_loser_i <- (1 - model$pfun(rt = t, pars = p_loser_i))
#         s_loser_i[s_loser_i > 1] <- 1
#         if(any((is.na(s_loser_i) | !is.finite(s_loser_i) | s_loser_i < 0))){
#           next # don't set ll to 0 (mucks up RACE)
#         }
#         survivor_losers <- survivor_losers * s_loser_i
#       }
#       out_val <- pdf_winner * survivor_losers
#     } else {
#       out_val <- pdf_winner
#     }
#     out_val[is.na(out_val) | !is.finite(out_val) | out_val < 0] <- 0
#     return(out_val)
#   }
#   
#   integrate_for_kth_winner <- function(k_winner_idx, p_all_acc, low, upp, model) {
#     if (low >= upp) return(0)
#     if (is.na(k_winner_idx) || k_winner_idx < 1 || k_winner_idx > nrow(p_all_acc) ) return(0)
#     
#     pars_ordered <- order_pars_for_winner(p_all_acc, k_winner_idx)
#     
#     # Basic check on integrand before calling integrate
#     # Test point can be tricky if low/upp are Inf.
#     test_t <- if (is.finite(low) && is.finite(upp)) (low + upp) / 2 else if (is.finite(low)) low + 1 else if (is.finite(upp)) upp -1 else 1
#     if (test_t < 0 && low == 0) test_t <- 1e-6 # Avoid negative if low is 0
#     if (test_t < low && low >= 0) test_t <- low + 1e-6 # Ensure test_t is within bounds or sensible
#     if (test_t > upp && upp >= 0) test_t <- upp - 1e-6
#     
#     integrand_val_at_test_t <- f_race_integrand(t = test_t, p_trial_this_winner_first = pars_ordered, model = model)
#     if (is.na(integrand_val_at_test_t) || !is.finite(integrand_val_at_test_t) || integrand_val_at_test_t < 0) {
#       # If integrand is already bad, no point integrating
#       # This might happen due to bad parameters leading to NA/Inf in dfun/pfun
#     }
#     
#     res <- suppressWarnings(try(stats::integrate(f_race_integrand,
#                                                  lower = low, upper = upp,
#                                                  p_trial_this_winner_first = pars_ordered,
#                                                  model = model,
#                                                  rel.tol = .Machine$double.eps^0.4
#     )$value, silent = TRUE))
#     if (inherits(res, "try-error") || is.na(res) || !is.finite(res) || res < 0) return(0)
#     return(res)
#   }
#   
#   get_trunc_corr_factor_for_kth_winner <- function(k_winner_idx, p_all_acc, model) {
#     if (!(LT > 0 || UT < Inf)) return(1.0)
#     if (is.na(k_winner_idx)) return(NA_real_)
#     
#     prob_untruncated <- integrate_for_kth_winner(k_winner_idx, p_all_acc, 0, Inf, model)
#     prob_truncated_interval <- integrate_for_kth_winner(k_winner_idx, p_all_acc, LT, UT, model)
#     
#     if (prob_truncated_interval > 1e-12) {
#       if(is.na(prob_untruncated) || prob_untruncated < 0) return(NA_real_)
#       return(prob_untruncated / prob_truncated_interval)
#     } else {
#       if (prob_untruncated > 1e-12) return(NA_real_)
#       return(1.0)
#     }
#   }
#   
#   # Separate unique trials into those that can be batched and those that need iteration
#   finite_rt_trial_indices <- integer(0)
#   other_rt_trial_indices <- integer(0)
#   
#   # Store data for batch processing
#   # These will be populated only for trials in observed_finite_rt_trial_indices
#   batch_rts_list <- list()
#   batch_pars_ordered_list <- list() # List of parameter matrices (winner first)
#   batch_pars_unord_list <- list()   # List of original parameter matrices (for truncation)
#   batch_R_j_idx_list <- list()      # List of winner indices
#   
#   for (j in 1:n_unique_trials) {
#     current_trial_par_indices <- ((j - 1) * n_acc + 1):(j * n_acc)
#     
#     rt_j <- dadm$rt[current_trial_par_indices[1]]
#     R_j_idx <- as.integer(dadm$R[current_trial_par_indices[1]]) # Ensure it's integer
#     
#     # Criteria for batching: finite RT, positive, known winner, within truncation bounds
#     if (is.finite(rt_j) && rt_j > 0 && !is.na(R_j_idx) && rt_j >= LT && rt_j <= UT) {
#       finite_rt_trial_indices <- c(finite_rt_trial_indices, current_trial_par_indices)
#       # Store necessary info for batch processing later
#       # Note: actual parameter extraction will happen once we know the size of the batch
#     } else {
#       other_rt_trial_indices <- c(other_rt_trial_indices, current_trial_par_indices)
#     }
#   }
#   
#   # --- Batch process observed finite RT trials ---
#   if (length(observed_finite_rt_trial_indices) > 0) {
#     num_batch_trials <- length(finite_rt_trial_indices)
# 
#     # Vectorized call for winner PDFs
#     # Assuming model$dfun can handle a vector of rts and a matrix of corresponding parameters
#     win <- model$dfun(rt = dadm$rt[finite_rt_trial_indices&dadm$winner==1], pars = pars[finite_rt_trial_indices,])
#     
#     # Calculate survivor probabilities for losers
#     survivor_losers_batch <- rep(1.0, num_batch_trials) # Initialize with 1s
#     if (n_acc > 1) {
#       for (k_loser in 1:(n_acc - 1)) {
#         # Vectorized call for loser CDFs
#         cdf_loser_k_batch <- model$pfun(rt = batch_rts, pars = loser_pars_batch_list[[k_loser]])
#         s_loser_k_batch <- (1 - cdf_loser_k_batch)
#         ok<-!is.na(s_loser_k_batch) & is.finite(s_loser_k_batch) & s_loser_k_batch >= 0 & s_loser_k_batch <= 1
#         s_loser_k_batch[s_loser_k_batch > 1] <- 1
#         survivor_losers_batch[ok] <- survivor_losers_batch[ok] * s_loser_k_batch[ok]
#       }
#     }
#     
#     # Combine for final probability density for each trial in the batch
#     prob_density_batch <- pdf_winner_batch * survivor_losers_batch
#     prob_density_batch[is.na(prob_density_batch) | !is.finite(prob_density_batch) | prob_density_batch < 0] <- 0
#     
#     # Apply truncation correction (iteratively for now) and calculate log-likelihood
#     for (i in 1:num_batch_trials) {
#       unique_trial_idx <- observed_finite_rt_trial_indices[i]
#       current_prob_density <- prob_density_batch[i]
#       
#       if (current_prob_density > .Machine$double.eps) {
#         # Use the stored un-ordered parameters for this specific trial for truncation
#         pars_for_this_trial_trunc <- pars_unord_for_trunc_batch[[i]]
#         R_j_for_this_trial_trunc <- R_j_indices_batch[i]
#         
#         trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_for_this_trial_trunc, pars_for_this_trial_trunc, model)
#         if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) {
#           ll_unique[unique_trial_idx] <- min_ll
#         } else {
#           final_prob_val <- current_prob_density * trunc_cf
#           if (final_prob_val > .Machine$double.eps) {
#             ll_unique[unique_trial_idx] <- log(final_prob_val)
#           } else {
#             ll_unique[unique_trial_idx] <- min_ll
#           }
#         }
#       } else {
#         ll_unique[unique_trial_idx] <- min_ll
#       }
#       ll_unique[unique_trial_idx] <- max(min_ll, ll_unique[unique_trial_idx])
#     }
#   }
#   
#   for (j in other_rt_trial_indices) {
#     current_trial_par_indices <- ((j - 1) * n_acc + 1):(j * n_acc)
#     # Make sure parameters are valid
#     if (!all(ok_params[current_trial_par_indices])) {
#       ll_unique[j] <- min_ll
#       next
#     }
#     
#     pars_condition_j_all_acc <- pars[current_trial_par_indices, , drop = FALSE]
#     
#     rt_j <- dadm$rt[current_trial_par_indices[1]]
#     R_j_idx <- as.integer(dadm$R[current_trial_par_indices[1]])
#     
#     current_ll_j <- min_ll
#     prob_density_or_mass <- 0
#     if (rt_j == -Inf) { # fast censoring
#       if (!is.na(R_j_idx)) {
#         prob_density_or_mass <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, LT, LC, model)
#         if (prob_density_or_mass > 0) {
#           trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model)
#           if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
#           else prob_density_or_mass <- prob_density_or_mass * trunc_cf
#         }
#       } else {
#         current_sum_prob <- 0
#         for (k_winner_loopvar in 1:n_acc) {
#           p_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, LT, LC, model)
#           if (p_k > 0) {
#             trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model)
#             if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k * trunc_cf_k)
#           }
#         }
#         prob_density_or_mass <- current_sum_prob
#       }
#     } else if (rt_j == Inf) { # Slow censoring
#       if (!is.na(R_j_idx)) {
#         prob_density_or_mass <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, UC, UT, model)
#         if (prob_density_or_mass > 0) {
#           trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model)
#           if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
#           else prob_density_or_mass <- prob_density_or_mass * trunc_cf
#         }
#       } else {
#         current_sum_prob <- 0
#         for (k_winner_loopvar in 1:n_acc) {
#           p_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, UC, UT, model)
#           if (p_k > 0) {
#             trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model)
#             if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k * trunc_cf_k)
#           }
#         }
#         prob_density_or_mass <- current_sum_prob
#       }
#     } else if (is.na(rt_j)) {
#       if (!is.na(R_j_idx)) {
#         prob_L <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, LT, LC, model)
#         prob_U <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, UC, UT, model)
#         prob_density_or_mass <- prob_L + prob_U
#         if (prob_density_or_mass > 0) {
#           trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model)
#           if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
#           else prob_density_or_mass <- prob_density_or_mass * trunc_cf
#         }
#       } else {
#         current_sum_prob <- 0
#         for (k_winner_loopvar in 1:n_acc) {
#           prob_L_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, LT, LC, model)
#           prob_U_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, UC, UT, model)
#           p_k_sum <- prob_L_k + prob_U_k
#           if (p_k_sum > 0) {
#             trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model)
#             if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k_sum * trunc_cf_k)
#           }
#         }
#         prob_density_or_mass <- current_sum_prob
#       }
#     }
#     
#     prob_density_or_mass <- max(0, prob_density_or_mass)
#     
#     if (prob_density_or_mass > .Machine$double.eps) {
#       current_ll_j <- log(prob_density_or_mass)
#     } else {
#       current_ll_j <- min_ll
#     }
#     
#     ll_unique[j] <- max(min_ll, current_ll_j)
#   }
#   
#   
#   if (!is.null(dadm$N) && length(dadm$N) == n_unique_trials) {
#     total_ll <- sum(ll_unique * dadm$N)
#   } else if (!is.null(attr(dadm, "expand"))) {
#     total_ll <- sum(ll_unique[attr(dadm, "expand")])
#   } else {
#     total_ll <- sum(ll_unique)
#   }
#   
#   return(total_ll)
# }

# ──────────────────────────────────────────────────────────────────────────────
# Helpers ---------------------------------------------------------------------

## create a string key that uniquely identifies the parameter set + bounds
.make_key <- function(lower, upper, pars_mat)
  paste0(
    formatC(lower, digits = 12, width = 0, format = "fg"), "_",
    formatC(upper, digits = 12, width = 0, format = "fg"), "_",
    paste(formatC(as.numeric(pars_mat), digits = 12, width = 0, format = "fg"),
          collapse = "_")
  )

## integrate the joint density for the *k-th winner* on (lower, upper)
.integrate_kth_winner <- function(k, pars_mat, lower, upper, model,
                                  rel.tol = 1e-7, ...) {
  n_acc <- nrow(pars_mat)
  if (any(is.na(pars_mat[k,]))) return(0)
  idx_los <- setdiff(seq_len(n_acc), k)
  idx_los <- idx_los[ rowSums(is.na(pars_mat[idx_los, , drop=FALSE])) == 0L ]
  integrand <- function(t) {
    dk <- model$dfun(t, pars = pars_mat[k, , drop = FALSE])
    if (length(idx_los) > 0) {
      Sl <- 1 - model$pfun(t, pars = pars_mat[idx_los, , drop = FALSE])
      dk = dk * prod(Sl,na.rm=TRUE)
    }
    dk
  }
  pracma::quadinf(integrand, lower, upper, tol = rel.tol)$Q
}

## probability of *any* observable response in (lower, upper)
.truncation_normaliser <- function(pars_mat, lower, upper, model, rel.tol = 1e-5, ...) {

  present <- which(rowSums(is.na(pars_mat)) == 0L)
  n_acc <- length(present)
  p_tot <- 0
  for (k in present) {
    if (any(is.na(pars_mat[k,]))) {next}
    p_tot <- p_tot +
    .integrate_kth_winner(k, pars_mat[present,], lower, upper, model,
                          rel.tol = rel.tol)
  }
  
  max(p_tot, .Machine$double.eps)             # guard against 0
}

# ──────────────────────────────────────────────────────────────────────────────
# Main likelihood -------------------------------------------------------------

log_likelihood_race_cens_trunc <- function(pars,dadm,model,min_ll=log(1e-10)) {
  cache_env <- new.env(parent = emptyenv())
  posdrift = ifelse(model$c_name=="LBAIO",FALSE,TRUE)
  ## ── basic dimensions ──────────────────────────────────────────────────────
  n_trials <- nrow(dadm)
  n_acc    <- length(levels(dadm$R))
  stopifnot(n_trials %% n_acc == 0L)          # one row per acc. per trial
  n_unique <- n_trials / n_acc
  dadm$RACE_num <- if("RACE" %in% names(dadm)) as.numeric(as.character(dadm$RACE)) else rep(n_acc, n_trials)
  if (any(names(dadm)=="RACE")){# Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA
  }
  #dadm = dadm[dadm$RACE_num<=as.numeric(dadm$lR),] # drop missing accumulator rows
  if (is.null(attr(pars,"ok"))){
    ok <- !logical(dim(pars)[1])
  } else ok <- attr(pars,"ok")
  
  ## ── convenience column refs ───────────────────────────────────────────────
  RT <- dadm$rt;  LT <- dadm$LT;  UT <- dadm$UT
  LC <- dadm$LC;  UC <- dadm$UC;  R_idx <- dadm$R
  
  ## ── initialise per-row log densities ──────────────────────────────────────
  lds <- rep(NA, n_trials)
  ll_unique <- rep(NA,n_unique)
  ## ── batch: finite RT, in-bounds, known winner ────────────────────────────
  finite_mask <- is.finite(RT) & RT > 0 &
    RT >= LT & RT <= UT & !is.na(R_idx) & (as.numeric(dadm$lR) <= dadm$RACE_num)
  trunc_cens_mask = !(is.finite(RT)) & (as.numeric(dadm$lR) <= dadm$RACE_num)
  if (any(finite_mask)) {
    
    # pdf for winners
    lds[dadm$winner & finite_mask] <-
      log(model$dfun(rt   = RT[dadm$winner & finite_mask],
                     pars = pars[dadm$winner & finite_mask, , drop = FALSE]))
    
    # survivor for losers
    if (n_acc > 1L) {
      lds[!dadm$winner & finite_mask] <- log(1 - model$pfun(rt   = RT[!dadm$winner & finite_mask],
                                                            pars = pars[!dadm$winner & finite_mask, , drop = FALSE]))
    }
    
    # truncate-window correction, cached by unique trial
    for (j in seq_len(n_unique)) {
      idx  <- ((j - 1L) * n_acc + 1L):(j * n_acc)
      if (finite_mask[idx[1]]) {
        n_acc_j = dadm$RACE_num[idx[1]]
        idx=idx[1:n_acc_j]
        if (!(LT[idx[1]] == 0 && UT[idx[1]] == Inf)){
          key <- .make_key(LT[idx[1]], UT[idx[1]], pars[idx, , drop = FALSE])
          invZ <- cache_env[[key]]
          if (is.null(invZ)) {
            invZ <- 1 / .truncation_normaliser(pars[idx, , drop = FALSE],
                                    LT[idx[1]], UT[idx[1]], model,
                                    rel.tol = 1e-7)
            cache_env[[key]] <- invZ
          }
          lds[idx] <- lds[idx] + log(invZ)
        }
        ll_j <- sum(lds[idx])
        ll_unique[j] <- max(min_ll, ll_j)
      }
    }
  }
  
  ## ── process other trials one-by-one (-Inf, +Inf, NA) ─────────────────────
  prob_win <- function(k, lo, hi) {
    .integrate_kth_winner(k, pars[idx, , drop = FALSE],
                          lo, hi, model, rel.tol = 1e-7)
  }

  for (j in seq_len(n_unique)) {
    idx  <- ((j - 1L) * n_acc + 1L):(j * n_acc)
    if (!trunc_cens_mask[idx[[1]]]) {next}
    rt   <- RT[idx[1]]          # all rows of a trial share RT
    Rj   <- as.numeric(R_idx[idx[1]])
    pval <- 0
    
    
    n_acc_j = dadm$RACE_num[idx[1]]
    idx=idx[1:n_acc_j]
    if(is.na(Rj)){ks = 1:n_acc_j}else{ks = Rj}
    if (identical(rt, -Inf)) {                       # fast censor
      lo <- LT[idx[1]]; hi <- LC[idx[1]]
      pval <- sum(vapply(ks, prob_win, numeric(1), lo, hi))
      if (!posdrift && is.na(Rj)) {                   # intrinsic omission
        v  <- pars[idx, "v"]; sv <- pars[idx, "sv"]
        pI <- prod(pnorm(0, v, sv),na.rm=TRUE)
        pval <- pval + pI
      }
      
    } else if (identical(rt, Inf)) {                 # slow censor
      lo <- UC[idx[1]]; hi <- UT[idx[1]]
      pval <- sum(vapply(ks, prob_win, numeric(1), lo, hi))
      if (!posdrift && is.na(Rj)) {
        v  <- pars[idx, "v"]; sv <- pars[idx, "sv"]
        pI <- prod(pnorm(0, v, sv),na.rm=TRUE)
        pval <- pval + pI
      }
      
    } else if (is.na(rt)) {                          # missing RT
      lo1 <- LT[idx[1]]; hi1 <- LC[idx[1]]
      lo2 <- UC[idx[1]]; hi2 <- UT[idx[1]]
      pval <- sum(vapply(ks, function(k)
        prob_win(k, lo1, hi1) + prob_win(k, lo2, hi2), numeric(1)))
      if (!posdrift && is.na(Rj)) {
        v  <- pars[idx, "v"]; sv <- pars[idx, "sv"]
        pI <- prod(pnorm(0, v, sv),na.rm=TRUE)
        pval <- pval + pI
      }
    }                                               # 0 or negative RT ⇒ prob 0
    
    ## truncation correction (if RT unobserved)
    if (!(LT[idx[1]] == 0 && UT[idx[1]] == Inf) && pval > 0) {
      key  <- .make_key(LT[idx[1]], UT[idx[1]], pars[idx, , drop = FALSE])
      invZ <- cache_env[[key]]
      if (is.null(invZ)) {
        invZ <- 1 / .truncation_normaliser(pars[idx, , drop = FALSE],
                                LT[idx[1]], UT[idx[1]], model,
                                rel.tol = 1e-5, posdrift = posdrift)
        cache_env[[key]] <- invZ
      }
      pval <- pval * invZ
    }
    
    ll_unique[j] <- if (pval > .Machine$double.eps) log(pval) else min_ll
  }
  ## ── aggregate over copies (expand) or N column ───────────────────────────
  return(sum(pmax(min_ll,ll_unique[attr(dadm,"expand")]),na.rm = TRUE))
}



log_likelihood_ddm <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood, with protection against numerical issues
{
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok")))
    like[attr(pars,"ok")] <- model$dfun(dadm$rt[attr(pars,"ok")],dadm$R[attr(pars,"ok")],
                                        pars[attr(pars,"ok"),,drop=FALSE])
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}


log_likelihood_ddmgng <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood for go/nogo model
{
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok"))) {
    isna <- is.na(dadm$rt)
    ok <- attr(pars,"ok") & !isna
    like[ok] <- model$dfun(dadm$rt[ok],dadm$R[ok],pars[ok,,drop=FALSE])
    ok <- attr(pars,"ok") & isna
    like[ok] <- # dont terminate on go boundary before timeout
      pmax(0,pmin(1,(1-model$pfun(dadm$TIMEOUT[ok],dadm$Rgo[ok],pars[ok,,drop=FALSE]))))
    
  }
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}



#### sdt choice likelihoods ----

log_likelihood_sdt <- function(pars,dadm, model,lb=-Inf, min_ll=log(1e-10))
  # probability of ordered discrete choices based on integrals of a continuous
  # distribution between thresholds, with fixed lower bound for first response
  # lb. Upper bound for last response is a fixed value in threshold vector
{
  first <- dadm$lR==levels(dadm$lR)[1]
  last <- dadm$lR==levels(dadm$lR)[length(levels(dadm$lR))]
  pars[last,"threshold"] <- Inf
  # upper threshold
  ut <- pars[dadm$winner,"threshold"]
  # lower threshold fixed at lb for first response
  pars[first &  dadm$winner,"threshold"] <- lb
  # otherwise threshold of response before one made
  notfirst <- !first &  dadm$winner
  pars[notfirst,"threshold"] <- pars[which(notfirst)-1,"threshold"]
  lt <- pars[dadm$winner,"threshold"]
  # log probability
  ll <- numeric(sum(dadm$winner))
  if (!is.null(attr(pars,"ok"))) { # Bad parameter region
    ok <- attr(pars,"ok")
    okw <- ok[dadm$winner]
    ll[ok] <- log(model$pfun(lt=lt[okw],ut=ut[okw],pars=pars[dadm$winner & ok,,drop=FALSE]))
  } else ll <- log(model$pfun(lt=lt,ut=ut,pars=pars[dadm$winner,,drop=FALSE]))
  ll <- ll[attr(dadm,"expand")]
  ll[is.na(ll)] <- 0
  sum(pmax(min_ll,ll))
}

# Two options:
# 1) component = NULL in which case we do all likelihoods in one block
# 2) component = integer, in which case we are blocking the ll and only want that one
log_likelihood_joint <- function(proposals, dadms, model_list, component = NULL){
  parPreFixs <- unique(gsub("[|].*", "", colnames(proposals)))
  i <- 0
  k <- 0
  total_ll <- 0
  for(dadm in dadms){
    i <- i + 1
    if(is.null(component) || component == i){
      if(is.data.frame(dadm)){
        k <- k + 1
        # Sometimes designs are the same across models (typically in fMRI)
        # Instead of storing multiple, we just store a pointer to the first original design
        if(is.numeric(attr(dadm, "designs"))){
          ref_idx <- attr(dadm, "designs")
          attr(dadm, "designs") <- attr(dadms[[ref_idx]], "designs")
        }
        parPrefix <- parPreFixs[k]
        # Unfortunately indexing can't get quicker than this as far as I can tell.
        columns_to_use <- sapply(strsplit(colnames(proposals), "|", fixed = TRUE), function(x) x == parPrefix)[1,]
        currentPars <- proposals[,columns_to_use, drop = F]
        colnames(currentPars) <- gsub(".*[|]", "", colnames(currentPars))
        total_ll <- total_ll +  calc_ll_manager(currentPars, dadm, model_list[[i]])
      }
    }
  }
  return(total_ll)
}

# Likelihood function for a redundant target race model.
# Vectorized version.
#
# Assumes dadm is structured with 4 rows per original trial, where each row
# corresponds to one of four conceptual accumulators.
# - dadm$lR: A factor column identifying the role of each row (accumulator role).
#   Expected levels are "S1D", "S2D", "S1A", "S2A". The order of these levels
#   in the factor is not critical, as the code explicitly extracts by name.
# - dadm$R: A column containing the actual observed response ("yes" or "no")
#   for the original trial, replicated across the 4 rows.
# - dadm$rt: The observed RT, replicated across the 4 rows.
# - pars: Matrix aligned with dadm, rows providing base parameters (e.g., "v", "b")
#   for the accumulator specified in dadm$lR.
log_likelihood_redundant_target_race <- function(pars, dadm, model, min_ll = log(1e-10))
{
  # --- Input Validations ---
  if (is.null(dadm$rt)) stop("dadm$rt is missing.")
  if (is.null(dadm$lR)) stop("dadm$lR (factor identifying accumulator roles S1D, S2D, S1A, S2A) is missing.")
  if (is.null(dadm$R)) stop("dadm$R (observed 'yes'/'no' response) is missing.")
  if (is.null(attr(dadm, "expand"))) stop("attr(dadm, 'expand') is missing.")
  if (nrow(dadm) %% 4 != 0) {
    stop("nrow(dadm) must be a multiple of 4 (4 accumulator rows per original trial).")
  }
  
  # Expected accumulator roles for internal mapping
  internal_role_order <- c("S1D", "S2D", "S1A", "S2A")
  if (!all(internal_role_order %in% levels(as.factor(dadm$lR)))) { # Check dadm$lR
    missing_roles <- internal_role_order[!internal_role_order %in% levels(as.factor(dadm$lR))]
    stop(paste("dadm$lR is missing expected levels for accumulator roles:", paste(missing_roles, collapse=", "),
               ". Expected levels are S1D, S2D, S1A, S2A."))
  }
  
  # --- PDF and CDF Calculation (once for all rows) ---
  f_all <- model$dfun(dadm$rt, pars)
  F_all <- model$pfun(dadm$rt, pars)
  
  n_blocks <- nrow(dadm) / 4 # Number of original trials
  
  # --- Extract f/F values for each role into vectors of length n_blocks ---
  # This assumes that dadm is structured such that filtering by dadm$lR
  # results in vectors of length n_blocks, correctly ordered by original trial.
  f1 <- f_all[dadm$lR == "S1D"]; F1 <- F_all[dadm$lR == "S1D"]
  f2 <- f_all[dadm$lR == "S2D"]; F2 <- F_all[dadm$lR == "S2D"]
  f3 <- f_all[dadm$lR == "S1A"]; F3 <- F_all[dadm$lR == "S1A"]
  f4 <- f_all[dadm$lR == "S2A"]; F4 <- F_all[dadm$lR == "S2A"]
  
  # Validation: Check if all extracted vectors have the correct length (n_blocks)
  expected_len <- n_blocks
  if (any(sapply(list(f1, F1, f2, F2, f3, F3, f4, F4), length) != expected_len)) {
    stop("Length mismatch after extracting role-specific f/F values. Check dadm structure and lR factor.")
  }
  
  # --- Get Observed Responses (one per original trial/block) ---
  # Taking from the first row of each conceptual 4-row block using dadm$R.
  observed_responses_per_block <- dadm$R[seq(1, nrow(dadm), by = 4)] # dadm$R holds "yes"/"no"
  if (length(observed_responses_per_block) != n_blocks) {
    stop("Could not correctly extract one observed response per trial block from dadm$R.")
  }
  
  ll_block_values <- numeric(n_blocks) # Stores one LL value per original trial
  is_yes_response <- observed_responses_per_block == "yes"
  is_no_response <- observed_responses_per_block == "no"
  
  # --- Apply Formulas Vectorially ---
  if (any(is_yes_response)) {
    term1_yes <- f1[is_yes_response] * (1 - F2[is_yes_response]) + f2[is_yes_response] * (1 - F1[is_yes_response])
    term2_yes <- 1 - (F3[is_yes_response] * F4[is_yes_response])
    
    term1_yes <- pmax(term1_yes, .Machine$double.eps)
    term2_yes <- pmax(term2_yes, .Machine$double.eps)
    ll_block_values[is_yes_response] <- log(term1_yes) + log(term2_yes)
  }
  
  if (any(is_no_response)) {
    term1_no <- f3[is_no_response] * F4[is_no_response] + f4[is_no_response] * F3[is_no_response]
    term2_no <- (1 - F1[is_no_response])
    term3_no <- (1 - F2[is_no_response])
    
    term1_no <- pmax(term1_no, .Machine$double.eps)
    term2_no <- pmax(term2_no, .Machine$double.eps)
    term3_no <- pmax(term3_no, .Machine$double.eps)
    ll_block_values[is_no_response] <- log(term1_no) + log(term2_no) + log(term3_no)
  }
  
  unhandled_responses <- !(is_yes_response | is_no_response)
  if (any(unhandled_responses)) {
    original_indices_unhandled <- which(unhandled_responses)
    # Get the problematic dadm$R values that correspond to these unhandled blocks
    # Need to index observed_responses_per_block or dadm$R at block level
    problematic_values <- unique(observed_responses_per_block[unhandled_responses])
    warning(paste("Unhandled values in dadm$R (observed response) found:", paste(problematic_values, collapse=", "),
                  ". Corresponding log-likelihoods will be NA, then min_ll."))
    ll_block_values[unhandled_responses] <- NA
  }
  
  # --- Final Summation (same logic as before) ---
  # Replicate the block's LL to its constituent rows.
  original_trial_idx_for_row <- rep(1:n_blocks, each = 4) # Assuming dadm is sorted by block
  ll_for_dadm_rows <- ll_block_values[original_trial_idx_for_row]
  
  final_ll_values <- pmax(min_ll, ll_for_dadm_rows[attr(dadm, "expand")])
  final_ll_values[is.na(final_ll_values)] <- min_ll
  
  return(sum(final_ll_values))
}


