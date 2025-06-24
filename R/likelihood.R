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
  lds[is.na(lds) | !ok] <- min_ll
  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
    ll[is.na(ll)] <- min_ll
    return(sum(pmax(min_ll,ll[attr(dadm,"expand")])))
  } else return(sum(pmax(min_ll,lds[attr(dadm,"expand")])))
}


log_likelihood_race_cens_trunc <- function(pars, dadm, model, min_ll = log(1e-10)) {
  # Hardcoded censoring and truncation values for now
  LT <- attr(dadm,"LT"); if (is.null(LT)) LT <- 0
  UT <- attr(dadm,"UT"); if (is.null(UT)) UT <- Inf
  LC <- attr(dadm,"LC"); if (is.null(LC)) LC <- 0
  UC <- attr(dadm,"UC"); if (is.null(UC)) UC <- Inf

  if (any(names(dadm) == "RACE")) {
    pars[as.numeric(dadm$lR) > as.numeric(as.character(dadm$RACE)), ] <- NA
  }

  ok_params <- if (!is.null(attr(pars, "ok"))) attr(pars, "ok") else rep(TRUE, dim(pars)[1])
  n_acc <- length(levels(dadm$R)) # Number of accumulators
  n_unique_trials <- dim(dadm)[1]/n_acc # Number of unique trial conditions

  ll_unique <- numeric(n_unique_trials) # Initialize vector to store log-likelihoods for each unique trial

  # Original single-trial integrand function (remains for use by stats::integrate)
  # t: time point
  # p_trial_this_winner_first: parameter matrix for a single trial, ordered with winner first
  # model: model object containing dfun (density) and pfun (CDF)
  f_race_integrand_single <- function(t, p_trial_this_winner_first, model) {
    p_winner <- p_trial_this_winner_first[1, , drop = FALSE] # Winner's parameters

    pdf_winner <- model$dfun(rt = t, pars = p_winner) # Density for winner
    pdf_winner[is.na(pdf_winner) | !is.finite(pdf_winner) | pdf_winner < 0] <- 0

    if (nrow(p_trial_this_winner_first) > 1) { # If there are losers
      survivor_losers <- 1
      for (i in 2:nrow(p_trial_this_winner_first)) { # Iterate over losers
        p_loser_i <- p_trial_this_winner_first[i, , drop = FALSE]
        s_loser_i <- (1 - model$pfun(rt = t, pars = p_loser_i)) # Survivor function for loser i
        s_loser_i[is.na(s_loser_i) | !is.finite(s_loser_i) | s_loser_i < 0] <- 0
        s_loser_i[s_loser_i > 1] <- 1 # Clamp survivor function (should not be >1 if CDF is valid)
        survivor_losers <- survivor_losers * s_loser_i
        if (survivor_losers == 0) break # Optimization
      }
      out_val <- pdf_winner * survivor_losers # Product of winner's PDF and losers' survivor functions
    } else { # Only one accumulator
      out_val <- pdf_winner
    }
    out_val[is.na(out_val) | !is.finite(out_val) | out_val < 0] <- 0 # Ensure non-negative finite output
    return(out_val)
  }

  # Batch integrand function for finite RTs
  # rts_vec: a vector of RTs for the batch
  # pars_list_ordered: a list of parameter matrices, each ordered (winner first), corresponding to rts_vec
  # model: model object
  # n_acc: number of accumulators
  f_race_integrand_R_batch <- function(rts_vec, pars_list_ordered, model, n_acc) {
    n_trials_batch <- length(rts_vec)
    if (n_trials_batch == 0) return(numeric(0))
    if (length(pars_list_ordered) != n_trials_batch) {
      stop("Length of rts_vec and pars_list_ordered must match in f_race_integrand_R_batch.")
    }

    out_vals <- numeric(n_trials_batch) # Initialize results vector

    # Loop through the batch (vectorization across trials with varying params is complex at this level)
    for (i in 1:n_trials_batch) {
      t_i <- rts_vec[i] # Current RT
      p_trial_i_ordered <- pars_list_ordered[[i]] # Current trial's ordered parameters

      p_winner_i <- p_trial_i_ordered[1, , drop = FALSE]
      pdf_winner_i <- model$dfun(rt = t_i, pars = p_winner_i) # Assumes model$dfun handles single t_i
      pdf_winner_i[is.na(pdf_winner_i) | !is.finite(pdf_winner_i) | pdf_winner_i < 0] <- 0

      if (n_acc > 1) { # If there are loser accumulators
        survivor_losers_i <- 1
        for (k_loser in 2:n_acc) { # Iterate over losers
          p_loser_k <- p_trial_i_ordered[k_loser, , drop = FALSE]
          s_loser_k <- (1 - model$pfun(rt = t_i, pars = p_loser_k)) # Assumes model$pfun handles single t_i
          s_loser_k[is.na(s_loser_k) | !is.finite(s_loser_k) | s_loser_k < 0] <- 0
          s_loser_k[s_loser_k > 1] <- 1
          survivor_losers_i <- survivor_losers_i * s_loser_k
          if (survivor_losers_i == 0) break # Optimization
        }
        out_vals[i] <- pdf_winner_i * survivor_losers_i
      } else { # Only one accumulator
        out_vals[i] <- pdf_winner_i
      }
    }
    out_vals[is.na(out_vals) | !is.finite(out_vals) | out_vals < 0] <- 0 # Ensure non-negative finite output
    return(out_vals)
  }
  
  # Helper to order parameters: winner first, then losers
  order_pars_for_winner <- function(p_all_acc, k_idx, unique_trial_id_for_error = NA) {
    if (is.na(k_idx) || k_idx < 1 || k_idx > nrow(p_all_acc)) {
      err_msg <- "Invalid k_idx in order_pars_for_winner"
      if (!is.na(unique_trial_id_for_error)) {
        err_msg <- paste0(err_msg, ": ", k_idx, " for unique trial ", unique_trial_id_for_error)
      } else {
        err_msg <- paste0(err_msg, ": ", k_idx)
      }
      stop(err_msg)
    }
    return(rbind(p_all_acc[k_idx, , drop=FALSE], p_all_acc[-k_idx, , drop=FALSE]))
  }
  
  # Helper for numerical integration for a k_th winner
  integrate_for_kth_winner <- function(k_winner_idx, p_all_acc, low, upp, model, unique_trial_id_for_error = NA) {
    if (low >= upp) return(0) # Integral over zero or negative range is 0
    if (is.na(k_winner_idx) || k_winner_idx < 1 || k_winner_idx > nrow(p_all_acc) ) return(0) # Invalid winner index
    
    pars_ordered <- order_pars_for_winner(p_all_acc, k_winner_idx, unique_trial_id_for_error)
    
    # Test integrand at a point to catch immediate issues (e.g., from bad parameters)
    test_t <- if (is.finite(low) && is.finite(upp)) (low + upp) / 2 else if (is.finite(low)) low + 1 else if (is.finite(upp)) upp -1 else 1
    if (test_t < 0 && low == 0) test_t <- 1e-6
    if (test_t < low && low >= 0) test_t <- low + 1e-6
    if (test_t > upp && upp >= 0) test_t <- upp - 1e-6
    
    integrand_val_at_test_t <- f_race_integrand_single(t = test_t, p_trial_this_winner_first = pars_ordered, model = model)
    if (is.na(integrand_val_at_test_t) || !is.finite(integrand_val_at_test_t) || integrand_val_at_test_t < 0) {
      # message(paste("Integrand is bad at test point for trial", unique_trial_id_for_error, "k_winner", k_winner_idx))
    }
    
    # Perform numerical integration
    res <- suppressWarnings(try(stats::integrate(f_race_integrand_single,
                                                 lower = low, upper = upp,
                                                 p_trial_this_winner_first = pars_ordered,
                                                 model = model,
                                                 rel.tol = .Machine$double.eps^0.4
    )$value, silent = TRUE))
    if (inherits(res, "try-error") || is.na(res) || !is.finite(res) || res < 0) return(0) # Handle integration errors
    return(res)
  }
  
  # Helper to calculate truncation correction factor
  get_trunc_corr_factor_for_kth_winner <- function(k_winner_idx, p_all_acc, model, unique_trial_id_for_error = NA) {
    if (!(LT > 0 || UT < Inf)) return(1.0) # No truncation if bounds are 0 to Inf
    if (is.na(k_winner_idx)) return(NA_real_) # Cannot calculate if winner is unknown
    
    # Integral over [0, Inf) - total probability for this winner without truncation
    prob_untruncated <- integrate_for_kth_winner(k_winner_idx, p_all_acc, 0, Inf, model, unique_trial_id_for_error)
    # Integral over [LT, UT] - probability within the truncation interval
    prob_truncated_interval <- integrate_for_kth_winner(k_winner_idx, p_all_acc, LT, UT, model, unique_trial_id_for_error)
    
    if (prob_truncated_interval > 1e-12) { # If probability in truncated interval is non-negligible
      if(is.na(prob_untruncated) || prob_untruncated < 0) return(NA_real_)
      return(prob_untruncated / prob_truncated_interval) # Correction factor
    } else { # Denominator (prob_truncated_interval) is effectively zero
      if (prob_untruncated > 1e-12) return(NA_real_) # Numerator is positive, but interval prob is zero -> issue
      return(1.0) # Both are zero, implies no density; effectively factor is 1 (or issue with params)
    }
  }

  # --- Main Likelihood Calculation ---

  # Prepare RTs and Response indices for all unique trials
  unique_trial_first_row_indices <- seq(1, nrow(dadm), by = n_acc)
  unique_trial_rts <- dadm$rt[unique_trial_first_row_indices]
  unique_trial_Rs <- as.integer(dadm$R[unique_trial_first_row_indices])

  ll_unique[] <- min_ll # Initialize all log-likelihoods to min_ll

  # Determine parameter validity for each unique trial block
  ok_params_unique_trial <- sapply(1:n_unique_trials, function(j) {
    all(ok_params[((j - 1) * n_acc + 1):(j * n_acc)])
  })

  # Identify indices for finite RT trials suitable for batching
  # These trials must have valid parameters, finite positive RT within truncation bounds, and a known winner.
  finite_rt_indices <- which(ok_params_unique_trial &
                             is.finite(unique_trial_rts) &
                             unique_trial_rts > 0 &
                             unique_trial_rts >= LT &
                             unique_trial_rts <= UT &
                             !is.na(unique_trial_Rs))

  # Identify other trials: those with bad parameters or not meeting finite_rt_indices criteria
  all_valid_param_indices <- which(ok_params_unique_trial)
  other_indices <- setdiff(all_valid_param_indices, finite_rt_indices)

  # --- Batch Processing for Finite RT Trials ---
  if (length(finite_rt_indices) > 0) {
    batch_rts_values <- unique_trial_rts[finite_rt_indices]
    batch_Rs_values <- unique_trial_Rs[finite_rt_indices]

    # Prepare lists of parameter matrices for the batch
    batch_pars_list_ordered <- vector("list", length(finite_rt_indices)) # For integrand calculation
    batch_pars_list_unord <- vector("list", length(finite_rt_indices))   # For truncation correction

    for (i in seq_along(finite_rt_indices)) {
      j_unique <- finite_rt_indices[i] # Original unique trial index
      current_trial_par_indices_map <- ((j_unique - 1) * n_acc + 1):(j_unique * n_acc)

      pars_j_all_acc <- pars[current_trial_par_indices_map, , drop = FALSE]
      batch_pars_list_unord[[i]] <- pars_j_all_acc
      batch_pars_list_ordered[[i]] <- order_pars_for_winner(pars_j_all_acc, batch_Rs_values[i], j_unique)
    }
    
    if (length(batch_rts_values) > 0) {
      # Calculate probability densities for the batch
      prob_densities_batch <- f_race_integrand_R_batch(batch_rts_values, batch_pars_list_ordered, model, n_acc)

      # Calculate truncation correction factors for the batch
      trunc_cfs_batch <- sapply(seq_along(finite_rt_indices), function(i) {
          j_unique <- finite_rt_indices[i]
          get_trunc_corr_factor_for_kth_winner(batch_Rs_values[i], batch_pars_list_unord[[i]], model, j_unique)
      })

      # Calculate log-likelihood for each trial in the batch
      for (k in seq_along(batch_rts_values)) {
        original_j_idx <- finite_rt_indices[k] # Map back to original unique trial index
        pd_k <- prob_densities_batch[k]
        tc_k <- trunc_cfs_batch[k]

        if (pd_k > .Machine$double.eps && !is.na(tc_k) && is.finite(tc_k) && tc_k > 0) {
          final_prob_k <- pd_k * tc_k
          ll_unique[original_j_idx] <- log(final_prob_k)
        } else {
          ll_unique[original_j_idx] <- min_ll
        }
        ll_unique[original_j_idx] <- max(min_ll, ll_unique[original_j_idx]) # Ensure not below min_ll
      }
    }
  }

  # --- Iterative Processing for Other Trials ---
  # (NA/Inf RTs, finite RTs outside truncation, or trials with issues not caught above but still valid params)
  if (length(other_indices) > 0) {
    for (j in other_indices) { # j is an original unique trial index with valid parameters
      # Parameters are already confirmed valid for trials in other_indices (which came from all_valid_param_indices)
      current_trial_par_indices_map <- ((j - 1) * n_acc + 1):(j * n_acc)
      pars_condition_j_all_acc <- pars[current_trial_par_indices_map, , drop = FALSE]
      rt_j <- unique_trial_rts[j]
      R_j_idx <- unique_trial_Rs[j]

      prob_density_or_mass <- 0 # Initialize probability for the current trial

      # Case 1: Finite RT but strictly outside truncation bounds (LT, UT).
      # These trials have zero probability density according to the truncated distribution.
      if (is.finite(rt_j) && rt_j > 0 && (rt_j < LT || rt_j > UT)) {
          prob_density_or_mass <- 0
      # Case 2: Fast censoring (RT = -Inf). Integral over [LT, LC].
      } else if (rt_j == -Inf) {
        if (!is.na(R_j_idx)) { # Winner known
          prob_density_or_mass <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, LT, LC, model, j)
          if (prob_density_or_mass > 0) { # Apply truncation correction if probability is non-zero
              trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model, j)
              if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
              else prob_density_or_mass <- prob_density_or_mass * trunc_cf
          }
        } else { # Winner unknown: sum over all possible winners
          current_sum_prob <- 0
          for (k_winner_loopvar in 1:n_acc) {
            p_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, LT, LC, model, j)
            if (p_k > 0) {
              trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model, j)
              if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k * trunc_cf_k)
            }
          }
          prob_density_or_mass <- current_sum_prob
        }
      # Case 3: Slow censoring (RT = Inf). Integral over [UC, UT].
      } else if (rt_j == Inf) {
        if (!is.na(R_j_idx)) { # Winner known
          prob_density_or_mass <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, UC, UT, model, j)
           if (prob_density_or_mass > 0) { # Apply truncation correction
              trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model, j)
              if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
              else prob_density_or_mass <- prob_density_or_mass * trunc_cf
          }
        } else { # Winner unknown: sum over all possible winners
          current_sum_prob <- 0
          for (k_winner_loopvar in 1:n_acc) {
            p_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, UC, UT, model, j)
            if (p_k > 0) {
              trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model, j)
              if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k * trunc_cf_k)
            }
          }
          prob_density_or_mass <- current_sum_prob
        }
      # Case 4: Missing RT (is.na(rt_j)). Sum of integrals over [LT, LC] and [UC, UT].
      } else if (is.na(rt_j)) {
        if (!is.na(R_j_idx)) { # Winner known
          prob_L <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, LT, LC, model, j)
          prob_U <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, UC, UT, model, j)
          prob_density_or_mass <- prob_L + prob_U
          if (prob_density_or_mass > 0) { # Apply truncation correction
              trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model, j)
              if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
              else prob_density_or_mass <- prob_density_or_mass * trunc_cf
          }
        } else { # Winner unknown: sum over all possible winners
          current_sum_prob <- 0
          for (k_winner_loopvar in 1:n_acc) {
            prob_L_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, LT, LC, model, j)
            prob_U_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, UC, UT, model, j)
            p_k_sum <- prob_L_k + prob_U_k
            if (p_k_sum > 0) {
               trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model, j)
               if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k_sum * trunc_cf_k)
            }
          }
          prob_density_or_mass <- current_sum_prob
        }
      }
      # All other conditions (e.g. rt_j is 0 or negative finite) result in prob_density_or_mass = 0.

      # Convert final probability to log-likelihood
      prob_density_or_mass <- max(0, prob_density_or_mass, na.rm = TRUE) # Ensure non-negative
      if (prob_density_or_mass > .Machine$double.eps) {
        ll_unique[j] <- log(prob_density_or_mass)
      } else {
        ll_unique[j] <- min_ll
      }
      ll_unique[j] <- max(min_ll, ll_unique[j]) # Ensure not less than min_ll
    }
  }

  # Final summation of log-likelihoods
  if (!is.null(dadm$N) && length(dadm$N) == n_unique_trials) { # If trial counts (N) are provided
    total_ll <- sum(ll_unique * dadm$N)
  } else if (!is.null(attr(dadm, "expand"))) { # If an expansion vector is provided
    total_ll <- sum(ll_unique[attr(dadm, "expand")])
  } else { # Default: sum unique likelihoods (each counted once)
    total_ll <- sum(ll_unique)
  }

  return(total_ll)
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
