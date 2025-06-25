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
  n_acc <- length(levels(dadm$R))
  n_unique_trials <- dim(dadm)[1]/n_acc
  
  ll_unique <- numeric(n_unique_trials)
  
  f_race_integrand <- function(t, p_trial_this_winner_first, model) {
    p_winner <- p_trial_this_winner_first[1, , drop = FALSE]
    
    pdf_winner <- model$dfun(rt = t, pars = p_winner)
    pdf_winner[is.na(pdf_winner) | !is.finite(pdf_winner) | pdf_winner < 0] <- 0
    
    if (nrow(p_trial_this_winner_first) > 1) {
      survivor_losers <- 1
      for (i in 2:nrow(p_trial_this_winner_first)) {
        p_loser_i <- p_trial_this_winner_first[i, , drop = FALSE]
        s_loser_i <- (1 - model$pfun(rt = t, pars = p_loser_i))
        s_loser_i[is.na(s_loser_i) | !is.finite(s_loser_i) | s_loser_i < 0] <- 0
        s_loser_i[s_loser_i > 1] <- 1
        survivor_losers <- survivor_losers * s_loser_i
      }
      out_val <- pdf_winner * survivor_losers
    } else {
      out_val <- pdf_winner
    }
    out_val[is.na(out_val) | !is.finite(out_val) | out_val < 0] <- 0
    return(out_val)
  }
  
  order_pars_for_winner <- function(p_all_acc, k_idx) {
    if (is.na(k_idx) || k_idx < 1 || k_idx > nrow(p_all_acc)) {
      # This can happen if R_j_idx is NA (response unknown) and we are iterating for observed RT
      # For observed RT, winner must be known.
      stop(paste("Invalid k_idx in order_pars_for_winner: ", k_idx, " for unique trial ", j))
    }
    return(rbind(p_all_acc[k_idx, , drop=FALSE], p_all_acc[-k_idx, , drop=FALSE]))
  }
  
  integrate_for_kth_winner <- function(k_winner_idx, p_all_acc, low, upp, model) {
    if (low >= upp) return(0)
    if (is.na(k_winner_idx) || k_winner_idx < 1 || k_winner_idx > nrow(p_all_acc) ) return(0)
    
    pars_ordered <- order_pars_for_winner(p_all_acc, k_winner_idx)
    
    # Basic check on integrand before calling integrate
    # Test point can be tricky if low/upp are Inf.
    test_t <- if (is.finite(low) && is.finite(upp)) (low + upp) / 2 else if (is.finite(low)) low + 1 else if (is.finite(upp)) upp -1 else 1
    if (test_t < 0 && low == 0) test_t <- 1e-6 # Avoid negative if low is 0
    if (test_t < low && low >= 0) test_t <- low + 1e-6 # Ensure test_t is within bounds or sensible
    if (test_t > upp && upp >= 0) test_t <- upp - 1e-6
    
    integrand_val_at_test_t <- f_race_integrand(t = test_t, p_trial_this_winner_first = pars_ordered, model = model)
    if (is.na(integrand_val_at_test_t) || !is.finite(integrand_val_at_test_t) || integrand_val_at_test_t < 0) {
      # If integrand is already bad, no point integrating
      # This might happen due to bad parameters leading to NA/Inf in dfun/pfun
    }
    
    res <- suppressWarnings(try(stats::integrate(f_race_integrand,
                                                 lower = low, upper = upp,
                                                 p_trial_this_winner_first = pars_ordered,
                                                 model = model,
                                                 rel.tol = .Machine$double.eps^0.4
    )$value, silent = TRUE))
    if (inherits(res, "try-error") || is.na(res) || !is.finite(res) || res < 0) return(0)
    return(res)
  }
  
  get_trunc_corr_factor_for_kth_winner <- function(k_winner_idx, p_all_acc, model) {
    if (!(LT > 0 || UT < Inf)) return(1.0)
    if (is.na(k_winner_idx)) return(NA_real_)
    
    prob_untruncated <- integrate_for_kth_winner(k_winner_idx, p_all_acc, 0, Inf, model)
    prob_truncated_interval <- integrate_for_kth_winner(k_winner_idx, p_all_acc, LT, UT, model)
    
    if (prob_truncated_interval > 1e-12) {
      if(is.na(prob_untruncated) || prob_untruncated < 0) return(NA_real_)
      return(prob_untruncated / prob_truncated_interval)
    } else {
      if (prob_untruncated > 1e-12) return(NA_real_)
      return(1.0)
    }
  }
  
  for (j in 1:n_unique_trials) {
    current_trial_par_indices <- ((j - 1) * n_acc + 1):(j * n_acc)
    
    # Make sure parameters are valid
    if (!all(ok_params[current_trial_par_indices])) {
      ll_unique[j] <- min_ll
      next
    }
    
    pars_condition_j_all_acc <- pars[current_trial_par_indices, , drop = FALSE]
    
    rt_j <- dadm$rt[current_trial_par_indices[1]]
    R_j_idx <- as.integer(dadm$R[current_trial_par_indices[1]])
    
    current_ll_j <- min_ll
    prob_density_or_mass <- 0
    
    # If RT is observed
    if (is.finite(rt_j) && rt_j > 0) {
      if (!is.na(R_j_idx)) {
        if (rt_j < LT || rt_j > UT) {
          prob_density_or_mass <- 0
        } else {
          pars_ordered_obs <- order_pars_for_winner(pars_condition_j_all_acc, R_j_idx)
          prob_density_or_mass <- f_race_integrand(t = rt_j, p_trial_this_winner_first = pars_ordered_obs, model = model)
        }
        if (prob_density_or_mass > 0) {
          trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model)
          if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
          else prob_density_or_mass <- prob_density_or_mass * trunc_cf
        }
      } else { prob_density_or_mass <- 0 } # Observed RT but R_j_idx is NA: should not happen
      
    } else if (rt_j == -Inf) { # fast censoring
      if (!is.na(R_j_idx)) {
        prob_density_or_mass <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, LT, LC, model)
        if (prob_density_or_mass > 0) {
          trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model)
          if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
          else prob_density_or_mass <- prob_density_or_mass * trunc_cf
        }
      } else {
        current_sum_prob <- 0
        for (k_winner_loopvar in 1:n_acc) {
          p_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, LT, LC, model)
          if (p_k > 0) {
            trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model)
            if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k * trunc_cf_k)
          }
        }
        prob_density_or_mass <- current_sum_prob
      }
    } else if (rt_j == Inf) { # Slow censoring
      if (!is.na(R_j_idx)) {
        prob_density_or_mass <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, UC, UT, model)
        if (prob_density_or_mass > 0) {
          trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model)
          if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
          else prob_density_or_mass <- prob_density_or_mass * trunc_cf
        }
      } else {
        current_sum_prob <- 0
        for (k_winner_loopvar in 1:n_acc) {
          p_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, UC, UT, model)
          if (p_k > 0) {
            trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model)
            if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k * trunc_cf_k)
          }
        }
        prob_density_or_mass <- current_sum_prob
      }
    } else if (is.na(rt_j)) {
      if (!is.na(R_j_idx)) {
        prob_L <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, LT, LC, model)
        prob_U <- integrate_for_kth_winner(R_j_idx, pars_condition_j_all_acc, UC, UT, model)
        prob_density_or_mass <- prob_L + prob_U
        if (prob_density_or_mass > 0) {
          trunc_cf <- get_trunc_corr_factor_for_kth_winner(R_j_idx, pars_condition_j_all_acc, model)
          if (is.na(trunc_cf) || !is.finite(trunc_cf) || trunc_cf < 0) prob_density_or_mass <- 0
          else prob_density_or_mass <- prob_density_or_mass * trunc_cf
        }
      } else {
        current_sum_prob <- 0
        for (k_winner_loopvar in 1:n_acc) {
          prob_L_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, LT, LC, model)
          prob_U_k <- integrate_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, UC, UT, model)
          p_k_sum <- prob_L_k + prob_U_k
          if (p_k_sum > 0) {
            trunc_cf_k <- get_trunc_corr_factor_for_kth_winner(k_winner_loopvar, pars_condition_j_all_acc, model)
            if (!is.na(trunc_cf_k) && is.finite(trunc_cf_k) && trunc_cf_k >=0) current_sum_prob <- current_sum_prob + (p_k_sum * trunc_cf_k)
          }
        }
        prob_density_or_mass <- current_sum_prob
      }
    }
    
    prob_density_or_mass <- max(0, prob_density_or_mass)
    
    if (prob_density_or_mass > .Machine$double.eps) {
      current_ll_j <- log(prob_density_or_mass)
    } else {
      current_ll_j <- min_ll
    }
    
    ll_unique[j] <- max(min_ll, current_ll_j)
  }
  
  if (!is.null(dadm$N) && length(dadm$N) == n_unique_trials) {
    total_ll <- sum(ll_unique * dadm$N)
  } else if (!is.null(attr(dadm, "expand"))) {
    total_ll <- sum(ll_unique[attr(dadm, "expand")])
  } else {
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
