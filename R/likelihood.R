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
# - dadm$accumulator_role: A factor column identifying the role of each row.
#   Expected levels are "S1D", "S2D", "S1A", "S2A". The order of these levels
#   in the factor is not critical, as the code explicitly extracts by name.
# - dadm$observed_response: A column containing the actual observed response
#   ("yes" or "no") for the original trial, replicated across the 4 rows.
# - dadm$rt: The observed RT, replicated across the 4 rows.
# - pars: Matrix aligned with dadm, rows providing base parameters (e.g., "v", "b")
#   for the accumulator specified in dadm$accumulator_role.
log_likelihood_redundant_target_race <- function(pars, dadm, model, min_ll = log(1e-10))
{
  # --- Input Validations ---
  if (is.null(dadm$rt)) stop("dadm$rt is missing.")
  if (is.null(dadm$accumulator_role)) stop("dadm$accumulator_role (factor identifying S1D, S2D, S1A, S2A) is missing.")
  if (is.null(dadm$observed_response)) stop("dadm$observed_response ('yes'/'no') is missing.")
  if (is.null(attr(dadm, "expand"))) stop("attr(dadm, 'expand') is missing.")
  if (nrow(dadm) %% 4 != 0) {
    stop("nrow(dadm) must be a multiple of 4 (4 accumulator rows per original trial).")
  }

  # Expected accumulator roles for internal mapping
  internal_role_order <- c("S1D", "S2D", "S1A", "S2A")
  if (!all(internal_role_order %in% levels(as.factor(dadm$accumulator_role)))) {
    missing_roles <- internal_role_order[!internal_role_order %in% levels(as.factor(dadm$accumulator_role))]
    stop(paste("dadm$accumulator_role is missing expected levels:", paste(missing_roles, collapse=", "),
               ". Expected levels are S1D, S2D, S1A, S2A."))
  }

  # --- PDF and CDF Calculation (once for all rows) ---
  f_all <- model$dfun(dadm$rt, pars)
  F_all <- model$pfun(dadm$rt, pars)

  n_blocks <- nrow(dadm) / 4 # Number of original trials

  # --- Extract f/F values for each role into vectors of length n_blocks ---
  # This assumes that dadm is structured such that filtering by accumulator_role
  # results in vectors of length n_blocks, correctly ordered by original trial.
  # This typically means dadm should be sorted by an original trial ID,
  # then by accumulator_role (or this extraction order matches such a sort).

  f1 <- f_all[dadm$accumulator_role == "S1D"]; F1 <- F_all[dadm$accumulator_role == "S1D"]
  f2 <- f_all[dadm$accumulator_role == "S2D"]; F2 <- F_all[dadm$accumulator_role == "S2D"]
  f3 <- f_all[dadm$accumulator_role == "S1A"]; F3 <- F_all[dadm$accumulator_role == "S1A"]
  f4 <- f_all[dadm$accumulator_role == "S2A"]; F4 <- F_all[dadm$accumulator_role == "S2A"]

  # Validation: Check if all extracted vectors have the correct length (n_blocks)
  expected_len <- n_blocks
  if (any(sapply(list(f1, F1, f2, F2, f3, F3, f4, F4), length) != expected_len)) {
    stop("Length mismatch after extracting role-specific f/F values. Check dadm structure and accumulator_role factor.")
  }

  # --- Get Observed Responses (one per original trial/block) ---
  # Taking from the first row of each conceptual 4-row block.
  # This assumes dadm is grouped by original trial, so seq(1, nrow(dadm), by = 4) picks these first rows.
  # If not, a more robust way might be needed if an original_trial_id column exists.
  # For now, this is consistent with the loop version's row_indices[1].
  observed_responses_per_block <- dadm$observed_response[seq(1, nrow(dadm), by = 4)]
  if (length(observed_responses_per_block) != n_blocks) {
      stop("Could not correctly extract one observed_response per trial block.")
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

  # Handle any responses that were not 'yes' or 'no'
  unhandled_responses <- !(is_yes_response | is_no_response)
  if (any(unhandled_responses)) {
      original_indices_unhandled <- which(unhandled_responses)
      first_row_in_dadm_unhandled <- (original_indices_unhandled - 1) * 4 + 1
      problematic_values <- unique(dadm$observed_response[first_row_in_dadm_unhandled])
      warning(paste("Unhandled dadm$observed_response values found:", paste(problematic_values, collapse=", "),
                  ". Corresponding log-likelihoods will be NA, then min_ll."))
      ll_block_values[unhandled_responses] <- NA # Will become min_ll later
  }


  # --- Final Summation (same logic as before) ---
  ll_for_dadm_rows <- numeric(nrow(dadm))
  # Replicate the block's LL to its constituent rows.
  # rep(ll_block_values, each=4) works if dadm is perfectly sorted by block, then accumulator within block.
  # A loop is safer if that order isn't strictly guaranteed for attr(dadm,"expand")'s sake.
  # However, the extraction of f1-f4 and observed_responses_per_block *does* assume a certain usable order.
  # Given that, rep should be fine if the upstream extraction is robust.
  # For max safety in this step, retain loop or use a more robust replication:
  original_trial_idx_for_row <- rep(1:n_blocks, each = 4) # Assuming dadm is sorted by block
  ll_for_dadm_rows <- ll_block_values[original_trial_idx_for_row]

  final_ll_values <- pmax(min_ll, ll_for_dadm_rows[attr(dadm, "expand")])
  final_ll_values[is.na(final_ll_values)] <- min_ll

  return(sum(final_ll_values))
}
