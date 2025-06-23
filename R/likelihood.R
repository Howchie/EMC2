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
#
# Assumes dadm is structured with 4 rows per original trial, where each row
# corresponds to one of four conceptual accumulators.
# - dadm$accumulator_role: A factor column identifying the role of each row
#   (e.g., levels "S1D", "S2D", "S1A", "S2A").
# - dadm$observed_response: A column containing the actual observed response
#   ("yes" or "no") for the original trial, replicated across the 4 rows.
# - dadm$rt: The observed RT, replicated across the 4 rows.
# - pars: Matrix aligned with dadm, rows providing base parameters (e.g., "v", "b")
#   for the accumulator specified in dadm$accumulator_role.
#
# The order S1D, S2D, S1A, S2A is used internally for f1, f2, f3, f4.
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

  # --- PDF and CDF Calculation ---
  # f_all and F_all will have one value per row in dadm (i.e., per accumulator-trial)
  f_all <- model$dfun(dadm$rt, pars)
  F_all <- model$pfun(dadm$rt, pars)

  n_blocks <- nrow(dadm) / 4 # Number of original trials represented in dadm
  ll_block_values <- numeric(n_blocks) # Stores one LL value per original trial

  # --- Define Fixed Internal Mapping for Formulas ---
  # The code will map dadm$accumulator_role levels to these fixed slots:
  # Slot 1: S1D (Stimulus 1 Detect)
  # Slot 2: S2D (Stimulus 2 Detect)
  # Slot 3: S1A (Stimulus 1 Absent/Non-Detect)
  # Slot 4: S2A (Stimulus 2 Absent/Non-Detect)
  internal_role_order <- c("S1D", "S2D", "S1A", "S2A")

  # Check if all required roles are present in dadm$accumulator_role levels
  # This ensures the user's factor levels align with what the function expects.
  if (!all(internal_role_order %in% levels(as.factor(dadm$accumulator_role)))) {
    missing_roles <- internal_role_order[!internal_role_order %in% levels(as.factor(dadm$accumulator_role))]
    stop(paste("dadm$accumulator_role is missing expected levels:", paste(missing_roles, collapse=", "),
               ". Expected levels are S1D, S2D, S1A, S2A."))
  }

  f_vec <- numeric(4) # To store f values for S1D, S2D, S1A, S2A for the current block
  F_vec <- numeric(4) # To store F values for S1D, S2D, S1A, S2A for the current block

  # --- Loop Through Each Original Trial (Block of 4 Rows) ---
  for (i in 1:n_blocks) {
    row_indices <- ((i - 1) * 4 + 1):(i * 4) # Get row indices for the current block

    # Extract f, F, and accumulator roles for the current block
    f_block_unsorted <- f_all[row_indices]
    F_block_unsorted <- F_all[row_indices]
    # Roles for the current block, ensuring it's character for matching
    roles_in_block_char <- as.character(dadm$accumulator_role[row_indices])

    # Populate f_vec and F_vec according to the fixed internal_role_order
    for (k in 1:4) {
      current_role_to_find <- internal_role_order[k]
      idx_in_block <- match(current_role_to_find, roles_in_block_char)

      if (is.na(idx_in_block)) {
        stop(paste0("Expected accumulator role '", current_role_to_find,
                    "' not found in block ", i, ". Roles found: ",
                    paste(roles_in_block_char, collapse=", ")))
      }
      f_vec[k] <- f_block_unsorted[idx_in_block]
      F_vec[k] <- F_block_unsorted[idx_in_block]
    }

    # Assign to f1-f4 and F1-F4 based on the fixed internal order
    f1 <- f_vec[1]; F1 <- F_vec[1] # S1D
    f2 <- f_vec[2]; F2 <- F_vec[2] # S2D
    f3 <- f_vec[3]; F3 <- F_vec[3] # S1A
    f4 <- f_vec[4]; F4 <- F_vec[4] # S2A

    # Get the observed response for this original trial (replicated across the 4 rows)
    observed_response_trial <- dadm$observed_response[row_indices[1]]

    # --- Apply Likelihood Formula Based on Observed Response ---
    if (observed_response_trial == "yes") {
      # Formula: (fS1D*(1-FS2D) + fS2D*(1-FS1D)) * (1 - FS1A*FS2A)
      # This is P(first of S1D,S2D finishes at t) * P(NOT(S1A finishes by t AND S2A finishes by t))
      term1_yes <- f1 * (1 - F2) + f2 * (1 - F1)
      term2_yes <- 1 - (F3 * F4)

      # Numerical stability: ensure terms are not negative before log
      term1_yes <- pmax(term1_yes, .Machine$double.eps)
      term2_yes <- pmax(term2_yes, .Machine$double.eps)

      ll_block_values[i] <- log(term1_yes) + log(term2_yes)
    } else if (observed_response_trial == "no") {
      # Formula: (fS1A*FS2A + fS2A*FS1A) * (1-FS1D) * (1-FS2D)
      # This is P(Max(S1A,S2A) finishes at t) * P(S1D not finished by t) * P(S2D not finished by t)
      term1_no <- f3 * F4 + f4 * F3
      term2_no <- (1 - F1)
      term3_no <- (1 - F2)

      # Numerical stability
      term1_no <- pmax(term1_no, .Machine$double.eps)
      term2_no <- pmax(term2_no, .Machine$double.eps)
      term3_no <- pmax(term3_no, .Machine$double.eps)

      ll_block_values[i] <- log(term1_no) + log(term2_no) + log(term3_no)
    } else {
      stop(paste("Invalid dadm$observed_response value:'", observed_response_trial,
                 "' in block ", i, ". Must be 'yes' or 'no'."))
    }
  } # End of loop through blocks

  # --- Final Summation ---
  # ll_block_values has one LL per original trial (block).
  # Replicate these values for each constituent row in dadm to use attr(dadm, "expand").
  ll_for_dadm_rows <- numeric(nrow(dadm))
  for (i in 1:n_blocks) {
    row_indices <- ((i - 1) * 4 + 1):(i * 4)
    ll_for_dadm_rows[row_indices] <- ll_block_values[i]
  }

  final_ll_values <- pmax(min_ll, ll_for_dadm_rows[attr(dadm, "expand")])
  # Ensure any NAs that might arise (e.g. if a log term was NA despite pmax) become min_ll
  final_ll_values[is.na(final_ll_values)] <- min_ll

  return(sum(final_ll_values))
}
