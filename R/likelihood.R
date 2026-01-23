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

# ZH Race likelihood with censoring and truncation
# Helpers ---------------------------------------------------------------------

## create a string key that uniquely identifies the parameter set + bounds
.make_key <- function(lower, upper, pars_mat)
  paste0(
    formatC(lower, digits = 12, width = 0, format = "fg"), "_",
    formatC(upper, digits = 12, width = 0, format = "fg"), "_",
    paste(formatC(as.numeric(pars_mat), digits = 12, width = 0, format = "fg"),
          collapse = "_")
  )

.log_diff_exp <- function(a, b) {
  if (!is.finite(a)) return(-Inf)
  if (!is.finite(b)) {
    if (is.infinite(b) && b < 0) return(a)
    return(-Inf)
  }
  if (a <= b) return(-Inf)
  m <- max(a, b)
  m + log(exp(a - m) - exp(b - m))
}

.log_survivor_prod <- function(t, pars_mat, model) {
  if (nrow(pars_mat) == 0) return(0)
  cdf <- model$pfun(rt = t, pars = pars_mat)
  if (any(!is.finite(cdf))) return(-Inf)
  cdf[cdf <= 0] <- 0
  if (any(cdf >= 1, na.rm = TRUE)) return(-Inf)
  cdf <- pmin(cdf, 1 - 1e-15)
  sum(log1p(-cdf))
}

.prob_min_in_interval <- function(lower, upper, pars_mat, model) {
  if (lower>=upper) return(0)
  logS_lower <- .log_survivor_prod(lower, pars_mat, model)
  if (!is.finite(logS_lower)) return(.Machine$double.eps)
  if (is.infinite(upper)) {
    return(max(exp(logS_lower), .Machine$double.eps))
  }
  logS_upper <- .log_survivor_prod(upper, pars_mat, model)
  log_prob <- .log_diff_exp(logS_lower, logS_upper)
  if (!is.finite(log_prob)) return(.Machine$double.eps)
  max(exp(log_prob), .Machine$double.eps)
}

## integrate the joint density for the *k-th winner* on (lower, upper)
.integrate_kth_winner <- function(k, pars_mat, lower, upper, model,
                                  rel.tol = 1e-7, ...) {
  if (lower>=upper) return(0)
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
  if (lower >= upper) return(.Machine$double.eps)
  # Prefer analytic survivor-product wherever possible (matches particle_ll.cpp)
  p_tot <- .prob_min_in_interval(lower, upper, pars_mat, model)
  if (is.finite(p_tot) && p_tot > .Machine$double.eps) {
    return(p_tot)
  }

  # Fallback to numerical integration if survivor path failed
  n_acc <- nrow(pars_mat)
  p_tot_int <- 0
  for (k in 1:n_acc) {
    p_tot_int <- p_tot_int +
      .integrate_kth_winner(k, pars_mat, lower, upper, model,
                            rel.tol = rel.tol)
  }
  max(p_tot_int, .Machine$double.eps)
}


# Main likelihood function -----------------------------------------------------
log_likelihood_race_cens_trunc <- function(pars,dadm,model,min_ll=log(1e-10)) {
  cache_env <- new.env(parent = emptyenv())
  posdrift = ifelse(model$c_name=="LBAIO",FALSE,TRUE)
  gng = any(grepl("GNG", model$c_name))
  has_pC <- !is.null(colnames(pars)) && ("pContaminant" %in% colnames(pars) && !(all(pars[,"pContaminant"]==0)))
  log_sum_exp2 <- function(a, b) {
    if (is.infinite(a) && a < 0) return(b)
    if (is.infinite(b) && b < 0) return(a)
    m <- max(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  ## basic dimensions
  n_trials <- nrow(dadm)
  n_acc    <- length(levels(dadm$R))
  stopifnot(n_trials %% n_acc == 0L)          # one row per acc. per trial
  n_unique <- n_trials / n_acc
  dadm$RACE_num <- if("RACE" %in% names(dadm))
    as.numeric(as.character(dadm$RACE)) else
    rep(n_acc, n_trials)
  if (any(names(dadm)=="RACE")){# Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA
  }
  #dadm = dadm[dadm$RACE_num<=as.numeric(dadm$lR),] # drop missing accumulator rows
  if (is.null(attr(pars,"ok"))){
    ok <- !logical(dim(pars)[1])
  } else ok <- attr(pars,"ok")

  ## convenience column refs
  RT <- dadm$rt;  LT <- dadm$LT;  UT <- dadm$UT
  LC <- dadm$LC;  UC <- dadm$UC;  R_idx <- dadm$R

  ## initialise per-row log densities
  lds <- rep(NA, n_trials)
  ll_unique <- rep(NA,n_unique)
  ## batch: finite RT, in-bounds, known winner
  finite_mask <- is.finite(RT) & RT > 0 &
    RT >= LT & RT <= UT & !is.na(R_idx) & (as.numeric(dadm$lR) <= dadm$RACE_num)
  trunc_cens_mask = ( (!is.finite(RT))) & (as.numeric(dadm$lR) <= dadm$RACE_num)
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
        ll_j <- sum(lds[idx])
        ll_unique[j] <- ll_j # max(min_ll, ll_j)
        if (!(LT[idx[1]] == 0 && UT[idx[1]] == Inf)){
          key <- .make_key(LT[idx[1]], UT[idx[1]], pars[idx, , drop = FALSE])
          invZ <- cache_env[[key]]
          if (is.null(invZ)) {
            invZ <- 1/.truncation_normaliser(pars[idx, , drop = FALSE],
                                               LT[idx[1]], UT[idx[1]], model,
                                               rel.tol = 1e-7)
            cache_env[[key]] <- invZ
          }
          ll_unique[j] <- ll_unique[j] + log(invZ)
        }
      }
    }
  }

  ## process other trials one-by-one (-Inf, +Inf, NA)
  prob_win <- function(k, lo, hi) {
    .integrate_kth_winner(k, pars[idx, , drop = FALSE],
                          lo, hi, model, rel.tol = 1e-7)
  }

  for (j in seq_len(n_unique)) {
    idx  <- ((j - 1L) * n_acc + 1L):(j * n_acc)
    winner = dadm$winner[idx]
    if (!trunc_cens_mask[idx[[1]]]) {next}
    rt   <- RT[idx[1]]          # all rows of a trial share RT
    Rj   <- as.numeric(R_idx[idx[1]])
    pval <- 0

    n_acc_j = dadm$RACE_num[idx[1]]
    idx_j=idx[1:n_acc_j]
    if(is.na(Rj)){ks = 1:n_acc_j}else{ks = Rj}
    if (identical(rt, -Inf)) {                       # fast censor
      lo <- LT[idx_j[1]]; hi <- LC[idx_j[1]]
      if (hi==0) {stop("LC must be non-zero if rt==-Inf")}
      pval <- sum(vapply(ks, prob_win, numeric(1), lo, hi))
    } else if (identical(rt, Inf)) {                 # slow censor
      if (gng) {
        lo <- 0; hi <- UC[idx_j[1]]
        k_nogo = which(winner)
        termA = prob_win(k_nogo,lo,hi)
        termB=1
        for (k in 1:n_acc_j) {
          termB = termB*(1 - model$pfun(rt = hi, pars = pars[idx[k], , drop = FALSE]))
        }
        pval = termA+termB
      } else {
        lo <- UC[idx_j[1]]; hi <- UT[idx_j[1]]
        if (lo>=hi) pval <- 0 else {
          if (length(idx_j)==1) {
            pval = 1 - (model$pfun(lo,pars[idx, , drop = FALSE])) # if a single acccumulator, design omissions are just 1-F(t)
          } else if (is.na(Rj)) {
            # Unknown winner: use analytic survivor product wherever possible
            prob_from_surv <- .prob_min_in_interval(lo, hi, pars[idx_j, , drop = FALSE], model)
            if (is.finite(prob_from_surv) && prob_from_surv > .Machine$double.eps) {
              pval <- prob_from_surv
            } else {
              pval <- sum(vapply(ks, prob_win, numeric(1), lo, hi))
            }
          } else {
            pval <- sum(vapply(ks, prob_win, numeric(1), lo, hi))
          }
        }
      }
    } else if (is.na(rt)) {
      # missing RT
      lo1 <- LT[idx_j[1]]; hi1 <- LC[idx_j[1]]
      lo2 <- UC[idx_j[1]]; hi2 <- UT[idx_j[1]]
      if (lo2==Inf) {stop("UC must be finite if rt==NA")}
      if (is.na(Rj)) {
        # Unknown winner: compute via survivor products, fall back to integrals if needed
        pval_lower <- .prob_min_in_interval(lo1, hi1, pars[idx_j, , drop = FALSE], model)
        pval_upper <- .prob_min_in_interval(lo2, hi2, pars[idx_j, , drop = FALSE], model)
        pval <- pval_lower + pval_upper
        if (!is.finite(pval) || pval <= .Machine$double.eps) {
          pval <- sum(vapply(ks, function(k)
            prob_win(k, lo1, hi1) + prob_win(k, lo2, hi2), numeric(1)))
        }
      } else {
        pval <- sum(vapply(ks, function(k)
          prob_win(k, lo1, hi1) + prob_win(k, lo2, hi2), numeric(1)))
      }
    }                                               # 0 or negative RT prob 0

    ## truncation correction (if RT unobserved)
    if (!(LT[idx_j[1]] == 0 && UT[idx_j[1]] == Inf) && pval > 0) {
      key  <- .make_key(LT[idx_j[1]], UT[idx_j[1]], pars[idx, , drop = FALSE])
      invZ <- cache_env[[key]]
      if (is.null(invZ)) {
        invZ <- 1 / .truncation_normaliser(pars[idx, , drop = FALSE],
                                           LT[idx_j[1]], UT[idx_j[1]], model,
                                           rel.tol = 1e-5, posdrift = posdrift)
        cache_env[[key]] <- invZ
      }
      pval <- pval * invZ
    }


    ll_unique[j] <- log(pval) # if (pval > .Machine$double.eps) log(pval) else min_ll
  }

  ## Mixture for contaminant omissions (pContaminant): with prob pC the trial is a contaminant omission,
  ## otherwise it follows the model likelihood. For finite RTs, contaminant omissions are impossible.
  if (has_pC) {
    pc_by_trial <- pars[seq(1L, n_trials, by = n_acc), "pContaminant"]
    pc_by_trial[is.na(pc_by_trial)] <- 0
    for (j in seq_len(n_unique)) {
      if (is.na(ll_unique[j])) next
      pC <- pc_by_trial[j]
      if (!is.finite(pC)) pC <- 0
      pC <- max(0, min(1, pC))
      log1m_pC <- if (pC < 1) log1p(-pC) else -Inf
      rt <- RT[(j - 1L) * n_acc + 1L]
      if (rt!=Inf) {
        ll_unique[j] <- ll_unique[j] + log1m_pC
      } else {
        ll_unique[j] <- log_sum_exp2(log(pC), log1m_pC + ll_unique[j])
      }
    }
  }
  ## aggregate over copies (expand) or N column
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
    isna <- is.na(dadm$rt)|is.infinite(dadm$rt) # allow Inf like other censored models
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
