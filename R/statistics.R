#' Information Criteria and Marginal Likelihoods
#'
#' Returns the BPIC/DIC or marginal deviance (-2*marginal likelihood) for a list of samples objects.
#'
#' @param sList List of samples objects
#' @param stage A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, defaults to `TRUE`, uses the minimal or mean likelihood (whichever is better) in the
#' calculation, otherwise always uses the mean likelihood.
#' @param BayesFactor Boolean, defaults to `TRUE`. Include marginal likelihoods as estimated using WARP-III bridge sampling.
#' Usually takes a minute per model added to calculate
#' @param WAIC Boolean, defaults to `TRUE`. Include WAIC (Watanabe-Akaike Information Criterion) computed
#' via the `loo` package. Requires computing per-trial log-likelihoods across all posterior samples.
#' @param LOO Boolean, defaults to `FALSE`. Include PSIS-LOO computed via the `loo` package.
#' Uses the same pointwise log-likelihood matrix as WAIC and is typically slower.
#' @param pointwise Character string, one of `"trial"` or `"subject"`, controlling
#' the unit of prediction for WAIC/LOO in pooled hierarchical comparisons.
#' `"trial"` (default) treats each trial as one pointwise unit. `"subject"` computes
#' the hierarchical leave-one-subject-out LOO by marginalising subject parameters:
#' for each posterior draw of the group-level parameters (theta_mu, theta_var),
#' `K` proposals are drawn from the group distribution and the marginal subject
#' log-likelihood is estimated via log-mean-exp. This gives well-conditioned PSIS
#' importance weights where summing conditional per-trial LLs would give Pareto-k >> 1.
#' Ignored unless `WAIC` or `LOO` is `TRUE`.
#' @param K Integer (default 200). Number of importance samples drawn from the group
#' prior per posterior iteration when `pointwise = "subject"`. Higher values give a
#' more accurate marginal likelihood estimate at proportionally higher cost.
#' @param cores_for_loo Integer (default 1). Cores for WAIC/LOO computation. Used in two places:
#' (1) parallelises the per-subject log-likelihood computation across subjects, and
#' (2) passed to \code{loo::loo()} as its \code{cores} argument for internal PSIS parallelism.
#' Only relevant when \code{WAIC} or \code{LOO} is \code{TRUE} and the pw_ll is not already cached.
#' @param cores_for_props Integer, how many cores to use for the Bayes factor calculation, here 4 is the default for the 4 different proposal densities to evaluate, only 1, 2 and 4 are sensible.
#' @param cores_per_prop Integer, how many cores to use for the Bayes factor calculation if you have more than 4 cores available. Cores used will be cores_for_props * cores_per_prop. Best to prioritize cores_for_props being 4 or 2
#' @param print_summary Boolean (default `TRUE`), print table of results
#' @param digits Integer, significant digits in printed table for information criteria
#' @param digits_p Integer, significant digits in printed table for model weights
#' @param ... Additional, optional arguments
#'
#' @return Matrix of effective number of parameters, mean deviance, deviance of
#' mean, DIC, BPIC, Marginal Deviance (if `BayesFactor=TRUE`) and associated weights.
#' @examples \donttest{
#' compare(list(samples_LNR), cores_for_props = 1)
#' # Typically we would define a list of two (or more) different models:
#' # # Here the full model is an emc object with the hypothesized effect
#' # # The null model is an emc object without the hypothesized effect
#' # design_full <- design(data = forstmann,model=DDM,
#' #                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' #                            constants=c(s=log(1)))
#' # # Now without a ~ E
#' # design_null <- design(data = forstmann,model=DDM,
#' #                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#' #                            constants=c(s=log(1)))
#' #
#' # full_model <- make_emc(forstmann, design_full)
#' # full_model <- fit(full_model)
#' #
#' # null_model <- make_emc(forstmann, design_null)
#' # null_model <- fit(null_model)
#' # sList <- list(full_model, null_model)
#' # # By default emc uses 4 cores to parallelize marginal likelihood estimation across proposals
#' # # So cores_per_prop = 3 results in 12 cores used.
#' # compare(sList, cores_per_prop = 3)
#' }
#' @export

compare <- function(sList,stage="sample",filter=NULL,use_best_fit=TRUE,
                        BayesFactor = TRUE, WAIC = TRUE, LOO = FALSE, pointwise = c("trial", "subject"),
                        K = 200, cores_for_loo = 1, cores_for_props = 4, cores_per_prop = 1,
                        print_summary=TRUE,digits=0,digits_p=3, ...) {
  if(is(sList, "emc")) sList <- list(sList)
  pointwise <- match.arg(pointwise)
  getp <- function(IC) {
    IC <- -(IC - min(IC))/2
    exp(IC)/sum(exp(IC))
  }
  if (is.numeric(filter)) defaultsf <- filter[1] else defaultsf <- 0
  sflist <- as.list(setNames(rep(defaultsf,length(sList)),names(sList)))
  if (is.list(filter)) for (i in names(filter))
    if (i %in% names(sflist)) sflist[[i]] <- filter[[i]]
  dots <- add_defaults(list(...), group_only = FALSE)
  ICs <- setNames(vector(mode="list",length=length(sList)),names(sList))
  for (i in 1:length(ICs)) ICs[[i]] <- IC(sList[[i]],stage=stage,
                                          filter=sflist[[i]],use_best_fit=use_best_fit,subject=list(...)$subject,print_summary=FALSE,
                                          group_only = dots$group_only)
  ICs <- data.frame(do.call(rbind,ICs))
  DICp <- getp(ICs$DIC)
  BPICp <- getp(ICs$BPIC)
  out <- cbind.data.frame(DIC=ICs$DIC,wDIC=DICp,BPIC=ICs$BPIC,wBPIC=BPICp,ICs[,-c(1:2)])

  if(WAIC || LOO){
    subj <- list(...)$subject
    WAICs <- if(WAIC) rep(NA_real_, length(sList)) else NULL
    LOOs <- if(LOO) rep(NA_real_, length(sList)) else NULL
    pw_ll_list <- list()
    for(i in seq_along(sList)){
      tryCatch({
        ll_mat <- if(!is.null(subj)) {
          .ll_matrix_subject(sList[[i]], stage=stage, filter=sflist[[i]], subject=subj)
        } else if (pointwise == "trial") {
          .ll_matrix_pooled(sList[[i]], stage=stage, filter=sflist[[i]], cores=cores_for_loo)
        } else {
          .marg_ll_matrix(sList[[i]], stage=stage, filter=sflist[[i]], K=K)
        }
        pw_ll_list[[i]] <- ll_mat
        if(WAIC) WAICs[i] <- waic_from_ll(ll_mat)
        if(LOO) LOOs[i] <- loo_from_ll(ll_mat, cores=cores_for_loo)
      }, error = function(e) {
        if(WAIC) warning("WAIC computation failed for model ", i, ": ", conditionMessage(e), call. = FALSE)
        if(LOO) warning("LOO computation failed for model ", i, ": ", conditionMessage(e), call. = FALSE)
      })
    }
    names(pw_ll_list) <- names(sList)
    if(WAIC){
      if(!all(is.na(WAICs))){
        WAICp <- getp(WAICs)
        out <- cbind.data.frame(WAIC=WAICs, wWAIC=WAICp, out)
      } else {
        WAIC <- FALSE
      }
    }
    if(LOO){
      if(!all(is.na(LOOs))){
        LOOp <- getp(LOOs)
        out <- cbind.data.frame(LOO=LOOs, wLOO=LOOp, out)
      } else {
        LOO <- FALSE
      }
    }
    attr(out, "pw_ll") <- pw_ll_list
  }

  if(BayesFactor){
    MLLs <- numeric(length(sList))
    for(i in 1:length(MLLs)){
      MLLs[i] <- run_bridge_sampling(sList[[i]], stage = stage, filter = sflist[[i]], both_splits = FALSE,
                                     cores_for_props = cores_for_props, cores_per_prop = cores_per_prop)
    }
    MD <- -2*MLLs
    modelProbability <- getp(MD)
    out <- cbind.data.frame(MD = MD, wMD = modelProbability, out)
  }
  if (print_summary) {
    tmp <- out
    weight_cols <- c("wDIC","wBPIC")
    if(WAIC) weight_cols <- c(weight_cols, "wWAIC")
    if(LOO) weight_cols <- c(weight_cols, "wLOO")
    if(BayesFactor) weight_cols <- c(weight_cols, "wMD")
    w_idx <- which(names(tmp) %in% weight_cols)
    for(wc in weight_cols) if(wc %in% names(tmp)) tmp[[wc]] <- round(tmp[[wc]], digits_p)
    if (length(w_idx) > 0) {
      tmp[,-w_idx] <- round(tmp[,-w_idx], digits=digits)
    } else {
      tmp[,] <- round(tmp[,], digits=digits)
    }
    print(tmp)
  }
  invisible(out)
}

std_error_IS2 <- function(IS_samples, n_bootstrap = 50000){
  log_marglik_boot= array(dim = n_bootstrap)
  for (i in 1:n_bootstrap){
    log_weight_boot = sample(IS_samples, length(IS_samples), replace = TRUE) #resample with replacement from the lw
    log_marglik_boot[i] <- median(log_weight_boot)
  }
  return(sd(log_marglik_boot))
}

#' Extract Pointwise Log-Likelihoods
#'
#' Extracts the pointwise log-likelihoods for each trial in the data and
#' appends them to the Data Augmented Design Matrix (DADM).
#'
#' @param emc An emc object.
#' @param stage A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`.
#' @param filter An integer or vector. iterations to exclude.
#' @param FUN A function to aggregate across iterations (e.g., \code{mean}, \code{median}, or \code{NULL} to return all iterations).
#'
#' @return A data frame (concatenated DADMs) with an additional `ll` column (or multiple columns if \code{FUN=NULL}).
#' @export
extract_pw_ll <- function(emc, stage = "sample", filter = 0, FUN = mean) {
  # Get the pointwise log-likelihood matrix [n_iter x total_trials]
  ll_mat <- .ll_matrix_pooled(emc, stage = stage, filter = filter)

  # Aggregation
  if (!is.null(FUN)) {
    ll_vals <- apply(ll_mat, 2, FUN)
  } else {
    ll_vals <- t(ll_mat)
    colnames(ll_vals) <- paste0("ll_iter", 1:nrow(ll_mat))
  }

  # Get the data
  # We use the standard get_data function to get trial-wise data
  all_data <- get_data(emc)
  # For joint models, get_data returns a list. Concatenate them.
  if (is.list(all_data) && !is.data.frame(all_data)) {
    all_data <- do.call(rbind, all_data)
  }

  if (nrow(all_data) != (if(is.null(FUN)) nrow(ll_vals) else length(ll_vals))) {
     stop("DADM rows do not match pointwise log-likelihood columns. Check if data was expanded or filtered.")
  }

  if (is.null(FUN)) {
    out <- cbind(all_data, ll_vals)
  } else {
    all_data$ll <- ll_vals
    out <- all_data
  }

  return(out)
}

#' Add Pointwise Log-Likelihoods to an emc Object
#'
#' Calculates the pointwise log-likelihoods for all samples in the emc object
#' and stores them in the \code{samples} list of each chain. This allows
#' for faster subsequent calls to \code{compare()} or \code{extract_pw_ll()}.
#'
#' @param emc An emc object.
#' @param cores_for_chains Integer. Number of cores to use for parallel calculation across chains.
#'
#' @return The same emc object with \code{pw_ll} added to each chain's samples.
#' @export
add_pw_ll <- function(emc, cores_for_chains = 1) {
  emc <- restore_duplicates(emc)
  results <- auto_mclapply(emc, add_pw_ll_chain, mc.cores = cores_for_chains)
  for (i in seq_along(results)) {
    if (inherits(results[[i]], "try-error")) {
      cond <- attr(results[[i]], "condition")
      msg <- if (!is.null(cond)) conditionMessage(cond) else as.character(results[[i]])
      warning("add_pw_ll: chain ", i, " failed (", msg, "); pw_ll not added for this chain.")
    } else {
      emc[[i]] <- results[[i]]
    }
  }
  class(emc) <- "emc"
  return(emc)
}

add_pw_ll_chain <- function(chain) {
  alpha <- chain$samples$alpha
  if (is.null(alpha) || length(dim(alpha)) < 3 || dim(alpha)[3] == 0) return(chain)
  n_subjects <- dim(alpha)[2]
  n_iter <- dim(alpha)[3]
  data <- chain$data
  model <- chain$model
  type <- if (is.null(chain$type)) "" else chain$type

  sub_pw_lls <- lapply(1:n_subjects, function(s) {
    sub <- chain$subjects[s]
    # alpha[, s, ] is matrix [n_pars x n_iter]
    proposals <- t(alpha[, s, ])
    ll_mat <- calc_ll_pw(proposals, data[[sub]], model)
    return(t(ll_mat))
  })
  chain$samples$pw_ll <- do.call(rbind, sub_pw_lls)
  return(chain)
}


# robust_diwish <- function (W, v, S) { #RJI_change: this function is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
#   if (!is.matrix(S)) S <- matrix(S)
#   if (!is.matrix(W)) W <- matrix(W)
#   p <- nrow(S)
#   gammapart <- sum(lgamma((v + 1 - 1:p)/2))
#   ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
#   if (corpcor::is.positive.definite(W, tol=1e-8)){
#     cholW<-base::chol(W)
#   }else{
#     return(1e-10)
#   }
#   if (corpcor::is.positive.definite(S, tol=1e-8)){
#     cholS <- base::chol(S)
#   }else{
#     return(1e-10)
#   }
#   halflogdetS <- sum(log(diag(cholS)))
#   halflogdetW <- sum(log(diag(cholW)))
#   invW <- chol2inv(cholW)
#   exptrace <- sum(S * invW)
#   lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
#   lpdf <- lnum - ldenom
#   out <- exp(lpdf)
#   if(!is.finite(out)) return(1e-100)
#   if(out < 1e-10) return(1e-100)
#   return(exp(lpdf))
# }

robust_diwish <- function (W, v, S) { #RJI_change: this function is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
  if (!is.matrix(S)) S <- matrix(S)
  if (!is.matrix(W)) W <- matrix(W)
  p <- nrow(S)
  gammapart <- sum(lgamma((v + 1 - 1:p)/2))
  ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
  cholW <- base::chol(nearPD(W)$mat)
  cholS <- base::chol(nearPD(S)$mat)
  halflogdetS <- sum(log(diag(cholS)))
  halflogdetW <- sum(log(diag(cholW)))
  invW <- chol2inv(cholW)
  exptrace <- sum(S * invW)
  lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
  lpdf <- lnum - ldenom
  out <- exp(lpdf)
  if(!is.finite(out)) return(1e-100)
  return(out)
}

dhalft <- function (x, scale = 25, nu = 1, log = FALSE)
{
  x <- as.vector(x)
  scale <- as.vector(scale)
  nu <- as.vector(nu)
  if (any(scale <= 0))
    stop("The scale parameter must be positive.")
  NN <- max(length(x), length(scale), length(nu))
  x <- rep(x, len = NN)
  scale <- rep(scale, len = NN)
  nu <- rep(nu, len = NN)
  dens <- log(2) - log(scale) + lgamma((nu + 1)/2) - lgamma(nu/2) -
    0.5 * log(pi * nu) - (nu + 1)/2 * log(1 + (1/nu) * (x/scale) *
                                            (x/scale))
  if (log == FALSE)
    dens <- exp(dens)
  return(dens)
}

rwish <- function(v, S){
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message = "S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * (p - 1)/2)
  }
  return(crossprod(Z %*% CC))
}


riwish <- function(v, S){
  return(solve(rwish(v, solve(S))))
}

logdinvGamma <- function(x, shape, rate){
  dgamma(1/x, shape, rate, log = TRUE) - 2 * log(x)
}


# #' Calculate information criteria (DIC, BPIC), effective number of parameters and
# #' constituent posterior deviance (D) summaries (meanD = mean of D, Dmean = D
# #' for mean of posterior parameters and minD = minimum of D).
# #'
# #' @param emc emc object or list of these
# #' @param stage A string. Specifies which stage you want to plot.
# #' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
# #' If a vector is supplied, only the iterations in the vector will be considered.
# #' @param use_best_fit Boolean, default TRUE use best of minD and Dmean in
# #' calculation otherwise always use Dmean
# #' @param print_summary Boolean (default TRUE) print table of results
# #' @param digits Integer, significant digits in printed table
# #' @param subject Integer or string selecting a single subject, default NULL
# #' returns sums over all subjects
# #' @param group_only Boolean. If `TRUE` will calculate the IC for the group-level only
# #'
# #' @return Table of DIC, BPIC, EffectiveN, meanD, Dmean, and minD

IC <- function(emc,stage="sample",filter=0,use_best_fit=TRUE,
               print_summary=TRUE,digits=0,subject=NULL,
               group_only = FALSE)
  # Gets DIC, BPIC, effective parameters, mean deviance, and deviance of mean
{
  # Mean log-likelihood for each subject
  if (length(subject)!=1) {
    ll <- get_pars(emc, stage = stage, filter = filter, selection = "LL", merge_chains = TRUE)
    alpha <- get_pars(emc,selection="alpha",stage=stage,filter=filter, by_subject = TRUE, merge_chains = TRUE)
  } else {
    ll <- get_pars(emc, stage = stage, filter = filter, selection = "LL", merge_chains = TRUE,subject=subject)
    alpha <- get_pars(emc,selection="alpha",stage=stage,filter=filter, by_subject = TRUE, merge_chains = TRUE,subject=subject)
  }
  # ZH optimization: Replace col-wise apply loop with vectorized C-level primitive
  ll_mat <- ll[[1]][[1]]
  minDs <- -2*ll_mat[cbind(max.col(t(ll_mat), ties.method="first"), 1:ncol(ll_mat))]
  mean_lls <- colMeans(ll_mat)
  mean_pars <- lapply(alpha, function(x) {
    rowSums(vapply(x, colSums, numeric(ncol(x[[1]])))) / sum(vapply(x, nrow, numeric(1)))
  })
  # log-likelihood for each subject using their mean parameter vector
  data <- emc[[1]]$data
  mean_pars_lls <- setNames(numeric(length(mean_pars)),names(mean_pars))
  for (sub in names(mean_pars)){
    mean_pars_lls[sub] <- calc_ll_manager(t(mean_pars[[sub]]),dadm = data[[sub]], emc[[1]]$model)
  }
  Dmeans <- -2*mean_pars_lls

  if (!is.null(subject)) {
    Dmeans <- Dmeans[subject[1]]
    mean_lls <- mean_lls[subject[1]]
    minDs <- minDs[subject[1]]
  } else{
    group_stats <- group_IC(emc, stage=stage,filter=filter, type = emc[[1]]$type)
    if(group_only){
      mean_lls <- group_stats$mean_ll
      minDs <- group_stats$minD
      Dmeans <- group_stats$Dmean
    } else{
      mean_lls <- c(mean_lls, group_stats$mean_ll)
      minDs <- c(minDs, group_stats$minD)
      Dmeans <- c(Dmeans, group_stats$Dmean)
    }
  }
  if (use_best_fit) minDs <- pmin(minDs,Dmeans)

  # mean deviance(-2*ll of all data)
  mD <- sum(-2 * mean_lls)
  # Deviance of mean
  Dmean <- sum(Dmeans)
  # mimimum Deviance
  minD <- sum(minDs)

  # Use deviance of mean as best fit or use actual best fit
  if (!use_best_fit) Dm <- Dmean else Dm <- minD

  # effective number of parameters
  pD <- mD - Dm
  # DIC = mean deviance + effective number of parameters
  DIC <- mD + pD
  # BPIC = mean deviance + 2*effective number of parameters
  # Note this is the "easy" BPIC, instead of the complex 2007 one
  BPIC <- mD + 2*pD
  out <- c(DIC = DIC, BPIC = BPIC, EffectiveN = pD,meanD=mD,Dmean=Dmean,minD=minD)
  names(out) <- c("DIC","BPIC","EffectiveN","meanD","Dmean","minD")
  if (print_summary) print(round(out,digits))
  invisible(out)
}


#' Bayes Factors
#'
#' returns the Bayes Factor for two models
#'
#' @param MLL1 Numeric. Marginal likelihood of model 1. Obtained with `run_bridge_sampling()`
#' @param MLL2 Numeric. Marginal likelihood of model 2. Obtained with `run_bridge_sampling()`
#'
#' @return The BayesFactor for model 1 over model 2
#' @examples \donttest{
#' # Normally one would compare two different models
#' # Here we use two times the same model:
#' M1 <- M0 <- run_bridge_sampling(samples_LNR, both_splits = FALSE, cores_for_props = 1)
#' get_BayesFactor(M1, M0)
#' }
#' @export
get_BayesFactor <- function(MLL1, MLL2){
  exp(MLL1 - MLL2)
}

#' Information Criteria For Each Participant
#'
#' Returns the BPIC/DIC based model weights for each participant in a list of samples objects
#'
#' @param sList List of samples objects
#' @param stage A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, defaults to `TRUE`, use minimal likelihood or mean likelihood
#' (whichever is better) in the calculation, otherwise always uses the mean likelihood.
#' @param print_summary Boolean (defaults to `TRUE`) print tables of model weight results
#' @param return_summary Return tables of model weight results
#' @param digits Integer, significant digits in printed table
#' @param n_cores Number of cores for parallel processing
#' @param subject Used to select subset of subjects (integer or character vector)
#' @param WAIC Boolean, defaults to `TRUE`. Include WAIC computed via the `loo` package for each subject.
#' @param LOO Boolean, defaults to `FALSE`. Include PSIS-LOO computed via the `loo` package for each subject.
#' @param pointwise Character string, one of `"trial"` or `"subject"`, passed to
#' `compare()` for pooled WAIC/LOO comparisons. If `subject` is specified for
#' `compare_subject()`, the pointwise unit remains the selected subject's trials.
#'
#' @return List of matrices for each subject of effective number of parameters,
#' mean deviance, deviance of mean, DIC, BPIC and associated weights.
#' @examples
#' # For a broader illustration see `compare`.
#' # Here we just take two times the same model, but normally one would compare
#' # different models
#' compare_subject(list(m0 = samples_LNR, m1 = samples_LNR))
#' @export
compare_subject <- function(sList,stage="sample",filter=0,use_best_fit=TRUE,
  print_summary=TRUE,digits=3,return_summary=FALSE,n_cores=1,subject=NULL,WAIC=TRUE,LOO=FALSE,
  pointwise = c("trial", "subject")) {
  pointwise <- match.arg(pointwise)

  compare_one <- function(subject,sList,stage,filter,use_best_fit)
    compare(sList,subject=subject,stage=stage,filter=filter,use_best_fit=use_best_fit,
            BayesFactor=FALSE,WAIC=WAIC,LOO=LOO,pointwise=pointwise,print_summary=FALSE)

  if(is(sList, "emc")) sList <- list(sList)
  subjects <- names(sList[[1]][[1]]$data)
  if (!is.null(subject)) if (is.integer(subject)) subjects <- subjects[subject] else
    if (is.character(subject))
      if (!all(subject %in% subjects)) stop("subject name(s) not in subjects") else
      subjects <- subject
  # is_single <- sapply(sList, function(x) return(x[[1]]$type == "single"))
  # if(any(!is_single)) warning("subject-by-subject comparison is best done with models of type `single`")
  if (n_cores>1) {
    out <- auto_mclapply(subjects,compare_one,sList=sList,stage=stage,filter=filter,
                   use_best_fit=use_best_fit,mc.cores=n_cores)
    names(out) <- subjects
  } else {
    out <- setNames(vector(mode="list",length=length(subjects)),subjects)
    for (i in subjects) {
      cat(".")
      out[[i]] <- compare(sList,subject=i,BayesFactor=FALSE,
        stage=stage,filter=filter,use_best_fit=use_best_fit,WAIC=WAIC,LOO=LOO,pointwise=pointwise,print_summary=FALSE)
     }
  }
  wDIC <- lapply(out,function(x)x["wDIC"])
  wBPIC <- lapply(out,function(x)x["wBPIC"])
  pDIC <- do.call(rbind,lapply(wDIC,function(x){
      setNames(data.frame(t(x)),paste("wDIC",rownames(x),sep="_"))}))
  pBPIC <- do.call(rbind,lapply(wBPIC,function(x){
      setNames(data.frame(t(x)),paste("wBPIC",rownames(x),sep="_"))}))
  has_waic <- WAIC && !is.null(out[[1]]$wWAIC)
  has_loo <- LOO && !is.null(out[[1]]$wLOO)
  if(has_waic){
    wWAIC <- lapply(out,function(x)x["wWAIC"])
    pWAIC <- do.call(rbind,lapply(wWAIC,function(x){
        setNames(data.frame(t(x)),paste("wWAIC",rownames(x),sep="_"))}))
  }
  if(has_loo){
    wLOO <- lapply(out,function(x)x["wLOO"])
    pLOO <- do.call(rbind,lapply(wLOO,function(x){
      setNames(data.frame(t(x)),paste("wLOO",rownames(x),sep="_"))}))
  }
  if (print_summary) {
    to_print <- list(DIC=round(pDIC,digits),BPIC=round(pBPIC,digits))
    if(has_waic) to_print <- c(list(WAIC=round(pWAIC,digits)), to_print)
    if(has_loo) to_print <- c(list(LOO=round(pLOO,digits)), to_print)
    print(to_print)
    mnams <- unlist(lapply(strsplit(dimnames(pDIC)[[2]],"_"),function(x){
      if (length(x)==2) x[[2]] else paste(x[-1],collapse="_")
    }))
    cat("\nWinners\n")
    winners <- list()
    if(has_waic) winners$WAIC <- table(mnams[max.col(pWAIC, ties.method = "first")])
    if(has_loo) winners$LOO <- table(mnams[max.col(pLOO, ties.method = "first")])
    winners$DIC <- table(mnams[max.col(pDIC, ties.method = "first")])
    winners$BPIC <- table(mnams[max.col(pBPIC, ties.method = "first")])
    print(do.call(rbind, winners))
  }
  if (return_summary) {
    summaries <- list(DIC=pDIC,BPIC=pBPIC)
    if(has_waic) summaries <- c(list(WAIC=pWAIC), summaries)
    if(has_loo) summaries <- c(list(LOO=pLOO), summaries)
    out <- summaries
  }
  invisible(out)
}


# #' Calculate a table of model probabilities based for a list of samples objects
# #' based on samples of marginal log-likelihood (MLL) added to these objects by
# #' run_IS2. Probabilities estimated by a bootstrap ath picks a vector of MLLs,
# #' one for each model in the list randomly with replacement nboot times,
# #' calculates model probabilities and averages
# #'
# #'
# #' @param mll List of samples objects with IS_samples attribute added by by run_IS2
# #' @param nboot Integer number of bootstrap samples, the default (1e5) usually
# #' gives stable results at 2 decimal places.
# #' @param print_summary Boolean (default TRUE) print table of results
# #' @param digits Integer, significant digits in printed table
# #'
# #' @return Vector of model probabilities with names from samples list.

compare_MLL <- function(mll,nboot=1e5,digits=2,print_summary=TRUE)
  # mll is a list of vectors of marginal log-likelihoods for a set of models
  # picks a vector of mlls for each model in the list randomly with replacement
  # nboot times, calculates model probabilities and averages.
{
  pmp <- function(x)
    # posterior model probability for a vector of marginal log-likelihoods
  {
    x <- exp(x-max(x))
    x/sum(x)
  }

  out <- sort(apply(apply(do.call(rbind,lapply(mll,function(x){
    attr(x,"IS_samples")[sample(length(attr(x,"IS_samples")),nboot,replace=TRUE)]
  })),2,pmp),1,mean),decreasing=TRUE)
  print(round(out,digits))
  invisible(out)
}





condMVN <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
{
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(given.ind) == 0)
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma))
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08)){
      sigma <- sigma + abs(diag(rnorm(nrow(sigma), sd = 1e-3)))
    }
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% chol2inv(chol(D))
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}


make_nice_summary <- function(object, stat = "max", stat_only = FALSE, stat_name = NULL, ...){
  if(is.null(stat_name)) stat_name <- stat
  row_names <- names(object)
  col_names <- unique(unlist(lapply(object, names)))
  if(all(row_names %in% col_names)){
    col_names <- row_names
  }
  out_mat <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
  for(i in 1:length(object)){
    idx <- col_names %in% names(object[[i]])
    out_mat[i,idx] <- object[[i]]
  }
  row_stat <- apply(out_mat, 1, FUN = get(stat), na.rm = T)
  out_mat <- cbind(out_mat, row_stat)

  if(nrow(out_mat) > 1){
    col_stat <- apply(out_mat, 2, FUN = get(stat), na.rm = T)
    col_stat[length(col_stat)] <- get(stat)(unlist(object))
    out_mat <- rbind(out_mat, c(col_stat))
    rownames(out_mat) <- c(row_names, stat_name)
  } else{
    rownames(out_mat) <- row_names
  }
  colnames(out_mat) <- c(col_names, stat_name)
  if(stat_only){
    out_mat <- out_mat[nrow(out_mat), ncol(out_mat)]
  }
  return(out_mat)
}


get_summary_stat <- function(emc, selection = "mu", fun, stat = NULL,
                             stat_only = FALSE, stat_name = NULL, digits = 3, ...){
  dots <- list(...)
  if(is.null(emc[[1]]$n_subjects) || length(dots$subject) == 1 || emc[[1]]$n_subjects == 1) dots$by_subject <- TRUE
  MCMC_samples <- do.call(get_pars, c(list(emc = emc, selection = selection), fix_dots(dots, get_pars)))
  out <- vector("list", length = length(MCMC_samples))
  for(i in 1:length(MCMC_samples)){
    # cat("\n", names(MCMC_samples)[[i]], "\n")
    if(length(fun) > 1){
      outputs <- list()
      for(j in 1:length(fun)){
        outputs[[j]] <- do.call(fun[[j]], c(list(MCMC_samples[[i]]), fix_dots(dots, fun[[j]])))
      }
      out[[i]] <- do.call(cbind, outputs)
      if(!is.null(stat_name)){
        if(ncol(out[[i]]) != length(stat_name)) stop("make sure stat_name is the same length as function output")
        colnames(out[[i]]) <- stat_name
      }
    } else{
      out[[i]] <- do.call(fun, c(list(MCMC_samples[[i]]), fix_dots(dots, fun)))#fun(MCMC_samples[[i]], ...)
    }
  }
  names(out) <- names(MCMC_samples)
  if(length(fun) == 1 & !is.matrix(out[[i]]) & !is.null(stat)){
    out <- make_nice_summary(out, stat, stat_only, stat_name)
    out <- round(out, digits)
  } else{
    out <- lapply(out, round, digits)
  }
  return(out)
}


get_posterior_quantiles <- function(x, probs = c(0.025, .5, .975)){
  summ <- summary(x, probs)
  return(summ$quantiles)
}

#' Model Averaging
#'
#' Computes model weights and a Bayes factor by comparing two groups of models based on their
#' Information Criterion (IC) values. The function works with either numeric vectors or data
#' frames containing multiple IC measures (e.g., MD, BPIC, DIC).
#'
#' When provided with numeric vectors, it computes the weights for the two groups by first
#' converting the IC values into relative weights and then normalizing them. When provided with
#' a data frame, it assumes that the data frame is the output of a call to `compare`
#' and applies averaging to each IC metric
#'
#' @param IC_for A numeric vector or the output of `compare`
#' @param IC_against A numeric vector or the output of `compare`
#'
#' @return A \code{data.frame} with the following columns:
#'   \describe{
#'     \item{\code{wFor}}{The aggregated weight of the models in favor.}
#'     \item{\code{wAgainst}}{The aggregated weight of the models against.}
#'     \item{\code{Factor}}{The Bayes factor (ratio of \code{wFor} to \code{wAgainst}).}
#'   }
#'   If \code{IC_for} is a data frame, a matrix with rows corresponding to each IC measure is returned.
#'
#' @examples
#' # First set up some example models (normally these would be alternative models)
#' samples_LNR2 <- subset(samples_LNR, length.out = 45)
#' samples_LNR3 <- subset(samples_LNR, length.out = 40)
#' samples_LNR4 <- subset(samples_LNR, length.out = 35)
#'
#' # Run compare on them, BayesFactor = F is set for speed.
#' ICs <- compare(list(S1 = samples_LNR, S2 = samples_LNR2,
#'                     S3 = samples_LNR3, S4 = samples_LNR4), BayesFactor = FALSE)
#'
#' # Model averaging can either be done with a vector of ICs:
#' model_averaging(ICs$BPIC[1:2], ICs$BPIC[2:4])
#'
#' # Or the output of compare:
#' model_averaging(ICs[1:2,], ICs[3:4,])
#'
#' @export
model_averaging <- function(IC_for, IC_against) {
  if(is.null(IC_for)) return(NULL)

  if(is.data.frame(IC_for)){
    # Recursive call to make it work with the output of compare
    MD <- model_averaging(IC_for$MD, IC_against$MD)
    BPIC <- model_averaging(IC_for$BPIC, IC_against$BPIC)
    DIC <- model_averaging(IC_for$DIC, IC_against$DIC)
    return(rbind(MD = MD, BPIC = BPIC, DIC = DIC))
  }
  # Combine the IC values from both groups
  all_IC <- c(IC_for, IC_against)

  # Find the smallest IC value (for numerical stability)
  min_IC <- min(all_IC)

  # Compute the unnormalized weights
  unnorm_weights <- exp(-0.5 * (all_IC - min_IC))

  # Normalize weights so they sum to 1
  weights <- unnorm_weights / sum(unnorm_weights)

  # Separate the weights for the two groups
  weight_for <- sum(weights[seq_along(IC_for)])
  weight_against <- sum(weights[(length(IC_for) + 1):length(all_IC)])

  # Compute the Bayes factor: evidence in favor relative to against
  bayes_factor <- weight_for / weight_against

  # Return the results as a data frame
  return(data.frame(
    wFor = weight_for,
    wAgainst = weight_against,
    Factor = bayes_factor
  ))
}
