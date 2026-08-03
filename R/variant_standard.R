sample_store_standard <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, ...){
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)

  # Get total parameters including regressors
  group_designs <- list(...)$group_design
  total_par_names <- add_group_par_names(par_names, group_designs)

  # Create samples structure
  samples <- list(
    theta_mu = array(NA_real_, dim = c(length(total_par_names), iters),
                     dimnames = list(total_par_names, NULL)),
    theta_var = array(NA_real_, dim = c(n_pars, n_pars, iters),
                     dimnames = list(par_names, par_names, NULL)),
    a_half = array(NA_real_, dim = c(n_pars, iters),
                  dimnames = list(par_names, NULL))
  )

  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

add_info_standard <- function(sampler, prior = NULL, ...){
  n_pars <- sum(!sampler$nuisance)
  group_design <- list(...)$group_design

  # Get all parameter names including regressors
  sampler$par_names_all <- add_group_par_names(sampler$par_names, group_design)

  sampler$par_group <- list(...)$par_groups
  if(is.null(sampler$par_group)) sampler$par_group <- rep(1, n_pars)
  sampler$is_blocked <- sampler$par_group %in% which(table(sampler$par_group) > 1)
  sampler$prior <- get_prior_standard(prior, n_pars, sample = F,
                                      group_design = group_design)
  sampler$group_designs <- group_design
  
  if (!is.null(group_design)) {
    gd <- add_group_design(sampler$par_names[!sampler$nuisance], group_design, sampler$n_subjects)
    X_all <- do.call(cbind, gd)
    sampler$XtX <- crossprod(X_all)
    n_cols <- vapply(gd, ncol, integer(1))
    sampler$design_row_idx <- rep(seq_along(n_cols), times = n_cols)
    sampler$design_n_cols <- n_cols
    sampler$gd <- gd
  } else {
    sampler$gd <- NULL
  }
  
  return(sampler)
}

calculate_mean_design <- function(group_designs, n_pars) {
  # Initialize the mean design matrix for each parameter
  mean_designs <- list()

  for (k in 1:n_pars) {
    if (is.null(group_designs) || is.null(group_designs[[k]])) {
      # If no design matrix, use a simple intercept
      mean_designs[[k]] <- matrix(1, nrow = 1, ncol = 1)
    } else {
      # Calculate the mean of the design matrix across subjects
      mean_designs[[k]] <- colMeans(group_designs[[k]], dims = 1)
    }
  }

  return(mean_designs)
}

calculate_implied_means <- function(mean_designs, beta_params, n_pars) {
  mu_implied <- numeric(n_pars)

  par_idx <- 0
  for (k in 1:n_pars) {
    x_k_mean <- mean_designs[[k]]
    mu_implied[k] <- x_k_mean %*% beta_params[par_idx + 1:length(x_k_mean)]
    par_idx <- par_idx + length(x_k_mean)
  }

  return(mu_implied)
}

get_prior_standard <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL, group_design = NULL,
                               par_groups = NULL){
  # Checking and default priors
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_pars(design, doMap = F))
  }
  if (!is.null(prior$A) & is.null(n_pars)) {
    n_pars <- length(prior$A)
  }

  # Number of additional parameters from design matrices
  if(!is.null(group_design)){
    n_additional <- n_additional_group_pars(group_design)
  } else{
    n_additional <- 0
  }

  # Set up combined theta_mu_mean to include both intercepts and slopes
  if(is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars + n_additional)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars + n_additional))
  }
  if(is.null(prior$v)){
    prior$v <- 2
  }
  if(is.null(prior$A)){
    prior$A <- rep(.3, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix
  attr(prior, "type") <- "standard"
  out <- prior
  if(sample){
    par_names <- names(sampled_pars(design, doMap = F))
    samples <- list()
    if(selection %in% c("mu", "beta", "alpha")){
      # Sample beta (all parameters including regressors)
      beta <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                             sigma = prior$theta_mu_var))
      rownames(beta) <- add_group_par_names(par_names, group_design)
      if(selection %in% c("beta", "mu", "alpha")){
        samples$theta_mu <- beta
        if(selection %in% c("mu", "alpha")){
          samples$par_names <- par_names
          samples$group_designs <- group_design
        }
      }
    }
    if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for(i in 1:N){
        a_half <- 1 / rgamma(n = n_pars,shape = 1/2,
                             rate = 1/(prior$A^2))
        attempt <- tryCatch({
          vars[,,i] <- riwish(prior$v + n_pars - 1, 2 * prior$v * diag(1 / a_half))
        },error=function(e) e, warning=function(w) w)
        if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
          sample_idx <- sample(1:(i-1),1)
          vars[,,i] <- vars[,,sample_idx]
        }
      }
      if(is.null(par_groups)) par_groups <- rep(1, n_pars)
      constraintMat <- matrix(0, n_pars, n_pars)
      for(i in 1:length(unique(par_groups))){
        idx <- par_groups == i
        constraintMat[idx, idx] <- Inf
      }
      vars <- constrain_lambda(vars, constraintMat)
      if(selection != "alpha") samples$theta_var <- vars
    }
    if(selection %in% "alpha"){
      samples <- list()
      samples$alpha <- get_alphas(beta, vars, design[[1]]$Ffactors$subjects,
                                  group_design = group_design)
    }
    out <- samples
  }
  return(out)
}

get_startpoints_standard <- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!pmwgs$nuisance)
  n_total_pars <- length(pmwgs$prior$theta_mu_mean) # Includes regressor parameters
  if (is.null(start_mu)) start_mu <- rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)[1,]
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3, diag(n_pars))
  start_a_half <- 1 / rgamma(n = n_pars, shape = 2, rate = 1)

  # Calculate subject-specific means using design matrices
  group_designs <- add_group_design(pmwgs$par_names[!pmwgs$nuisance], pmwgs$group_designs, pmwgs$n_subjects)
  subj_mu <- calculate_subject_means(group_designs, start_mu)
  return(list(tmu = start_mu, tvar = start_var, tvinv = ginv(start_var),
              a_half = start_a_half, subj_mu = subj_mu))
}

get_group_level_standard <- function(parameters, s){
  # This function is modified for other versions
  mu <- parameters$subj_mu[,s]
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

fill_samples_standard <- function(samples, group_level, proposals, j = 1, n_pars){
  samples$a_half[, j] <- group_level$a_half
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}

# Data-based precision and mean for the joint group-level update.
# The per-subject map M_i (p x M) is block-sparse: row k occupies the k-th
# design's columns, so Sum_i M_i' tvinv M_i and Sum_i M_i' tvinv alpha_i
# assemble blockwise from crossproducts over subjects:
#   prec[rows_k, cols_l] = tvinv[k, l] * X_k' X_l
#   mean[rows_k]         = X_k' (tvinv alpha)[k, ]
# This replaces an O(n * p * M^2) subject loop with p^2 small crossproducts.
group_design_moments <- function(group_designs, tvinv, alpha, M) {
  p <- length(group_designs)
  prec_data <- matrix(0, M, M)
  mean_data <- numeric(M)
  n_cols  <- vapply(group_designs, ncol, integer(1))
  col_off <- cumsum(c(0L, n_cols))
  ta <- tvinv %*% alpha  # p x n
  for (k in seq_len(p)) {
    Xk   <- group_designs[[k]]
    rows <- col_off[k] + seq_len(n_cols[k])
    mean_data[rows] <- crossprod(Xk, ta[k, ])
    for (l in seq_len(p)) {
      tv_kl <- tvinv[k, l]
      if (tv_kl == 0) next
      cols <- col_off[l] + seq_len(n_cols[l])
      prec_data[rows, cols] <- tv_kl * crossprod(Xk, group_designs[[l]])
    }
  }
  list(prec_data = prec_data, mean_data = mean_data)
}

gibbs_step_standard <- function(sampler, alpha) {
  #
  # alpha:  (p x n) subject-level parameter matrix
  #
  # sampler$group_designs:   list of length p, each (n x m_k)
  # sampler$is_blocked:    length p logical; is_blocked[k] = TRUE => param k is in some covariance block
  # sampler$par_group:     length p integer group labels (only used if is_blocked=TRUE)
  #
  # sampler$prior contains at least:
  #   - theta_mu_mean (length total_params) and theta_mu_invar (total_params x total_params)
  #   - v, A => hyperparams for the IW/Gamma updates
  #
  # last_sample_standard(...) must provide:
  #   - last$tmu (total_params-vector) of group means (includes both intercepts and slopes)
  #   - last$tvar, last$tvinv (p x p)
  #   - last$a_half (p)
  #
  # The function returns an updated (tmu, tvar, tvinv, a_half, alpha),
  # doing partial block updates for the covariance tvar, and a joint
  # update for all mu parameters (including former beta parameters)


  group_designs <- sampler$group_designs
  is_blocked    <- sampler$is_blocked[!sampler$nuisance]
  par_group     <- sampler$par_group[!sampler$nuisance]
  prior         <- sampler$prior

  last    <- last_sample_standard(sampler$samples)
  tmu     <- last$tmu        # group means, dimension total_params (includes former beta)
  tvar    <- last$tvar       # (p x p)
  tvinv   <- last$tvinv
  a_half  <- last$a_half     # length p

  p <- nrow(alpha)           # number of subject-level parameters
  n <- ncol(alpha)           # number of subjects

  marginal_idx <- sampler$marginalised_idx
  if (is.null(marginal_idx)) marginal_idx <- rep(FALSE, p)
  if (length(marginal_idx) != p && !is.null(sampler$nuisance)) {
    marginal_idx <- marginal_idx[!sampler$nuisance]
  }
  if (length(marginal_idx) != p) marginal_idx <- rep(FALSE, p)
  marginal_spec <- sampler$marginalise
  marginal_idx <- as.logical(marginal_idx)

  # Some backwards compatibility
  if(is.null(is_blocked) || length(is_blocked) == 0){
    is_blocked <- rep(T, p)
  }
  if(is.null(par_group)){
    par_group <- rep(1, p)
  }
  if (any(marginal_idx) && !is.null(marginal_spec)) {
    # Never let a marginalized coordinate define or enlarge an IW block.
    is_blocked[marginal_idx] <- FALSE
  }

  has_group_designs <- !is.null(sampler$gd) || (!is.null(group_designs) && length(group_designs) > 0)
  if (has_group_designs) {
    if (!is.null(sampler$gd)) {
      group_designs <- sampler$gd
    } else {
      group_designs <- add_group_design(sampler$par_names[!sampler$nuisance], group_designs, n)
    }
  } else {
    group_designs <- NULL
  }

  # Stage 1 marginalization uses a fixed prior eta for the integrated
  # coordinate. Remove its precision cross-terms from the group update and
  # replace its residual by zero so the reconstructed draws cannot alter the
  # covariance of the ordinary parameters.
  if (any(marginal_idx) && !is.null(marginal_spec)) {
    tvinv[marginal_idx, ] <- 0
    tvinv[, marginal_idx] <- 0
    tvinv[cbind(which(marginal_idx), which(marginal_idx))] <- 1 / marginal_spec$sigma^2
  }

  # Calculate total parameters (including regressors)
  M <- length(tmu)  # Total parameters (former mu + beta)

  ##--------------------------------------------------
  ## 1) Build data-based precision & mean
  ##--------------------------------------------------
  if (has_group_designs) {
    if (!is.null(sampler$XtX)) {
      Tvinv_expanded <- tvinv[sampler$design_row_idx, sampler$design_row_idx]
      prec_data <- Tvinv_expanded * sampler$XtX
      
      mean_data <- numeric(M)
      ta <- tvinv %*% alpha
      col_off <- cumsum(c(0L, sampler$design_n_cols))
      for (k in seq_along(sampler$design_n_cols)) {
        rows <- col_off[k] + seq_len(sampler$design_n_cols[k])
        mean_data[rows] <- crossprod(group_designs[[k]], ta[k, ])
      }
    } else {
      moments   <- group_design_moments(group_designs, tvinv, alpha, M)
      prec_data <- moments$prec_data
      mean_data <- moments$mean_data
    }
  } else {
    prec_data <- n * tvinv
    mean_data <- as.vector(tvinv %*% rowSums(alpha))
  }

  prec_post <- prior$theta_mu_invar + prec_data
  b_post <- as.vector(prior$theta_mu_invar %*% prior$theta_mu_mean + mean_data)

  # Draw new parameter vector
  R_post   <- chol(prec_post)
  w        <- forwardsolve(t(R_post), b_post)
  tmu_new  <- backsolve(R_post, w + rnorm(M))

  ##--------------------------------------------------
  ## 2) Compute residuals = alpha - (X * tmu)
  ##--------------------------------------------------
  resid <- matrix(0, nrow=p, ncol=n)

  # Calculate subject-level means using the new parameters
  if (has_group_designs) {
    subj_mu <- calculate_subject_means(group_designs, tmu_new)
  } else {
    subj_mu <- matrix(tmu_new, nrow = p, ncol = n)
  }

  # Compute residuals
  for (i in seq_len(n)) {
    resid[, i] <- alpha[, i] - subj_mu[, i]
  }
  if (any(marginal_idx) && !is.null(marginal_spec)) {
    resid[marginal_idx, ] <- 0
  }

  ##--------------------------------------------------
  ## 3) Partial-block update of tvar_new and tvinv_new
  ##--------------------------------------------------
  tvar_new <- matrix(0, p, p)
  tvinv_new <- matrix(0, p, p)

  blocked_idx   <- which(is_blocked)
  unblocked_idx <- which(!is_blocked)

  # 3a) For each block among the "blocked" subset, do an inverse-Wishart
  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[ par_group[blocked_idx] == g ]
      d <- length(group_idx)
      resid_block <- resid[group_idx, , drop=FALSE]
      cov_block   <- tcrossprod(resid_block)

      B_half_block <- 2 * prior$v * diag(1 / a_half[group_idx], d) + cov_block
      df_block     <- prior$v + d - 1 + n

      Sigma_block  <- riwish(df_block, B_half_block)
      tvar_new[group_idx, group_idx] <- Sigma_block
      tvinv_new[group_idx, group_idx] <- chol2inv(chol(Sigma_block))
    }
  }

  # 3b) For each "unblocked" param, do a vectorized Gamma for diagonal-only
  if (length(unblocked_idx) > 0) {
    sum_sq_vec <- rowSums(resid[unblocked_idx, , drop=FALSE]^2)
    # shape = v/2 + n/2
    shape_vec <- prior$v/2 + n/2
    # rate = v/a_half + 0.5 * sum_sq
    rate_vec  <- prior$v / a_half[unblocked_idx] + 0.5 * sum_sq_vec

    tvinv_diag <- rgamma(length(unblocked_idx), shape=shape_vec, rate=rate_vec)
    tvar_diag  <- 1 / tvinv_diag

    tvar_new[cbind(unblocked_idx, unblocked_idx)] <- tvar_diag
    tvinv_new[cbind(unblocked_idx, unblocked_idx)] <- tvinv_diag
  }

  if (any(marginal_idx) && !is.null(marginal_spec)) {
    # Locate the stacked group-design coefficients for each marginalized
    # subject parameter. An intercept-only t0 design has one coefficient; if a
    # caller supplied a richer group design, keep eta in the intercept and
    # hold the additional coefficients at zero.
    if (has_group_designs) {
      col_off <- cumsum(c(0L, vapply(group_designs, ncol, integer(1L))))
    } else {
      col_off <- 0:p
    }
    for (k in which(marginal_idx)) {
      cols <- (col_off[k] + 1L):col_off[k + 1L]
      tmu_new[cols] <- c(marginal_spec$mu, rep(0, length(cols) - 1L))
      tvar_new[k, ] <- 0
      tvar_new[, k] <- 0
      tvar_new[k, k] <- marginal_spec$sigma^2
      
      tvinv_new[k, ] <- 0
      tvinv_new[, k] <- 0
      tvinv_new[k, k] <- 1 / marginal_spec$sigma^2
    }
    if (has_group_designs) {
      subj_mu <- calculate_subject_means(group_designs, tmu_new)
    } else {
      subj_mu <- matrix(tmu_new, nrow = p, ncol = n)
    }
    subj_mu[marginal_idx, ] <- marginal_spec$mu
  }
  ##--------------------------------------------------
  ## 4) Update a_half
  ##--------------------------------------------------
  # shape depends on block dimension: (v + block_dim[k]) / 2
  block_dim <- integer(p)
  block_dim[unblocked_idx] <- 1
  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[ par_group[blocked_idx] == g ]
      block_dim[group_idx] <- length(group_idx)
    }
  }

  shape_vec <- (prior$v + block_dim) / 2
  rate_vec  <- prior$v * diag(tvinv_new) + 1/(prior$A^2)
  a_half_new <- 1 / rgamma(p, shape=shape_vec, rate=rate_vec)
  if (any(marginal_idx) && !is.null(marginal_spec)) {
    a_half_new[marginal_idx] <- a_half[marginal_idx]
  }

  ##--------------------------------------------------
  ## 5) Return updated values
  ##--------------------------------------------------
  out <- list(
    tmu    = tmu_new,          # updated parameters
    tvar   = tvar_new,
    tvinv  = tvinv_new,
    a_half = a_half_new,
    subj_mu = subj_mu,
    alpha  = alpha
  )
  return(out)
}


last_sample_standard <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],  # Now includes beta
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
}

get_conditionals_standard <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  pts2_unwound <- apply(samples$theta_var[idx,idx,],3,unwind)
  all_samples <- rbind(samples$alpha[idx, s,],samples$theta_mu[idx,],pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration], unwind(samples$theta_var[idx,idx,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

filtered_samples_standard <- function(sampler, filter, ...){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter, drop = F],
    theta_var = sampler$samples$theta_var[, , filter, drop = F],
    alpha = sampler$samples$alpha[, , filter, drop = F],
    iteration = length(filter)
  )
}


calc_log_jac_chol <- function(x) {
  ## work out dimension -------------------------------------------------
  n <- floor(sqrt(2 * length(x) + 0.25) - 0.5)   # solves n(n+1)/2 = length(x)

  ## rebuild the raw lower‑triangular matrix ----------------------------
  mat <- matrix(NA_real_, n, n)
  mat[lower.tri(mat, diag = TRUE)] <- x          # raw params (diag on log‑scale)

  ## step 1: create L ---------------------------------------------------
  L <- mat
  diag(L) <- exp(diag(mat))                      # positive diagonal
  # (off‑diagonals stay on natural scale)

  ## log‑Jacobian -------------------------------------------------------
  logL <- diag(mat)                              # log of L_ii
  log_jac <- n * log(2) + sum((n - seq_len(n) + 2) * logL)

  return(log_jac)
}

unwind_chol <- function(x,reverse=FALSE) {

  if (reverse) {
    n=sqrt(2*length(x)+0.25)-0.5 ## Dim of matrix.
    out=array(0,dim=c(n,n))
    out[lower.tri(out,diag=TRUE)]=x
    diag(out)=exp(diag(out))
    out=out%*%t(out)
  } else {
    y=t(base::chol(x))
    diag(y)=log(diag(y))
    out=y[lower.tri(y,diag=TRUE)]
  }
  return(out)
}

# bridge_sampling ---------------------------------------------------------
bridge_group_and_prior_and_jac_standard <- function(
    proposals_group,   # (n_iter x [theta_mu, theta_a, var1, var2])
    proposals_list,    # list of length n_iter, each item is the subject-level alpha draws
    info
) {
  # 'info' should contain:
  #   info$n_pars        = p        (# of rows in alpha)
  #   info$n_subjects    = n_subj   (# of subjects)
  #   info$is_blocked    = logical p  (which rows are blocked)
  #   info$par_group     = integer p  (group IDs for blocked rows)
  #   info$group_designs   = list of length p, each (n_subj x m_k)
  #   info$prior:
  #       - $theta_mu_mean (p+B), $theta_mu_var ((p+B)x(p+B))
  #       - $v, $A (IW & gamma hyperparams)
  #
  # The 'proposals_group' columns are arranged:
  #   1) theta_mu (p+B columns) - includes both intercepts and slopes
  #   2) theta_a  (p columns) -> log(a_half)
  #   3) theta_var1 (sum(!has_cov) columns) -> log of unblocked diagonal
  #   4) theta_var2 ((d*(d+1))/2 columns)    -> cholesky for the blocked subset
  #
  # 'proposals_list' is a list of length n_iter, each item is an (n_subj * p)-vector or
  #   something that we can reshape into (n_subj x p), i.e. the subject-level alphas.

  par_group <- info$par_group
  has_cov   <- info$is_blocked
  block_groups <- unique(par_group[has_cov])
  prior     <- info$prior

  p         <- info$n_pars        # dimension of each alpha_i
  n_subj    <- info$n_subjects
  group_designs <- add_group_design(info$par_names, info$group_designs, n_subj)

  # Total parameters (including regressors)
  total_pars <- length(prior$theta_mu_mean)
  # Extract columns:
  theta_mu   <- proposals_group[, seq_len(total_pars), drop=FALSE]  # (n_iter x (p+B))
  theta_a    <- proposals_group[, total_pars + seq_len(p), drop=FALSE]  # (n_iter x p)

  if(any(!has_cov)){
    theta_var1 <- proposals_group[, (total_pars + p + 1) : (total_pars + p + sum(!has_cov)), drop=FALSE]  # (n_iter x sum(!has_cov))
  }

  if(any(has_cov)){
    min_idx <- (total_pars + p + sum(!has_cov))
    theta_var2_list <- list()
    for(block in block_groups){
      cur_idx <- has_cov[par_group == block]
      max_idx <- min_idx + (sum(cur_idx)*(sum(cur_idx)+1)) / 2
      theta_var2_list[[block]] <- proposals_group[, (min_idx + 1) :
                                      max_idx, drop=FALSE]  # (n_iter x (d*(d+1))/2)
      min_idx <- max_idx
    }
  }

  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)

  # Precompute the prior on theta_mu outside the loop (vectorized)
  prior_mu_log <- dmvnorm(theta_mu,
                          mean  = prior$theta_mu_mean,
                          sigma = prior$theta_mu_var,
                          log   = TRUE)

  # For the Jacobian of var1 => rowSums(theta_var1)
  # For the Jacobian of a    => rowSums(theta_a)
  if(any(!has_cov)){
    jac_var1 <- rowSums(theta_var1)
  } else{
    jac_var1 <- 0
  }
  jac_a_vec    <- rowSums(theta_a)

  #--- Main loop over MCMC draws ---
  # We reconstruct var_curr & compute the group-likelihood, plus prior on var1/var2/a
  # then we store in sum_out[i].  After the loop, we add prior_mu_log + jac_var1 + jac_a.
  for (i in seq_len(n_iter)) {
    # 1) Reconstruct partial-block var_curr
    var_curr <- matrix(0, p, p)

    # Unblocked diagonal => exp(theta_var1[i,])
    if (any(!has_cov)) {
      var_curr[!has_cov, !has_cov] <- diag(exp(theta_var1[i, ]), sum(!has_cov))
    }

    # Jacobian for var2
    jac_var2 <- 0
    prior_var2 <- 0
    # Blocked => unwind_chol
    if (any(has_cov)) {
      # Fill in theta_var2
      for(block in block_groups){
        cur_idx <- par_group == block
        var_block <- unwind_chol(theta_var2_list[[block]][i, ], reverse=TRUE)
        var_curr[cur_idx, cur_idx] <- var_block
        # Prior influence
        d_block   <- sum(cur_idx)
        prior_var2 <- prior_var2 + log( robust_diwish(
          W = var_block,
          v = prior$v + d_block - 1,
          S = 2 * prior$v * diag( 1 / exp(theta_a[i,cur_idx]) )
        ))
        # Jacobian influence
        jac_var2 <- jac_var2 + calc_log_jac_chol(theta_var2_list[[block]][i, ])
      }
    }

    # 3) Group likelihood = sum_{s=1..n_subj} dmvnorm(alpha_s, mu_s, var_curr).
    #    Project the stacked beta vector to per-subject means with the shared C++
    #    routine (same as the sampler start/gibbs steps), then center each
    #    subject's alpha and evaluate the zero-mean MVN density in one call.
    subj_mu <- calculate_subject_means(group_designs, theta_mu[i, ])  # p x n_subj
    alpha_i <- matrix(vapply(proposals_list, function(pr) pr[i, ], numeric(p)), nrow = p)  # p x n_subj
    
    U <- tryCatch(chol(var_curr), error = function(e) NULL)
    if (is.null(U)) {
      U <- tryCatch(chol(var_curr + diag(1e-6, p)), error = function(e) NULL)
    }
    if (is.null(U)) {
      group_ll <- -Inf
    } else {
      rooti <- backsolve(U, diag(p))
      log_const <- sum(log(diag(rooti))) - 0.5 * p * log(2*pi)
      z <- t(alpha_i - subj_mu) %*% rooti
      group_ll <- -0.5 * sum(z^2) + n_subj * log_const
    }

    # 4) Prior on var1, var2, and a => same partial-block logic
    # "var1" => Inverse-Gamma with shape=v/2, rate=v/exp(theta_a[i, !has_cov])
    # "var2" => robust_diwish( var_block, v= v + sum(has_cov)-1, S= 2 v diag(1/ a_half ) ), but a_half=exp(theta_a)
    # "a" => a_half = exp(theta_a), shape=1/2, rate=1/(A^2)
    prior_var1 <- 0
    if (any(!has_cov)) {
      # sum of logdinvGamma for each unblocked row
      # x = exp(theta_var1), shape= v/2, rate= v / exp(theta_a)
      # note that 'theta_a[i, !has_cov]' picks the relevant subset
      prior_var1 <- sum( logdinvGamma(
        x     = exp(theta_var1[i, ]),
        shape = prior$v/2,
        rate  = prior$v / exp( theta_a[i, !has_cov, drop=FALSE] )
      ))
    }

    prior_a <- sum(logdinvGamma(
      x     = exp(theta_a[i, ]),
      shape = 1/2,
      rate  = 1/(prior$A^2)
    ))
    # Combine
    sum_out[i] <- group_ll + prior_var1 + prior_var2 + jac_var2 + prior_a
  } # end loop i

  #--- 5) Add the vectorized prior on theta_mu, plus Jacobians var1, a
  sum_out <- sum_out + prior_mu_log + jac_var1 + jac_a_vec

  return(sum_out)
}


bridge_add_info_standard <- function(info, samples){
  if(is.null(samples$is_blocked)){
    has_cov <- rep(T, samples$n_pars)
  } else{
    has_cov <- samples$is_blocked
  }
  info$par_group <- samples$par_group
  if(is.null(info$par_group)){
    info$par_group <- rep(1, samples$n_pars)
  }
  info$is_blocked <- has_cov
  info$group_designs <- samples$group_designs

  # Calculate total parameters (including former beta)
  n_total_pars <- nrow(samples$samples$theta_mu)

  # Calculate group index including all parameters (theta_mu now includes beta)
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):
    (samples$n_pars*samples$n_subjects + n_total_pars + samples$n_pars +
       sum(!has_cov) + (sum(has_cov) * (sum(has_cov) +1))/2)

  return(info)
}

bridge_add_group_standard <- function(all_samples, samples, idx){
  # Add theta_mu (now includes all parameters)
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(log(samples$samples$a_half[,idx])))

  # Handle variance matrices
  par_group <- samples$par_group
  has_cov <- samples$is_blocked
  if(is.null(has_cov)){
    has_cov <- rep(T, samples$n_pars)
  }
  if(is.null(par_group)){
    par_group <- rep(1, samples$n_pars)
  }

  if(any(!has_cov)){
    all_samples <- cbind(all_samples, t(log(matrix(apply(samples$samples$theta_var[!has_cov,!has_cov,idx, drop = F], 3, diag), ncol = nrow(all_samples)))))
  }
  if(any(has_cov)){
    for(block in unique(par_group[has_cov])){
      cur_idx <- par_group == block
      all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$theta_var[cur_idx,cur_idx,idx, drop = F], 3, unwind_chol), ncol = nrow(all_samples))))
    }
  }
  return(all_samples)
}

# for IC ------------------------------------------------------------------

group__IC_standard <- function(emc, stage="sample", filter=NULL) {
  # 1) Retrieve draws
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    return_mcmc = FALSE, merge_chains = TRUE)
  # With group designs the subject-level means are built from the full stacked
  # coefficient vector (M = intercepts + slopes), so fetch "beta" (theta_beta),
  # not the mapped p-vector "mu": calculate_subject_means and standard_subj_ll
  # need the M-vector, and "mu" would abort inside the design projection.
  mu_selection <- if (is.null(emc[[1]]$group_designs)) "mu" else "beta"
  theta_mu <- get_pars(emc, selection = mu_selection, stage = stage, filter = filter,
                       return_mcmc = FALSE, merge_chains = TRUE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                        return_mcmc = FALSE, merge_chains = TRUE, remove_constants = F)

  p <- dim(alpha)[1]           # number of subject-level parameters
  n_subj <- dim(alpha)[2]      # number of subjects
  N <- dim(alpha)[3]           # number of samples

  # 2) Averages. Under group designs theta_mu holds the stacked betas, so
  # mean_mu is the mean beta vector (length M); with no group design it is the
  # mean mu vector (length p). Both feed the matching branch below.
  mean_alpha <- rowMeans(alpha, dims = 2)
  mean_mu <- rowMeans(theta_mu)
  mean_var <- rowMeans(theta_var, dims = 2)

  if(is.null(emc[[1]]$group_designs)){
    lls <- numeric(N)
    for(i in 1:N){
      lls[i] <- sum(dmvnorm(t(alpha[,,i]), theta_mu[,i], theta_var[,,i], log = T))
    }
    mean_pars_ll <-  sum(dmvnorm(t(mean_alpha), mean_mu, mean_var, log = TRUE))
  } else{
    # 3) Build design matrices
    group_designs <- add_group_design(emc[[1]]$par_names, emc[[1]]$group_designs, n_subj)
    # 4) Per-draw log-likelihood (argument order must match the signature:
    #    theta_var, theta_mu, alpha, n_subj, group_designs)
    lls <- standard_subj_ll(theta_var, theta_mu, alpha, n_subj, group_designs)
    # 5) Likelihood at posterior mean
    # Calculate subject-level means using the mean parameters
    mean_subj_means <- calculate_subject_means(group_designs, mean_mu)
    mean_pars_ll <- 0
    for(s in seq_len(n_subj)) {
      mean_pars_ll <- mean_pars_ll +
        dmvnorm(mean_alpha[, s], mean_subj_means[, s], mean_var, log=TRUE)
    }

  }

  minD <- -2*max(lls)
  mean_ll <- mean(lls)
  Dmean <- -2*mean_pars_ll
  list(mean_ll = mean_ll, Dmean = Dmean, minD = minD)
}

standard_subj_ll <- function(theta_var, theta_mu, alpha, n_subj, group_designs)
{
  # p is the subject-parameter dimension (rows of alpha / Sigma). theta_mu holds
  # the stacked betas (M rows, M >= p under group designs), so p must come from
  # alpha, not theta_mu, or the backsolve/quadratic form get the wrong dimension.
  N <- dim(theta_var)[3];  p <- nrow(alpha)
  log2pi <- log(2*pi)

  # pre‑allocate
  ll <- numeric(N)

  for (i in seq_len(N)) {

    Sigma <- theta_var[,,i]
    U     <- chol(Sigma)              # will error if non‑PD
    rooti <- backsolve(U, diag(p))
    log_const <- sum(log(diag(rooti))) - 0.5 * p * log2pi

    mu  <- calculate_subject_means(group_designs, theta_mu[,i])
    z   <- t(alpha[,,i] - mu) %*% rooti      # n_subj × p
    ll[i] <- -0.5 * sum(z^2) + n_subj*log_const
  }
  ll
}
