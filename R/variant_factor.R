
sample_store_factor <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, ...) {
  n_factors <- list(...)$n_factors
  Lambda_mat <- list(...)$Lambda_mat
  if(is.null(n_factors)){
    n_factors <- ncol(Lambda_mat)
  }
  if(!is.null(Lambda_mat)){
    f_names <- colnames(Lambda_mat)
  } else{
    f_names <- paste0("F", 1:n_factors)
  }
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, f_names, NULL)),
    lambda_untransf = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, f_names, NULL)),
    epsilon_inv = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL)),
    psi_inv = array(NA_real_, dim = c(n_factors, iters), dimnames = list(NULL, NULL)),
    eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, NULL, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}


add_info_factor <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  n_factors <- args$n_factors
  Lambda_mat <- args$Lambda_mat
  if(is.null(n_factors)) n_factors <- ncol(Lambda_mat)
  n_pars <- sum(!sampler$nuisance)
  if(is.null(Lambda_mat)){
    Lambda_mat <- matrix(Inf, nrow = n_pars, ncol = n_factors)
    diag(Lambda_mat) <- 1
    Lambda_mat[upper.tri(Lambda_mat, diag = F)] <- 0
  }
  if(!is.null(args$signFix)){
    signFix <- args$signFix
  } else{
    signFix <- F
  }

  attr(sampler, "signFix") <- signFix
  attr(sampler, "Lambda_mat") <- Lambda_mat

  sampler$prior <- get_prior_factor(prior, sum(!sampler$nuisance), sample = F, n_factors = n_factors)
  sampler$n_factors <- n_factors
  return(sampler)
}

get_prior_factor <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL,
                             Lambda_mat = NULL, n_factors = NULL){


  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_pars(design, doMap = F))
  }
  if(is.null(Lambda_mat)) Lambda_mat <- attr(prior, "Lambda_mat")
  if(is.null(Lambda_mat)){
    Lambda_mat <- matrix(Inf, nrow = n_pars, ncol = n_factors)
    diag(Lambda_mat) <- 1
    Lambda_mat[upper.tri(Lambda_mat, diag = F)] <- 0
  }
  if(is.null(n_factors)) n_factors <- ncol(Lambda_mat)

  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- rep(1, n_pars)
  }
  if(is.null(prior$theta_lambda_var)){
    prior$theta_lambda_var <- rep(.2, n_pars)
  }
  if(is.null(prior$ap)){
    prior$ap <- 2
  }
  if(is.null(prior$bp)){
    prior$bp <- .5
  }
  if(is.null(prior$as)){
    prior$as <- rep(2.5, n_pars)
  }
  if(is.null(prior$bs)){
    prior$bs <- rep(.1, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- diag(1/prior$theta_mu_var)
  prior$theta_lambda_invar <-1/prior$theta_lambda_var
  # Things I save rather than re-compute inside the loops.
  attr(prior, "type") <- "factor"
  attr(prior, "Lambda_mat") <- Lambda_mat
  out <- prior
  if(sample){
    samples <- list()
    par_names <- names(sampled_pars(design, doMap = F))
    if(!selection %in% c("mu", "sigma2", "covariance", "alpha", "correlation", "Sigma", "loadings", "residuals")){
      stop("for variant factor, you can only specify the prior on the mean, variance, covariance, loadings, residuals, or the correlation of the parameters")
    }
    if(selection %in% c("mu", "alpha")){
      mu <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                               sigma = diag(prior$theta_mu_var)))
      rownames(mu) <- par_names
      if(selection %in% c("mu")){
        samples$theta_mu <- mu
      }
    }
    if(selection %in% c("loadings", "std_loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      lambda <- array(0, dim = c(n_pars, n_factors, N))
      for(i in 1:n_factors){
        lambda[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$theta_lambda_var)))
      }
      lambda <- constrain_lambda(lambda, Lambda_mat)
      rownames(lambda) <- par_names
      if(is.null(colnames(Lambda_mat))){
        colnames(lambda) <- paste0("F", 1:n_factors)
      } else{
        colnames(lambda) <- colnames(Lambda_mat)
      }
      if(selection %in% c("loadings", "std_loadings")){
        samples$lambda <- lambda
      }
    }
    if(selection %in% c("residuals", "std_loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      residuals <- t(matrix(rgamma(n_pars*N, shape = prior$as, rate = prior$bs),
                          ncol = n_pars, byrow = T))
      rownames(residuals) <- par_names
      if(selection %in% c("residuals", "std_loadings")){
        samples$epsilon_inv <- residuals
      }
    }
    if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for(i in 1:N){
        sigma <- 1/residuals[,i]
        psi <- 1/rgamma(n_factors, prior$ap, prior$bp)
        loadings <- lambda[,,i]
        vars[,,i] <- loadings %*% diag(psi, n_factors) %*% t(loadings) + diag(sigma)
      }
      if(selection != "alpha") samples$theta_var <- vars
    }
    if(selection %in% "alpha"){
      samples$alpha <- get_alphas(mu, vars, design[[1]]$Ffactors$subjects)
    }
    out <- samples
  }
  return(out)
}

get_startpoints_factor<- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!pmwgs$nuisance)
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$prior$theta_mu_mean, sd = sqrt(pmwgs$prior$theta_mu_var))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3,diag(n_pars))
  start_psi_inv <- rep(1, pmwgs$n_factors)
  start_sig_err_inv <- rep(1, n_pars)
  start_lambda <- matrix(0, nrow = n_pars, ncol = pmwgs$n_factors)
  Lambda_mat <- attr(pmwgs, "Lambda_mat")
  start_lambda[Lambda_mat != Inf] <- Lambda_mat[Lambda_mat != Inf]
  start_eta <- matrix(0, nrow = pmwgs$n_subjects, ncol = pmwgs$n_factors)
  return(list(tmu = start_mu, tvar = start_var, lambda = start_lambda, lambda_untransf = start_lambda,
              sig_err_inv = start_sig_err_inv, psi_inv = start_psi_inv,
              eta = start_eta))
}

fill_samples_factor <- function(samples, group_level, proposals, j = 1, n_pars){
  emc_set_last_slice(samples$lambda, j, as.numeric(group_level$lambda))
  emc_set_last_slice(samples$lambda_untransf, j, as.numeric(group_level$lambda_untransf))
  emc_set_last_slice(samples$epsilon_inv, j, as.numeric(group_level$sig_err_inv))
  emc_set_last_slice(samples$psi_inv, j, as.numeric(group_level$psi_inv))
  emc_set_last_slice(samples$eta, j, as.numeric(group_level$eta))
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}

# mvtnorm::rmvnorm(n, mean, sigma) for the group draws: the native form
# (src/gibbs_factor.cpp) of its default method, the same draw on the same
# random stream without the R function's overhead, else that function from
# the stream's state before the call.
.emc_rmvnorm <- function(n, mean = NULL, sigma) {
  seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  out <- if (!is.null(seed) && is.matrix(sigma) && is.double(sigma) &&
             length(n) == 1L && (is.null(mean) || is.double(mean)))
    emc_rmvnorm_cpp(n, if (!is.null(mean)) as.numeric(mean), sigma)
  if (!is.null(out)) {
    if (!is.null(names(mean))) colnames(out) <- names(mean)
    return(out)
  }
  if (!is.null(seed)) assign(".Random.seed", seed, envir = globalenv())
  if (is.null(mean)) rmvnorm(n, sigma = sigma) else rmvnorm(n, mean, sigma)
}

# solve(a) for the group draws: the same inverse without the R function's
# overhead for a real square matrix without dimnames, else that function.
.emc_solve <- function(a) {
  out <- if (is.matrix(a) && is.double(a) && is.null(dimnames(a)) && nrow(a) == ncol(a))
    emc_solve_cpp(a)
  if (is.null(out)) solve(a) else out
}

# The group draw: the native form (src/gibbs_factor.cpp), which gives the R
# form's draw on the same random stream, else the R form from the stream's
# state before the call.
gibbs_step_factor <- function(sampler, alpha){
  seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  out <- if (!is.null(seed) && isTRUE(getOption("emc2.native_gibbs_factor", TRUE)))
    .emc_gibbs_factor_native(sampler, alpha)
  if (!is.null(out)) return(out)
  if (!is.null(seed)) assign(".Random.seed", seed, envir = globalenv())
  gibbs_step_factor_r(sampler, alpha)
}

# What the native draw takes: the state (`last`, as last_sample_factor()
# gives it), which loadings are free and the prior's vectors at full length.
# NULL where the native form does not cover the hierarchy.
.emc_gibbs_factor_inputs <- function(sampler, last) {
  prior <- sampler$prior
  free <- attr(sampler, "Lambda_mat") == Inf
  n <- sampler$n_subjects
  p <- sum(!sampler$nuisance)
  k <- sampler$n_factors
  full <- function(x, len) is.numeric(x) && length(x) %in% c(1L, len)
  if (p < 2L || !is.matrix(free) || !all(dim(free) == c(p, k)) || anyNA(free) ||
      !is.matrix(prior$theta_mu_invar) || !all(dim(prior$theta_mu_invar) == c(p, p)) ||
      length(prior$theta_mu_mean) != p || length(prior$theta_lambda_invar) != p ||
      !full(prior$as, p) || !full(prior$bs, p) || !full(prior$ap, k) || !full(prior$bp, k))
    return(NULL)
  list(state = list(eta = matrix(as.numeric(last$eta), n, k),
                    lambda = matrix(as.numeric(last$lambda), p, k),
                    mu = as.numeric(last$mu), sig_err_inv = as.numeric(last$sig_err_inv),
                    psi_inv = as.numeric(last$psi_inv)),
       free = free,
       prior = list(mu_invar = prior$theta_mu_invar,
                    mu_mean = as.numeric(prior$theta_mu_mean),
                    shape_err = rep_len(as.numeric(prior$as + n / 2), p),
                    rate_err = rep_len(as.numeric(prior$bs), p),
                    lambda_invar = as.numeric(prior$theta_lambda_invar),
                    shape_psi = rep_len(as.numeric(prior$ap + n / 2), k),
                    rate_psi = rep_len(as.numeric(prior$bp), k)))
}

# The native draw's state in the form gibbs_step_factor_r() returns, given
# the subjects' states (subjects x parameters).
.emc_gibbs_factor_pars <- function(res, alpha_t) {
  mu <- matrix(res$mu, 1L, dimnames = list(NULL, colnames(alpha_t)))
  list(tmu = mu, tvar = res$var, lambda_untransf = res$lambda,
       lambda = res$lambda_scaled, eta = res$eta,
       sig_err_inv = res$sig_err_inv, psi_inv = res$psi_inv, alpha = t(alpha_t))
}

# The native draw, or NULL where the native form does not cover the
# hierarchy or declines the draw.
.emc_gibbs_factor_native <- function(sampler, alpha) {
  n <- sampler$n_subjects
  p <- sum(!sampler$nuisance)
  if (!is.matrix(alpha) || !all(dim(alpha) == c(p, n))) return(NULL)
  inp <- .emc_gibbs_factor_inputs(sampler, last_sample_factor(sampler$samples))
  if (is.null(inp)) return(NULL)
  alpha_t <- t(alpha)
  res <- emc_gibbs_factor_cpp(alpha_t, inp$state, inp$free, inp$prior)
  if (is.null(res)) return(NULL)
  .emc_gibbs_factor_pars(res, alpha_t)
}

gibbs_step_factor_r <- function(sampler, alpha){
  # Gibbs step for group means with parameter expanded factor analysis from Ghosh & Dunson 2009
  # mu = theta_mu, var = theta_var
  last <- last_sample_factor(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # extract previous values (for ease of reading)
  alpha <- t(alpha)
  n_subjects <- sampler$n_subjects
  n_pars <- sum(!sampler$nuisance)
  n_factors <- sampler$n_factors

  # Save the original constraint matrix (with fixed values) for later rescaling.
  Lambda_constraints <- hyper$Lambda_mat

  # Create a binary indicator: TRUE for free parameters (Inf) and FALSE for fixed ones.
  Lambda_mat <- (Lambda_constraints == Inf)

  eta <- matrix(last$eta, n_subjects, n_factors)
  psi_inv <- diag(last$psi_inv, n_factors)
  sig_err_inv <- diag(last$sig_err_inv)
  lambda <- matrix(last$lambda, n_pars, n_factors)
  mu <- last$mu

  # The particle step draws alpha from the eta-marginal subject prior and
  # leaves eta untouched, so eta must be redrawn from its full conditional
  # before anything conditions on it (partially collapsed Gibbs; van Dyk &
  # Park 2008). Without this, mu below conditions on scores drawn for the
  # previous alpha.
  eta_sig <- solve(psi_inv + t(lambda) %*% sig_err_inv %*% lambda)
  eta[,] <- t(eta_sig %*% t(lambda) %*% sig_err_inv %*% t(sweep(alpha, 2, mu))) +
    rmvnorm(n_subjects, sigma = eta_sig)

  # Update mu
  mu_sig <- solve(n_subjects * sig_err_inv + prior$theta_mu_invar)
  mu_mu <- mu_sig %*% (sig_err_inv %*% colSums(alpha - eta %*% t(lambda)) +
                         prior$theta_mu_invar %*% prior$theta_mu_mean)
  mu <- rmvnorm(1, mu_mu, mu_sig)
  colnames(mu) <- colnames(alpha)

  # calculate mean-centered observations
  alphatilde <- sweep(alpha, 2, mu)

  # Update eta (latent factors)
  eta_sig <- solve(psi_inv + t(lambda) %*% sig_err_inv %*% lambda)
  eta_mu <- eta_sig %*% t(lambda) %*% sig_err_inv %*% t(alphatilde)
  eta[,] <- t(eta_mu) + rmvnorm(n_subjects, sigma = eta_sig)

  # Update sig_err (error precisions)
  sig_err_inv <- diag(rgamma(n_pars, shape = prior$as + n_subjects/2,
                             rate = prior$bs + colSums((alphatilde - eta %*% t(lambda))^2)/2))

  # Update lambda (factor loadings) for free entries only
  for (j in 1:n_pars) {
    constraint <- Lambda_mat[j,]  # TRUE for free parameters
    if(any(constraint)){  # Only update if there are free entries in row j
      etaS <- eta[, constraint, drop = FALSE]
      # Fixed loadings (e.g. an anchor of 1) are part of the row's mean, so
      # the free-coefficient regression uses the response net of them.
      y_j <- alphatilde[, j]
      if (any(!constraint))
        y_j <- y_j - eta[, !constraint, drop = FALSE] %*% lambda[j, !constraint]
      lambda_sig <- solve(sig_err_inv[j,j] * t(etaS) %*% etaS +
                            diag(prior$theta_lambda_invar[j], sum(constraint)))
      lambda_mu <- (lambda_sig * sig_err_inv[j,j]) %*% (t(etaS) %*% y_j)
      lambda[j, constraint] <- rmvnorm(1, lambda_mu, lambda_sig)
    }
  }

  # Update psi_inv (latent factor precisions)
  psi_inv[,] <- diag(rgamma(n_factors, shape = prior$ap + n_subjects/2,
                            rate = prior$bp + colSums(eta^2)/2), n_factors)

  var <- lambda %*% solve(psi_inv) %*% t(lambda) + diag(1/diag(sig_err_inv))
  # lambda is the working (parameter-expanded) matrix that the Gibbs step samples
  # and that carries the Lambda_mat constraints; lambda_scaled folds the factor
  # scales in for reporting. Only the working matrix may be fed back to the next
  # iteration -- rescaling it repeatedly drives the fixed entries to zero.
  lambda_scaled <- lambda %*% matrix(diag(sqrt(1/diag(psi_inv)), n_factors), nrow = n_factors)

  return(list(tmu = mu, tvar = var, lambda_untransf = lambda,
              lambda = lambda_scaled, eta = eta,
              sig_err_inv = diag(sig_err_inv), psi_inv = diag(psi_inv), alpha = t(alpha)))
}

last_sample_factor <- function(store) {
  list(
    mu = store$theta_mu[, store$idx],
    eta = matrix(store$eta[, , store$idx], nrow = dim(store$eta)[1L], ncol = dim(store$eta)[2L],
                 dimnames = dimnames(store$eta)[1:2]),
    lambda = matrix(store$lambda_untransf[, , store$idx], nrow = dim(store$lambda_untransf)[1L],
                    ncol = dim(store$lambda_untransf)[2L],
                    dimnames = dimnames(store$lambda_untransf)[1:2]),
    psi_inv = store$psi_inv[,store$idx],
    sig_err_inv = store$epsilon_inv[,store$idx]
  )
}

get_conditionals_factor <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  .emc_mean_conditional(s, samples, n_pars, iteration, idx)
}

filtered_samples_factor <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    lambda_untransf = sampler$samples$lambda_untransf[, , filter, drop = F],
    psi_inv = sampler$samples$psi_inv[, filter, drop = F],
    epsilon_inv = sampler$samples$epsilon_inv[, filter],
    eta = sampler$samples$eta[, , filter, drop = F],
    theta_var = sampler$samples$theta_var[,,filter],
    alpha = sampler$samples$alpha[, , filter],
    Lambda_mat = attributes(sampler)$Lambda_mat,
    n_factors = sampler$n_factors,
    iteration = length(filter)
  )
}
unwind_lambda <- function(lambda, constraintMat, reverse = F){
  if(reverse){
    out <- constraintMat
    out[constraintMat == Inf] <- lambda
  } else{
    out <- as.numeric(lambda[constraintMat == Inf])
  }
  return(out)
}

constrain_lambda <- function(lambda, constraintMat){
  for(i in 1:dim(lambda)[3]){
    tmp <- lambda[,,i]
    tmp[constraintMat != Inf] <- 0
    lambda[,,i] <- tmp
  }
  return(lambda)
}

# bridge_sampling ---------------------------------------------------------

bridge_add_info_factor <- function(info, samples){
  info$n_factors <- samples$n_factors
  info$Lambda_mat <- attr(samples, "Lambda_mat")
  # Free loadings + group-level mean + factor variances + residual variances
  # note not factor scores, since we use the marginal model as adviced by Merkle et al. 2023
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):(samples$n_pars*samples$n_subjects + sum(info$Lambda_mat == Inf) + samples$n_pars + samples$n_factors + samples$n_pars)
  return(info)
}

bridge_add_group_factor <- function(all_samples, samples, idx){
  Lambda_mat <- attr(samples, "Lambda_mat")
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[, idx, drop = FALSE]))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$lambda_untransf[,,idx,drop = F], 3, unwind_lambda, Lambda_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(log(samples$samples$epsilon_inv[, idx, drop = FALSE])))
  all_samples <- cbind(all_samples, t(log(samples$samples$psi_inv[,idx, drop = F])))
  return(all_samples)
}

bridge_group_and_prior_and_jac_factor <- function(proposals_group, proposals_list, info){
  prior <- info$prior
  proposals <- do.call(cbind, proposals_list)
  theta_mu <- proposals_group[, seq_len(info$n_pars), drop = FALSE]
  lambda <- proposals_group[, (info$n_pars + 1):(info$n_pars + sum(info$Lambda_mat == Inf)), drop = FALSE]
  theta_epsilon_inv <- proposals_group[,(1 + info$n_pars + sum(info$Lambda_mat == Inf)): (info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars), drop = FALSE]
  psi_inv <- proposals_group[,(1 + info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars):
                                     (info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars + info$n_factors), drop = F]

  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)
  for(i in 1:n_iter){ # Reconstruct constrained matrices for each draw.
    lambda_curr <- unwind_lambda(lambda[i,], info$Lambda_mat, reverse = T)
    epsilon_curr <- diag(1/exp(theta_epsilon_inv[i,]))
    psi_curr <- diag(1/exp(psi_inv[i,]), info$n_factors)
    theta_var_curr <- lambda_curr %*% psi_curr %*% t(lambda_curr) + epsilon_curr
    proposals_curr <- matrix(proposals[i,], ncol = info$n_pars, byrow = T)
    group_ll <- sum(dmvnorm(proposals_curr, theta_mu[i,], theta_var_curr, log = T))
    # Gamma priors on the precisions, which are sampled on the log scale.
    prior_epsilon <- sum(dgamma(exp(theta_epsilon_inv[i,]), shape = prior$as, rate = prior$bs, log = TRUE))
    prior_psi <- sum(dgamma(exp(psi_inv[i,]), shape = prior$ap, rate = prior$bp, log = TRUE))
    sum_out[i] <- group_ll + prior_epsilon + prior_psi
  }

  prior_lambda <- dmvnorm(lambda, mean = rep(0, ncol(lambda)),
                          sigma = diag(prior$theta_lambda_var, ncol(lambda)), log = T)
  prior_mu <- dmvnorm(theta_mu, mean = prior$theta_mu_mean, sigma = diag(prior$theta_mu_var), log =T)

  jac_epsilon <- rowSums(theta_epsilon_inv)
  jac_psi <- rowSums(psi_inv)
  return(sum_out + prior_mu + prior_lambda + jac_epsilon + jac_psi) # Output is of length nrow(proposals)
}
