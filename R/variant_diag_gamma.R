sample_store_diag_gamma <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, ...) {
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

fill_samples_diag_gamma <- function(samples, group_level, proposals, j = 1, n_pars){
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}


add_info_diag_gamma <- function(sampler, prior = NULL, ...){
  sampler$prior <- get_prior_diag_gamma(prior, sum(!sampler$nuisance), sample = F)
  return(sampler)
}

get_prior_diag_gamma <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL){
  # Checking and default priors
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_pars(design, doMap = F))
  }
  if (!is.null(prior$theta_mu_mean)) {
    n_pars <- length(prior$theta_mu_mean)
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars))
  }
  if(is.null(prior$shape)){
    prior$shape <- rep(2, n_pars)
  }
  if(is.null(prior$rate)){
    prior$rate <- rep(.3, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix
  attr(prior, "type") <- "diagonal-gamma"
  out <- prior
  if(sample){
    par_names <- names(sampled_pars(design, doMap = F))
    samples <- list()
    if(selection %in% c("mu", "alpha")){
      mu <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                               sigma = prior$theta_mu_var))
      rownames(mu) <- par_names
      if(selection %in% c("mu")){
        samples$theta_mu <- mu
      }
    }
    if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for(i in 1:N){
        vars[,,i] <- diag(1/ rgamma(n = n_pars, shape = prior$shape, rate = prior$rate))
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

get_startpoints_diag_gamma <- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!pmwgs$nuisance)
  if (is.null(start_mu)) start_mu <- rnorm(n_pars, mean = pmwgs$prior$theta_mu_mean, sd = sqrt(diag(pmwgs$prior$theta_mu_var)))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- diag(1/rgamma(n_pars, 10, 5), n_pars) #Bit stupid maybe as startpoint
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var)))
}

last_sample_diag_gamma <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv
  )
}


# Means from their normal conditional (independent across parameters when the
# prior precision is diagonal, joint otherwise), then each precision from its
# gamma conditional, in C++ (src/ensemble_move.cpp, the kernel the ensemble
# rounds share).
gibbs_step_diag_gamma <- function(sampler, alpha){
  last <- last_sample_diag_gamma(sampler$samples)
  n_pars <- sum(!sampler$nuisance)
  alpha <- as.matrix(alpha)
  prior <- sampler$prior
  P0 <- prior$theta_mu_invar
  if (!is.matrix(P0)) {
    prior$theta_mu_invar <- diag(as.numeric(P0), n_pars)
  }
  res <- emc_gibbs_diag_gamma_cpp(alpha, diag(last$tvinv), prior)
  tmu <- res$tmu
  names(tmu) <- sampler$par_names[!sampler$nuisance]
  return(list(tmu = tmu, tvar = diag(1 / res$tvinv, n_pars), tvinv = diag(res$tvinv, n_pars),
              alpha = alpha))
}


