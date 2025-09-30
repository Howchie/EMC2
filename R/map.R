do_transform <- function(pars, model)
{
  pars=to_matrix(pars)
  transform=model$transform
  ptypes <- get_p_types(colnames(pars))
  isexp    <- transform$func[ptypes] == "exp"
  isprobit <- transform$func[ptypes] == "pnorm"
  islogis <- transform$func[ptypes] == "plogit"

  ## exp link:  lower + exp(real)
  pars[, isexp] <- sweep(
    exp(pars[, isexp, drop = FALSE]), 2,
    transform$lower[ptypes[isexp]], "+")
  
  
  ## probit link: lower + (upper lower) * pnorm(real)
  pars[, isprobit] <- sweep(
    sweep(pnorm(pars[, isprobit, drop = FALSE]), 2,
          transform$upper[ptypes[isprobit]] -
            transform$lower[ptypes[isprobit]], "*"),
    2, transform$lower[ptypes[isprobit]], "+")
  
  ## logis link: lower + (upper lower) * plogis(real)
  pars[, islogis] <- sweep(
    sweep(plogis(pars[, islogis, drop = FALSE]), 2,
          transform$upper[ptypes[islogis]] -
            transform$lower[ptypes[islogis]], "*"),
    2, transform$lower[ptypes[islogis]], "+")
  pars
}
#' @export
#' 
do_reverse_transform <- function(pars, model)
{
  pars=to_matrix(pars)
  transform=model$transform
  ptypes <- get_p_types(colnames(pars))
  islog    <- transform$func[ptypes] == "exp"
  isqnorm <- transform$func[ptypes] == "pnorm"
  isqlogit <- transform$func[ptypes] == "plogit"
  ## Log link:
  # residual on the natural scale
  pars[, islog] <- mapply(
    function(x, b) ifelse(x <= b | is.na(x), log(b), log(x)),
    as.data.frame(pars[, islog, drop = FALSE]),
    as.numeric(model$bound$minmax[1, ptypes[islog]])
  )
  
  # Probit link:
  # Get lower and upper bounds per column
  lower <- model$bound$minmax[1, ptypes[isqnorm]]
  upper <- model$bound$minmax[2, ptypes[isqnorm]]
  
  # Clamp the values in pars to [lower, upper] per column
  clamped <- sweep(
    pmin(pmax(pars[, isqnorm, drop = FALSE], rep(lower, each = nrow(pars))),
         rep(upper, each = nrow(pars))),
    2, lower, "-"
  )
  scaled <- sweep(clamped, 2, upper - lower, "/")
  pars[, isqnorm] <- qnorm(scaled)
  
  # Logis link:
  # Get lower and upper bounds per column
  lower <- model$bound$minmax[1, ptypes[isqlogit]]
  upper <- model$bound$minmax[2, ptypes[isqlogit]]
  
  # Clamp the values in pars to [lower, upper] per column
  clamped <- sweep(
    pmin(pmax(pars[, isqlogit, drop = FALSE], rep(lower, each = nrow(pars))),
         rep(upper, each = nrow(pars))),
    2, lower, "-"
  )
  scaled <- sweep(clamped, 2, upper - lower, "/")
  pars[, isqlogit] <- qlogis(scaled)
  
  pars

}
#' @export
#' 
do_reverse_transform_variance <- function(mu_nat, var, model, prop=TRUE)
{ # input both mu and variance on the natural scale
  # doesn't handle a covariance matrix yet
  if (!prop) {var_prop = var/mu_nat} else {var_prop=var} # defaults to specifying variance as a proportion of the mean on natural scale but also allows for natural scale mu/var
  vec_input <- is.null(dim(mu_nat))            # remember original shape
  mu_nat    <- to_matrix(mu_nat)
  if(is.null(dim(var_prop)) & length(var_prop==ncol(mu_nat))){
    var_prop = diag(var, ncol(mu_nat))
  } else {var_prop  <- to_matrix(var_prop)}
  if(ncol(var_prop)==1) {var_prop=rep(var_prop,ncol(mu_nat))}
  ## -- pre-compute flags ------------------------------------------------------
  ptypes   <- get_p_types(colnames(mu_nat))
  trf      <- model$transform$func[ptypes]
  is_log   <- trf == "exp"
  is_qnorm <- trf == "pnorm"
  is_qlogis<- trf == "plogis"
  is_nat   <- trf == "identity"
  ## -- allocate outputs -------------------------------------------------------
  var_tr  <- matrix(NA,  ncol = ncol(mu_nat), nrow = nrow(var_prop))
  par_tr  <- matrix(NA, ncol = ncol(mu_nat), nrow = nrow(mu_nat))
  diag(var_tr)[is_nat] <- diag(var_prop)[is_nat]
  par_tr[, is_nat] <- mu_nat[, is_nat, drop = FALSE]
  ## exp link:  lower + exp(real)
  # residual on the natural scale
  diag(var_tr)[is_log] <- log1p(diag(var_prop)[is_log] /
                              mu_nat[1:nrow(mu_nat), is_log])
  par_tr[1:nrow(mu_nat), is_log] <- ifelse(
    mu_nat[, is_log, drop = FALSE] <= model$bound$minmax[1, ptypes[is_log]] |
      is.na(mu_nat[, is_log, drop = FALSE]),
    log(model$bound$minmax[1, ptypes[is_log]]),
    log(mu_nat[, is_log, drop = FALSE]) - 0.5 * diag(var_tr)[is_log]
  )
  
  ## probit link - to get variance we need a root finder (this probably needs to be thoroughly checked by someone not ZH)
  var_nat=diag(as.vector(mu_nat),ncol=ncol(mu_nat))*var_prop
  make_root <- function(mu, target_var, lower, upper) {
    #target_var <- pmin(target_var, 0.999 * mu * (1 - mu))
    W <- upper - lower
    p <- (mu - lower) / W
    if (any(p <= 0) || any(p >= 1))          stop("mean outside (lower, upper)")
    vmax <- W^2 * p * (1 - p)
    if (target_var > vmax)          stop(paste("target variance for some parameters unattainable"))
    
    a <- qnorm(p)                   # STANDARD-NORMAL CUT-POINT, nothing fancy
    function(sig2) {
      rho <- sig2 / (1 + sig2)      
      W^2 * (pbivnorm::pbivnorm(a, a, rho) - p^2) - target_var
    }
  }
  for (j in 1:ncol(var_nat)) {
    if(is_qnorm[j]) {
      lower=model$bound$minmax[1,ptypes[j]]
      upper=model$bound$minmax[2,ptypes[j]]
      W <- upper - lower
      p <- (mu_nat[j] - lower) / W
      if (any(p <= 0) || any(p >= 1))          stop("mean outside (lower, upper)")
      a <- qnorm(p)                   # STANDARD-NORMAL CUT-POINT, nothing fancy
      

      if(!var_nat[j,j]==0) {
        var_tr[j, j] <- uniroot(make_root(mu_nat[j], var_nat[j,j], 
                                                     model$bound$minmax[1,ptypes[j]], model$bound$minmax[2,ptypes[j]]),
                                           c(-1e-6, 50))$root
        par_tr[1, j] = a*sqrt(1+var_tr[j,j])
      }
    }
  }
  
  # Logis has no analytic solution so combine integration and root finding
  # Mean of plogis(N(mu, sigma^2))
  logit_mean <- function(mu, sigma) {
    integrand <- function(x) plogis(x) * dnorm(x, mean = mu, sd = sigma)
    integrate(integrand, -20, 20)$value
  }
  
  # Variance of plogis(N(mu, sigma^2))
  logit_var <- function(mu, sigma) {
    m <- logit_mean(mu, sigma)
    integrand <- function(x) plogis(x)^2 * dnorm(x, mean = mu, sd = sigma)
    integrate(integrand, -20, 20)$value - m^2
  }
  find_mu_logit <- function(target_mean, sigma) {
    uniroot(function(mu) logit_mean(mu, sigma) - target_mean, c(-10, 10))$root
  }
  
  find_sigma_logit <- function(target_var, mu) {
    uniroot(function(sigma) logit_var(mu, sigma) - target_var, c(1e-6, 10))$root
  }
  for (j in 1:ncol(var_nat)) {
    if (is_qlogis[j]) {
      lower = model$bound$minmax[1, ptypes[j]]
      upper = model$bound$minmax[2, ptypes[j]]
      W <- upper - lower
      target_mean = (mu_nat[j] - lower) / W
      target_var = var_nat[1, j] / W^2
      # Step 1: Fix sigma (start with 1), find mu such that mean matches
      sigma0 = 1
      mu0 = find_mu_logit(target_mean, sigma0)
      
      # Step 2: Now root-find sigma for the target variance
      sigma = find_sigma_logit(target_var, mu0)
      
      # Step 3: Re-solve for mu at the found sigma
      mu = find_mu_logit(target_mean, sigma)
      
      var_tr[j, j] = sigma^2
      par_tr[1, j] = mu
      }
  }
  
out = data.frame(pars=par_tr[1,],var=diag(var_tr),row.names = colnames(mu_nat))
}

do_pre_transform <- function(p_vector, model)
{
  transform=model$pre_transform
  ptypes <- names(p_vector)
  isexp    <- transform$func[ptypes] == "exp"
  isprobit <- transform$func[ptypes] == "pnorm"
  
  ## exp link
  p_vector[isexp] <- transform$lower[ptypes[isexp]] + exp(p_vector[isexp])
  
  ## probit link
  p_vector[isprobit] <- transform$lower[ptypes[isprobit]] +
    (transform$upper[ptypes[isprobit]] -
       transform$lower[ptypes[isprobit]]) *
    pnorm(p_vector[isprobit])
  p_vector
}


# This form used in random number generation
do_bound <- function(pars,bound, lR = NULL) {
  pars=to_matrix(pars)
  tpars <- t(pars[,colnames(bound$minmax),drop=FALSE])
  ok <- tpars > bound$minmax[1,] & tpars < bound$minmax[2,]
  if (!is.null(bound$exception)) ok[names(bound$exception),] <-
    ok[names(bound$exception),] |
    (tpars[names(bound$exception),] == bound$exception)
  bound <- colSums(ok) == nrow(ok)
  if(!is.null(lR)){
    lvl <- length(unique(lR))
    bound <- rep(colSums(matrix(bound, lvl)) == lvl, each = lvl)
  }
  return(bound)
}

# This form used in get_pars
add_bound <- function(pars,bound, lR = NULL) {
  pars=to_matrix(pars)
  attr(pars, "ok") <- do_bound(pars,bound, lR = lR)
  pars
}


#### Functions to look at parameters ----

#### Functions to look at parameters ----

map_p <- function(p,dadm,model)
  # Map p to dadm and returns matrix of mapped parameters
  # p is either a vector or a matrix (ncol = number of subjects) of p_vectors
  # dadm is a design matrix with attributes containing model information
{
  
  # Check if p is a matrix and validate column names match parameter names
  if ( is.matrix(p) ) {
    if (!all(sort(dimnames(p)[[2]])==sort(attr(dadm,"p_names"))))
      stop("p col.names must be: ",paste(attr(dadm,"p_names"),collapse=", "))
    if (!all(levels(dadm$subjects) %in% dimnames(p)[[1]]))
      stop("p must have rows named for every subject in dadm")
    p <- p[dadm$subjects,,drop=FALSE]
  } else if (!all(sort(names(p))==sort(attr(dadm,"p_names")))) # If p is vector, check names
    stop("p names must be: ",paste(attr(dadm,"p_names"),collapse=", "))
  
  # Get parameter names from model and create output matrix
  do_p <- names(model$p_types)
  pars <- matrix(nrow=nrow(dadm),ncol=length(do_p),dimnames=list(NULL,do_p))
  
  # If there are any trends do these first, they might be used later in mapping
  # Otherwise we're not applying the trend premap, but we are doing it pre-transform
  # So these trend parameters are post-map, pre-transform and have to be included in the pars output
  premap_idx <- rep(F, length(do_p))
  if(!is.null(model$trend) &&
     (attr(model$trend, "premap") || attr(model$trend, "pretransform"))){
    trend_names <- get_trend_pnames(model$trend)
    pretrend_idx <- do_p %in% trend_names
    if((attr(model$trend, "premap"))){
      # These can be removed from the pars matrix at the end
      # Since they are already used before the mapping
      premap_idx <- pretrend_idx
    }
    # Reorder parameters to make design matrix for trends first
    do_p <- c(do_p[pretrend_idx], do_p[!pretrend_idx])
  } else{
    pretrend_idx <- rep(F, length(do_p))
  }
  k <- 1
  # Loop through each parameter
  for (i in do_p) {
    cur_design <- attr(dadm,"designs")[[i]]
    # Handle vector vs matrix input differently
    if ( !is.matrix(p) ) {
      pm <- t(as.matrix(p[colnames(cur_design)]))
      pm <- pm[rep(1,nrow(pars)),,drop=FALSE]
    } else pm <- p[,colnames(cur_design),drop=FALSE]
    
    # Apply pre-mapped trends if they exist
    if (!is.null(model$trend) && attr(model$trend, "premap")) {
      trend <- model$trend
      isin <- names(trend) %in% colnames(pm)
      if (any(isin)){ # At this point the trend has already been mapped and transformed
        for (j in names(trend)[isin]) {
          cur_trend <- trend[[j]]
          # We can select the trend pars from the already update pars matrix
          trend_pars <- pars[,cur_trend$trend_pnames]
          pm[,j] <- run_trend(dadm, cur_trend, pm[,j], trend_pars)
        }
      }
    }
    
    # Apply design matrix and sum parameter effects
    tmp <- pm*cur_design[attr(cur_design,"expand"),,drop=FALSE]
    tmp[is.nan(tmp)] <- 0 # Handle 0 weight x Inf parameter cases
    tmp <- apply(tmp,1,sum)
    # If this is a premap trend parameter, transform it here already
    # We'll need it transformed later in this loop (for trending other parameters)
    if(k <= sum(pretrend_idx)){
      tmp <- as.matrix(tmp)
      colnames(tmp) <- i
      tmp <- do_transform(tmp, model)
    }
    k <- k + 1
    pars[,i] <- tmp
  }
  # Return only non-trend parameters
  return(pars[,!premap_idx,drop=FALSE])
}


get_pars_matrix <- function(p_vector,dadm,model) {
  # Order:
  # 1 pretransform
  # 2 add constants
  # 3 map
  # # - if premap trend:
  # #   First make trend pars matrix
  # #   Apply trends premap and remove trend pars from pars matrix
  # # - map
  # 4 if pretransform trend:
  # # - apply trends and remove trend pars from pars matrix
  # 5 transform
  # 6 if posttransform trend:
  # # - apply trends and remove trend pars from pars matrix
  # 7 trial-wise transform
  # 8 bound
  
  # Niek should constants be included in pre_transform? I think not?
  if (is.null(dim(p_vector))) { # Check if pars is a vector
    p_vector <- do_pre_transform(p_vector, model)
  } else{p_vector <- t(apply(p_vector, 1, do_pre_transform, model))}
  # If there's any premap trends, they're done in map_p
  pars <- map_p(add_constants(p_vector,attr(dadm,"constants")),dadm, model)
  if(!is.null(model$trend) && attr(model$trend, "pretransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(dadm, model$trend, pars)
  }
  pars <- do_transform(pars, model)
  if(!is.null(model$trend) && attr(model$trend, "posttransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(dadm, model$trend, pars)
  }
  pars <- model$Ttransform(pars, dadm)
  pars <- add_bound(pars, model$bound, dadm$lR)
  return(pars)
}


make_pmat <- function(p_vector,design)
  # puts vector form of p_vector into matrix form
{
  ss <- design$Ffactors$subjects
  out <- matrix(rep(p_vector,each=length(ss)),nrow=length(ss),
                dimnames=list(ss,names(p_vector)))
  if(is.null(colnames(out))){
    colnames(out) <- names(sampled_pars(design))
  }
  return(out)
}

add_constants <- function(p,constants)
{
  if (is.null(constants)) return(p)
  if (is.matrix(p)) {
    nams <- c(dimnames(p)[[2]],names(constants))
    p <- cbind(p,matrix(rep(constants,each=dim(p)[1]),nrow=dim(p)[1]))
    dimnames(p)[[2]] <- nams
    return(p)
  } else{
    return(c(p,constants))
  }
  
}

add_recalculated_pars <- function(pmat, model, cnams){
  modifiers <- unlist(lapply(strsplit(cnams,"_"),function(x){paste0(x[-1], collapse = "_")}))
  par_names <- colnames(pmat)
  unq_pars <- unique(par_names)
  par_table <- table(par_names)
  counts <- lapply(par_table, function(x) return(1:x))
  combn <- do.call(expand.grid, counts)
  colnames(combn) <- names(par_table)
  
  out <- list()
  modfs <- list()
  for(r in 1:nrow(combn)){
    pmat_in <- matrix(NA, nrow = nrow(pmat), ncol = length(unq_pars))
    colnames(pmat_in) <- unq_pars
    cur_modifiers <- setNames(numeric(length(unq_pars)), unq_pars)
    for(par in unq_pars){
      pmat_in[,par] <- pmat[,which(colnames(pmat) == par)[combn[r,par]]]
      cur_modifiers[par] <- modifiers[which(colnames(pmat) == par)[combn[r,par]]]
    }
    added <- model()$Ttransform(pmat_in)
    added <- added[, !(colnames(added) %in% colnames(pmat_in)), drop = F]
    attr(added, "ok") <- NULL
    out[[r]] <- added
    modfs[[r]] <- cur_modifiers
  }
  m_out <- matrix(0, nrow = nrow(pmat), ncol = 0)
  if(ncol(out[[1]]) == 0) return(NULL)
  for(i in 1:ncol(out[[1]])){
    cur_par <- lapply(out, function(x) x[,i, drop = F])
    not_dups <- !duplicated(cur_par)
    cur_combn <- combn[not_dups,]
    pars_vary <- colnames(cur_combn)[colMeans(cur_combn) != 1]
    to_add <- do.call(cbind, cur_par[not_dups])
    cur_modfs <- modfs[not_dups]
    fnams <- sapply(cur_modfs, function(x) paste0(unique(unlist(strsplit(x[pars_vary], split = "_"))), collapse = "_"))
    if(!all(fnams == "")) colnames(to_add) <- paste0(colnames(to_add)[1], "_", fnams)
    m_out <- cbind(m_out, to_add)
  }
  return(m_out)
}

get_p_types <- function(nams, reverse = FALSE){
  if(reverse){
    out <- unlist(lapply(strsplit(nams,"_"),function(x){x[[-1]]}))
  } else{
    out <- unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
  }
  return(out)
}

fill_transform <- function(transform, model, p_vector,
                           supported=c("identity","exp","pnorm","plogis"),
                           has_lower=c("exp","pnorm","plogis"),has_upper=c("pnorm","plogis"),
                           is_pre = FALSE){
  if(!is.null(transform)){
    if (!all(transform$func %in% supported)){
      stop("Only ", paste(supported, collapse = ", "), " transforms supported")
    }
    if(!is_pre){
      if (!all(names(transform$func) %in% names(model()$p_types)))stop("transform on parameter not in the model p_types")
      if (!all(names(transform$lower) %in% names(model()$p_types)))stop("transform on parameter not in the model p_types")
      if (!all(names(transform$upper) %in% names(model()$p_types)))stop("transform on parameter not in the model p_types")
    } else{
      if (!all(names(transform$func) %in% names(p_vector))) stop("pre_transform on parameter not in the sampled_pars")
      if (!all(names(transform$lower) %in% names(p_vector))) stop("transform on parameter not in the model p_types")
      if (!all(names(transform$upper) %in% names(p_vector))) stop("transform on parameter not in the model p_types")
    }
  }
  model_list <- model()
  p_names <- names(model_list$p_types)
  if(is_pre){
    p_names <- names(p_vector)
  }
  filled_func <- filled_lower <- filled_upper <- setNames(rep(NA, length(p_names)), p_names)
  # First update from the model
  if(is_pre){
    filled_func[names(model_list$pre_transform$func)] <- model_list$pre_transform$func
    filled_lower[names(model_list$pre_transform$lower)] <- model_list$pre_transform$lower
    filled_upper[names(model_list$pre_transform$upper)] <- model_list$pre_transform$upper
  } else{
    filled_func[names(model_list$transform$func)] <- model_list$transform$func
    filled_lower[names(model_list$transform$lower)] <- model_list$transform$lower
    filled_upper[names(model_list$transform$upper)] <- model_list$transform$upper
  }
  # Now updates from the user transform (they override the model)
  filled_func[names(transform$func)] <- transform$func
  filled_lower[names(transform$lower)] <- transform$lower
  filled_upper[names(transform$upper)] <- transform$upper
  
  # Now fill in remainers
  filled_func[is.na(filled_func)] <- "identity"
  filled_lower[is.na(filled_lower) & (filled_func %in% has_lower)] <- 0
  filled_upper[is.na(filled_upper) & (filled_func %in% has_upper)] <- 1
  # Remainers must be identity and have no bounds
  filled_lower[is.na(filled_lower)] <- -Inf
  filled_upper[is.na(filled_upper)] <- Inf
  return(list(func=filled_func,lower=filled_lower,upper=filled_upper))
}

fill_bound <- function(bound, model) {
  filled_bound <- model()$bound
  if (!is.null(bound)) {
    if (names(bound)[1] != "minmax")
      stop("first entry of bound must be named minmax")
    if (!all(colnames(bound$minmax) %in% names(model()$p_types)))
      stop("minmax column names must correspond to parameter types")
    if (!is.null(bound$exception) &&
        (!all(names(bound$exception) %in% names(model()$p_types))))
      stop("exception names must correspond to parameter types")
    filled_bound$minmax[,colnames(bound$minmax)] <- bound$minmax
    if (!is.null(bound$exception)) {
      filled_bound$exception <- c(bound$exception,filled_bound$exception)
      filled_bound$exception <- filled_bound$exception[!duplicated(names(filled_bound$exception))]
    }
  }
  filled_bound
}



generate_design_equations <- function(design_matrix,
                                      factor_cols = NULL,
                                      numeric_cols = NULL,
                                      transform,
                                      pre_transform) {
  # 1. If user hasn't specified which columns are factors or numeric, guess:
  if (is.null(factor_cols)) {
    # We'll assume anything that is a factor or character is a "factor column"
    factor_cols <- names(design_matrix)[
      sapply(design_matrix, function(x) is.factor(x) || is.character(x))
    ]
  }
  if (is.null(numeric_cols)) {
    # We'll assume everything that is numeric is for the design (contrast) columns
    numeric_cols <- names(design_matrix)[
      sapply(design_matrix, is.numeric)
    ]
  }
  
  # 2. Build the "equation" strings from numeric cols
  make_equation_string <- function(row_i) {
    eq_terms <- c()
    for (colname in numeric_cols) {
      new_name <- colname
      if(colname %in% names(pre_transform)){
        cur_trans <- pre_transform[colname]
        if(cur_trans != "identity"){
          new_name <- paste0(cur_trans, "(", colname, ")")
        }
      }
      val <- as.numeric(row_i[[colname]])
      # Skip zeros
      if (abs(val) < 1e-15) next
      
      # If val is exactly +1 or -1, skip numeric part:
      if (abs(val - 1) < 1e-15) {
        # val == +1
        term_str <- paste0("+ ", new_name)
      } else if (abs(val + 1) < 1e-15) {
        # val == -1
        term_str <- paste0("- ", new_name)
      } else {
        # Some other numeric coefficient => e.g. "+ 0.5 * col"
        sign_str <- ifelse(val >= 0, "+ ", "- ")
        term_str <- paste0(sign_str, format(abs(val), digits = 3), " * ", new_name)
      }
      eq_terms <- c(eq_terms, term_str)
    }
    
    if (length(eq_terms) == 0) {
      # If everything was zero
      return("0")
    }
    
    # Combine eq_terms
    eq_string <- paste(eq_terms, collapse = " ")
    # Remove the leading '+ ' if present
    sub("^\\+ ", "", eq_string)
  }
  # 3. Compute the equation string for each row
  eq_strings <- apply(design_matrix, 1, make_equation_string)
  
  # 4. Determine the widths needed for each factor column
  #    so everything lines up in columns nicely.
  col_widths <- sapply(factor_cols, function(fc) {
    # The width must accommodate either the column name or the largest factor level
    max(nchar(fc), max(nchar(as.character(design_matrix[[fc]])), na.rm = TRUE))
  })
  
  # We'll label the final column as "Equation"
  eq_header <- ""
  # We don't necessarily need to align the entire equation column itself
  # (since it may vary in length), but let's keep a short label for the header.
  
  # 5. Print the header row
  #    E.g. if factor_cols = c("S", "E"), we do:
  #    " S     E    : Equation"
  factor_header_str <- paste(mapply(function(fc, w) {
    sprintf("%-*s", w, fc)   # left-justify the column name
  }, factor_cols, col_widths),
  collapse = "  ")
  
  cat("  ", factor_header_str, "  ", eq_header, "\n", sep="")
  
  # Optionally, add a simple "rule" line or blank line:
  # (Uncomment if you like to separate header from rows)
  # cat(" ", paste(rep("-", nchar(factor_header_str) + nchar(eq_header) + 5),
  #               collapse=""), "\n", sep="")
  
  # 6. Print each row: factor columns in their spaces, then ": <equation>"
  for (i in seq_len(nrow(design_matrix))) {
    row_factor_str <- paste(mapply(function(fc, w) {
      # Each factor value is left-justified in the same width as the header
      sprintf("%-*s", w, as.character(design_matrix[[fc]][i]))
    }, factor_cols, col_widths),
    collapse = "  ")
    if(transform == "identity"){
      cat(" ", row_factor_str, "  : ", eq_strings[i], "\n", sep="")
    } else{
      cat(" ", row_factor_str, "  : ", transform, "(", eq_strings[i], ")", "\n", sep="")
    }
    
  }
}


verbal_dm <- function(design){
  map <- attr(sampled_pars(design, add_da = TRUE, doMap = TRUE), "map")
  map_no_da <- attr(sampled_pars(design, doMap = TRUE), "map")
  transforms <- design$model()$transform$func
  pre_transforms <- design$model()$pre_transform$func
  for(i in 1:length(map)){
    m <- map[[i]]
    if(ncol(m) == 1) next
    cat(paste0("$", names(map)[i]), "\n")
    mnd <- map_no_da[[i]]
    par_idx <- colnames(m) %in% colnames(mnd)
    generate_design_equations(m, colnames(m)[!par_idx], colnames(m)[par_idx],
                              transform = transforms[names(map)[i]],
                              pre_transform = pre_transforms)
    cat("\n")
  }
}

#' Call the C++ version of `get_pars_matrix`
#'
#' This helper allows inspection of the parameter matrix generated on the C++
#' side using the same steps as `calc_ll()`.
#' @param p_vector Numeric vector of parameters.
#' @param dadm Data frame used when mapping parameters.
#' @param model Model object as produced by `design$model()`.
#' @return Numeric matrix of parameters after mapping and transformations.
#' @export
get_pars_matrix_cpp <- function(p_vector, dadm, model) {
  get_pars_matrix_rcpp(p_vector,
                       attr(dadm, "constants"),
                       model$transform,
                       model$pre_transform,
                       model$p_types,
                       model$design,
                       nrow(dadm),
                       dadm,
                       model$trend)
}


