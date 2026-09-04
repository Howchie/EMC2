normalize_marginalise <- function(marginalise) {
  if (is.null(marginalise)) return(NULL)
  param <- NULL
  n_nodes <- NULL

  if (is.list(marginalise)) {
    param <- marginalise$param
    if (is.null(param) && length(marginalise) > 0) param <- marginalise[[1]]
    n_nodes <- marginalise$n_nodes
    if (is.null(n_nodes) && !is.null(marginalise$nodes)) n_nodes <- marginalise$nodes
  } else if (is.character(marginalise) || is.atomic(marginalise)) {
    nm <- names(marginalise)
    if (!is.null(nm) && "param" %in% nm) {
      param <- marginalise[["param"]]
    } else if (!is.null(nm)) {
      idx_param <- which(!nm %in% c("n_nodes", "nodes"))
      if (length(idx_param) > 0) param <- marginalise[idx_param[1]] else param <- marginalise[1]
    } else {
      param <- marginalise[1]
    }
    if (!is.null(nm) && "n_nodes" %in% nm) {
      n_nodes <- as.integer(marginalise[["n_nodes"]])
    } else if (!is.null(nm) && "nodes" %in% nm) {
      n_nodes <- as.integer(marginalise[["nodes"]])
    } else if (length(marginalise) > 1 && is.null(nm)) {
      if (suppressWarnings(!is.na(as.integer(marginalise[2])))) {
        n_nodes <- as.integer(marginalise[2])
      }
    }
  } else {
    stop("marginalise must be a character vector or a list identifying a parameter name")
  }

  param <- unname(as.character(unlist(param)))[1L]
  if (is.na(param) || !nzchar(param)) {
    stop("marginalise must name at least one parameter")
  }

  res <- list(param = param)
  if (!is.null(n_nodes) && !is.na(n_nodes)) res$n_nodes <- as.integer(n_nodes)
  res
}

validate_marginalise_design <- function(marginalise, design, model) {
  norm <- normalize_marginalise(marginalise)
  if (is.null(norm)) return(NULL)
  p <- norm$param
  if (length(p) != 1L) {
    stop("marginalise currently supports exactly one shared parameter")
  }
  spec <- model()
  if (is.null(spec$type) || !grepl("RACE", spec$type, fixed = TRUE) ||
      grepl("DDM", spec$type, fixed = TRUE)) {
    stop("marginalise is currently supported for race models only (not DDM)")
  }
  if (is.null(spec$p_types) || !p %in% names(spec$p_types)) {
    stop("marginalise parameter '", p, "' is not a model p_type")
  }
  f <- design$Flist[[p]]
  if (is.null(f)) {
    matches <- vapply(design$Flist, function(z) {
      as.character(stats::terms(z))[[2L]] == p
    }, logical(1L))
    if (any(matches)) f <- design$Flist[[which(matches)[1L]]]
  }
  if (is.null(f)) {
    stop("marginalise parameter '", p, "' must be specified in the design")
  }
  rhs <- tryCatch(as.character(stats::terms(f))[[3L]], error = function(e) NULL)
  if (!identical(rhs, "1")) {
    stop("marginalise parameter '", p, "' must use an intercept-only design (~ 1)")
  }
  mm <- spec$bound$minmax
  if (is.null(mm) || is.null(colnames(mm)) || !p %in% colnames(mm) ||
      !is.finite(mm[1L, p])) {
    stop("marginalise parameter '", p, "' must have a finite lower bound")
  }
  if (!is.null(design$constants) && p %in% names(design$constants)) {
    stop("marginalise parameter '", p, "' cannot be fixed as a constant")
  }
  if (!is.null(norm$n_nodes)) {
    if (length(norm$n_nodes) != 1L || is.na(norm$n_nodes) || norm$n_nodes < 2L) {
      stop("marginalise n_nodes must be an integer >= 2")
    }
  }
  if (is.null(norm$n_nodes)) {
    return(p)
  }
  norm
}

#' Specify a Design and Model
#'
#' This function combines information regarding the data, type of model, and
#' the model specification.
#'
#' @param formula A list. Contains the design formulae in the
#' format `list(y ~ x, a ~ z)`.
#' @param factors A named list containing all the factor variables that span
#' the design cells and that should be taken into account by the model.
#' The name `subjects` must be used to indicate the participant factor variable,
#' also in the data.
#'
#' Example: `list(subjects=levels(dat$subjects), condition=levels(dat$condition))`
#'
#' @param Rlevels A character vector. Contains the response factor levels.
#' Example: `c("right", "left")`
#' @param model A function, specifies the model type.
#' Choose from the drift diffusion model (`DDM()`, `DDMt0natural()`),
#' the log-normal race model (`LNR()`), the linear ballistic model (`LBA()`),
#' the racing diffusion model (`RDM()`, `RDMt0natural()`), or define your own
#' model functions.
#' @param data A data frame. `data` can be used to automatically detect
#'  `factors`, `Rlevels` and `covariates` in a dataset. The variable `R` needs
#'  to be a factor variable indicating the response variable. Any numeric column
#'  except `trials` and `rt` are treated as covariates, and all remaining factor
#'  variables are internally used in `factors`.
#' @param contrasts Optional. A named list specifying a design matrix.
#' Example for supplying a customized design matrix:
#' `list(lM = matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))))`
#' @param matchfun A function. Only needed for race models. Specifies whether a
#' response was correct or not. Example: `function(d)d$S==d$lR` where lR refers
#' to the latent response factor.
#' @param constants A named vector that sets constants. Any parameter in
#' `sampled_pars` can be set constant.
#' @param covariates Names of numeric covariates.
#' @param functions List of functions to create new factors based on those in
#' the factors argument. These new factors can then be used in `formula`.
#' @param report_p_vector Boolean. If TRUE (default), it returns the vector of
#' parameters to be estimated.
#' @param custom_p_vector A character vector. If specified, a custom likelihood
#' function can be supplied.
#' @param trend A trend list, as made by \code{\link{make_trend}}
#' @param pre_transform_terms Optional named list keyed by model parameter. Each
#' element is a character vector naming design-matrix coefficients that should be
#' summed on the transformed scale before the model transform is applied, with
#' the remaining coefficients added linearly on the natural scale.
#' @param transform A list with custom transformations to be applied to the parameters of the model,
#' if the conventional transformations aren't desired.
#' See `DDM()` for an example of such transformations
#' @param bound A list with custom bounds to be applied to the parameters of the model,
#' if the conventional bound aren't desired.
#' see `DDM()` for an example of such bounds. Bounds are used to set limits to
#' the likelihood landscape that cannot reasonable be achieved with `transform`
#' @param marginalise Optional. Names one shared race-model parameter (e.g.
#' `"t0"`) to be *marginalised* out of the sampler by numerical quadrature
#' rather than sampled directly. This is available for group-level samplers and
#' for `type = "single"` alpha-only samplers. Either a character scalar
#' (`marginalise = "t0"`), or a list giving the parameter plus quadrature
#' controls (`marginalise = list(param = "t0", n_nodes = 12)`). `n_nodes` is the
#' number of Gauss-Legendre quadrature nodes per proposal and defaults to `12`
#' (the rule is centred on each proposal's own conditional mode, so few nodes
#' are needed). The named
#' parameter must (i) belong to a race model (not `DDM()`), (ii) be a sampled
#' model parameter (not a `constant`, and not a `custom_p_vector` design),
#' (iii) use an intercept-only design (`~ 1`), and (iv) have a finite lower
#' bound. See Details. This is a sampling/efficiency feature and does not change
#' the model's likelihood.
#' @param ... Additional, optional arguments
#'
#' @details
#' # Marginalising a shared parameter (`marginalise`)
#'
#' Some race-model parameters that are shared across a subject's accumulators —
#' most notably the non-decision time `t0` — trade off strongly against the
#' race-speed parameters (thresholds/rates). In go/no-go race designs this
#' trade-off produces a ridge in the posterior that makes the MCMC chains mix
#' poorly (a highly correlated, ill-conditioned proposal geometry). `marginalise`
#' addresses this by integrating the shared parameter out of the per-subject
#' likelihood instead of sampling it directly.
#'
#' When `marginalise` is set, that parameter is held out of the MCMC proposal
#' vector. Its per-subject likelihood contribution is computed by Gauss-Legendre
#' quadrature over the parameter (on its sampled/log scale) using `n_nodes`
#' nodes, weighted by a fixed Gaussian `eta` taken from the parameter's existing
#' prior mean and variance (its group-level distribution is pinned to that prior,
#' not estimated). For `type = "single"`, the same fixed prior is used directly;
#' no `theta_mu` or `theta_var` samples are created. A representative value is
#' drawn back (reconstructed) from the quadrature grid and stored in each
#' posterior sample, so downstream summaries still report the parameter. The
#' integration interval is clipped to the parameter's finite lower bound and to
#' just below the subject's fastest response time, which is why a finite lower
#' bound and the `~ 1` design are required.
#'
#' Marginalising does not alter the likelihood: turning it off yields an
#' identical model, only sampled with the parameter left in. It typically
#' improves mixing/ESS per iteration at the cost of extra quadrature work per
#' likelihood evaluation. Currently exactly one parameter may be marginalised,
#' and only for race models.
#'
#' @return A design list.
#' @examples
#'
#' # load example dataset
#' dat <- forstmann
#'
#' # create a function that takes the latent response (lR) factor (d) and returns a logical
#' # defining the correct response for each stimulus. Here the match is simply
#' # such that the S factor equals the latent response factor
#' matchfun <- function(d)d$S==d$lR
#'
#' # When working with lM and lR, it can be useful to design  an
#' # "average and difference" contrast matrix. For binary responses, it has a
#' # simple canonical form
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
#'
#' # Create a design for a linear ballistic accumulator model (LBA) that allows
#' # thresholds to be a function of E and lR. The final result is a 9 parameter model.
#' design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
#'                             formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#'                             contrasts=list(v=list(lM=ADmat)),
#'                             constants=c(sv=log(1)))
#' @export
#'
#'
design <- function(formula = NULL,factors = NULL,Rlevels = NULL,model,data=NULL,
                   contrasts=NULL,matchfun=NULL,constants=NULL,covariates=NULL,
                   functions=NULL,report_p_vector=TRUE, custom_p_vector = NULL,
                   trend=NULL,
                   pre_transform_terms = NULL,
                   transform = NULL, bound = NULL, TC = NULL,
                   LT=NULL,LC=NULL,UC=NULL,UT=NULL,
                   fixed_accumulator_roles = NULL, marginalise = NULL,...){

  TC <- check_missing(TC, data = data)
  if (!is.null(LT)) TC$LT <- LT
  if (!is.null(LC)) TC$LC <- LC
  if (!is.null(UC)) TC$UC <- UC
  if (!is.null(UT)) TC$UT <- UT

  optionals <- list(...)
  if (is.list(model) && !is.function(model)) {
    model_list <- model
    model <- function(){return(model_list)}
  }
  if (!is.function(model)) {
    stop("model must be a function or a model list (e.g., LBA())")
  }
  model_spec <- tryCatch(model(), error = function(e) NULL)
  if (is.list(model_spec) && isTRUE(model_spec$correlated) &&
      !identical(model_spec$correlation_type, "rdmswtn_gaussian_copula") &&
      is.null(matchfun)) {
    stop("BAwLcorr requires matchfun so it can construct the lM role indicator.")
  }
  if(!is.null(optionals$pre_transform)){
    pre_transform <- optionals$pre_transform
  } else {
    pre_transform <- NULL
  }

  if(any(names(factors) %in% c("trial", "R", "rt", "lR", "lM","UC","LC","UT","LT","winner"))){
    stop("Please do not use any of these factor names: winner, trial, R, rt, lR, lM, UC, LC, UT, LT")
  }
  if (!"subjects" %in% names(factors)) factors$subjects <- 1 # Add subjects when data are not supplied.
  factors <- setNames(
    lapply(factors, function(x)
      if (!is.factor(x)) factor(x, levels = unique(x)) else x
    ),
    names(factors)
  ) # factors were being passed even if not factor types, which is inconsistent with the case where data is passed
  if(any(grepl("_", names(factors)))){
    stop("_ in variable names detected. Please refrain from using any underscores.")
  }

  if(!is.null(custom_p_vector)){

    if (!is.null(marginalise)) {
      stop("marginalise is not supported with custom_p_vector designs")
    }

    model_list <- function(){list(log_likelihood = model)}
    if(!is.null(list(...)$rfun)){
      model_list$rfun <- list(...)$rfun
    }
    design <- list(Flist = formula, model = model_list, Ffactors = factors)

    attr(design, "sampled_p_names") <-custom_p_vector
    attr(design, "custom_ll") <- TRUE
    class(design) <- "emc.design"
    return(design)
  }
  if (!is.null(data)) {
    if(!"R" %in% colnames(data)) stop("make sure R is specified in data")
    if(!"subjects" %in% colnames(data)) stop("make sure subjects identifier is present in data")
    data$subjects <- factor(data$subjects)
    facs <- lapply(data,levels)
    nfacs <- facs[unlist(lapply(facs,is.null))]
    facs <- facs[!unlist(lapply(facs,is.null))]
    Rlevels <- facs[["R"]]
    factors <- facs[names(facs)!="R"]
    nfacs <- nfacs[!(names(nfacs) %in% c("trials","rt"))]
    all_preds <- unlist(lapply(lapply(formula, `[[`, 3L), all.vars))
    
    # Identify variables needed from matchfun
    match_vars <- if (!is.null(matchfun)) all.vars(body(matchfun)) else NULL
    
    # Identify variables needed from user-defined factor functions
    function_vars <- if (!is.null(functions)) {
      unique(unlist(lapply(functions, function(fun) {
        setdiff(all.vars(body(fun)), names(formals(fun)))
      })))
    } else NULL

    # Identify variables needed from trends
    trend_vars <- if (!is.null(trend)) {
      unique(c(
        unlist(lapply(trend, function(x) x$covariate)),
        unlist(lapply(trend, function(x) x$at))
      ))
    } else NULL

    # Map functions receive the full accumulator data frame and may use
    # columns that are not themselves trend covariates (for example,
    # `cov_left` and `cov_right` in the documented accumulator-map example).
    # Keep any data columns referenced by a map when rebuilding the minimal
    # design for sampled_pars()/mapped_pars().
    map_vars <- if (!is.null(trend)) {
      unique(unlist(lapply(trend, function(x) {
        if (is.null(x$map)) return(character())
        unique(unlist(lapply(x$map, function(fun) {
          intersect(all.names(body(fun), functions = TRUE), names(data))
        })))
      })))
    } else NULL

    # Required factors: used in formula, matchfun, functions, or trend, plus 'subjects'
    needed_factors <- unique(c(all_preds, match_vars, function_vars, trend_vars,
                               map_vars, "subjects"))
    factors <- factors[names(factors) %in% needed_factors]

    if (length(nfacs)>0){
      covariates <- names(nfacs)
      # Covariates must be in the set of all predictors, specifically used in
      # trends, or required by an accumulator map. Map inputs can be character
      # columns (for example, `cov_left` and `cov_right`), so they are not
      # represented in `factors` above.
      covariates <- covariates[covariates %in%
                               c(all_preds, function_vars, trend_vars, map_vars)]
      if(length(covariates) == 0) covariates <- NULL
    }
  } else {if(is.null(Rlevels)) stop("make sure Rlevels is specified")}
  if (!is.null(trend)) {
    formula <- check_trend(trend,c(names(functions), covariates), model, formula)
  }
  
  ## Handle GNG models silently
  # Select the GNG likelihood variant when "nogo" is among the response levels.
  if ("nogo" %in% Rlevels) {
    m_list <- model()
    if(!("GNG"%in%m_list$type)){
      if(!grepl("GNG", m_list$type)) m_list$type <- paste0(m_list$type, "GNG")
      if (!is.null(m_list$c_name)) {if(!grepl("GNG", m_list$c_name)) m_list$c_name <- paste0(m_list$c_name, "GNG")}
      model <- function() m_list
    }
  }
  # Check if all parameters in the model are specified in the formula
  nams <- unlist(lapply(formula,function(x) as.character(stats::terms(x)[[2]])))
  if (!all(sort(names(model()$p_types)) %in% sort(nams)) & is.null(custom_p_vector)){
    p_types <- model()$p_types
    not_specified <- sort(names(p_types))[!sort(names(p_types)) %in% sort(nams)]
    canonical <- model()$p_types_canonical
    warn_pars <- if (!is.null(canonical)) intersect(not_specified, canonical) else not_specified
    if (length(warn_pars) > 0)
      message(paste0("Parameter(s) ", paste0(warn_pars, collapse = ", "), " not specified in formula and assumed constant."))
    additional_constants <- p_types[not_specified]
    names(additional_constants) <- not_specified
    constants <- c(constants, additional_constants[!names(additional_constants) %in% names(constants)])
    for(add_constant in not_specified) formula[[length(formula)+ 1]] <- as.formula(paste0(add_constant, "~ 1"))
  }


  # pGuess is applied only in the compiled likelihoods (decision 8 of the
  # guess-contaminant design): the R reference paths deliberately still
  # implement pContaminant alone.  A model with no c_name dispatches to the R
  # path, where a free pGuess would be sampled and silently ignored -- a wrong
  # answer with no symptom.  Turn that into a loud one.
  if ("pGuess" %in% names(model()$p_types) && is.null(model()$c_name)) {
    pg_const <- if ("pGuess" %in% names(constants)) constants[["pGuess"]] else NULL
    pg_free <- ("pGuess" %in% nams) && is.null(pg_const)
    if (pg_free || (!is.null(pg_const) && is.finite(pg_const)))
      stop("pGuess requires a compiled likelihood, but this model has no c_name.")
  }

  design <- list(Flist=formula,Ffactors=factors,Rlevels=Rlevels,
                 Clist=contrasts,matchfun=matchfun,constants=constants,
                 Fcovariates=covariates,Ffunctions=functions,model=model,
                 TC=TC,LT=TC$LT,LC=TC$LC,UC=TC$UC,UT=TC$UT,
                 fixed_accumulator_roles = fixed_accumulator_roles)
  class(design) <- "emc.design"
  # Preserve the data context when the design is constructed from data.  This
  # allows sampled_pars(design) and mapped_pars(design) to use observed cells
  # without requiring the caller to pass the same data again.
  if (!is.null(data)) attr(design, "data") <- data
  if (!is.null(trend)) {
    # check for at = 'lR'
    if(any(sapply(trend, function(x) x$at)=='lR') & model()$type!='RACE') {
      warning('A trend has `at="lR"`, but this is not a race model. Setting `at` to NULL')
      for(i in 1:length(trend)) if(trend[[i]]$at=='lR') trend[[i]]$at <- NULL
    }
    model <- update_model_trend(trend, model)
    model_list <- model()
    model <- function(){return(model_list)}
  }
  p_vector <- sampled_pars(design)
  lhs_terms <- unlist(lapply(formula, function(x) as.character(stats::terms(x)[[2]])))

  # Check if any terms are not in model parameters
  if (!is.null(formula) && !all(lhs_terms %in% names(model()$p_types))) {
    invalid_terms <- lhs_terms[!lhs_terms %in% names(model()$p_types)]
    stop(paste0("Parameter(s) ", paste0(invalid_terms, collapse=", "),
                " in formula not found in model p_types"))
  }
  model_list <- model()
  model_list$transform <- fill_transform(transform,model)
  model_list$transform$pre_sum_terms <- validate_pre_transform_terms(
    pre_transform_terms = pre_transform_terms,
    design = design,
    model = model
  )
  model_list$bound <- fill_bound(bound,model)
  model_list$pre_transform <- fill_transform(pre_transform, model = model, p_vector = p_vector, is_pre = TRUE)
  model <- function(){return(model_list)}
  design$model <- model
  attr(design,"p_vector") <- p_vector
  if (!is.null(marginalise)) {
    attr(design, "marginalise") <- validate_marginalise_design(marginalise, design, model)
  }
  if (report_p_vector) {
    summary(design, data = data)
  }
  return(design)
}

validate_pre_transform_terms <- function(pre_transform_terms, design, model) {
  if (is.null(pre_transform_terms)) return(NULL)
  if (!is.list(pre_transform_terms) || is.null(names(pre_transform_terms)) ||
      any(names(pre_transform_terms) == "")) {
    stop("pre_transform_terms must be a named list keyed by model parameter")
  }

  model_pars <- names(model()$p_types)
  if (!all(names(pre_transform_terms) %in% model_pars)) {
    bad <- setdiff(names(pre_transform_terms), model_pars)
    stop("pre_transform_terms has parameter(s) not in model p_types: ",
         paste(bad, collapse = ", "))
  }

  min_design <- minimal_design(
    design,
    drop_subjects = FALSE,
    drop_R = FALSE,
    add_acc = TRUE,
    do_functions = FALSE,
    verbose = FALSE
  )
  dadm <- design_model(
    min_design,
    design,
    model,
    add_acc = FALSE,
    compress = FALSE,
    verbose = FALSE,
    rt_check = FALSE
  )
  dm_names <- names(attr(dadm, "designs"))

  out <- vector("list", length(pre_transform_terms))
  names(out) <- names(pre_transform_terms)
  for (param in names(pre_transform_terms)) {
    terms <- pre_transform_terms[[param]]
    if (!is.character(terms) || length(terms) < 1 || anyNA(terms)) {
      stop("pre_transform_terms$", param, " must be a non-empty character vector")
    }
    if (!param %in% dm_names) {
      stop("pre_transform_terms$", param, " refers to a parameter with no design matrix")
    }
    param_terms <- colnames(attr(dadm, "designs")[[param]])
    if (!all(terms %in% param_terms)) {
      bad <- setdiff(terms, param_terms)
      stop("pre_transform_terms$", param, " contains coefficient(s) not in the design: ",
           paste(bad, collapse = ", "))
    }
    out[[param]] <- unique(terms)
  }
  out
}

#' Contrast Enforcing Equal Prior Variance on each Level
#'
#' Typical contrasts impose different levels of marginal prior variance for the different levels.
#' This contrast can be used to ensure that each level has equal marginal priors (Rouder, Morey, Speckman, & Province; 2012).
#'
#' @param n An integer. The number of items for which to create the contrast
#'
#' @return A contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.bayes),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }
contr.bayes <- function(n) {
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
      levels <- seq_len(n)
    else stop("not enough degrees of freedom to define contrasts")
  }
  else levels <- n
  levels <- as.character(levels)
  n <- length(levels)
  cont <- diag(n)
  a <- n
  I_a <- diag(a)
  J_a <- matrix(1, nrow = a, ncol = a)
  Sigma_a <- I_a - J_a/a
  cont <- eigen(Sigma_a)$vectors[,seq_len(a-1), drop = FALSE]
  return(cont)
}


#' Contrast Enforcing Increasing Estimates
#'
#' Each level will be estimated additively from the previous level
#'
#' @param n an integer. The number of items for which to create the contrast.
#'
#' @return a contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.increasing),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }
contr.increasing <- function(n)
{
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
      levels <- seq_len(n)
    else stop("not enough degrees of freedom to define contrasts")
  }
  else levels <- n
  levels <- as.character(levels)
  n <- length(levels)
  contr <- matrix(0,nrow=n,ncol=n-1,dimnames=list(NULL,2:n))
  contr[lower.tri(contr)] <- 1
  contr
}

#' Contrast Enforcing Decreasing Estimates
#'
#' Each level will be estimated as a reduction from the previous level
#'
#' @param n an integer. The number of items for which to create the contrast.
#'
#' @return a contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.decreasing),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }
contr.decreasing <- function(n) {
  out <- contr.increasing(n)
  out[dim(out)[1]:1,]
}

#' Anova Style Contrast Matrix
#'
#' Similar to `contr.helmert`, but then scaled to estimate differences between conditions. Use in `design()`.
#'
#' @param n An integer. The number of items for which to create the contrast
#'
#' @return A contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.anova),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }

contr.anova <- function(n) {
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
       levels <- seq_len(n) else
         stop("not enough degrees of freedom to define contrasts")
  } else levels <- n
  levels <- as.character(levels)
  n <- length(levels)
  contr <- stats::contr.helmert(n)
  contr/rep(2*matrixStats::colMaxs(abs(contr)),each=dim(contr)[1])
}

add_accumulators <- function(data,matchfun=NULL,simulate=FALSE, type = "RACE", Fcovariates=NULL, fixed_accumulator_roles = NULL) {
  if(is.null(type) || !type %in% c("RACE", "RACEGNG", "SDT", "MT", "TC", "timed")) return(data)
  if (!is.factor(data$R)) stop("data must have a factor R")
  factors <- names(data)[!names(data) %in% c("R","rt","trials",Fcovariates)]
  if (!is.null(fixed_accumulator_roles)) {
    if (!is.factor(fixed_accumulator_roles) || length(fixed_accumulator_roles) == 0) {
      stop("fixed_accumulator_roles must be a factor with at least one level.")
    }
    role_levels <- levels(fixed_accumulator_roles)
    nacc <- length(role_levels)
    datar <- data[rep(seq_len(nrow(data)), times = nacc), , drop = FALSE]
    datar$lR <- factor(rep(role_levels, each = nrow(data)), levels = role_levels)
    datar <- datar[order(rep(seq_len(nrow(data)), nacc), datar$lR), , drop = FALSE]
    if (!is.null(matchfun)) {
      lM <- matchfun(datar)
      if (!is.factor(lM)){
        datar$lM <- factor(lM)
      } else{
        datar$lM <- factor(lM,levels=levels(lM))
      }
    }
    
  } else if (type %in% c("RACE","SDT", "RACEGNG", "timed")) {
    r_levels <- levels(data$R)
    if (type == "timed") {
      if ("time" %in% r_levels) stop("Data already has 'time' response level")
      r_levels <- c(r_levels, "time")
    }
    nacc <- length(r_levels)
    datar <- data[rep(seq_len(nrow(data)), times = nacc), , drop = FALSE]
    datar$lR <- factor(rep(r_levels, each = nrow(data)), levels = r_levels)
    datar <- datar[order(rep(seq_len(nrow(data)), nacc), datar$lR), , drop = FALSE]
    if (!is.null(matchfun)) {
      lM <- matchfun(datar)
      if (!is.factor(lM)){
        datar$lM <- factor(lM)
      } else{
        datar$lM <- factor(lM,levels=levels(lM))
      }
    }
    # Advantage NAFC
    nam <- unlist(lapply(strsplit(dimnames(datar)[[2]],"lS"),function(x)x[[1]]))
    islS <- nam ==""
    if (any(islS)) {
      if (sum(islS) != length(levels(data$R)))
        stop("The number of lS columns in the data must equal the length of Rlevels")
      lR <- unlist(lapply(strsplit(dimnames(datar)[[2]],"lS"),function(x){
        if (length(x)==2) x[[2]] else NULL}))
      if (!all(lR %in% levels(data$R)))
        stop("x  in lSx must be in Rlevels")
      if (any(names(datar=="lSmagnitude")))
        stop("Do not use lSmagnitude as a factor name")
      lSmagnitude <- as.character(datar$lR)
      for (i in levels(datar$lR)) {
        isin <- datar$lR==i
        lSmagnitude[isin] <- as.character(datar[isin,paste0("lS",i)])
      }
      factors <- factors[!(factors %in% dimnames(datar)[[2]][islS])]
      datar$lSmagnitude <- as.numeric(lSmagnitude)
    }
  }
  if (type %in% c("MT","TC")) {
    datar <- cbind(do.call(rbind,lapply(1:2,function(x){data})),
      lR=factor(rep(1:2,each=dim(data)[1]),levels=1:2))
    if (!is.null(matchfun)) {
      lM <- matchfun(datar)
      if (any(is.na(lM)) || !(is.logical(lM)))
        stop("matchfun not scoring properly")
      datar$lM <- factor(lM)
    }
  }
  row.names(datar) <- NULL
  if (simulate) {
    if(!'rt' %in% Fcovariates) {
      datar$rt <- NA
    }
  } else {
    if (!is.null(fixed_accumulator_roles)) {
      datar$winner <- FALSE
    } else if (type %in% c("MT","TC")) {
      datar$winner <- NA
    } else {
      R <- datar$R
      R[is.na(R)] <- levels(datar$lR)[1]
      datar$winner <- as.character(datar$lR) == as.character(R)
    }
    if (is.null(fixed_accumulator_roles) && "RACEGNG"%in%type) {
      is_inf <- is.infinite(datar$rt)
      if (any(is_inf)) {
        datar$winner[is_inf] <- datar$lR[is_inf] == "nogo" # here we code all missing responses as "nogo" winner because the make_missing trick doesn't handle that
      }
    }
  }
  datar
}

design_model_custom_ll <- function(data, design, model){
  if (!is.factor(data$subjects)) {
    data$subjects <- factor(data$subjects)
    warning("subjects column was converted to a factor")
  }
  dadm <- data
  model_input <- model
  attr(dadm, "model") <- function(){
    return(list(log_likelihood = model_input))
  }
  attr(dadm,"sampled_p_names") <- attr(design, "sampled_p_names")
  attr(dadm, "custom_ll") <- TRUE
  return(dadm)

}


#' Resolve the uniform-guess window for `pGuess`
#'
#' `pGuess` mixes a uniform "guess" density into observed RTs (Ratcliff &
#' Tuerlinckx, 2002; HDDM's `w_outlier`).  That density needs a window, and
#' EMC2 resolves ONE scalar window per dadm rather than a per-row column -- the
#' guess kernel is a single double the C++ side reads once per likelihood call.
#'
#' Resolution order:
#' \enumerate{
#'   \item `TC$guess_window`, if the user supplied it -- wins outright.  This is
#'     the exposed knob and normally would not be touched.  `TC$w_outlier` is
#'     accepted as the HDDM spelling and converted to `c(0, 1/(n_resp * w))`.
#'   \item Otherwise `LG = max(LT, LC)`, `UG = min(UC, UT)` from the existing
#'     bound columns.
#'   \item If `UG` is still infinite, `UG = max(5, floor(max finite rt) + 1)`.
#' }
#'
#' The fallback in (3) is deliberate.  HDDM's fixed 5 s window is improper by
#' its authors' own admission ("in practice, the outlier model is applied to all
#' RTs, even those larger than 5").  Keeping the density exactly 0 above `UG`
#' while defaulting `UG` to 5 would leave a 7 s outlier with no mixture
#' protection at all -- precisely the trial the model exists to catch.  Scaling
#' up to the next whole second above `max(rt)` stays proper and covers every
#' observed RT; for data that fit inside 5 s the behaviour is identical to HDDM.
#'
#' Because the window is `[max(LT, LC), min(UC, UT)]` by construction, a guess
#' can never be censored or truncated away, which is what removes all
#' interval-mass machinery from the likelihood.
#'
#' @param dadm A design-augmented data model (or any data frame carrying the
#'   `LT`/`LC`/`UT`/`UC` bounds as columns).
#' @param TC Optional truncation/censoring list; `TC$guess_window` or
#'   `TC$w_outlier` override the derived window.
#' @param verbose Whether to report the resolved window and the implied HDDM
#'   `w_outlier`.
#' @return A list with `window` (length-2 numeric) and `n_resp`, or `NULL` when
#'   no proper window exists.
#' @keywords internal
resolve_guess_window <- function(dadm, TC = NULL, verbose = FALSE) {
  bound <- function(nm, default) {
    v <- if (nm %in% colnames(dadm)) dadm[[nm]] else NULL
    if (is.null(v) || length(v) == 0) default else v
  }

  # A guess is an OVERT response, so it can be neither a withheld (nogo) trial
  # nor a timeout; those pseudo-levels do not count towards the response set.
  resp_levels <- levels(dadm$R)
  if (is.null(resp_levels)) resp_levels <- levels(dadm$lR)
  n_resp <- sum(!(resp_levels %in% c("nogo", "time")))
  if (n_resp < 1) return(NULL)

  win <- TC$guess_window
  if (is.null(win) && !is.null(TC$w_outlier)) {
    # HDDM spelling.  w_outlier = 0.1 only means "5 s" because HDDM has exactly
    # two responses; on a four-accumulator race the same number silently implies
    # a 2.5 s window.  Convert it here and store the window, which is what the
    # likelihood actually uses.
    w <- as.numeric(TC$w_outlier)[1]
    if (!is.finite(w) || w <= 0) stop("TC$w_outlier must be a positive number")
    win <- c(0, 1 / (n_resp * w))
  }

  if (!is.null(win)) {
    win <- as.numeric(win)
    if (length(win) != 2 || !all(is.finite(win)) || !(win[2] > win[1]))
      stop("TC$guess_window must be a length-2 numeric c(lower, upper) with upper > lower")
    LG <- win[1]; UG <- win[2]
    # A supplied window must cover the data.  The derived window below cannot
    # fail this -- it is built from the truncation/censoring bounds the RTs were
    # already checked against -- but an explicit `guess_window`, and especially
    # the `w_outlier` spelling, can be narrower than the observed RTs without
    # the user noticing.  The guess density is uniform on [LG, UG] and zero
    # outside it, so a trial beyond the edge has no guess component at all: the
    # mixture would silently stop being a mixture for exactly the slow outliers
    # pGuess exists to catch.  Refuse rather than quietly re-deriving a wider
    # window, so the user picks the number themselves.
    rt_obs <- dadm$rt
    rt_obs <- rt_obs[!is.na(rt_obs) & is.finite(rt_obs)]
    if (length(rt_obs)) {
      tol <- sqrt(.Machine$double.eps) * pmax(1, abs(UG), abs(LG), abs(rt_obs))
      if (any(rt_obs > UG + tol) || any(rt_obs < LG - tol)) {
        src <- if (!is.null(TC$guess_window)) "TC$guess_window" else "TC$w_outlier"
        stop("Guess window [", signif(LG, 4), ", ", signif(UG, 4),
             "] set by ", src, " does not cover the observed RTs [",
             signif(min(rt_obs), 4), ", ", signif(max(rt_obs), 4),
             "]. The uniform guess density is zero outside its window, so ",
             "trials beyond the edge would get no guess component. Widen the ",
             "window",
             if (src == "TC$w_outlier")
               paste0(" (w_outlier <= ", signif(1 / (n_resp * (max(rt_obs) - LG)), 4),
                      " for these data)")
             else "",
             ", or drop it and let it be derived from the LT/LC/UT/UC bounds.")
      }
    }
  } else {
    LT <- bound("LT", 0); LC <- bound("LC", 0)
    UT <- bound("UT", Inf); UC <- bound("UC", Inf)
    lows <- pmax(LT, LC)
    highs <- pmin(UT, UC)
    # Subject-wise bounds (from make_missing) can vary by row; take the union so
    # the window still covers every retained trial.
    if (length(unique(lows)) > 1 || length(unique(highs)) > 1)
      message("Guess window: bound columns are not constant across rows; using their union.")
    LG <- min(lows)
    UG <- max(highs)
    if (!is.finite(UG)) {
      rt <- dadm$rt
      max_rt <- suppressWarnings(max(rt[is.finite(rt)]))
      UG <- if (is.finite(max_rt)) max(5, floor(max_rt) + 1) else 5
    }
  }

  if (!is.finite(LG) || !is.finite(UG) || !(UG > LG)) return(NULL)
  if (verbose)
    message("Guess window: [", signif(LG, 4), ", ", signif(UG, 4), "] over ",
            n_resp, " responses (implied HDDM w_outlier = ",
            signif(1 / (n_resp * (UG - LG)), 4), ")")
  list(window = c(LG, UG), n_resp = n_resp)
}


compress_dadm <- function(da,designs,Fcov,Ffun)
    # out keeps only unique rows in terms of all parameters design matrices
    # R, lR and rt (at given resolution) from full data set
  {
  LT <- if ("LT" %in% colnames(da)) da$LT else 0
  UT <- if ("UT" %in% colnames(da)) da$UT else Inf
  LC <- if ("LC" %in% colnames(da)) da$LC else 0
  UC <- if ("UC" %in% colnames(da)) da$UC else Inf
    nacc <- length(unique(da$lR))
    # Covariate maps are part of the trial design just like the ordinary
    # covariate columns.  Normalise them before building the contraction key so
    # rows with different map values cannot be merged and subsequently receive
    # the map from whichever row happened to be retained.
    covariate_maps <- attr(da, "covariate_maps")
    if (!is.null(covariate_maps)) {
      covariate_maps <- lapply(covariate_maps, function(map) {
        if (is.null(dim(map))) {
          if (length(map) != nrow(da)) {
            stop("Each covariate map must have one value per expanded data row")
          }
          map <- matrix(map, ncol = 1L)
        }
        if (!is.matrix(map) || !is.numeric(map) || nrow(map) != nrow(da) ||
            ncol(map) < 1L) {
          stop("Each covariate map must be a numeric matrix with one row per expanded data row")
        }
        map
      })
    }
    # contract output
    design_cells_list <- lapply(designs, function(x) {
      do.call(paste, c(unname(as.data.frame(x[attr(x, "expand"), , drop = FALSE])), sep = "_"))
    })
    design_cells <- do.call(paste, c(design_cells_list, sep = "+"))

    cells <- paste(design_cells,
                da$subjects, da$R, da$lR, da$rt,
                LT, UT, LC, UC,  # Include truncation and censoring bounds in the cell key.
                sep="+"
              )
    # Make sure that if row is included for a trial so are other rows
    if (!is.null(Fcov)) {
      if (is.null(names(Fcov))) nFcov <- Fcov else nFcov <- names(Fcov)
      cells <- paste(cells, do.call(paste, c(unname(as.data.frame(da[, nFcov, drop = FALSE])), sep = "+")), sep = "+")
    }
    if (!is.null(Ffun))
      cells <- paste(cells, do.call(paste, c(unname(as.data.frame(da[, Ffun, drop = FALSE])), sep = "+")), sep = "+")

    covariate_map_cells <- NULL
    if (length(covariate_maps)) {
      map_cells <- lapply(covariate_maps, function(map) {
        do.call(paste, c(unname(as.data.frame(map)), sep = "+"))
      })
      covariate_map_cells <- do.call(paste, c(unname(map_cells), sep = "+"))
      cells <- paste(cells, covariate_map_cells, sep = "+")
    }

    if (nacc>1) {
      cells_mat <- matrix(cells, nrow = nacc)
      base_cells <- do.call(paste, c(unname(as.data.frame(t(cells_mat))), sep = "_"))
      cells <- paste0(rep(base_cells, each = nacc), "_", rep(1:nacc, times = length(base_cells)))
    }

    contract <- !duplicated(cells)
    out <- da[contract,,drop=FALSE]
    # Covariate maps are stored as attributes on the expanded data frame.
    # Subsetting a data frame does not subset arbitrary attributes, so without
    # explicitly carrying the maps through the contraction they retain the
    # expanded row count and no longer align with `out`.  Preserve matrix
    # dimensions as well: a one-column map must remain an n x 1 matrix rather
    # than being simplified to a vector.
    if (!is.null(covariate_maps)) {
      covariate_maps <- lapply(covariate_maps, function(map) map[contract, , drop = FALSE])
      attr(out, "covariate_maps") <- covariate_maps
    }
    attr(out,"contract") <- contract
    attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
    lR1 <- da$lR==levels(da$lR)[[1]]
    attr(out,"expand_winner") <- as.numeric(factor(cells[lR1],levels=unique(cells[lR1])))
    attr(out,"s_expand") <- da$subjects
    attr(out,"designs") <- lapply(designs,function(x){
      attr(x,"expand") <- attr(x,"expand")[contract]; x})

    # indices to use to contract further ignoring rt then expand back
    cells_nort <- paste(
      design_cells, da$subjects, da$R, da$lR, LT, UT, LC, UC, sep = "+"
    )
    if (!is.null(covariate_map_cells))
      cells_nort <- paste(cells_nort, covariate_map_cells, sep = "+")
    cells_nort <- cells_nort[contract]
    attr(out,"unique_nort") <- !duplicated(cells_nort)
    attr(out,"expand_nort") <- as.numeric(factor(cells_nort,levels=unique(cells_nort)))

    # indices to use to contract ignoring rt and response (R), then expand back
    cells_nortR <- paste(
      design_cells, da$subjects, LT, UT, LC, UC, sep = "+"
    )
    if (!is.null(covariate_map_cells))
      cells_nortR <- paste(cells_nortR, covariate_map_cells, sep = "+")
    cells_nortR <- cells_nortR[contract] #  ,da$lR
    attr(out,"unique_nortR") <- !duplicated(cells_nortR)
    attr(out,"expand_nortR") <- as.numeric(factor(cells_nortR,levels=unique(cells_nortR)))

    out
}

check_rt <- function(b, d, upper = TRUE)
  # Check bounds respected if present; b and d are numeric vectors
{
  idx <- !is.na(d) & is.finite(d)
  d <- d[idx]
  b <- b[idx]
  if (!length(d)) return(invisible(TRUE))
  b_num <- as.numeric(as.character(b))
  tol <- sqrt(.Machine$double.eps) * pmax(1, abs(b_num), abs(d))
  if (upper) {
    ok <- d <= (b_num + tol)
  } else {
    ok <- d >= (b_num - tol)
  }
  if (!all(ok)) stop("Bound not respected in data")
  invisible(TRUE)
}

rt_check_function <- function(data){
  # Truncation
  if ("UT"%in%colnames(data)) {
    # Use upper truncation bound, not censoring bound
    if (any(is.finite(data$UT))) {
      check_rt(data$UT, data$rt)
    }
  }
  if ("LT"%in%colnames(data)) {
    if (any(data$LT<0)) stop("Lower truncation cannot be negative")
    if (any(data$LT>0)) {
      check_rt(data$LT,data$rt,upper=FALSE)
    }
  }
  if ("UT"%in%colnames(data) & "LT"%in%colnames(data)) {
    DT <- data$UT - data$LT
    if (!is.null(DT)) {
      tol <- sqrt(.Machine$double.eps) * pmax(1, abs(data$UT), abs(data$LT))
      if (any(DT < -tol, na.rm = TRUE)) stop("UT must be greater than LT")
    }
  }

  # Censoring
  if ("UC"%in%colnames(data)) {
    if (any(is.finite(data$UC))) {
      check_rt(data$UC,data$rt)
    }
    if ("UT"%in%colnames(data) && any(is.finite(data$UC) & (data$UT < data$UC)) )
      stop("Upper truncation must not be less than upper censor")
  }
  if ("LC"%in%colnames(data)) {
    if (any(data$LC<0)) stop("Lower censor cannot be negative")
    if (any(data$LC>0)) { 
      check_rt(data$LC,data$rt,upper=FALSE)
    }
    if ("LT"%in%colnames(data) && any(data$LC!=0 & (data$LT>data$LC)))
      stop("Lower censor must not be less than lower truncation")
  }
  if (any(data$rt[!is.na(data$rt)]==-Inf) & !("LC"%in%colnames(data)))
    stop("Data must have an LC column if any rt = -Inf")
  if (any(data$rt[!is.na(data$rt)]==Inf) & !("UC"%in%colnames(data)))
    stop("Data must have a UC column if any rt = Inf")
  if ("UC"%in%colnames(data) & "LC"%in%colnames(data)) {
    DC <- data$UC - data$LC
    if (!is.null(DC)) {
      tol <- sqrt(.Machine$double.eps) * pmax(1, abs(data$UC), abs(data$LC))
      if (any(DC < -tol, na.rm = TRUE)) stop("UC must be greater than LC")
    }
  }
}

check_dm_identifiability <- function(designs, constants)
  # Warns if the design matrix for any parameter type is rank deficient over
  # its sampled (non-constant) columns, i.e., some sampled parameters trade off
  # exactly and only their linear combinations are identified. Rank is computed
  # on the (possibly compressed) design matrices; unique rows preserve rank.
{
  for (ptype in names(designs)) {
    dm <- designs[[ptype]]
    cols <- setdiff(colnames(dm), names(constants))
    if (length(cols) == 0) next
    X <- as.matrix(as.data.frame(dm)[, cols, drop = FALSE])
    if (!is.numeric(X) || any(!is.finite(X))) next
    ev <- eigen(crossprod(X), symmetric = TRUE)
    if (max(ev$values) <= 0) {
      involved <- cols
      n_dependent <- length(cols)
    } else {
      dependent <- ev$values < max(ev$values) * 1e-10
      if (!any(dependent)) next
      n_dependent <- sum(dependent)
      involved <- cols[rowSums(abs(ev$vectors[, dependent, drop = FALSE]) > 1e-8) > 0]
    }
    warning("The design matrix for '", ptype, "' is rank deficient (rank ",
            length(cols) - n_dependent, " < ", length(cols), " sampled parameters), ",
            "so these parameters are not individually identified and their ",
            "estimates can diverge even when the fit looks good. Parameters ",
            "involved: ", paste(involved, collapse = ", "), ". Fix ",
            n_dependent, " of them with constants (e.g., constants = c(`",
            involved[1], "` = 0)) or reparameterize the formula.",
            call. = FALSE)
  }
}

design_model <- function(data,design,model=NULL,
                         add_acc=TRUE,rt_resolution=1/60,verbose=TRUE,
                         compress=TRUE,rt_check=TRUE, add_da = FALSE, all_cells_dm = FALSE,
                         compress_dms = TRUE, drop_unobserved = FALSE)
{
  add_bound_column_if_needed <- function(df, col, value, default) {
    if (col %in% names(df)) return(df)
    if (is.null(value)) return(df)
    if (length(value) == 1 && isTRUE(all.equal(as.numeric(value), as.numeric(default)))) return(df)
    df[[col]] <- rep_len(value, nrow(df))
    df
  }

  if (is.null(model)) {
    if (is.null(design$model))
      stop("Model must be supplied if it has not been added to design")
    model <- design$model
  }
  if (model()$type=="SDT") rt_check <- FALSE
  # design_model() defaults rt_resolution to 1/60 and most internal callers take
  # that default rather than passing one, so make_emc()'s override would not
  # reach them.  Enforce the model's own answer here, where every path funnels.
  if (!is.null(rt_resolution) && !model_compress_ok(model)) rt_resolution <- NULL
  fixed_accumulator_roles <- design$fixed_accumulator_roles
  if(grepl("MRI", model()$type)){
    dadm <- data
    attr(dadm, "design_matrix") <- attr(design, "design_matrix")
    p_names <- names(model()$p_types)
    attr(dadm,"p_names") <- p_names
    sampled_p_names <- p_names[!(p_names %in% names(design$constants))]
    attr(dadm,"sampled_p_names") <- sampled_p_names
    return(dadm)
  }
  if (any(names(model()$p_types) %in% names(data)))
    stop("Data cannot have columns with the same names as model parameters")
  if (!is.factor(data$subjects)) {
    data$subjects <- factor(data$subjects)
    warning("subjects column was converted to a factor")
  }

  if (!any(names(data)=="trials")) data$trials <- 1:dim(data)[1]
  data <- add_bound_column_if_needed(data, "LT", design$LT, 0)
  data <- add_bound_column_if_needed(data, "LC", design$LC, 0)
  data <- add_bound_column_if_needed(data, "UT", design$UT, Inf)
  data <- add_bound_column_if_needed(data, "UC", design$UC, Inf)
  if(rt_check){rt_check_function(data)}
  if (!add_acc) da <- data else
    da <- add_accumulators(data,design$matchfun,type=model()$type,Fcovariates=design$Fcovariates,fixed_accumulator_roles=fixed_accumulator_roles)
  order_idx <- order(da$subjects)
  da <- da[order_idx,] # fixes different sort in add_accumulators depending on subject type

  # Only create Ffunction columns when missing.
  # Many designs use stochastic Ffunctions (e.g., SSD assignment), so overwriting an existing
  # column would silently change the design encoded in the data and break likelihood checks.
  if (!is.null(design$Ffunctions)) for (i in names(design$Ffunctions)) {
    if (i %in% names(da)) next
    newF <- stats::setNames(data.frame(design$Ffunctions[[i]](da)), i)
    da[, i] <- newF
  }

  # Add covariate_map as attribute to da
  if(!is.null(model()$trend)) {
    trend_list <- model()$trend
    for(i in 1:length(trend_list)) {
      if(!is.null(trend_list[[i]]$map)) {
        if(!'covariate_maps' %in% names(attributes(da))) attr(da, 'covariate_maps') <- list()
        covariate_map_names <- names(trend_list[[i]]$map)
        covariate_map_functions <- trend_list[[i]]$map
        for(map_n in 1:length(covariate_map_names)) {
          covs <- trend_list[[i]]$covariate
          map_name <- covariate_map_names[map_n]
          map <- covariate_map_functions[[map_n]](dadm=da, covs)

          # A single covariate is allowed to be returned as a vector, but the
          # downstream trend engine requires a numeric matrix. Normalize that
          # case here and validate all maps before attaching them to `da`.
          if (is.null(dim(map))) {
            if (length(covs) != 1L || length(map) != nrow(da)) {
              stop("Covariate map '", map_name,
                   "' must have one value per expanded data row")
            }
            map <- matrix(map, ncol = 1L,
                          dimnames = list(NULL, covs))
          }
          if (!is.matrix(map) || nrow(map) != nrow(da) ||
              ncol(map) != length(covs) || !is.numeric(map)) {
            stop("Covariate map '", map_name,
                 "' must be a numeric matrix with one row per expanded data row",
                 " and one column per covariate")
          }
          if (is.null(colnames(map))) colnames(map) <- covs
          attr(da, 'covariate_maps')[[map_name]] <- map
        }
      }
    }
  }

  if (is.null(model()$p_types) | is.null(model()$Ttransform))
    stop("p_types and Ttransform must be supplied")
  if (!all(unlist(lapply(design$Flist,class))=="formula"))
    stop("Flist must contain formulas")
  nams <- unlist(lapply(design$Flist,function(x)as.character(stats::terms(x)[[2]])))
  names(design$Flist) <- nams
  if (is.null(design$Clist)) design$Clist=list(stats::contr.treatment)
  if (!is.list(design$Clist)) stop("Clist must be a list")
  pnames <- names(model()$p_types)
  if (!is.list(design$Clist[[1]])[1]){
    design$Clist <- stats::setNames(lapply(1:length(pnames),
                                           function(x)design$Clist),pnames)
  } else {
   missing_p_types <- pnames[!(pnames %in% names(design$Clist))]
   if (length(missing_p_types)>0) {
     nok <- length(design$Clist)
      for (i in 1:length(missing_p_types)) {
        design$Clist[[missing_p_types[i]]] <- list(stats::contr.treatment)
        names(design$Clist)[nok+i] <- missing_p_types[i]
      }
    }
  }
  for (i in pnames) if (!is.null(design$Flist[[i]])) attr(design$Flist[[i]],"Clist") <- design$Clist[[i]]

  out <- lapply(design$Flist, make_dm, da = da, Fcovariates = design$Fcovariates,
                add_da = add_da, all_cells_dm = all_cells_dm, compress_dms = compress_dms,
                drop_unobserved = drop_unobserved)

  # sampled_pars()/mapped_pars() deliberately build a compact parameter map
  # with add_acc = FALSE; the invariant is checked on the actual expanded
  # likelihood data in make_emc/design_model(add_acc = TRUE).
  model_spec <- tryCatch(model(), error = function(e) NULL)
  if (is.list(model_spec) && isTRUE(model_spec$correlated) && isTRUE(add_acc) &&
      !identical(model_spec$correlation_type, "rdmswtn_gaussian_copula")) {
    rho_dm <- out[["rho"]]
    if (is.null(rho_dm) || is.null(da$lM)) {
      stop("BAwLcorr requires a matchfun-generated lM role indicator in the expanded data.")
    }
    rho_expand <- attr(rho_dm, "expand")
    rho_full <- if (!is.null(rho_expand) && length(rho_expand) == nrow(da)) {
      rho_dm[rho_expand, , drop = FALSE]
    } else {
      rho_dm
    }
    n_lR <- length(levels(da$lR))
    if (n_lR < 1L || nrow(rho_full) != nrow(da) || nrow(da) %% n_lR != 0L) {
      stop("BAwLcorr could not verify the within-trial rho design.")
    }
    # Constants do not create a free row-specific effect.  Removing their
    # columns makes the check accept rho ~ coupled with the PM/reference
    # intercept fixed to zero, while still rejecting rho ~ 0 + lR.
    free_cols <- !(colnames(rho_full) %in% names(design$constants))
    rho_free <- rho_full[, free_cols, drop = FALSE]
    if (ncol(rho_free) > 0L) {
      for (j in seq_len(nrow(da) / n_lR)) {
        rows <- ((j - 1L) * n_lR + 1L):(j * n_lR)
        block <- rho_free[rows, , drop = FALSE]
        active <- rowSums(abs(block) > 1e-12) > 0
        if (sum(active) > 1L) {
          ref <- block[which(active)[1L], , drop = FALSE]
          same <- apply(block[active, , drop = FALSE], 1L, function(x)
            isTRUE(all.equal(as.numeric(x), as.numeric(ref), tolerance = 1e-12)))
          if (!all(same)) {
            stop("BAwLcorr rho must be shared within each trial; use a trial-level formula (for example rho ~ 1 or rho ~ TrialType). Row-level formulas such as rho ~ 0 + lR are not allowed.")
          }
        }
      }
    }
  }
  if (is.list(model_spec) &&
      identical(model_spec$correlation_type, "rdmswtn_gaussian_copula") &&
      isTRUE(add_acc)) {
    rho_dm <- out[["rho"]]
    if (is.null(rho_dm)) {
      stop("RDMSWTNcorr requires a rho design.")
    }
    rho_expand <- attr(rho_dm, "expand")
    rho_full <- if (!is.null(rho_expand) && length(rho_expand) == nrow(da)) {
      rho_dm[rho_expand, , drop = FALSE]
    } else {
      rho_dm
    }
    n_lR <- length(levels(da$lR))
    if (n_lR < 1L || nrow(rho_full) != nrow(da) ||
        nrow(da) %% n_lR != 0L) {
      stop("RDMSWTNcorr could not verify the direct-pair rho design.")
    }
    free_cols <- !(colnames(rho_full) %in% names(design$constants))
    rho_free <- rho_full[, free_cols, drop = FALSE]
    for (j in seq_len(nrow(da) / n_lR)) {
      rows <- ((j - 1L) * n_lR + 1L):(j * n_lR)
      block <- rho_free[rows, , drop = FALSE]
      active <- if (ncol(block)) rowSums(abs(block) > 1e-12) > 0 else
        rep(FALSE, n_lR)
      if ("RACE" %in% names(da)) {
        n_acc <- suppressWarnings(as.integer(as.character(da$RACE[rows[1L]])))
        if (is.finite(n_acc)) {
          active <- active & (seq_len(n_lR) <= n_acc)
        }
      }
      if (sum(active) > 2L) {
        stop("RDMSWTNcorr rho design may select at most two accumulator rows per trial; use a participation factor and fix opted-out coefficients to zero.")
      }
      if (sum(active) == 2L &&
          !isTRUE(all.equal(as.numeric(block[active, , drop = FALSE][1L, ]),
                            as.numeric(block[active, , drop = FALSE][2L, ]),
                            tolerance = 1e-12))) {
        stop("RDMSWTNcorr's two participating rows must share the same rho design.")
      }
    }
  }
  if (!is.null(rt_resolution) & !is.null(da$rt))
    da$rt <- .floor_to_rt_resolution(da$rt, rt_resolution)
  if (compress){
    dadm <- compress_dadm(da,designs=out, Fcov=design$Fcovariates,Ffun=names(design$Ffunctions))
    # Change expansion names
    if(!is.null(dadm$lR)){
      attr(dadm,"expand") <- attr(dadm,"expand_winner")
      attr(dadm,"expand_winner") <- NULL
    }
  }  else {
    dadm <- da
    attr(dadm,"designs") <- out
    attr(dadm,"s_expand") <- da$subjects
    if(is.null(dadm$lR)){
      attr(dadm,"expand") <- 1:nrow(dadm)
    } else{
      attr(dadm,"expand") <- 1:(nrow(dadm)/length(unique(dadm$lR)))
    }
  }
  p_names <-  unlist(lapply(out,function(x){dimnames(x)[[2]]}),use.names=FALSE)
  dropped_constants <- character(0)
  if (isTRUE(drop_unobserved)) {
    dropped <- unique(unlist(lapply(out, attr, which = "dropped_columns"),
                             use.names = FALSE))
    dropped_constants <- intersect(names(design$constants), dropped)
    if (length(dropped_constants) > 0) {
      design$constants <- design$constants[!names(design$constants) %in% dropped_constants]
    }
  }
  bad_constants <- names(design$constants)[!(names(design$constants) %in% p_names)]
  if (length(bad_constants) > 0)
    stop("Constant(s) ",paste(bad_constants,collapse=" ")," not in design")

  # Pick out constants
  sampled_p_names <- p_names[!(p_names %in% names(design$constants))]
  attr(dadm,"p_names") <- p_names
  attr(dadm,"sampled_p_names") <- sampled_p_names
  if (verbose) check_dm_identifiability(out, design$constants)
 
  # `type` is a length-1 string (e.g. "DDM", "DDMGNG"); `%in%` would fail for "DDMGNG".
  if (grepl("DDM", model()$type)) nunique <- dim(dadm)[1] else
    nunique <- dim(dadm)[1]/length(levels(dadm$lR))
  if (verbose & compress) {
    if (all(c(dadm$LC,dadm$LT)==0) & all(is.infinite(c(dadm$UC,dadm$UT))))
      mismes <- NULL else {
        okDADM <- !is.na(dadm$R) | !is.infinite(dadm$rt)
        okDA <- !is.na(da$R) | !is.infinite(da$rt)
        mismes <- paste0(" (with no missing ",round(sum(okDA)/sum(okDADM),1),"x)")
      }
    message("Likelihood speedup factor: ",
    round(dim(da)[1]/dim(dadm)[1],1),mismes,", ",nunique," unique trials")
  }
  attr(dadm,"model") <- model
  attr(dadm,"constants") <- design$constants
  attr(dadm,"ok_trials") <- is.finite(data$rt)
  attr(dadm,"s_data") <- data$subjects
  # One scalar guess window per dadm, attached AFTER compression so it needs no
  # place in the compression key and cannot be lost by contraction.  Only
  # resolved when the model actually declares pGuess; see resolve_guess_window().
  if ("pGuess" %in% names(model()$p_types)) {
    pg_const <- if ("pGuess" %in% names(design$constants)) design$constants[["pGuess"]] else NULL
    pg_disabled <- !is.null(pg_const) && !is.finite(pg_const)
    gw <- resolve_guess_window(dadm, design$TC, verbose = verbose && !pg_disabled)
    if (!is.null(gw)) {
      attr(dadm,"guess_window") <- gw$window
      attr(dadm,"guess_n_resp") <- gw$n_resp
    }
  }
  dadm
}


make_full_dm <- function(form, Clist, da) {
  if (is.null(Clist)) Clist <- attr(form, "Clist")
  pnam <- stats::terms(form)[[2]]
  da[[pnam]] <- 1
  # Check if there are any nested CList entries to only contrast for this parameter
  if(any(names(Clist) == pnam)){
    replac <- Clist[[pnam]]
    for(i in 1:length(replac)){
      Clist[[names(replac)[i]]] <- replac[[i]]
    }
    Clist[[pnam]] <- NULL
  }
  for (i in names(Clist)) {
    if (i %in% names(da)) {
      if (!is.factor(da[[i]])) {
        stop(i, " must be a factor (design factors has a parameter name?)")
      }

      levs <- levels(da[[i]])
      nl <- length(levs)

      if (class(Clist[[i]])[1] == "function") {
        stats::contrasts(da[[i]]) <- do.call(Clist[[i]], list(n = levs))
      } else {
        if (!is.matrix(Clist[[i]]) || nrow(Clist[[i]]) != nl) {
          if (all(levs %in% row.names(Clist[[i]]))) {
            Clist[[i]] <- Clist[[i]][levs, ]
          } else {
            stop("Clist for ", i, " not a ", nl, " row matrix")
          }
        } else {
          dimnames(Clist[[i]])[[1]] <- levs
        }
        stats::contrasts(da[[i]], how.many = ncol(Clist[[i]])) <- Clist[[i]]
      }
    }
  }

  out <- stats::model.matrix(form, da)

  if (dim(out)[2] == 1) {
    dimnames(out)[[2]] <- as.character(pnam)
  } else {
    if (attr(stats::terms(form), "intercept") != 0) {
      cnams <- paste(pnam, dimnames(out)[[2]][-1], sep = "_")
      dimnames(out)[[2]] <- c(pnam, cnams)
    } else {
      dimnames(out)[[2]] <- paste(pnam, dimnames(out)[[2]], sep = "_")
    }
  }

  return(out)
}

make_dm <- function(form,da,Clist=NULL,Fcovariates=NULL, add_da = FALSE, all_cells_dm = FALSE,
                    compress_dms = TRUE, drop_unobserved = FALSE)
  # Makes a design matrix based on formula form from augmented data frame da
{

  compress_dm <- function(dm, da = NULL, all_cells_dm = FALSE)
    # out keeps only unique rows, out[attr(out,"expand"),] gets back original.
  {
    cells <- do.call(paste, c(unname(as.data.frame(dm)), sep = "_"))
    ass <- attr(dm,"assign")
    contr <- attr(dm,"contrasts")
    if(!is.null(da)){
      da_cells <- do.call(paste, c(unname(as.data.frame(da)), sep = "_"))
      dups <- duplicated(paste0(cells, da_cells))
    } else{
      dups <- duplicated(cells)
    }
    out <- dm[!dups,,drop=FALSE]
    if(!is.null(da) & !all_cells_dm){
      if(nrow(da) != 0){
        out <- cbind(da[!dups,colnames(da) != "subjects",drop=FALSE], out)
      }
    }
    attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
    attr(out,"assign") <- ass
    attr(out,"contrasts") <- contr
    out
  }
  out <- make_full_dm(form, Clist, da)

  dropped_columns <- character(0)
  if (isTRUE(drop_unobserved)) {
    keep <- colSums(abs(out), na.rm = TRUE) > 0
    dropped_columns <- colnames(out)[!keep]
    if (any(!keep)) {
      out <- out[, keep, drop = FALSE]
      attr(out, "assign") <- attr(out, "assign")[keep]
    }
  }

  if (!compress_dms) {
    if (isTRUE(drop_unobserved)) attr(out, "dropped_columns") <- dropped_columns
    return(out)
  }

  if(add_da){
    da <- da[,all.vars(form)[-1], drop = F]
    out <- compress_dm(out, da, all_cells_dm)
  } else{
    out <- compress_dm(out)
  }
  if (isTRUE(drop_unobserved)) attr(out, "dropped_columns") <- dropped_columns
  return(out)
}

# data generation

# Used in make_data and make_emc
add_trials <- function(dat)
  # Add trials column, 1:n for each subject
{
  n <- table(dat$subjects)
  if (!any(names(dat)=="trials")) dat <- cbind.data.frame(dat,trials=NA)
  for (i in names(n)) dat$trials[dat$subjects==i] <- 1:n[i]
  dat
}

dm_list <- function(dadm)
  # Makes data model into subjects list for use by likelihood
  # Assumes each subject has the same design.
{

  sub_design <- function(designs,isin)
    lapply(designs,function(x) {
      attr(x,"expand") <- attr(x,"expand")[isin]
      x
    })

  LT <- if ("LT" %in% colnames(dadm)) dadm$LT else 0
  UT <- if ("UT" %in% colnames(dadm)) dadm$UT else Inf
  LC <- if ("LC" %in% colnames(dadm)) dadm$LC else 0
  UC <- if ("UC" %in% colnames(dadm)) dadm$UC else Inf
  if(length(LT)==1) {dadm$LT = rep(LT,nrow(dadm))}
  else{dadm$LT=LT}
  if(length(UT)==1) {dadm$UT = rep(UT,nrow(dadm))}
  else{dadm$UT=UT}
  if(length(LC)==1) {dadm$LC = rep(LC,nrow(dadm))}
  else{dadm$LC=LC}
  if(length(UC)==1) {dadm$UC = rep(UC,nrow(dadm))}
  else{dadm$UC=UC}
  model <- attr(dadm,"model")
  p_names <- attr(dadm,"p_names")
  sampled_p_names <- attr(dadm,"sampled_p_names")
  designs <- attr(dadm,"designs")
  expand <- attr(dadm,"expand")
  s_expand <- attr(dadm,"s_expand")
  unique_nort <- attr(dadm,"unique_nort")
  expand_nort <- attr(dadm,"expand_nort")
  unique_nortR <- attr(dadm,"unique_nortR")
  expand_nortR <- attr(dadm,"expand_nortR")
  dms_mri <- attr(dadm, "design_matrix")

  # winner on expanded dadm
  expand_winner <- attr(dadm,"expand")
  # subjects for first level of lR in expanded dadm
  slR1=dadm$subjects[expand][dadm$lR[expand]==levels(dadm$lR)[[1]]]

  dl <- stats::setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                        levels(dadm$subjects))
  for (i in levels(dadm$subjects)) {
    isin <- dadm$subjects == i # dadm
    dl[[i]] <- dadm[isin, ]
    dl[[i]]$subjects <- factor(as.character(dl[[i]]$subjects))

    if(!is.null(attr(dadm, 'covariate_maps'))) {
      covariate_maps <- attr(dadm, 'covariate_maps')
      for(ii in 1:length(covariate_maps)) {
        map <- covariate_maps[[ii]]
        if (is.null(dim(map))) {
          if (length(map) != nrow(dadm)) {
            stop("Each covariate map must have one value per data row")
          }
          map <- matrix(map, ncol = 1L)
        }
        if (!is.matrix(map) || nrow(map) != nrow(dadm)) {
          stop("Each covariate map must be a matrix with one row per data row")
        }
        covariate_maps[[ii]] <- map[isin, , drop = FALSE]
      }
      attr(dl[[i]], 'covariate_maps') <- covariate_maps
    }

    if(is.null(attr(dadm, "custom_ll"))){

      isin1 <- s_expand==i             # da
      isin2 <- attr(dadm,"s_data")==i  # data
      if(length(isin2) > 0){
        attr(dl[[i]],"expand") <- expand_winner[isin2]-min(expand_winner[isin2]) + 1
      }
      attr(dl[[i]], "model") <- NULL
      attr(dl[[i]], "p_names") <- p_names
      attr(dl[[i]], "sampled_p_names") <- sampled_p_names
      attr(dl[[i]], "designs") <- sub_design(designs, isin)
      attr(dl[[i]], "contract") <- NULL
      attr(dl[[i]], "expand_winner") <- NULL
      attr(dl[[i]], "ok_trials") <- NULL
      attr(dl[[i]], "s_data") <- NULL
      attr(dl[[i]], "s_expand") <- NULL
      attr(dl[[i]], "prior") <- NULL
      if (!is.null(dms_mri)) {
        attr(dl[[i]], "designs") <- make_mri_sampling_design(dms_mri[[i]], sampled_p_names)
        attr(dl[[i]], "design_matrix") <- NULL
      }

      attr(dl[[i]], "unique_nort") <- NULL
      attr(dl[[i]], "expand_nort") <- NULL
      # LL cache attrs are data-shape specific; drop any inherited cache from the
      # full dadm so per-subject caching is always rebuilt safely.
      attr(dl[[i]], "emc2_ll_cache_version") <- NULL
      attr(dl[[i]], "emc2_all_finite_trials") <- NULL
      attr(dl[[i]], "finite_rt_mask") <- NULL
      attr(dl[[i]], "finite_rt_unique_trial_indices") <- NULL
      attr(dl[[i]], "other_unique_trial_indices") <- NULL
      attr(dl[[i]], "RACE_nacc_by_row") <- NULL
      attr(dl[[i]], "RACE_mask") <- NULL
      dl[[i]] <- .cache_ll_data_attrs(dl[[i]], force_rebuild = TRUE)
    }
  }


  return(dl)
}

#' Update EMC Objects to the Current Version
#'
#' This function updates EMC objects created with older versions of the package to be compatible with the current version.
#'
#' @param emc An EMC object to update
#' @return An updated EMC object compatible with the current version
#' @examples
#' # Update the model to current version
#' updated_model <- update2version(samples_LNR)
#'
#' @export
update2version <- function(emc){
  # For older versions, ensure that the class is emc:
  class(emc) <- "emc"
  get_new_model <- function(old_model, pars){
    if(old_model()$c_name == "LBA"){
      model <- LBA
    } else if(old_model()$c_name == "DDM"){
      model <- DDM
    } else if(old_model()$c_name == "RDM"){
      model <- RDM
    } else if(old_model()$c_name == "LNR"){
      model <- LNR
    } else{
      stop("current model not supported for updating, sorry!!")
    }
    model_list <- model()
    model_list$transform <- fill_transform(transform = old_model()$transform,model)
    model_list$pre_transform <- fill_transform(transform = old_model()$pre_transform, model = model, p_vector = pars, is_pre = TRUE)
    model_list$bound <- fill_bound(bound = NULL,model)
    model <- function(){return(model_list)}
    return(model)
  }

  new_expand <- function(x){
    if(!is.null(x$winner)){
      old_exp <- attr(x, "expand")
      if(length(unique(old_exp) != length(unique(x$winner)))){
        # In older versions we were working with a different expand version
        new_x <- x[old_exp,]
        new_x <- new_x[new_x$winner,]

        reduced <- unique(new_x)       # keeps the first appearance of every row

        ## ——— 2. create the "expand" index ———
        key_full    <- do.call(paste, c(new_x,      sep = "\r"))   # one string per row
        key_reduced <- do.call(paste, c(reduced, sep = "\r"))   # the same for reduced
        attr(x, "expand") <- match(key_full, key_reduced)
      }
    }
    return(x)
  }
  update_expand <- function(emc){
    first_data <- emc[[1]]$data[[1]]
    if(is.data.frame(first_data)){
      emc[[1]]$data <- lapply(emc[[1]]$data, new_expand)
    } else{
      emc[[1]]$data <- lapply(emc[[1]]$data, function(y){
        y <- lapply(y, new_expand)})
    }
    return(emc)
  }

  emc <- update_expand(emc)
  design_list <- get_design(emc)
  design_list <- lapply(design_list, function(x){
    if(!is.null(x$Fcovariates)){
      x$Fcovariates <- names(x$Fcovariates)
    }
    return(x)
  })

  if(is.null(emc[[1]]$type)){
    type <- attr(emc[[1]], "variant_funs")$type
    emc <- lapply(emc, FUN = function(x){
      x$type <- type
      return(x)
    })
  } else{
    type <- emc[[1]]$type
  }
  # Restore model metadata when it is absent from the emc object.
  first_data <- emc[[1]]$data[[1]]
  if(is.null(emc[[1]]$model)){
    if(is.data.frame(first_data)){
      old_model <- attr(first_data, "model")
      new_model <- get_new_model(old_model, sampled_pars(design_list[[1]]))
      design_list[[1]]$model <- new_model
      emc[[1]]$model <- new_model
    } else{
      old_model <- lapply(first_data, function(x) attr(x, "model"))
      new_model <- mapply(get_new_model, old_model, lapply(design_list, sampled_pars))
      design_list <- mapply(function(x, y){
        x$model <- y
        return(list(x))
      }, design_list, new_model)
      emc[[1]]$model <- new_model
    }

  } else{
    if(is.data.frame(first_data)){
      emc[[1]]$model <- get_new_model(emc[[1]]$model, sampled_pars(design_list[[1]]))
      design_list[[1]]$model <- emc[[1]]$model
    } else{
      old_model <- emc[[1]]$model
      new_model <- mapply(get_new_model, old_model, lapply(design_list, sampled_pars))
      design_list <- mapply(function(x, y){
        x$model <- y
        return(list(x))
      }, design_list, new_model)
      emc[[1]]$model <- new_model
    }
  }

  group_design <- attr(emc[[1]]$prior, "group_design")

  prior_new <- emc[[1]]$prior
  attr(prior_new, "type") <- type
  prior_new <- prior(design_list, type, update = prior_new, group_design = group_design)
  class(prior_new) <- "emc.prior"
  emc <- lapply(emc, function(x){
    x$prior <- prior_new
    return(x)
  })
  class(emc) <- "emc"
  attr(emc, "design_list") <- NULL
  return(emc)
}


# Some s3 classes for design objects ---------------------------------------------------------

# Restrict a factorial mapping backbone to the factor combinations that occur
# in the supplied data.  This is kept separate from `mapped_pars()` so that
# data-backed parameter discovery in `sampled_pars()` and `make_emc()` uses
# exactly the same definition of an observed cell.
.restrict_to_observed_factors <- function(mapping_data, design, data,
                                          remove_subjects = TRUE,
                                          extra_factors = NULL) {
  if (is.null(data)) return(mapping_data)
  if (!is.data.frame(data)) stop("data must be a data frame")

  observed_factors <- names(design$Ffactors)
  if (isTRUE(remove_subjects)) observed_factors <- setdiff(observed_factors, "subjects")
  # Group-design predictors are part of the displayed mapping data even when
  # they are not in the subject-level formulae.  Restrict them to observed
  # combinations when the caller supplied those columns in `data`.
  extra_factors <- intersect(extra_factors, names(data))
  observed_factors <- unique(c(observed_factors, extra_factors))
  missing_factors <- setdiff(observed_factors, names(data))
  if (length(missing_factors) > 0) {
    stop("data is missing design factor(s): ", paste(missing_factors, collapse = ", "))
  }

  observed_factors <- intersect(observed_factors, names(mapping_data))
  if (length(observed_factors) == 0) return(mapping_data)

  mapping_key <- function(df) {
    values <- lapply(df[, observed_factors, drop = FALSE], as.character)
    values <- lapply(values, function(value) {
      value[is.na(value)] <- "<NA>"
      value
    })
    do.call(paste, c(values, sep = "\r"))
  }

  observed_keys <- unique(mapping_key(data))
  mapping_data[mapping_key(mapping_data) %in% observed_keys, , drop = FALSE]
}

# Give trial-level covariates a cell structure.
#
# `minimal_design()` builds a factorial backbone out of the design *factors*.  A
# covariate is not a factor, so it has no cell of its own: minimal_design either
# fills it with random values (nothing was supplied for it) or recycles a
# supplied vector positionally across whatever rows happen to line up.  Both are
# meaningless for a covariate that is confounded with the factors -- and that is
# exactly the case whenever a trend is worth inspecting, since the covariate is
# what drives the trend.
#
# So resolve covariates the same way the factors are resolved:
#   * values passed explicitly in `covariates` are crossed into every cell, so
#     the caller sees the mapping at each value they asked for;
#   * otherwise the values observed for that covariate *within each design cell*
#     are used, expanding to one row per (cell x observed value combination).
# Cells with many distinct values (a continuous covariate) are thinned to an
# evenly spaced slice of at most `n_covariates` of them so the output stays
# readable.
.expand_covariate_values <- function(mapping_data, design, data, supplied = NULL,
                                     remove_subjects = TRUE, n_covariates = 10) {
  covariates <- design$Fcovariates
  if (is.null(covariates) || length(covariates) == 0) return(mapping_data)
  covariates <- intersect(covariates, names(mapping_data))
  if (length(covariates) == 0) return(mapping_data)
  if (!is.list(supplied)) supplied <- NULL

  cross <- intersect(covariates, names(supplied))
  from_data <- setdiff(covariates, cross)
  if (!is.null(data)) from_data <- intersect(from_data, names(data)) else from_data <- character(0)

  thin <- function(idx) {
    if (length(idx) <= n_covariates) return(idx)
    idx[unique(round(seq(1, length(idx), length.out = n_covariates)))]
  }

  # 1. Explicitly supplied values: cross them into every cell.
  for (cv in cross) {
    values <- unique(supplied[[cv]])
    values <- values[order(values)]
    values <- values[thin(seq_along(values))]
    if (length(values) == 0) next
    mapping_data <- mapping_data[rep(seq_len(nrow(mapping_data)), each = length(values)), ,
                                 drop = FALSE]
    mapping_data[[cv]] <- rep(values, length.out = nrow(mapping_data))
  }
  rownames(mapping_data) <- NULL
  if (length(from_data) == 0) return(mapping_data)

  # 2. Remaining covariates: take the values observed in each design cell.
  cells <- names(design$Ffactors)
  if (isTRUE(remove_subjects)) cells <- setdiff(cells, "subjects")
  cells <- intersect(intersect(cells, names(mapping_data)), names(data))

  cell_key <- function(df) {
    if (length(cells) == 0) return(rep("", nrow(df)))
    values <- lapply(df[, cells, drop = FALSE], as.character)
    values <- lapply(values, function(value) {
      value[is.na(value)] <- "<NA>"
      value
    })
    do.call(paste, c(values, sep = "\r"))
  }

  observed <- unique(data[, c(cells, from_data), drop = FALSE])
  observed <- observed[do.call(order, unname(as.list(observed[, from_data, drop = FALSE]))), ,
                       drop = FALSE]
  by_cell <- lapply(split(seq_len(nrow(observed)), cell_key(observed)), thin)

  matches <- by_cell[cell_key(mapping_data)]
  reps <- vapply(matches, function(m) max(1L, length(m)), integer(1L))
  out <- mapping_data[rep(seq_len(nrow(mapping_data)), reps), , drop = FALSE]
  # Cells with nothing observed keep whatever minimal_design() put there.
  fill <- unlist(lapply(matches, function(m) if (length(m)) m else NA_integer_),
                 use.names = FALSE)
  keep <- !is.na(fill)
  out[keep, from_data] <- observed[fill[keep], from_data, drop = FALSE]
  rownames(out) <- NULL
  out
}

# Evaluate a stacked group-level beta vector for the subjects represented in a
# mapping data frame.  Group-design matrices are stored for the subjects used
# to fit the model, but mapped_pars() may construct a different set of
# (pseudo-)subjects to display.  Rebuilding each block from its original
# formula keeps the calculation aligned with the stored contrast coding while
# allowing those display rows to carry new factor combinations.
.group_design_subject_pars <- function(p_vector, mapping_data, group_design,
                                       base_par_names, sampled_par_names = NULL) {
  if (is.null(group_design) || length(base_par_names) == 0) return(p_vector)

  if (is.data.frame(p_vector)) p_vector <- as.matrix(p_vector)
  matrix_names <- NULL
  if (length(dim(p_vector)) > 0) {
    if (ncol(p_vector) == 1L) {
      matrix_names <- rownames(p_vector)
      p_vector <- p_vector[, 1L]
    } else if (nrow(p_vector) == 1L) {
      matrix_names <- colnames(p_vector)
      p_vector <- p_vector[1L, ]
    } else {
      stop("p_vector must be a vector (or a one-row/one-column matrix) when used with group_design")
    }
  }
  supplied_names <- names(p_vector)
  if (is.null(supplied_names) && !is.null(matrix_names) &&
      !identical(matrix_names, as.character(seq_along(matrix_names)))) {
    supplied_names <- matrix_names
  }
  p_vector <- as.numeric(p_vector)

  expanded_names <- add_group_par_names(base_par_names, group_design)
  if (is.null(supplied_names)) {
    if (length(p_vector) == length(expanded_names)) {
      supplied_names <- expanded_names
    } else if (length(p_vector) == length(base_par_names)) {
      supplied_names <- base_par_names
    } else if (!is.null(sampled_par_names) && length(p_vector) == length(sampled_par_names)) {
      supplied_names <- sampled_par_names
    } else {
      stop("Unnamed p_vector must have length ", length(expanded_names),
           " (or ", length(base_par_names), ") when used with group_design")
    }
  }
  names(p_vector) <- supplied_names

  if (is.null(mapping_data$subjects)) {
    stop("mapping data must contain a subjects column when used with group_design")
  }
  subject_ids <- unique(as.character(mapping_data$subjects))
  out <- matrix(NA_real_, nrow = length(subject_ids), ncol = length(base_par_names),
                dimnames = list(subject_ids, base_par_names))

  formulas <- attr(group_design, "Flist")
  lhs_names <- if (length(formulas)) {
    vapply(formulas, function(f) as.character(stats::terms(f)[[2L]]), character(1L))
  } else character(0)

  for (j in seq_along(base_par_names)) {
    par_name <- base_par_names[[j]]
    block <- group_design[[par_name]]
    if (is.null(block)) {
      if (!par_name %in% names(p_vector)) {
        stop("p_vector is missing subject-level parameter '", par_name, "'")
      }
      values <- rep(unname(p_vector[[par_name]]), nrow(mapping_data))
    } else {
      if (!length(formulas) || !(par_name %in% lhs_names)) {
        stop("Could not find the group-design formula for parameter '", par_name, "'")
      }
      f <- formulas[[match(par_name, lhs_names)]]
      rhs_vars <- all.vars(stats::terms(f)[[3L]])
      group_data <- attr(group_design, "data")

      # The group-design object contains the matrix evaluated on the original
      # subjects.  Reuse those rows whenever the display data contains the
      # same predictor values.  Besides being cheaper, this preserves custom
      # contrasts for group designs created before contrast metadata was
      # stored on the object.
      X <- NULL
      if (!is.null(group_data) && nrow(block) == nrow(group_data) &&
          all(rhs_vars %in% names(group_data))) {
        key <- function(df, vars) {
          if (!length(vars)) return(rep("", nrow(df)))
          values <- lapply(df[, vars, drop = FALSE], as.character)
          values <- lapply(values, function(value) {
            value[is.na(value)] <- "<NA>"
            value
          })
          do.call(paste, c(values, sep = "\r"))
        }
        rows <- match(key(mapping_data, rhs_vars), key(group_data, rhs_vars))
        if (all(!is.na(rows))) X <- block[rows, , drop = FALSE]
      }
      if (is.null(X)) {
        group_contrasts <- attr(group_design, "contrasts")
        if (is.null(group_contrasts)) group_contrasts <- attr(block, "contrasts")
        X <- build_design(f, mapping_data,
                          contrasts.arg = group_contrasts)
      }
      target_names <- colnames(block)
      if (ncol(X) != length(target_names)) {
        stop("Group-design matrix for '", par_name, "' has ", ncol(X),
             " columns for the mapping data, but the fitted design has ",
             length(target_names), ".")
      }
      colnames(X) <- target_names

      coefficients <- numeric(length(target_names))
      names(coefficients) <- target_names
      present <- intersect(target_names, names(p_vector))
      coefficients[present] <- p_vector[present]
      # A subject-level vector is a useful shorthand for an intercept-only
      # group effect (and makes an unexpanded sampled_pars() vector map
      # sensibly): use it for the intercept and set omitted slopes to zero.
      intercept <- intersect(par_name, target_names)
      if (length(intercept) && !(intercept[[1L]] %in% names(p_vector)) &&
          par_name %in% names(p_vector)) {
        coefficients[intercept[[1L]]] <- p_vector[[par_name]]
      }
      if (!length(present) && !length(intercept) && !(par_name %in% names(p_vector))) {
        stop("p_vector is missing group-level coefficients for '", par_name, "'")
      }
      values <- as.vector(X %*% coefficients)
    }
    # A group-level mean is constant for all trial rows of a subject.  Keep
    # one value per subject for get_pars_matrix_oo(), which uses subject names
    # to expand it back to the trial rows.
    out[, j] <- values[match(subject_ids, as.character(mapping_data$subjects))]
  }
  out
}

#' Parameter Mapping Back to the Design Factors
#'
#' Maps parameters of the cognitive model back to the experimental design. If p_vector
#' is left unspecified will print a textual description of the mapping.
#' Otherwise the p_vector can be created using ``sampled_pars()``.
#' The returned matrix shows whether/how parameters
#' differ across the experimental factors.
#'
#' @param x an `emc`, `emc.prior` or `emc.design` object
#' @param p_vector Optional. Specify parameter vector to get numeric mappings.
#' Must be in the form of ``sampled_pars(design)`` (or
#' ``sampled_pars(design, group_design = group_design)`` when a group design
#' is supplied).
#' @param group_design Optional `emc.group_design` object. When supplied, group-
#'   level coefficients are evaluated for the group-design covariates before
#'   mapping the resulting subject-level parameters to the experimental design.
#' @param model Optional model type (if not already specified in ``design``)
#' @param digits Integer. Will round the output parameter values to this many decimals
#' @param ... optional arguments
#' @param remove_subjects Boolean. Whether to include subjects as a factor in the design
#' @param covariates Values for the covariates specified in the design, as a
#'   named list (e.g. `list(delta = 0:5)`). Each covariate's values are crossed
#'   into every design cell. Covariates not named here are taken from `data`,
#'   using the values each one actually takes within each design cell, so that a
#'   trend driven by a covariate is visible in the output. A covariate that can
#'   be resolved from neither is filled with random values, with a warning.
#' @param data Optional data frame. If supplied, only factor combinations observed
#'   in the data are returned.
#' @param use_data Logical. Whether to restrict mappings to combinations observed
#'   in `data`. Defaults to `TRUE`.
#' @param n_covariates Integer. Maximum number of distinct values of each
#'   covariate to report per design cell. A continuous covariate can take
#'   thousands of values; the reported ones are an evenly spaced slice of those
#'   observed (or supplied). Defaults to 10.
#' @return Matrix with a column for each factor in the design and for   each model parameter type (``p_type``).
#' @examples
#' # First define a design:
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' mapped_pars(design_DDMaE)
#' # Then create a p_vector:
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'           t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # This will map the parameters of the p_vector back to the design
#' mapped_pars(design_DDMaE, p_vector)
#'
#' @export
mapped_pars <- function(x, p_vector = NULL, model=NULL,
                        digits=3,remove_subjects=TRUE,
                        covariates=NULL, data = NULL, use_data = TRUE,
                        n_covariates = 10, group_design = NULL, ...)
  # Show augmented data and corresponding mapped parameter
{
  UseMethod("mapped_pars")
}

#' @rdname mapped_pars
#' @export
mapped_pars.emc.design <- function(x, p_vector = NULL, model=NULL,
                                   digits=3,remove_subjects=TRUE,
                                   covariates=NULL, data = NULL, use_data = TRUE,
                                   n_covariates = 10, group_design = NULL, ...){
  if(is.null(x)) return(NULL)
  if(is.null(x$Ffactors)){
    x <- x[[1]]
  }
  if(!is.null(attr(x, "custom_ll"))){
    stop("Mapped_pars not available for this design type")
  }
  design <- x
  design_data <- attr(design, "data")
  if(is.null(p_vector)){
    return(verbal_dm(design))
  }
  remove_RACE <- TRUE
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if (is.null(covariates))
    Fcovariates <- design$Fcovariates else
      Fcovariates <- covariates
  if (is.null(model)) if (is.null(design$model))
    stop("Must specify model as not in design") else model <- design$model
  if (remove_subjects) design$Ffactors$subjects <- design$Ffactors$subjects[1]
  if (is.null(data)) data <- design_data
  mapping_data <- minimal_design(
    design, covariates = Fcovariates, verbose = F, drop_R = F,
    add_acc = F, drop_subjects = F, do_functions = F,
    group_design = group_design
  )
  group_factors <- if (!is.null(group_design)) {
    unique(unlist(lapply(attr(group_design, "Flist"), function(f) {
      all.vars(stats::terms(f)[[3L]])
    })))
  } else character(0)
  if (isTRUE(use_data)) {
    mapping_data <- .restrict_to_observed_factors(mapping_data, design, data,
                                                  remove_subjects = remove_subjects,
                                                  extra_factors = group_factors)
  }
  # A covariate has no cell of its own in the factorial backbone, so give it one
  # from the supplied values or from the data.  Without this a trend driven by a
  # covariate is invisible here: the covariate column is noise or is recycled
  # into cells it never occurs in.
  mapping_data <- .expand_covariate_values(
    mapping_data, design, if (isTRUE(use_data)) data else NULL,
    supplied = if (is.list(covariates)) covariates else NULL,
    remove_subjects = remove_subjects, n_covariates = n_covariates
  )

  # A group design is defined over subjects, not over the trial-level factors
  # in `design`.  `minimal_design(..., group_design = ...)` adds those
  # subject-level predictors to the mapping rows (and creates one pseudo
  # subject per group-design row when subjects were removed).  Repeated
  # subjects are only an implementation detail in that case; collapse them
  # to the distinct displayed cells before constructing the model data.
  if (!is.null(group_design) && isTRUE(remove_subjects) && nrow(mapping_data) > 0) {
    keep <- setdiff(names(mapping_data), c("subjects", "trials"))
    if (length(keep) > 0) {
      mapping_data <- mapping_data[!duplicated(mapping_data[, keep, drop = FALSE]), ,
                                   drop = FALSE]
      rownames(mapping_data) <- NULL
    }
  }
  if (!is.null(group_design) && !is.factor(mapping_data$subjects)) {
    mapping_data$subjects <- factor(mapping_data$subjects)
  }
  unresolved <- setdiff(design$Fcovariates,
                        c(if (is.list(covariates)) names(covariates),
                          if (isTRUE(use_data)) names(data),
                          c("LT", "LC", "UT", "UC")))
  if (length(unresolved) > 0) {
    warning("covariate(s) ", paste(unresolved, collapse = ", "), " could not be ",
            "resolved from `data` or `covariates`, so they were filled with random ",
            "values; any trend on them is meaningless here. Supply `data`, or pass ",
            "covariates = list(", unresolved[1], " = <values>).", call. = FALSE)
  }
  dadm <- design_model(mapping_data, design,model,rt_check=FALSE,compress=FALSE,
                       verbose = FALSE,
                       drop_unobserved = isTRUE(use_data) && !is.null(data))
  ok <- !(names(dadm) %in% c("subjects","trials","R","rt","winner"))

  # `get_pars_matrix_oo()` consumes subject-level design coefficients.  For a
  # group design, first evaluate each group-design block at the group
  # predictors attached to the mapping rows and multiply it by the supplied
  # beta coefficients.  This is the deterministic counterpart of the
  # subject-mean construction used by the standard sampler.
  map_p <- p_vector
  if (!is.null(group_design)) {
    map_p <- .group_design_subject_pars(p_vector, mapping_data, group_design,
                                        base_par_names = names(suppressMessages(sampled_pars(design))),
                                        sampled_par_names = attr(dadm, "sampled_p_names"))
  }
  out <- cbind(dadm[,ok, drop = F],round(get_pars_matrix_oo(map_p,dadm, design$model()),digits))
  if (model()$type=="SDT")  out <- out[dadm$lR!=levels(dadm$lR)[length(levels(dadm$lR))],]
  if (model()$type=="DDM")  out <- out[,!(names(out) %in% c("lR","lM"))]
  if (any(names(out)=="RACE") && remove_RACE)
    out <- out[as.numeric(out$lR) <= as.numeric(as.character(out$RACE)),,drop=FALSE]

  return(out)
}



#' Get Model Parameters from a Design
#'
#' Makes a vector with zeroes, with names and length corresponding to the
#' model parameters of the design.
#'
#' @param x an `emc.design` object made with `design()` or an `emc` object.
#' @param group_design an `emc.group_design` object made with `group_design()`
#' @param doMap logical. If `TRUE` will also include an attribute `map`
#' with the design matrices that perform the mapping back to the design
#' @param add_da Boolean. Whether to include the relevant data columns in the map attribute
#' @param all_cells_dm Boolean. Whether to include all levels of a factor in the mapping attribute,
#' even when one is dropped in the design
#' @param data Optional data frame (or list of data frames for a joint model).
#'   When supplied, only observed factor combinations are used.
#' @param use_data Logical. Whether to restrict the parameter design to
#'   combinations observed in `data`. Defaults to `TRUE`.
#'
#'
#' @return Named vector.
#' @examples
#' # First define a design
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then for this design get which cognitive model parameters are sampled:
#' sampled_pars(design_DDMaE)
#'
#' @export
sampled_pars <- function(x,group_design=NULL,doMap=FALSE, add_da = FALSE,
                         all_cells_dm = FALSE, data = NULL, use_data = TRUE)
{
  UseMethod("sampled_pars")
}

#' @rdname sampled_pars
#' @export
sampled_pars.emc.design <- function(x,group_design=NULL,doMap=FALSE, add_da = FALSE,
                                    all_cells_dm = FALSE, data = NULL,
                                    use_data = TRUE){
  design <- x
  if(is.null(design)) return(NULL)
  if("Flist" %in% names(design)){
    design <- list(design)
  }
  out <- c()
  map_list <- list()
  if(is.null(names(design))){
    names(design) <- as.character(1:length(design))
  }
  for(j in 1:length(design)){
    cur_name <- names(design)[j]
    cur_design <- design[[j]]
    if(!is.null(attr(cur_design, "custom_ll"))){
      pars <- numeric(length(attr(cur_design,"sampled_p_names")))
      if(length(design) != 1){
        map_list[[cur_name]] <- NA
        names(pars) <- paste(cur_name,  attr(cur_design,"sampled_p_names"), sep = "|")
      } else{
        names(pars) <- attr(cur_design,"sampled_p_names")
      }
      out <- c(out, pars)
      next
    }
    model <- cur_design$model
    if (is.null(model)) stop("Must supply model as not in design")

    observed_data <- data
    use_covariate_data <- !is.null(data)
    if (!is.null(observed_data) && is.list(observed_data) &&
        !is.data.frame(observed_data)) {
      observed_data <- observed_data[[j]]
    }
    if (is.null(observed_data)) observed_data <- attr(cur_design, "data")
    if (!isTRUE(use_data)) observed_data <- NULL

    if(grepl("MRI", model()$type)){
      pars <- model()$p_types
      if(length(design) != 1){
        names(pars) <- paste(cur_name,  names(pars), sep = "|")
        map_list[[cur_name]] <- NA
      }
      pars[1:length(pars)] <- 0
      out <- c(out, pars)
      next
    }
    cur_design$Ffactors$subjects <- 1
    min_design <- minimal_design(cur_design, drop_subjects = F, drop_R = F, verbose = F,
                                 emc = if (use_covariate_data) observed_data else NULL,
                                 do_functions = F, add_acc = T)
    if (isTRUE(use_data) && !is.null(observed_data)) {
      min_design <- .restrict_to_observed_factors(min_design, cur_design,
                                                  observed_data,
                                                  remove_subjects = TRUE)
    }
    dadm <- design_model(
      min_design,
      cur_design,model,add_acc=FALSE,verbose=FALSE,rt_check=FALSE,compress=FALSE, add_da = add_da,
      all_cells_dm = all_cells_dm,
      drop_unobserved = isTRUE(use_data) && !is.null(observed_data))
    sampled_p_names <- attr(dadm,"sampled_p_names")
    if(length(design) != 1){
      map_list[[cur_name]] <- lapply(attributes(dadm)$designs,function(x){x[,,drop=FALSE]})
      sampled_p_names <- paste(cur_name, sampled_p_names, sep = "|")
    }
    out <- c(out, stats::setNames(numeric(length(sampled_p_names)),sampled_p_names))
    if(length(design) == 1){
      if (doMap) attr(out,"map") <- lapply(attributes(dadm)$designs,function(x){x[,,drop=FALSE]})
    }
  }
  if(length(design) != 1) attr(out, "map") <- map_list
  if(!add_da & any(duplicated(names(out)))) stop("duplicate parameter names found! Usually this happens when joint designs share indicator names")
  if(!is.null(group_design)){
    map <- attr(out, "map")
    par_names <- names(out)
    par_names <- add_group_par_names(par_names, group_design)
    out <- setNames(numeric(length(par_names)), par_names)
    attr(out, "map") <- map
  }
  return(out)
}

unique_rows_by_index <- function(df, columns) {
  # Keep only valid columns that are in the dataframe
  valid_columns <- intersect(columns, names(df))

  if (length(valid_columns) == 0) {
    # No valid columns to check uniqueness — return first row
    return(df[1, , drop = FALSE])
  }

  # Compute uniqueness based on valid columns
  unique_idx <- !duplicated(df[, valid_columns, drop = FALSE])

  # Return full rows from original dataframe
  df[unique_idx, , drop = FALSE]
}

#' Summary method for emc.design objects
#'
#' Prints a summary of the design object, including sampled parameters and design matrices.
#' For continuous covariates just prints one row, instead of all covariates.
#'
#' @param object An object of class `emc.design` containing the design to summarize
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the design matrices
#' @export
summary.emc.design <- function(object, ...){
  p_vector <- sampled_pars(object, doMap = TRUE)
  cat("\n Sampled Parameters: \n")
  print(names(p_vector))
  cat("\n Design Matrices: \n")
  map_out <- sampled_pars(object,add_da = TRUE, doMap = TRUE, data = list(...)$data)
  print_map <-
  print(lapply(attr(map_out, "map"), unique_rows_by_index,
               c(names(object$Ffactors), names(object$Ffunctions), 'lM', 'lR')), row.names = FALSE)
  return(invisible(map_out))
}

#' @export
print.emc.design <- function(x, ...){
  if("Ffactors" %in% names(x)){
    x <- list(x)
  }
  for(i in 1:length(x)){
    for(j in 1:length(x[[i]]$Flist)){
      cat(deparse(x[[i]]$Flist[j][[1]]), "\n")
    }
  }
}

#' @rdname plot_design
#' @export
plot_design.emc.design <- function(x, data = NULL, factors = NULL, plot_factor = NULL, n_data_sim = 10,
                            p_vector = NULL, functions = NULL, ...){
  if(is.null(p_vector)) stop("p_vector must be supplied if only the design is given")
  plot(x, p_vector, data = data, factors = factors, plot_factor = plot_factor, n_data_sim = n_data_sim,
       functions = functions, ...)
}


#' Plot method for emc.design objects
#'
#' Makes design illustration by plotting simulated data based on the design
#'
#' @param x An object of class `emc.design` containing the design to plot
#' @param p_vector A named vector of parameter values to use for data generation
#' @param data Optional data frame to overlay on the design plot. If NULL, data will be simulated.
#' @param factors Character vector. Factors to use for varying parameters in the plot
#' @param plot_factor Optional character. Make separate plots for each level of this factor
#' @param n_data_sim Integer. If data is NULL, number of simulated datasets to generate for the plot. Default is 10.
#' @param functions Optional named list of functions that create additional columns in the data
#' @param ... Additional arguments passed to `make_design_plot`
#' @return No return value, called for side effect of plotting
#' @export
plot.emc.design <- function(x, p_vector, data = NULL, factors = NULL, plot_factor = NULL, n_data_sim = 10,
                            functions = NULL, ...){
  if(!"Ffactors" %in% names(x)){
    if(length(x) != 1){
      stop("Current design type not supported for plotting")
    } else{
      x <- x[[1]]
    }
  }
  x$Ffunctions <- c(x$Ffunctions, functions)
  # Get a mapped parameter for each cell of the design
  pars <- mapped_pars(x, p_vector)
  if(is.null(data)){
    data <- vector("list", n_data_sim)
    # If no data is supplied generate some data sets
    for(i in 1:n_data_sim){
      data[[i]] <- make_data(p_vector, design = x, n_trials = 50)
    } # and bind them back together
    data <- do.call(rbind, data)
  }
  data <- data[!is.na(data$rt) & !is.infinite(data$rt),]
  data <- design_model(data, x, compress = FALSE, rt_resolution = 1e-15)

  if(is.null(x$model()$c_name)) stop("Current design type not supported for plotting")
  type <- ifelse(x$model()$c_name == "DDM", "DDM", ifelse(x$model()$c_name == "LNR", "LNR", "race"))
  within_noise <- ifelse(x$model()$c_name == "LBA", FALSE, TRUE)
  # Split only relevant for DDM
  dots <- add_defaults(list(...), split = "R", within_noise = within_noise, plot_legend = TRUE)
  if(type != "DDM"){
    dots$split = NULL
    data <- data[data$winner,]
  }
  do.call(make_design_plot, c(list(data = data, pars = pars, factors = factors, main = dots$main,
                   plot_factor = plot_factor,
                   type = type), fix_dots(dots, make_design_plot)))
}


#' @exportS3Method
sampled_pars.default <- function(x,group_design=NULL,doMap=FALSE, add_da = FALSE,
                                 all_cells_dm = FALSE, data = NULL,
                                 use_data = TRUE){
  if(is.null(x)) return(NULL)
  if(!is.null(attr(x, "custom_ll"))){
    pars <- numeric(length(attr(x,"sampled_p_names")))
    names(pars) <- attr(x,"sampled_p_names")
    return(pars)
  }
  if(!is.null(x$Ffactors)){
    x <- list(x)
    class(x) <- "emc.design"
  }
  out <- sampled_pars.emc.design(x, group_design = group_design, doMap = doMap,
                                 add_da = add_da, all_cells_dm = all_cells_dm,
                                 data = data, use_data = use_data)
  return(out)
}
