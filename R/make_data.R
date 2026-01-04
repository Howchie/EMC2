resolve_bound <- function(bound_name, supplied, default) {
    if (!is.null(supplied)) return(supplied)
    if (!is.null(data) && bound_name %in% colnames(data)) return(data[[bound_name]])
    if (!is.null(data) && !is.null(attr(data, bound_name))) return(attr(data, bound_name))
    default
}


bound_to_numeric <- function(x, bound_name = "bound") {
  if (is.null(x)) return(NULL)
  if (is.numeric(x)) return(x)
  out <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(out)) && any(!is.na(x))) {
    stop(paste0(bound_name, " must be numeric (or coercible to numeric)."))
  }
  names(out) <- names(x)
  out
}

expand_bound_rowwise <- function(b, data, bound_name = "bound") {
  if (any(names(data)==bound_name)) data$bound_name
  b <- bound_to_numeric(b, bound_name)
  n <- nrow(data)

  if (length(b) == 1) {
    return(rep(b, n))
  }
  if (length(b) == n) {
    return(b)
  }

  if (!is.null(names(b)) && "subjects" %in% names(data)) {
    present <- unique(as.character(data$subjects))
    if (!all(present %in% names(b))) {
      stop(paste0("Subject-wise ", bound_name, " must be named for all subjects present in data$subjects."))
    }
    return(unname(b[as.character(data$subjects)]))
  }

  stop(paste0(bound_name, " must be scalar, length nrow(data), or a subject-named vector."))
}


make_missing <- function(data, LT = 0, UT = Inf, LC = 0, UC = Inf,
                         LCresponse = FALSE, UCresponse = FALSE,
                         LCdirection = TRUE, UCdirection = TRUE,
                         pContaminant=NULL) {

#force_direction=FALSE,force_response=FALSE

  LT <- resolve_bound("LT", LT, 0)
  UT <- resolve_bound("UT", UT, Inf)
  LC <- resolve_bound("LC", LC, 0)
  UC <- resolve_bound("UC", UC, Inf)

  # Ensure truncation/censoring bounds are row-wise numeric columns (scalar,
  # subject-wise, or trial-wise)
  data$LT <- expand_bound_rowwise(LT, data, "LT")
  data$UT <- expand_bound_rowwise(UT, data, "UT")
  data$LC <- expand_bound_rowwise(LC, data, "LC")
  data$UC <- expand_bound_rowwise(UC, data, "UC")

  isgng <- data$R == "nogo"
  if (any(isgng)) {
    if (!UCdirection) {
      stop("UCdirction cannot be FALSE for a go/nogo model")
    }
    data[isgng,"UC"] <- 0
    data[isgng,"LC"] <- 0
    data[isgng,"UT"] <- Inf
    data[isgng,"LT"] <- 0
    # Wont work if go trials don't correspond
    UCdirection <- TRUE
    UCresponse <- TRUE
  }

  # Only keep trials in (LT, UT) or infinite
  pick <- is.infinite(data$rt) | (data$rt > data$LT & data$rt < data$UT)
  pick[is.na(pick)] <- TRUE
  data <- data[pick, ]

  # Lower censoring
  pickL <- data$rt < data$LC
  pickL[is.na(pickL)] <- FALSE
  if (LCdirection) data$rt[pickL] <- -Inf else data$rt[pickL] <- NA
  if (!LCresponse) data$R[pickL] <- NA

  # Upper censoring
  pickU <- data$rt > data$UC
  pickU[is.na(pickU)] <- FALSE
  if (UCdirection) data$rt[pickU] <- Inf else data$rt[pickU] <- NA
  if (!UCresponse) data$R[pickU] <- NA

  if (!is.null(pContaminant)) {
    contam <- rbinom(nrow(data), 1, pContaminant) == 1
    data[contam, "rt"] <- Inf
    data[contam, "R"] <- NA
  }

  data
}

#' Simulate Data
#'
#' Simulates data based on a model design and a parameter vector (`p_vector`) by one of two methods:
#' 1) Creating a fully crossed and balanced design specified by the design,
#' with number of trials per cell specified by the `n_trials` argument
#' 2) Using the design of a data frame supplied, which allows creation
#' of unbalanced and other irregular designs, and replacing previous data with
#' simulated data
#'
#' To create data for multiple subjects see ``?make_random_effects()``.
#'
#' @param parameters parameter vector used to simulate data.
#' Can also be a matrix with one row per subject (with corresponding row names)
#' or an emc object with sampled parameters
#' (in which case posterior medians of `alpha` are used to simulate data)
#' @param design Design list created by ``design()``
#' @param n_trials Integer. If ``data`` is not supplied, number of trials to create per design cell
#' @param data Data frame. If supplied, the factors are taken from the data. Determines the number of trials per level of the design factors and can thus allow for unbalanced designs
#' @param expand Integer. Replicates the ``data`` (if supplied) expand times to increase number of trials per cell.
#' @param staircase Default NULL, used with stop-signal paradigm simulation to specify a staircase
#' algorithm. If non-null and a list then passed through as is, if not it is assigned the
#' default list structure: list(p=.25,SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)
#' @param functions List of functions you want to apply to the data generation.
#' @param LT Lower truncation bound below which data are removed (scalar, subject-wise vector, or trial-wise vector).
#' @param UT Upper truncation bound above which data are removed (scalar, subject-wise vector, or trial-wise vector).
#' @param LC Lower censoring bound (scalar, subject-wise vector, or trial-wise vector).
#' @param UC Upper censoring bound (scalar, subject-wise vector, or trial-wise vector).
#' @param LCresponse Logical. If FALSE, set responses to NA on lower-censored trials.
#' @param UCresponse Logical. If FALSE, set responses to NA on upper-censored trials.
#' @param LCdirection Logical. If TRUE, lower-censored RTs are coded as -Inf; if FALSE, as NA.
#' @param UCdirection Logical. If TRUE, upper-censored RTs are coded as Inf; if FALSE, as NA.
#' @param ... Additional optional arguments
#' @return A data frame with simulated data
#' @examples
#' # First create a design
#' design_DDMaE <- design(factors = list(S = c("left", "right"),
#'                                            E = c("SPD", "ACC"),
#'                                            subjects = 1:30),
#'                             Rlevels = c("left", "right"), model = DDM,
#'                             formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                             constants=c(s=log(1)))
#' # Then create a p_vector:
#' parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_EACC=log(2), t0=log(.2),
#'               Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#'
#' # Now we can simulate data
#' data <- make_data(parameters, design_DDMaE, n_trials = 30)
#'
#' # We can also simulate data based on a specific dataset
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                             formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                             constants=c(s=log(1)))
#' parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'               t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#'
#' data <- make_data(parameters, design_DDMaE, data = forstmann)
#' @export

make_data <- function(parameters,design = NULL,n_trials=NULL,data=NULL,expand=1, staircase = NULL,
                      functions = NULL,
                      LT=NULL,LC=NULL,UC=NULL,UT=NULL,
                      LCresponse = FALSE, UCresponse = FALSE,
                      LCdirection = TRUE, UCdirection = TRUE,...)

#' @param force_direction Logical. If TRUE (default), do not infer censor direction from `data$rt`.
#' @param force_response Logical. If TRUE (default), do not infer whether censored trials keep responses from `data$R`.

{
  if (!is.null(staircase)){
    staircase <- check_staircase(staircase)
  }
  check_bounds <- FALSE

  return_Ffunctions <- FALSE
  post_functions <- NULL
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(is(parameters, "emc")){
    if(is.null(design)) design <- get_design(parameters)[[1]] # Currently not supported for multiple designs
    if(is.null(data)) data <- get_data(parameters)
    parameters <- do.call(rbind, credint(parameters, probs = 0.5, selection = "alpha", by_subject = TRUE))
  }

  # Make sure parameters are in the right format, either matrix or vector
  sampled_p_names <- names(sampled_pars(design))
  if(is.null(dim(parameters))){
    if(is.null(names(parameters))) names(parameters) <- sampled_p_names
  } else {
    if(!is.null(data)){
      if(nrow(parameters) == length(unique(data$subjects))){
        design$Ffactors$subjects <- unique(data$subjects)
      }
    }
    if(length(rownames(parameters)) != length(design$Ffactors$subjects)){
      stop("input parameter matrix must have number of rows equal to number of subjects specified in design")
    }
    if(is.null(colnames(parameters))) colnames(parameters) <- sampled_p_names
    rownames(parameters) <- design$Ffactors$subjects
  }

  if(!is.null(attr(design, "custom_ll"))){
    data <- list()
    for(i in 1:nrow(parameters)){
      data[[i]] <- design$model()$rfun(parameters[i,], n_trials = n_trials, subject = i)
    }
    return(do.call(rbind, data))
  }

  model <- design$model

  if(grepl("MRI", model()$type)){
    return(make_data_wrapper_MRI(parameters, data, design))
  }
  if(is.data.frame(parameters)) parameters <- as.matrix(parameters)
  if (!is.matrix(parameters)) parameters <- make_pmat(parameters,design)
  if ( is.null(data) ) {
    design$Ffactors$subjects <- rownames(parameters)
    if ( is.null(n_trials) )
      stop("If data is not provided need to specify number of trials")
    design_in <- design
    design_in$Fcovariates <- design_in$Fcovariates[!design$Fcovariates %in% names(functions)]
    # `make_data()` expands accumulators later via `add_accumulators()`/`design_model()`.
    # Avoid evaluating design$Ffunctions here because they may depend on accumulator
    # columns (e.g., `lR`) that are not present until after accumulator expansion.
    data <- minimal_design(design_in, covariates = list(...)$covariates,
                              drop_subjects = F, n_trials = n_trials,
                            add_acc = FALSE, do_functions = FALSE,drop_R = FALSE)
  } else {
		data <- add_trials(data[order(data$subjects),])
  }

  if(!is.null(functions)){
    for(i in 1:length(functions)){
      data[[names(functions)[i]]] <- functions[[i]](data)
    }
  }
  if (!is.factor(data$subjects)) data$subjects <- factor(data$subjects)
  if (!is.null(model)) {
    if (!is.function(model)) stop("model argument must  be a function")
    if ( is.null(model()$p_types) ) stop("model()$p_types must be specified")
    if ( is.null(model()$Ttransform) ) stop("model()$Ttransform must be specified")
  }

  simulate_unconditional_on_data <- return_trialwise_parameters <- FALSE
  dots_local <- list(...)
  if (isFALSE(dots_local$conditional_on_data)) {
    simulate_unconditional_on_data <- TRUE
  } else if (!is.null(dots_local$conditional_on_data)) {
    simulate_unconditional_on_data <- !isTRUE(dots_local$conditional_on_data)
  }
  return_trialwise_parameters <- isTRUE(dots_local$return_trialwise_parameters)

  ## For both conditional and unconditional simulations...
  pars <- t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
  pars <- add_constants(pars,design$constants)
  if (!is.null(staircase)) {
    attr(pars, "staircase") <- staircase
  }
  if (simulate_unconditional_on_data) {
    # unconditional simulation must produce R/rt already; do not fall through to rfun below
    if (expand > 1) {
      data_list <- vector("list", expand)
      twp_last <- NULL
      for (rep_i in seq_len(expand)) {
        res_i <- make_data_unconditional(
          data = data, pars = pars, design = design, model = model,
          return_trialwise_parameters = (rep_i == expand && return_trialwise_parameters)
        )
        data_list[[rep_i]] <- cbind(rep = rep_i, res_i$data)
        if (rep_i == expand) twp_last <- res_i$trialwise_parameters
      }
      data <- do.call(rbind, data_list)
      rownames(data) <- NULL
      trialwise_parameters <- twp_last
    } else {
      res <- make_data_unconditional(
        data = data, pars = pars, design = design, model = model,
        return_trialwise_parameters = return_trialwise_parameters
      )
      data <- res$data
      trialwise_parameters <- res$trialwise_parameters
    }
    attr(data, "p_vector") <- parameters
    if (return_trialwise_parameters) attr(data, "trialwise_parameters") <- trialwise_parameters
    return(data)
  } else {
    data <- design_model(
      add_accumulators(data, design$matchfun, simulate = TRUE, type = model()$type, Fcovariates = design$Fcovariates),
      design, model,
      add_acc = FALSE, compress = FALSE, verbose = FALSE,
      rt_check = FALSE
    )
    pars <- map_p(pars, data, model(), return_trialwise_parameters)

    if (!is.null(model()$trend)) {
      phases <- vapply(model()$trend, function(x) x$phase, character(1))
      if (any(phases == "pretransform")) {
        # apply only pretransform trends and remove their trend parameters
        pars <- prep_trend_phase(
          data, model()$trend, pars, "pretransform",
          return_trialwise_parameters
        )
      }
    }
    pars <- do_transform(pars, model()$transform)
    if (!is.null(model()$trend)) {
      phases <- vapply(model()$trend, function(x) x$phase, character(1))
      if (any(phases == "posttransform")) {
        # apply only posttransform trends and remove their trend parameters
        pars <- prep_trend_phase(
          data, model()$trend, pars, "posttransform",
          return_trialwise_parameters
        )
      }
    }
    if (return_trialwise_parameters)
      trialwise_parameters <- cbind(pars, attr(pars, "trialwise_parameters"))

    pars <- model()$Ttransform(pars, data)
    if (is.null(optionals$nobound)) {
      pars <- add_bound(pars, model()$bound, data$lR)
    } else {
      attr(pars, "ok") <- rep(TRUE, nrow(pars))
    }

    pars_ok <- attr(pars, "ok")
    if (mean(!pars_ok) > .1) {
      warning("More than 10% of parameter values fall out of model bounds, see <model_name>$bounds()")
      return(FALSE)
    }
  }
  if (expand>1) {
      data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
                    data.frame(lapply(data,rep,times=expand)))
      pars <- apply(pars,2,rep,times=expand)
  }
  if (!is.null(staircase)) {
      attr(pars, "staircase") <- staircase
  }
  if (any(names(data)=="RACE")) {
      Rrt <- RACE_rfun(data, pars, model)
  } else Rrt <- model()$rfun(data,pars)
  if (!any(dimnames(pars)[[2]]=="pContaminant")) pContaminant <- NULL else
      pContaminant <- pars[,"pContaminant"][data$lR==levels(data$lR)[1]]

  dropNames <- c("lR","lM")
  if (!return_Ffunctions && !is.null(design$Ffunctions))
    dropNames <- c(dropNames,names(design$Ffunctions))
  if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
  data <- data[,!(names(data) %in% dropNames)]
  for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[, i]

  data <- make_missing(data,LT, UT, LC, UC,
    LCresponse, UCresponse, LCdirection, UCdirection,pContaminant)

  attr(data,"p_vector") <- parameters;
  if(!is.null(post_functions)){
    for(i in 1:length(post_functions)){
      data[[names(post_functions)[i]]] <- post_functions[[i]](data)
    }
  }
  if(return_trialwise_parameters) attr(data, 'trialwise_parameters') <- trialwise_parameters
  data
}

# make_data <- function(parameters,design = NULL,n_trials=NULL,data=NULL,expand=1, staircase = NULL,
#                       functions = NULL,
#                       LT=NULL,LC=NULL,UC=NULL,UT=NULL,
#                       LCresponse = FALSE, UCresponse = FALSE,
#                       LCdirection = TRUE, UCdirection = TRUE,
#                       force_direction = TRUE, force_response = TRUE,
#                       rtContaminantNA = FALSE,...)
# {
#   if (!is.null(staircase)){
#     staircase <- check_staircase(staircase)
#   }
#   check_bounds <- FALSE
#
#   resolve_bound <- function(bound_name, supplied, default) {
#     if (!is.null(supplied)) return(supplied)
#     if (!is.null(data) && bound_name %in% colnames(data)) return(data[[bound_name]])
#     if (!is.null(data) && !is.null(attr(data, bound_name))) return(attr(data, bound_name))
#     # if (!is.null(design) && !is.null(design[[bound_name]])) return(design[[bound_name]])
#     default
#   }
#
#   LT <- resolve_bound("LT", LT, 0)
#   UT <- resolve_bound("UT", UT, Inf)
#   LC <- resolve_bound("LC", LC, 0)
#   UC <- resolve_bound("UC", UC, Inf)
#
#
#   return_Ffunctions <- FALSE
#   post_functions <- NULL
#   optionals <- list(...)
#   for (name in names(optionals) ) {
#     assign(name, optionals[[name]])
#   }
#   if(is(parameters, "emc")){
#     if(is.null(design)) design <- get_design(parameters)[[1]] # Currently not supported for multiple designs
#     if(is.null(data)) data <- get_data(parameters)
#     parameters <- do.call(rbind, credint(parameters, probs = 0.5, selection = "alpha", by_subject = TRUE))
#   }
#
#   # Make sure parameters are in the right format, either matrix or vector
#   sampled_p_names <- names(sampled_pars(design))
#   if(is.null(dim(parameters))){
#     if(is.null(names(parameters))) names(parameters) <- sampled_p_names
#   } else {
#     if(!is.null(data)){
#       if(nrow(parameters) == length(unique(data$subjects))){
#         design$Ffactors$subjects <- unique(data$subjects)
#       }
#     }
#     if(length(rownames(parameters)) != length(design$Ffactors$subjects)){
#       stop("input parameter matrix must have number of rows equal to number of subjects specified in design")
#     }
#     if(is.null(colnames(parameters))) colnames(parameters) <- sampled_p_names
#     rownames(parameters) <- design$Ffactors$subjects
#   }
#
#   if(!is.null(attr(design, "custom_ll"))){
#     data <- list()
#     for(i in 1:nrow(parameters)){
#       data[[i]] <- design$model()$rfun(parameters[i,], n_trials = n_trials, subject = i)
#     }
#     return(do.call(rbind, data))
#   }
#
#   model <- design$model
#
#   if(grepl("MRI", model()$type)){
#     return(make_data_wrapper_MRI(parameters, data, design))
#   }
#   if(is.data.frame(parameters)) parameters <- as.matrix(parameters)
#   if (!is.matrix(parameters)) parameters <- make_pmat(parameters,design)
#   if ( is.null(data) ) {
#     design$Ffactors$subjects <- rownames(parameters)
#     if ( is.null(n_trials) )
#       stop("If data is not provided need to specify number of trials")
#     design_in <- design
#     design_in$Fcovariates <- design_in$Fcovariates[!design$Fcovariates %in% names(functions)]
#     # `make_data()` expands accumulators later via `add_accumulators()`/`design_model()`.
#     # Avoid evaluating design$Ffunctions here because they may depend on accumulator
#     # columns (e.g., `lR`) that are not present until after accumulator expansion.
#     data <- minimal_design(design_in, covariates = list(...)$covariates,
#                               drop_subjects = F, n_trials = n_trials,
#                             add_acc = FALSE, do_functions = FALSE,
#                             drop_R = FALSE, UC = UC, UT = UT, LC = LC, LT = LT)
#   } else {
#     # bounds handled below after add_trials ordering
# 		if (!force_direction) {
# 		  ok <- data$rt==-Inf; ok[is.na(ok)] <- FALSE
# 		  if (missing(LCdirection)) LCdirection <- any(ok)
# 		  ok <- data$rt==Inf; ok[is.na(ok)] <- FALSE
# 		  if (missing(UCdirection)) UCdirection <- any(ok)
# 		}
# 		if (!force_response) {
# 		  if (!any(is.infinite(data$rt)) & any(is.na(data$R))) {
# 			if (missing(LCresponse)) LCresponse <- FALSE
# 			if (missing(UCresponse)) UCresponse <- FALSE
# 		  } else {
# 			if (missing(LCresponse)) {
# 			  ok <- data$rt == -Inf
# 			  bad <- is.na(ok)
# 			  LCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
# 			}
# 			if (missing(UCresponse)) {
# 			  ok <- data$rt == Inf
# 			  bad <- is.na(ok)
# 			  UCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
# 			}
# 		  }
# 		}
# 		data <- add_trials(data[order(data$subjects),])
#   }
#
#   # Ensure truncation/censoring bounds are row-wise numeric columns (scalar,
#   # subject-wise, or trial-wise)
#   data$LT <- expand_bound_rowwise(LT, data, "LT")
#   data$UT <- expand_bound_rowwise(UT, data, "UT")
#   data$LC <- expand_bound_rowwise(LC, data, "LC")
#   data$UC <- expand_bound_rowwise(UC, data, "UC")
#
#   if(!is.null(functions)){
#     for(i in 1:length(functions)){
#       data[[names(functions)[i]]] <- functions[[i]](data)
#     }
#   }
#   if (!is.factor(data$subjects)) data$subjects <- factor(data$subjects)
#   if (!is.null(model)) {
#     if (!is.function(model)) stop("model argument must  be a function")
#     if ( is.null(model()$p_types) ) stop("model()$p_types must be specified")
#     if ( is.null(model()$Ttransform) ) stop("model()$Ttransform must be specified")
#   }
#
#   simulate_unconditional_on_data <- return_trialwise_parameters <- FALSE
#   dots_local <- list(...)
#   if (isFALSE(dots_local$conditional_on_data)) {
#     simulate_unconditional_on_data <- TRUE
#   } else if (!is.null(dots_local$conditional_on_data)) {
#     simulate_unconditional_on_data <- !isTRUE(dots_local$conditional_on_data)
#   }
#   return_trialwise_parameters <- isTRUE(dots_local$return_trialwise_parameters)
#
#   ## For both conditional and unconditional simulations...
#   pars <- t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
#   pars <- add_constants(pars,design$constants)
#   if (!is.null(staircase)) {
#     attr(pars, "staircase") <- staircase
#   }
#   if (simulate_unconditional_on_data) {
#     # unconditional simulation must produce R/rt already; do not fall through to rfun below
#     if (expand > 1) {
#       data_list <- vector("list", expand)
#       twp_last <- NULL
#       for (rep_i in seq_len(expand)) {
#         res_i <- make_data_unconditional(
#           data = data, pars = pars, design = design, model = model,
#           return_trialwise_parameters = (rep_i == expand && return_trialwise_parameters)
#         )
#         data_list[[rep_i]] <- cbind(rep = rep_i, res_i$data)
#         if (rep_i == expand) twp_last <- res_i$trialwise_parameters
#       }
#       data <- do.call(rbind, data_list)
#       rownames(data) <- NULL
#       trialwise_parameters <- twp_last
#     } else {
#       res <- make_data_unconditional(
#         data = data, pars = pars, design = design, model = model,
#         return_trialwise_parameters = return_trialwise_parameters
#       )
#       data <- res$data
#       trialwise_parameters <- res$trialwise_parameters
#     }
#
#     # # Add trial-level pContaminant if present so it stays aligned through truncation
#     # if (!is.null(colnames(pars)) && "pContaminant" %in% colnames(pars)) {
#     #   dm_pc <- design_model(
#     #     add_accumulators(data, design$matchfun, simulate = TRUE, type = model()$type, Fcovariates = design$Fcovariates),
#     #     design, model,
#     #     add_acc = FALSE, compress = FALSE, verbose = FALSE,
#     #     rt_check = FALSE
#     #   )
#     #   pm_pc <- map_p(pars, dm_pc, model(), FALSE)
#     #   if (!is.null(model()$trend)) {
#     #     phases <- vapply(model()$trend, function(x) x$phase, character(1))
#     #     if (any(phases == "pretransform")) {
#     #       pm_pc <- prep_trend_phase(dm_pc, model()$trend, pm_pc, "pretransform", FALSE)
#     #     }
#     #   }
#     #   pm_pc <- do_transform(pm_pc, model()$transform)
#     #   if (!is.null(model()$trend)) {
#     #     phases <- vapply(model()$trend, function(x) x$phase, character(1))
#     #     if (any(phases == "posttransform")) {
#     #       pm_pc <- prep_trend_phase(dm_pc, model()$trend, pm_pc, "posttransform", FALSE)
#     #     }
#     #   }
#     #   pm_pc <- model()$Ttransform(pm_pc, dm_pc)
#     #   if (is.null(optionals$nobound)) pm_pc <- add_bound(pm_pc, model()$bound, dm_pc$lR)
#     #   if (any(dimnames(pm_pc)[[2]] == "pContaminant") && any(pm_pc[, "pContaminant"] > 0)) {
#     #     pc <- pm_pc[dm_pc$lR == levels(dm_pc$lR)[1], "pContaminant"]
#     #     data$.pContaminant <- pc
#     #   }
#     # }
#     #
#     # # Apply truncation/censoring to unconditional simulation output
#     # data <- make_missing(data, data$LT, data$UT, data$LC, data$UC,
#     #                      LCresponse, UCresponse, LCdirection, UCdirection)
#
#     # Apply contamination after truncation so indices stay aligned
#     if (!is.null(data$.pContaminant)) {
#       pc2 <- data$.pContaminant
#       data$.pContaminant <- NULL
#       if (!any(is.infinite(data$rt)) & any(is.na(data$R)))
#         stop("Cannot have contamination and censoring with no direction and response")
#       contam <- rbinom(nrow(data), 1, pc2) == 1
#       data[contam,"R"] <- NA
#
#       has_LC <- any(is.finite(data$LC) & data$LC != 0, na.rm = TRUE)
#       has_UC <- any(is.finite(data$UC), na.rm = TRUE)
#       if (has_LC | has_UC) { # censoring
#         if (has_LC & has_UC & (LCdirection & UCdirection) & !rtContaminantNA) {
#           stop("Cannot have contamination with a mixture of censor directions")
#         }
#         if (
#           rtContaminantNA &&
#           ((has_LC && !LCresponse && !LCdirection) || (has_UC && !UCresponse && !UCdirection))
#         ) {
#           stop("Cannot have contamination and censoring with no direction and response")
#         }
#
#         if (rtContaminantNA || (!LCdirection && !UCdirection)) {
#           data[contam, "rt"] <- NA
#         } else if (has_LC && !has_UC) {
#           data[contam, "rt"] <- if (LCdirection) -Inf else NA
#         } else if (has_UC && !has_LC) {
#           data[contam, "rt"] <- if (UCdirection) Inf else NA
#         } else {
#           if (LCdirection) data[contam, "rt"] <- -Inf else
#             if (UCdirection) data[contam, "rt"] <- Inf else data[contam, "rt"] <- NA
#         }
#       } else {
#         data[contam, "rt"] <- NA
#       }
#     }
#     attr(data, "p_vector") <- parameters
#     if (return_trialwise_parameters) attr(data, "trialwise_parameters") <- trialwise_parameters
#     return(data)
#   } else {
#     data <- design_model(
#       add_accumulators(data, design$matchfun, simulate = TRUE, type = model()$type, Fcovariates = design$Fcovariates),
#       design, model,
#       add_acc = FALSE, compress = FALSE, verbose = FALSE,
#       rt_check = FALSE
#     )
#     pars <- map_p(pars, data, model(), return_trialwise_parameters)
#
#     if (!is.null(model()$trend)) {
#       phases <- vapply(model()$trend, function(x) x$phase, character(1))
#       if (any(phases == "pretransform")) {
#         # apply only pretransform trends and remove their trend parameters
#         pars <- prep_trend_phase(
#           data, model()$trend, pars, "pretransform",
#           return_trialwise_parameters
#         )
#       }
#     }
#     pars <- do_transform(pars, model()$transform)
#     if (!is.null(model()$trend)) {
#       phases <- vapply(model()$trend, function(x) x$phase, character(1))
#       if (any(phases == "posttransform")) {
#         # apply only posttransform trends and remove their trend parameters
#         pars <- prep_trend_phase(
#           data, model()$trend, pars, "posttransform",
#           return_trialwise_parameters
#         )
#       }
#     }
#     if (return_trialwise_parameters) trialwise_parameters <- cbind(pars, attr(pars, "trialwise_parameters"))
#
#     pars <- model()$Ttransform(pars, data)
#     if (is.null(optionals$nobound)) {
#       pars <- add_bound(pars, model()$bound, data$lR)
#     } else {
#       attr(pars, "ok") <- rep(TRUE, nrow(pars))
#     }
#
#     pars_ok <- attr(pars, "ok")
#     if (mean(!pars_ok) > .1) {
#       warning("More than 10% of parameter values fall out of model bounds, see <model_name>$bounds()")
#       return(FALSE)
#     }
#   }
#   if (expand>1) {
#       data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
#                     data.frame(lapply(data,rep,times=expand)))
#       pars <- apply(pars,2,rep,times=expand)
#   }
#   if (!is.null(staircase)) {
#       attr(pars, "staircase") <- staircase
#   }
#   # ZH added for contaminant miss handling (compute after `expand` so it stays aligned)
#   pc <- NULL
#   if (any(dimnames(pars)[[2]] == "pContaminant") && any(pars[, "pContaminant"] > 0)) {
#     pc <- pars[data$lR == levels(data$lR)[1], "pContaminant"]
#   }
#   # if (model()$type=="GNG") {
#   #   if (any(names(data)=="RACE")) {
#   #     Rrt <- GNG_rfun(data, pars, model)
#   #   } else Rrt <- model()$rfun(data,pars)
#   # } else
#   if (any(names(data)=="RACE")) {
#       Rrt <- RACE_rfun(data, pars, model)
#   } else Rrt <- model()$rfun(data,pars)
#   dropNames <- c("lR","lM")
#   if (!return_Ffunctions && !is.null(design$Ffunctions))
#     dropNames <- c(dropNames,names(design$Ffunctions))
#   if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
#   data <- data[,!(names(data) %in% dropNames)]
#   for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[, i]
#   if (!is.null(pc)) data$.pContaminant <- pc
#   # ZH Added for censoring/truncation
#   # For censoring/truncation we want row-wise bounds: use the
#   # current LT/UT/LC/UC columns, which have already been set up.
#   data <- make_missing(
#     data[, names(data) != "winner"],
#     data$LT, data$UT, data$LC, data$UC,
#     LCresponse, UCresponse, LCdirection, UCdirection
#   )
#   # ZH added for contaminant miss handling
#   if ( !is.null(pc) ) {
#     pc2 <- data$.pContaminant
#     data$.pContaminant <- NULL
#     if (!any(is.infinite(data$rt)) & any(is.na(data$R)))
#       stop("Cannot have contamination and censoring with no direction and response")
#     contam <- rbinom(nrow(data), 1, pc2) == 1
#     data[contam,"R"] <- NA
#
#     has_LC <- any(is.finite(data$LC) & data$LC != 0, na.rm = TRUE)
#     has_UC <- any(is.finite(data$UC), na.rm = TRUE)
#     if (has_LC | has_UC) { # censoring
#       # if (has_LC & has_UC & (LCdirection & UCdirection) & !rtContaminantNA) {
#       #   stop("Cannot have contamination with a mixture of censor directions")
#       # }
#       # if (
#       #   rtContaminantNA &&
#       #   ((has_LC && !LCresponse && !LCdirection) || (has_UC && !UCresponse && !UCdirection))
#       # ) {
#       #   stop("Cannot have contamination and censoring with no direction and response")
#       # }
#       if (rtContaminantNA || (!LCdirection && !UCdirection)) {
#         data[contam, "rt"] <- NA
#       } else if (has_LC && !has_UC) {
#         data[contam, "rt"] <- if (LCdirection) -Inf else NA
#       } else if (has_UC && !has_LC) {
#         data[contam, "rt"] <- if (UCdirection) Inf else NA
#       } else {
#         if (LCdirection) data[contam, "rt"] <- -Inf else
#           if (UCdirection) data[contam, "rt"] <- Inf else data[contam, "rt"] <- NA
#       }
#     } else {
#       data[contam, "rt"] <- NA
#     }
#   }
#
#   attr(data,"p_vector") <- parameters;
#   if(!is.null(post_functions)){
#     for(i in 1:length(post_functions)){
#       data[[names(post_functions)[i]]] <- post_functions[[i]](data)
#     }
#   }
#   if(return_trialwise_parameters) attr(data, 'trialwise_parameters') <- trialwise_parameters
#   data
# }

RACE_rfun <- function(data, pars, model){
  Rrt <- matrix(ncol=2,nrow=dim(data)[1]/length(levels(data$lR)),
         dimnames=list(NULL,c("R","rt")))
  RACE <- data[data$lR==levels(data$lR)[1],"RACE"]
  ok <- as.numeric(data$lR) <= as.numeric(as.character(data$RACE))
  for (i in levels(RACE)) {
    pick <- data$RACE==i
    data_in <- data[pick & ok,]
    data_in$lR <- factor(data$lR[pick & ok])
    tmp <- pars[pick & ok,]
    attr(tmp, "ok") <- rep(T, nrow(tmp))
    if (!is.null(attr(pars, "staircase"))) attr(tmp, "staircase") <- attr(pars, "staircase")
    Rrti <- model()$rfun(data_in,tmp)
    Rrti$R <- as.numeric(Rrti$R)
    Rrt[RACE==i,] <- as.matrix(Rrti)
  }
  Rrt <- data.frame(Rrt)
  Rrt$R <- factor(Rrt$R, labels = levels(data$lR), levels = 1:length(levels(data$lR)))
  return(Rrt)
}

# GNG_rfun <- function(data, pars, model){
# 	if (!"nogo"%in%levels(data$lR)) {
# 		stop("Go/No-Go models must have a nogo level for Response")
# 	}
# 	if (any(is.infinite(data$UC))) {
# 		stop("UC must be finite for Go/No-Go Models")
# 	}
#   Rrt <- matrix(ncol=2,nrow=dim(data)[1]/length(levels(data$lR)),
#          dimnames=list(NULL,c("R","rt")))
#   RACE <- data[data$lR==levels(data$lR)[1],"RACE"]
#   D = data[data$lR==levels(data$lR)[1],"UC"]
#   ok <- as.numeric(data$lR) <= as.numeric(as.character(data$RACE))
#   for (i in levels(RACE)) {
#     pick <- data$RACE==i
#     data_in <- data[pick & ok,]
#     data_in$lR <- factor(data$lR[pick & ok])
#     tmp <- pars[pick & ok,]
#     attr(tmp, "ok") <- rep(T, nrow(tmp))
#     if (!is.null(attr(pars, "staircase"))) attr(tmp, "staircase") <- attr(pars, "staircase")
#     Rrti <- model()$rfun(data_in,tmp)
#     Rrti$R <- as.numeric(Rrti$R)
#     Rrt[RACE==i,] <- as.matrix(Rrti)
#   }
#   Rrt <- data.frame(Rrt)
#   Rrt$R <- factor(Rrt$R, labels = levels(data$lR), levels = 1:length(levels(data$lR)))
#   Rrt$R[Rrt$rt>D]="nogo"
#   Rrt$rt[Rrt$R=="nogo"]=Inf
#   return(Rrt)
# }

add_Ffunctions <- function(data,design)
  # Adds columns created by Ffunctions (if not already there)
{
  Fdf <- data.frame(lapply(design$Ffunctions,function(f){f(data)}))
  ok <- !(names(Fdf) %in% names(data))
  if (!any(ok)) data else
    data <-  cbind.data.frame(data,Fdf[,ok,drop=FALSE])
}

#' Generate Subject-Level Parameters
#'
#' Simulates subject-level parameters in the format required by ``make_data()``.
#'
#' @param design A design list. The design as specified by `design()`
#' @param group_means A numeric vector. The group level means for each parameter, in the same order as `sampled_pars(design)`
#' @param n_subj An integer. The number of subjects to generate parameters for. If `NULL` will be inferred from design
#' @param variance_proportion A double. Optional. If ``covariances`` are not specified, the variances will be created by multiplying the means by this number. The covariances will be 0.
#' @param covariances A covariance matrix. Optional. Specify the intended covariance matrix.
#'
#' @return A matrix of subject-level parameters.
#' @examples
#' # First create a design
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                             formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                             constants=c(s=log(1)))
#' # Then create a group-level means vector:
#' group_means =c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'                t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # Now we can create subject-level parameters
#' subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 19)
#'
#' # We can also define a covariance matrix to simulate from
#' subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 19,
#'              covariances = diag(.1, length(group_means)))
#'
#' # The subject level parameters can be used to generate data
#' make_data(subj_pars, design_DDMaE, n_trials = 10)
#' @export

make_random_effects <- function(design, group_means, n_subj = NULL, variance_proportion = .2, covariances = NULL){
  if(is.null(n_subj)){
    n_subj <- length(design$Ffactors$subjects)
    subnames <- design$Ffactors$subjects
  } else{
    subnames <- as.character(1:n_subj)
  }
  if(length(group_means) != length(sampled_pars(design))) stop("You must specify as many means as parameters in your design")
  if(is.null(covariances)) covariances <- diag(abs(group_means)*variance_proportion)
  random_effects <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=covariances)
  colnames(random_effects) <- names(sampled_pars(design))
  rownames(random_effects) <- subnames
  return(random_effects)
}
