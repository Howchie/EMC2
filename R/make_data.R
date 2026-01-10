# resolve_bound <- function(data, bound_name, supplied, default) {
#     if (!is.null(supplied)) return(supplied)
#     if (!is.null(data) && bound_name %in% colnames(data)) return(data[[bound_name]])
#     if (!is.null(data) && !is.null(attr(data, bound_name))) return(attr(data, bound_name))
#     default
# }
#
#
# bound_to_numeric <- function(x, bound_name = "bound") {
#   if (is.null(x)) return(NULL)
#   if (is.numeric(x)) return(x)
#   out <- suppressWarnings(as.numeric(as.character(x)))
#   if (all(is.na(out)) && any(!is.na(x))) {
#     stop(paste0(bound_name, " must be numeric (or coercible to numeric)."))
#   }
#   names(out) <- names(x)
#   out
# }
#
# expand_bound_rowwise <- function(b, data, bound_name = "bound") {
#   if (any(names(data)==bound_name)) data$bound_name
#   b <- bound_to_numeric(b, bound_name)
#   n <- nrow(data)
#
#   if (length(b) == 1) {
#     return(rep(b, n))
#   }
#   if (length(b) == n) {
#     return(b)
#   }
#
#   if (!is.null(names(b)) && "subjects" %in% names(data)) {
#     present <- unique(as.character(data$subjects))
#     if (!all(present %in% names(b))) {
#       stop(paste0("Subject-wise ", bound_name, " must be named for all subjects present in data$subjects."))
#     }
#     return(unname(b[as.character(data$subjects)]))
#   }
#
#   stop(paste0(bound_name, " must be scalar, length nrow(data), or a subject-named vector."))
# }
#
#
# make_missing <- function(data, LT = 0, UT = Inf, LC = 0, UC = Inf,
#                          LCresponse = FALSE, UCresponse = FALSE,
#                          LCdirection = TRUE, UCdirection = TRUE,
#                          pContaminant=NULL) {
#
# #force_direction=FALSE,force_response=FALSE
#
#   LT <- resolve_bound("LT", LT, 0)
#   UT <- resolve_bound("UT", UT, Inf)
#   LC <- resolve_bound("LC", LC, 0)
#   UC <- resolve_bound("UC", UC, Inf)
#
#   # Ensure truncation/censoring bounds are row-wise numeric columns (scalar,
#   # subject-wise, or trial-wise)
#   data$LT <- expand_bound_rowwise(LT, data, "LT")
#   data$UT <- expand_bound_rowwise(UT, data, "UT")
#   data$LC <- expand_bound_rowwise(LC, data, "LC")
#   data$UC <- expand_bound_rowwise(UC, data, "UC")
#
#   isgng <- data$R == "nogo"
#   if (any(isgng)) {
#     if (!UCdirection) {
#       stop("UCdirction cannot be FALSE for a go/nogo model")
#     }
#     data[isgng,"UC"] <- 0
#     data[isgng,"LC"] <- 0
#     data[isgng,"UT"] <- Inf
#     data[isgng,"LT"] <- 0
#     # Wont work if go trials don't correspond
#     UCdirection <- TRUE
#     UCresponse <- TRUE
#   }
#
#   # Only keep trials in (LT, UT) or infinite
#   pick <- is.infinite(data$rt) | (data$rt > data$LT & data$rt < data$UT)
#   pick[is.na(pick)] <- TRUE
#   data <- data[pick, ]
#
#   # Lower censoring
#   pickL <- data$rt < data$LC
#   pickL[is.na(pickL)] <- FALSE
#   if (LCdirection) data$rt[pickL] <- -Inf else data$rt[pickL] <- NA
#   if (!LCresponse) data$R[pickL] <- NA
#
#   # Upper censoring
#   pickU <- data$rt > data$UC
#   pickU[is.na(pickU)] <- FALSE
#   if (UCdirection) data$rt[pickU] <- Inf else data$rt[pickU] <- NA
#   if (!UCresponse) data$R[pickU] <- NA
#
#   if (!is.null(pContaminant)) {
#     contam <- rbinom(nrow(data), 1, pContaminant) == 1
#     data[contam, "rt"] <- Inf
#     data[contam, "R"] <- NA
#   }
#
#   data
# }



get_missing <- function(supplied, data, bound_name, default,type) {

  subjectwise <- FALSE
  if (is.function(supplied)) supplied <- supplied(data)
  if (is.null(supplied)) {
    if (bound_name %in% colnames(data)) bound <- data[[bound_name]] else
      bound <- default
  } else {
    if (type=="logical" & !is.logical(supplied))
        stop(bound_name," must be logical")
    if (type=="numeric" & !is.numeric(supplied))
        stop(bound_name," must be numeric")
    subjectwise <- all(hasName(supplied,levels(data$subjects)))
    if (!subjectwise) bound <- supplied else
      bound <- supplied[as.character(data$subjects)]
  }
  out <- rep(bound,length.out=nrow(data))
  attr(out,"subjectwise") <- subjectwise
  out
}

#' Add information about missing values to data and modify/filter accordingly.
#'
#' Columns corresponding to LC, UT, LC, UC, and pContaminant arguments are added
#' to the return, specifying, respectively, if a row is subject to lower or upper
#' truncation, lower or upper censoring, or contamination with the given probability.
#' Contamination makes R = NA and the rt = Inf, is applied before any other
#' missing operation and any row subject to contamination is not affected by
#' truncation or censoring. Truncation removes rows. Lower/upper  censoring sets
#'   R = NA if LCresponse/UCresponse = TRUE and
#'   rt = -Inf/Inf if LCdireciton/UCdireciton = TRUE (else NA).
#' All of these arguments can be: 1) NULL, in which case if the data frame has a column
#' of that name it is used, or if not the default is used (see argument definition),
#' 2) a scalar/logical (same value for every data row), 3) a subject-named  vector
#' (same value for each subject), 4) a vector with length matching data rows, or
#' 5) a function taking data as it argument that creates a column of appropriate
#' values. Note that if this function returns a subject named vector it will be
#' expanded as in case (3). For rows with response (column R) "nogo" the following
#' are enforced: UC=LC=LT=0, UT=Inf, UCdirection=UCresponse=TRUE
#'
#' @param data Data frame to be modified
#' @param LT Lower truncation bound below which data are removed, default 0.
#' @param UT Upper truncation bound above which data are removed, default Inf.
#' @param LC Lower censoring bound, default 0.
#' @param UC Upper censoring bound, default Inf.
#' @param LCresponse Logical. Default FALSE, set responses to NA on lower-censored trials.
#' @param UCresponse Logical. Default FALSE, set responses to NA on upper-censored trials.
#' @param LCdirection Logical. Default TRUE, lower-censored RTs are coded as -Inf; if FALSE, as NA.
#' @param UCdirection Logical. Default TRUE, upper-censored RTs are coded as Inf; if FALSE, as NA.
#' @param pContaminant Probability of contamination, default 0.
#' @param no_truncate Logical, default FALSE, for TRUE don't apply truncation to row (except if GO/NOGO).
#' @param no_censor Logical, default FALSE, for TRUE don't apply censor to row (except if GO/NOGO).
#' @param verbose Logical. Default FALSE, if TRUE report effects of filtering.
#' @param rt_resolution A double, see make_emc, specified here so binning of rt and LC/UC/LT/UT is consistent.
#'        The default is 1/60 as in make_emc, but when make_missing is called by make_data the default is to
#'        do nothing unless an explicit is value passed in the missing list.
#' @return A filtered and modified data frame with added/updated LC, UC, LT and UT columns
#' @examples
#' First make some data
#'   designRDM <- design(model = RDM,
#'   factors = list(subjects = 1:2, S = c("left", "right")),Rlevels = c("left", "right"),
#'   matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
#'   formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
#'   constants = c(s = log(1)))
#' p_vector <-log(c(B=2,A=.5,t0=0.2,v=1,v_lMTRUE=2,s_lMTRUE=.8)
#' dat <- make_data(p_vector, designRDM,n_trials = 10)
#'
#' Filter data frame without LT/UC/LT/UT columns (as in most real data files)
#' data <- dat
#' mdata <- make_missing(dat,LT=.7,LC=.75,UC=1.5,UT=1.6,verbose=TRUE)
#' @export

make_missing <- function(data, LT = NULL, UT = NULL, LC = NULL, UC = NULL,
  LCresponse = NULL, UCresponse = NULL,LCdirection = NULL, UCdirection = NULL,
  no_truncate=FALSE,no_censor=FALSE,
  pContaminant=NULL,verbose=FALSE,rt_resolution=1/60,digits = 2)
{

  pContaminant <- get_missing(pContaminant, data, "pContaminant",0,"numeric")
  if (!all(pContaminant==0)) {
    contam <- rbinom(nrow(data), 1, pContaminant) == 1
    data[contam, "rt"] <- Inf
    data[contam, "R"] <- NA
    if (verbose) {
      if (!attr(pContaminant,"subjectwise")) stat <- mean(contam) else
      stat <- tapply(contam,data$subjects,mean)
      message("% contaminated")
      print(round(100*stat,digits))
    }
  }

  no_truncate <- get_missing(no_truncate, data, "no_truncate",FALSE,"logical")
  LT <- get_missing(LT, data, "LT",0,"numeric")
  data$LT[!no_truncate] <- as.numeric(LT[!no_truncate])
  UT <- get_missing(UT, data, "UT",Inf,"numeric")
  data$UT[!no_truncate] <- as.numeric(UT[!no_truncate])


  no_censor <- get_missing(no_censor, data, "no_censor",FALSE,"logical")
  LC <- get_missing(LC, data, "LC",0,"numeric")
  data$LC[!no_censor] <- as.numeric(LC[!no_censor])
  UC <- get_missing(UC, data, "UC",Inf,"numeric")
  data$UC[!no_censor] <- as.numeric(UC[!no_censor])

  if (!is.null(rt_resolution)) {
    data$rt <- floor(data$rt/rt_resolution)*rt_resolution
    data$LC <- floor(data$LC/rt_resolution)*rt_resolution
    data$UC <- floor(data$UC/rt_resolution)*rt_resolution
    data$LT <- floor(data$LT/rt_resolution)*rt_resolution
    data$UT <- floor(data$UT/rt_resolution)*rt_resolution
  }

  LCresponse <- get_missing(LCresponse, data, "LCresponse",FALSE,"logical")
  UCresponse <- get_missing(UCresponse, data, "UCresponse",FALSE,"logical")
  LCdirection <- get_missing(LCdirection, data, "LCdirection",TRUE,"logical")
  UCdirection <- get_missing(UCdirection, data, "UCdirection",TRUE,"logical")

  isgng <- data$R == "nogo"
  isgng[is.na(isgng)] <- FALSE
  if (any(isgng)) {
    data[isgng,"UC"] <- 0
    data[isgng,"LC"] <- 0
    data[isgng,"UT"] <- Inf
    data[isgng,"LT"] <- 0
    UCdirection[isgng] <- TRUE
    UCresponse[isgng] <- TRUE
  }

  if (any(data$LT>data$LC & data$LC!=0)) stop("LT > LC not allowed")
  if (any(data$UC>data$UT & data$UC!=Inf)) stop("UC > UT not allowed")

  # Only keep trials in LT-UT (inclusive) or infinite or NA
  cutL <- is.finite(data$rt) & (data$rt < data$LT)
  cutL[is.na(cutL)] <- FALSE; cutL[no_truncate] <- FALSE
  cutU <- is.finite(data$rt) & (data$rt > data$UT)
  cutU[is.na(cutU)] <- FALSE; cutU[no_truncate] <- FALSE
  if (verbose) {
    if (!all(data$LT==0)) {
      if (!attr(LT,"subjectwise")) stat <- mean(cutL) else
        stat <- tapply(cutL,data$subjects,mean)
      message("% lower truncation")
      print(round(100*stat,digits))
    }
    if (!all(data$UT==Inf)) {
      if (!attr(UT,"subjectwise")) stat <- mean(cutU) else
        stat <- tapply(cutU,data$subjects,mean)
      message("% upper truncation")
      print(round(100*stat,digits))
    }
  }

  # Truncate
  data <- data[!cutL & !cutU, ]
  LCresponse <- LCresponse[!cutL & !cutU]
  UCresponse <- UCresponse[!cutL & !cutU]
  LCdirection <- LCdirection[!cutL & !cutU]
  UCdirection <- UCdirection[!cutL & !cutU]
  no_censor <- no_censor[!cutL & !cutU]

  # Censoring proportions (like truncation dont censor if equal to LC or UC)
  cutL <- is.finite(data$rt) & (data$rt < data$LC)
  cutL[is.na(cutL)] <- TRUE; cutL[no_censor] <- FALSE
  cutU <- is.finite(data$rt) & (data$rt > data$UC)
  cutU[is.na(cutU)] <- TRUE; cutU[no_censor] <- FALSE
  if (verbose) {
    if (!all(data$LC==0)) {
      if (!attr(LT,"subjectwise")) stat <- mean(cutL) else
        stat <- tapply(cutL,data$subjects,mean)
      message("% lower censoring (after truncation)")
      print(round(100*stat,digits))
    }
    if (!all(data$UC==Inf)) {
      if (!attr(UT,"subjectwise")) stat <- mean(cutU) else
        stat <- tapply(cutU,data$subjects,mean)
      message("% upper censoring (after truncation)")
      print(round(100*stat,digits))
    }
  }

  # Lower censoring
  data$rt[cutL &  LCdirection] <- -Inf
  data$rt[cutL & !LCdirection] <- NA
  data$R[cutL & !LCresponse] <- NA

  # Upper censoring
  data$rt[cutU &  LCdirection] <- Inf
  data$rt[cutU & !LCdirection] <- NA
  data$R[cutU & !LCresponse] <- NA

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
#' @param functions List of functions to create factors used in data generation.
#' @param return_functions Logical, should factors created by functions be returned, default FALSE.
#' @param missing List named with the arguments of the make_missing function (except data).
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
                      functions = NULL,return_functions = FALSE,
                      TC = NULL,...)

{
  if (!is.null(staircase)){
    staircase <- check_staircase(staircase)
  }

  TC <- add_defaults(TC,no_truncate=FALSE,no_censor=FALSE,verbose=FALSE,digits=2)

  # check_bounds <- FALSE

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

  if (is.null(TC$pContaminant) & any(dimnames(pars)[[2]]=="pContaminant"))
    TC$pContaminant <- pars[,"pContaminant"][data$lR==levels(data$lR)[1]]

  dropNames <- c("lR","lM")
  if (!return_functions && !is.null(design$Ffunctions))
    dropNames <- c(dropNames,names(design$Ffunctions))
  if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
  data <- data[,!(names(data) %in% dropNames)]
  for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[, i]

  data <- make_missing(data,LT=TC$LT,LC=TC$LC,UC=TC$UC,UT=TC$UT,
    LCresponse = TC$LCresponse, UCresponse = TC$UCresponse,
    LCdirection = TC$LCdirection, UCdirection = TC$UCdirection,
    pContaminant=TC$pContaminant,
    no_truncate=TC$no_truncate,no_censor=TC$no_censor,
    verbose=TC$verbose,rt_resolution=TC$rt_resolution,digits=TC$digits)

  attr(data,"p_vector") <- parameters;
  if(!is.null(post_functions)){
    for(i in 1:length(post_functions)){
      data[[names(post_functions)[i]]] <- post_functions[[i]](data)
    }
  }
  if(return_trialwise_parameters) attr(data, 'trialwise_parameters') <- trialwise_parameters
  data
}


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
