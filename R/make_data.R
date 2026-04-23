
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
#' Truncation removes rows. Contamination makes R = NA and the rt = Inf, and is
#' applied after truncation but before censoring. Lower/upper censoring sets
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

  no_truncate <- get_missing(no_truncate, data, "no_truncate",FALSE,"logical")
  LT <- get_missing(LT, data, "LT",0,"numeric")
  data$LT[!no_truncate] <- as.numeric(LT[!no_truncate])
  UT <- get_missing(UT, data, "UT",Inf,"numeric")
  data$UT[!no_truncate] <- as.numeric(UT[!no_truncate])

  # Go/no-go data are defined by explicit non-responses. Truncation is not a
  # coherent missing-data mechanism in this setting, so disable LT/UT globally.
  is_gng_data <- FALSE
  if ("R" %in% names(data)) {
    if (is.factor(data$R)) is_gng_data <- "nogo" %in% levels(data$R)
    if (!is_gng_data) is_gng_data <- any(data$R == "nogo", na.rm = TRUE)
  }
  if (is_gng_data) {
    requested_trunc <- any(((data$LT != 0) | is.finite(data$UT)) & !no_truncate, na.rm = TRUE)
    if (requested_trunc) {
      warning("Ignoring LT/UT truncation for go/no-go data (R includes 'nogo').")
    }
    no_truncate[] <- TRUE
    data$LT[] <- 0
    data$UT[] <- Inf
  }


  no_censor <- get_missing(no_censor, data, "no_censor",FALSE,"logical")
  LC <- get_missing(LC, data, "LC",0,"numeric")
  data$LC[!no_censor] <- as.numeric(LC[!no_censor])
  UC <- get_missing(UC, data, "UC",Inf,"numeric")
  data$UC[!no_censor] <- as.numeric(UC[!no_censor])

  if (!is.null(rt_resolution)) {
    data$rt <- .floor_to_rt_resolution(data$rt, rt_resolution)
    data$LC <- .floor_to_rt_resolution(data$LC, rt_resolution)
    data$UC <- .floor_to_rt_resolution(data$UC, rt_resolution)
    data$LT <- .floor_to_rt_resolution(data$LT, rt_resolution)
    data$UT <- .floor_to_rt_resolution(data$UT, rt_resolution)
  }

  LCresponse <- get_missing(LCresponse, data, "LCresponse",FALSE,"logical")
  UCresponse <- get_missing(UCresponse, data, "UCresponse",FALSE,"logical")
  LCdirection <- get_missing(LCdirection, data, "LCdirection",TRUE,"logical")
  UCdirection <- get_missing(UCdirection, data, "UCdirection",TRUE,"logical")

  LT_eff <- data$LT
  LC_eff <- data$LC
  UT_eff <- data$UT
  UC_eff <- data$UC

  isgng <- data$R == "nogo"
  isgng[is.na(isgng)] <- FALSE
  if (any(isgng)) {
    UC_eff[isgng] <- 0
    LC_eff[isgng] <- 0
    UT_eff[isgng] <- Inf
    LT_eff[isgng] <- 0
    UCdirection[isgng] <- TRUE
    UCresponse[isgng] <- FALSE
  }

  tol_l <- sqrt(.Machine$double.eps) * pmax(1, abs(LT_eff), abs(LC_eff))
  tol_u <- sqrt(.Machine$double.eps) * pmax(1, abs(UC_eff), abs(UT_eff))
  if (any((LT_eff - LC_eff) > tol_l & LC_eff != 0, na.rm = TRUE)) stop("LT > LC not allowed")
  if (any((UC_eff - UT_eff) > tol_u & UC_eff != Inf, na.rm = TRUE)) stop("UC > UT not allowed")

  # Only keep trials in LT-UT (inclusive) or infinite or NA
  cutL <- is.finite(data$rt) & (data$rt < LT_eff)
  cutL[is.na(cutL)] <- FALSE; cutL[no_truncate] <- FALSE
  cutU <- (data$rt > UT_eff)
  cutU[is.na(cutU)] <- FALSE; cutU[no_truncate] <- FALSE
  if (verbose) {
    if (!all(LT_eff==0)) {
      if (!attr(LT,"subjectwise")) stat <- mean(cutL) else
        stat <- tapply(cutL,data$subjects,mean)
      message("% lower truncation")
      print(round(100*stat,digits))
    }
    if (!all(UT_eff==Inf)) {
      if (!attr(UT,"subjectwise")) stat <- mean(cutU) else
        stat <- tapply(cutU,data$subjects,mean)
      message("% upper truncation")
      print(round(100*stat,digits))
    }
  }
  
  # Truncate
  data <- data[!cutL & !cutU, ]
  LT_eff <- LT_eff[!cutL & !cutU]
  LC_eff <- LC_eff[!cutL & !cutU]
  UT_eff <- UT_eff[!cutL & !cutU]
  UC_eff <- UC_eff[!cutL & !cutU]
  LCresponse <- LCresponse[!cutL & !cutU]
  UCresponse <- UCresponse[!cutL & !cutU]
  LCdirection <- LCdirection[!cutL & !cutU]
  UCdirection <- UCdirection[!cutL & !cutU]
  no_censor <- no_censor[!cutL & !cutU]

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

  # Censoring proportions (like truncation dont censor if equal to LC or UC)
  cutL <- is.finite(data$rt) & (data$rt < LC_eff)
  cutL[is.na(cutL)] <- TRUE; cutL[no_censor] <- FALSE
  cutU <- (data$rt > UC_eff)
  cutU[is.na(cutU)] <- TRUE; cutU[no_censor] <- FALSE
  if (verbose) {
    if (!all(LC_eff==0)) {
      if (!attr(LT,"subjectwise")) stat <- mean(cutL) else
        stat <- tapply(cutL,data$subjects,mean)
      message("% lower censoring (after truncation)")
      print(round(100*stat,digits))
    }
    if (!all(UC_eff==Inf)) {
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
  data$rt[cutU &  UCdirection] <- Inf
  data$rt[cutU & !UCdirection] <- NA
  data$R[cutU & !UCresponse] <- NA

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
  # This handles censoring and truncation where TC is not specified -- first check data, then design as a fallback (need to agree on the accepted order)
  if (is.null(TC)) {
    TC <- list()
    TC <- add_defaults(TC,
      no_truncate=FALSE,no_censor=FALSE,verbose=FALSE,digits=2,
      LCresponse=FALSE,UCresponse=FALSE,LCdirection=TRUE,UCdirection=TRUE,
      pContaminant=NULL,rt_resolution=NULL
    )
    if (is.null(data)) {
      missing_cols <- c("LT","LC","UC","UT")
    } else {
      missing_cols <- c("LT","LC","UC","UT")[!(c("LT","LC","UC","UT") %in% names(data))]
    }
    if (length(missing_cols) > 0) {
      for (nm in missing_cols) {
        if (!is.null(design[[nm]])) TC[[nm]] <- design[[nm]]
      }
    }
  } else {
    if (!is.list(TC)) stop("TC must be a list")
    TC <- add_defaults(TC,
      no_truncate=FALSE,no_censor=FALSE,verbose=FALSE,digits=2,
      LCresponse=FALSE,UCresponse=FALSE,LCdirection=TRUE,UCdirection=TRUE,
      pContaminant=NULL,rt_resolution=NULL
    )
  }

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
  ssd_meta <- NULL
  if(!is.null(functions)){
    for(i in 1:length(functions)){
      fun <- functions[[i]]
      value <- fun(data)
      meta <- attr(value, "emc_ssd")
      if (!is.null(meta)) {
        attr(value, "emc_ssd") <- NULL
        ssd_meta <- meta$staircase
      }
      data[[names(functions)[i]]] <- value
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
  use_vectorised <- isTRUE(dots_local$use_vectorised)
  if ("kernel_output_codes" %in% names(dots_local)) {
    kernel_output_codes <- dots_local$kernel_output_codes
  } else {
    kernel_output_codes <- c(1L)
  }
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
        if (use_vectorised) {
          res_i <- make_data_unconditional_vectorised(
            data = data, pars = pars, design = design, model = model,
            return_trialwise_parameters = (rep_i == expand && return_trialwise_parameters),
            kernel_output_codes = kernel_output_codes
          )
        } else {
          res_i <- make_data_unconditional(
            data = data, pars = pars, design = design, model = model,
            return_trialwise_parameters = (rep_i == expand && return_trialwise_parameters),
            kernel_output_codes = kernel_output_codes
          )
        }
        data_list[[rep_i]] <- cbind(rep = rep_i, res_i$data)
        if (rep_i == expand) twp_last <- res_i$trialwise_parameters
      }
      data <- do.call(rbind, data_list)
      rownames(data) <- NULL
      trialwise_parameters <- twp_last
    } else {
      if (use_vectorised) {
        res <- make_data_unconditional_vectorised(
          data = data, pars = pars, design = design, model = model,
          return_trialwise_parameters = return_trialwise_parameters,
          kernel_output_codes = kernel_output_codes
        )
      } else {
        res <- make_data_unconditional(
          data = data, pars = pars, design = design, model = model,
          return_trialwise_parameters = return_trialwise_parameters,
          kernel_output_codes = kernel_output_codes
        )
      }
      data <- res$data
      trialwise_parameters <- res$trialwise_parameters
    }
    attr(data, "p_vector") <- parameters
    if (return_trialwise_parameters) attr(data, "trialwise_parameters") <- trialwise_parameters
    return(data)
  } else {
    data <- design_model(
      add_accumulators(
        data,
        design$matchfun,
        simulate = TRUE,
        type = model()$type,
        Fcovariates = design$Fcovariates,
        fixed_accumulator_roles = design$fixed_accumulator_roles
      ),
      design, model,
      add_acc = FALSE, compress = FALSE, verbose = FALSE,
      rt_check = FALSE
    )
    pars <- get_pars_oo(parameters, data, model())
    if (return_trialwise_parameters)
      trialwise_parameters <- pars

    pars <- model()$Ttransform(pars, data)
    if (!is.null(optionals$nobound)) attr(pars,"ok") <- rep(TRUE,nrow(pars)) else
      pars <- fix_bound(pars, model()$bound, data$lR,fix=!is.null(optionals$shrink2bound))
    pars_ok <- attr(pars, 'ok')
    if(mean(!pars_ok) > .1){
      warning("More than 10% of parameter values fall out of model bounds, see <model_name>$bounds()")
      return(FALSE)
    }
  }
  if (expand>1) {
      data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
                    data.frame(lapply(data,rep,times=expand)))
      pars <- apply(pars,2,rep,times=expand)
  }
  lR_levels <- if (is.null(data$lR)) character(0) else levels(data$lR)
  if (!is.null(ssd_meta)) {
    ssd_meta$labels <- lR_levels
    if (!is.null(ssd_meta$specs)) {
      for (nm in names(ssd_meta$specs)) {
        if (is.list(ssd_meta$specs[[nm]])) {
          ssd_meta$specs[[nm]]$labels <- lR_levels
        }
      }
    }
    attr(data, "staircase") <- ssd_meta
    attr(pars, "staircase") <- ssd_meta
  }
  c_name <- model()$c_name
  if (any(names(data)=="RACE")) {
      Rrt <- RACE_rfun(data, pars, model)
  } else if (!is.null(c_name) && grepl("RedundantTarget", c_name, fixed = TRUE)) {
    Rrt <- RedundantTarget_rfun(data, pars, model)
  } else if (any(names(data)=="LogicalRule") && !is.null(c_name) && grepl("LogicalRules", c_name)) {
    Rrt <- LogicalRules_rfun(data, pars, model)
  } else Rrt <- model()$rfun(data,pars)
  
  if (is.null(TC$pContaminant) & any(dimnames(pars)[[2]]=="pContaminant"))
    TC$pContaminant <- pars[,"pContaminant"][data$lR==levels(data$lR)[1]]

  dropNames <- c("lR","lM","winner")
  if (!return_functions && !is.null(design$Ffunctions))
    dropNames <- c(dropNames,names(design$Ffunctions))
  # ZH added to protect SSD simulation
  if ("SSD"%in%dropNames) {dropNames=dropNames[!grepl("SSD",dropNames)]}
  if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
  data <- data[,!(names(data) %in% dropNames)]
  if (!is.null(ssd_meta)) {
    attr(data, "staircase") <- ssd_meta
    attr(pars, "staircase") <- ssd_meta
  }
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

apply_logical_rules <- function(LogicalRule, A_t, nA_t, B_t, nB_t) {
  A_yes <- A_t < nA_t
  B_yes <- B_t < nB_t

  tieA <- A_t == nA_t
  tieB <- B_t == nB_t
  if (any(tieA)) A_yes[tieA] <- runif(sum(tieA)) < 0.5
  if (any(tieB)) B_yes[tieB] <- runif(sum(tieB)) < 0.5

  tA <- pmin(A_t, nA_t)
  tB <- pmin(B_t, nB_t)

  R <- rep(NA_character_, length(A_t))
  RT <- rep(NA_real_, length(A_t))

  is_or <- LogicalRule == "OR"
  is_and <- LogicalRule == "AND"
  is_xor <- LogicalRule == "XOR"
  is_id <- LogicalRule == "ID"

  if (any(is_or)) {
    R[is_or] <- ifelse(A_yes[is_or] | B_yes[is_or], "yes", "no")
    RT_yes_or <- pmin(ifelse(A_yes, A_t, Inf), ifelse(B_yes, B_t, Inf))
    RT_no_or <- pmax(nA_t, nB_t)
    RT[is_or] <- ifelse(R[is_or] == "yes", RT_yes_or[is_or], RT_no_or[is_or])
  }

  if (any(is_and)) {
    R[is_and] <- ifelse(A_yes[is_and] & B_yes[is_and], "yes", "no")
    RT_yes_and <- pmax(A_t, B_t)
    RT_no_and <- pmin(ifelse(!A_yes, nA_t, Inf), ifelse(!B_yes, nB_t, Inf))
    RT[is_and] <- ifelse(R[is_and] == "yes", RT_yes_and[is_and], RT_no_and[is_and])
  }

  if (any(is_xor)) {
    R[is_xor] <- ifelse(xor(A_yes[is_xor], B_yes[is_xor]), "yes", "no")
    RT[is_xor] <- pmax(tA[is_xor], tB[is_xor])
  }

  if (any(is_id)) {
    R[is_id] <- ifelse(!A_yes[is_id] & !B_yes[is_id], "NN",
      ifelse(A_yes[is_id] & !B_yes[is_id], "AN",
        ifelse(!A_yes[is_id] & B_yes[is_id], "NB", "AB")
      )
    )
    RT[is_id] <- pmax(tA[is_id], tB[is_id])
  }

  data.frame(R = R, rt = RT)
}

LogicalRules_rfun <- function(data, pars, model) {
  n_trials <- dim(data)[1] / length(levels(data$lR))
  races <- levels(data$lR)
  required_racers <- c("A", "n_A", "B", "n_B")
  if (!all(required_racers %in% races)) stop("RACE must have: A, n_A, B, n_B")

  logical_rule <- as.character(data$LogicalRule[data$lR == races[1]])
  allowed_rules <- c("OR", "AND", "XOR", "ID")
  if (!all(logical_rule %in% allowed_rules)) {
    stop("LogicalRule must be one of: OR, AND, XOR, ID")
  }

  # Simulate finishing times for each accumulator role.
  Rrti <- matrix(NA_real_, nrow = n_trials, ncol = length(races), dimnames = list(NULL, races))
  for (i in races) {
    pick <- data$lR == i
    data_in <- data[pick, , drop = FALSE]
    data_in$lR <- factor(data$lR[pick])
    p <- pars[pick, , drop = FALSE]
    attr(p, "ok") <- rep(TRUE, nrow(p))
    if (!is.null(attr(pars, "staircase"))) attr(p, "staircase") <- attr(pars, "staircase")
    Rrti[, i] <- model()$rfun(data_in, p)$rt
  }

  out <- apply_logical_rules(
    LogicalRule = logical_rule,
    A_t = Rrti[, "A"],
    nA_t = Rrti[, "n_A"],
    B_t = Rrti[, "B"],
    nB_t = Rrti[, "n_B"]
  )

  has_id <- any(logical_rule == "ID")
  has_binary <- any(logical_rule %in% c("OR", "AND", "XOR"))
  if (has_id && has_binary) {
    out$R <- factor(out$R, levels = c("no", "yes", "NN", "AN", "NB", "AB"))
  } else if (has_id) {
    out$R <- factor(out$R, levels = c("NN", "AN", "NB", "AB"))
  } else {
    out$R <- factor(out$R, levels = c("no", "yes"))
  }

  out
}

RedundantTarget_rfun <- function(data, pars, model) {
  races <- levels(data$lR)
  if (!all(c("A", "B") %in% races)) stop("RedundantTarget model requires accumulator roles A and B.")
  has_nogo <- "nogo" %in% races
  n_trials <- nrow(data) / length(races)
  if (n_trials <= 0) return(data.frame(R = factor(character(0), levels = levels(data$R)), rt = numeric(0)))

  stim_col <- if ("S" %in% names(data)) {
    "S"
  } else if ("stimulus" %in% names(data)) {
    "stimulus"
  } else if ("condition" %in% names(data)) {
    "condition"
  } else {
    stop("RedundantTarget model requires stimulus column `S` (or `stimulus`/`condition`).")
  }
  cond <- as.character(data[data$lR == races[1], stim_col])
  cond[cond %in% c("BA", "A+B", "B+A")] <- "AB"
  if (!all(cond %in% c("A", "B", "AB"))) {
    stop("RedundantTarget stimulus must be one of A, B, AB.")
  }

  Rrti <- matrix(Inf, nrow = n_trials, ncol = length(races), dimnames = list(NULL, races))
  for (i in races) {
    pick <- data$lR == i
    data_in <- data[pick, , drop = FALSE]
    data_in$lR <- factor(data$lR[pick], levels = i)
    p <- pars[pick, , drop = FALSE]
    attr(p, "ok") <- rep(TRUE, nrow(p))
    if (!is.null(attr(pars, "staircase"))) attr(p, "staircase") <- attr(pars, "staircase")
    Rrti[, i] <- model()$rfun(data_in, p)$rt
  }

  go_levels <- levels(data$R)
  go_levels <- go_levels[go_levels != "nogo"]
  go_resp <- if ("yes" %in% go_levels) "yes" else if (length(go_levels) > 0) go_levels[1] else "yes"

  rt <- rep(Inf, n_trials)
  R <- rep(NA_character_, n_trials)

  tA <- Rrti[, "A"]
  tB <- Rrti[, "B"]
  tN <- if (has_nogo) Rrti[, "nogo"] else rep(Inf, n_trials)

  go_t <- rep(Inf, n_trials)
  go_t[cond == "A"] <- tA[cond == "A"]
  go_t[cond == "B"] <- tB[cond == "B"]
  go_t[cond == "AB"] <- pmin(tA[cond == "AB"], tB[cond == "AB"])

  emitted <- is.finite(go_t) & (go_t < tN)
  rt[emitted] <- go_t[emitted]
  R[emitted] <- go_resp

  data.frame(R = factor(R, levels = levels(data$R)), rt = rt)
}

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
  if(is.null(covariances)) {
    model <- design$model()
    p_template <- matrix(0, nrow = 1, ncol = length(group_means))
    colnames(p_template) <- get_p_types(names(group_means))
    # Use a baseline in sampled space for consistency
    vars <- (abs(group_means) + 1) * variance_proportion
    vars[vars == 0] <- variance_proportion
    covariances <- diag(vars)
  }
  random_effects <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=covariances)
  colnames(random_effects) <- names(sampled_pars(design))
  rownames(random_effects) <- subnames
  return(random_effects)
}
