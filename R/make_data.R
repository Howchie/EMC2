
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
    # A one-per-trial vector must never be read as a subject lookup.  Trial-level
    # vectors carry data row names ("1","2",...), so with numerically labelled
    # subjects every subject level is also a row name and hasName() is TRUE for
    # a vector that is not subjectwise at all: each subject would then be given
    # the value of the row whose name matches its label (all subject 1's trials).
    # A genuine subject lookup is named by exactly the subject levels, which
    # also settles the case where the two lengths happen to coincide.
    subj_levels <- levels(data$subjects)
    subjectwise <- all(hasName(supplied,subj_levels)) &&
      (length(supplied) != nrow(data) || setequal(names(supplied),subj_levels))
    if (!subjectwise) bound <- supplied else
      bound <- supplied[as.character(data$subjects)]
  }
  out <- rep(bound,length.out=nrow(data))
  attr(out,"subjectwise") <- subjectwise
  out
}

.active_nogo_rows <- function(data) {
  if (!("R" %in% names(data)) || !is.factor(data$R)) {
    return(rep(FALSE, nrow(data)))
  }
  nogo_code <- match("nogo", levels(data$R))
  if (is.na(nogo_code)) {
    return(rep(FALSE, nrow(data)))
  }
  active_nogo <- rep(TRUE, nrow(data))
  if ("RACE" %in% names(data) && is.factor(data$RACE)) {
    nacc_by_level <- suppressWarnings(as.integer(levels(data$RACE)))
    if (!anyNA(nacc_by_level)) {
      race_codes <- as.integer(data$RACE)
      active_nogo[] <- FALSE
      ok <- !is.na(race_codes)
      active_nogo[ok] <- nacc_by_level[race_codes[ok]] >= nogo_code
    }
  }
  active_nogo
}

#' Add information about missing values to data and modify/filter accordingly.
#'
#' @details
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
#' expanded as in case (3). For rows where a nogo response is active, the following
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
#' @param pContaminant Probability of an *omission* contaminant (`rt = Inf`,
#' `R = NA`), default 0.
#' @param pGuess Probability of a uniform *guess* contaminant, drawn among the
#' trials the omission did not take, so `P(guess) = (1 - pContaminant) * pGuess`.
#' A guess replaces `rt` with a draw from `guess_window` and `R` with a uniform
#' draw over the overt response levels.  Default 0.
#' @param guess_window Length-2 numeric `c(lower, upper)` for the guess RT.
#' Defaults to the effective truncation/censoring window, widened to
#' `max(5, floor(max(rt)) + 1)` when it has no finite upper edge -- the same rule
#' the likelihood uses (see `resolve_guess_window`).
#' @param no_truncate Logical, default FALSE, for TRUE don't apply truncation to row (except if GO/NOGO).
#' @param no_censor Logical, default FALSE, for TRUE don't apply censor to row (except if GO/NOGO).
#' @param verbose Logical. Default FALSE, if TRUE report effects of filtering.
#' @param rt_resolution A double, see make_emc, specified here so binning of rt and LC/UC/LT/UT is consistent.
#'        The default is 1/60 as in make_emc, but when make_missing is called by make_data the default is to
#'        do nothing unless an explicit is value passed in the missing list.
#' @param digits Integer. Number of decimal places for verbose output, default 2.
#' @return A filtered and modified data frame with added/updated LC, UC, LT and UT columns
#' @examples
#' \dontrun{
#' # First make some data
#'   designRDM <- design(model = RDM,
#'   factors = list(subjects = 1:2, S = c("left", "right")),Rlevels = c("left", "right"),
#'   matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
#'   formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
#'   constants = c(s = log(1)))
#' p_vector <- log(c(B=2,A=.5,t0=0.2,v=1,v_lMTRUE=2,s_lMTRUE=.8))
#' dat <- make_data(p_vector, designRDM,n_trials = 10)
#'
#' # Filter data frame without LT/UC/LT/UT columns (as in most real data files)
#' data <- dat
#' mdata <- make_missing(dat,LT=.7,LC=.75,UC=1.5,UT=1.6,verbose=TRUE)
#' }
#' @export

make_missing <- function(data, LT = NULL, UT = NULL, LC = NULL, UC = NULL,
  LCresponse = NULL, UCresponse = NULL,LCdirection = NULL, UCdirection = NULL,
  no_truncate=FALSE,no_censor=FALSE,
  pContaminant=NULL,pGuess=NULL,guess_window=NULL,
  verbose=FALSE,rt_resolution=1/60,digits = 2)
{

  no_truncate <- get_missing(no_truncate, data, "no_truncate",FALSE,"logical")
  LT <- get_missing(LT, data, "LT",0,"numeric")
  data$LT[!no_truncate] <- as.numeric(LT[!no_truncate])
  UT <- get_missing(UT, data, "UT",Inf,"numeric")
  data$UT[!no_truncate] <- as.numeric(UT[!no_truncate])

  active_nogo <- .active_nogo_rows(data)
  if (any(active_nogo)) {
    requested_trunc <- any(((data$LT != 0) | is.finite(data$UT)) & !no_truncate & active_nogo, na.rm = TRUE)
    if (requested_trunc) {
      warning("Ignoring LT/UT truncation for rows with an active nogo response.")
    }
    no_truncate[active_nogo] <- TRUE
    data$LT[active_nogo] <- 0
    data$UT[active_nogo] <- Inf
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

  isgng <- active_nogo & data$R == "nogo"
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

  # Uniform guess contaminant (pGuess).  A guess is part of the generative
  # process, so it is drawn BEFORE the truncation cut below and is exposed to
  # it like any other response -- the previous placement (after the cut) made
  # guessing an act of the experimenter rather than the participant.  In
  # practice this changes nothing: resolve_guess_window() builds the window as
  # [max(LT, LC), min(UC, UT)] (R/design.R), so a guess lands inside the
  # retention window by construction and can never be truncated away.  That is
  # also what keeps the ordering consistent with the likelihood, whose guess
  # term carries no truncation normaliser (src/contaminant_mixture.h).
  #
  # Nesting with the omission -- P(guess) = (1 - pC) * pG -- is preserved by
  # order alone: pContaminant is applied after the cut and overwrites a guessed
  # rt with +Inf, so an omission still wins whenever both fire.
  pGuess <- get_missing(pGuess, data, "pGuess",0,"numeric")
  if (!all(pGuess==0)) {
    gw <- .resolve_sim_guess_window(data, guess_window)
    # A guess is an overt response: never withheld, never a timeout.
    resp_levels <- levels(data$R)
    guess_levels <- resp_levels[!(resp_levels %in% c("nogo","time"))]
    if (length(guess_levels) == 0L)
      stop("pGuess needs at least one overt response level (excluding nogo/time)")
    guess <- rbinom(nrow(data), 1, pGuess) == 1
    if (any(guess)) {
      data$rt[guess] <- runif(sum(guess), gw[1], gw[2])
      data$R[guess] <- factor(sample(guess_levels, sum(guess), replace = TRUE),
                              levels = resp_levels)
    }
    if (verbose) {
      if (!attr(pGuess,"subjectwise")) stat <- mean(guess) else
        stat <- tapply(guess,data$subjects,mean)
      message("% guesses (window [",signif(gw[1],4),", ",signif(gw[2],4),"])")
      print(round(100*stat,digits))
    }
  }

  # Only keep trials in LT-UT (inclusive) or infinite or NA
  cutL <- is.finite(data$rt) & (data$rt < LT_eff & is.finite(data$rt))
  cutL[is.na(cutL)] <- FALSE; cutL[no_truncate] <- FALSE
  cutU <- (data$rt > UT_eff & is.finite(data$rt))
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
  # Applied after the truncation cut: an omission is an rt = +Inf atom, which
  # no finite truncation bound can remove, so pC is the omission proportion
  # among retained trials either way.  It overwrites any guess drawn above,
  # which is what enforces P(guess) = (1 - pC) * pG.
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


# Guess window used when SIMULATING.  Deliberately the same rule as
# resolve_guess_window() in R/design.R, which is what the likelihood uses: the
# window is the effective truncation/censoring window, so a simulated guess can
# never be censored or truncated away, and an unbounded upper edge falls back to
# max(5, floor(max rt) + 1) rather than a hard 5 s.
.resolve_sim_guess_window <- function(data, guess_window = NULL) {
  if (!is.null(guess_window)) {
    gw <- as.numeric(guess_window)
    if (length(gw) != 2 || !all(is.finite(gw)) || !(gw[2] > gw[1]))
      stop("guess_window must be a length-2 numeric c(lower, upper) with upper > lower")
    return(gw)
  }
  get1 <- function(nm, default) {
    v <- data[[nm]]
    if (is.null(v) || length(v) == 0) default else v
  }
  LG <- min(pmax(get1("LT", 0), get1("LC", 0)))
  UG <- max(pmin(get1("UT", Inf), get1("UC", Inf)))
  if (!is.finite(UG)) {
    rt <- data$rt
    max_rt <- suppressWarnings(max(rt[is.finite(rt)]))
    UG <- if (is.finite(max_rt)) max(5, floor(max_rt) + 1) else 5
  }
  if (!is.finite(LG) || !(UG > LG))
    stop("Cannot resolve a proper guess window; supply guess_window explicitly.")
  c(LG, UG)
}


check_missing <- function(TC,data=NULL,design=NULL) {
  # This handles censoring and truncation where TC is not specified.
  # First check data, then design
  if (is.null(TC)) {
    TC <- list()
    TC <- add_defaults(TC,LT=0,LC=0,UT=Inf,UC=Inf,
      no_truncate=FALSE,no_censor=FALSE,verbose=FALSE,digits=2,
      LCresponse=FALSE,UCresponse=FALSE,LCdirection=TRUE,UCdirection=TRUE,
      pContaminant=NULL,pGuess=NULL,guess_window=NULL,rt_resolution=NULL
    )
    if (!is.null(data)) {
      for (i in c("LT","LC","UC","UT")) {
        if (!is.null(data[[i]])) TC[[i]] <- data[[i]]
      }
    } else if (!is.null(design) && !is.null(design$TC)) {
        for (i in names(TC)) TC[[i]] <- design$TC[[i]]
    }
  } else {
    if (!is.list(TC)) stop("TC must be a list")
    TC <- add_defaults(TC,LT=0,LC=0,UT=Inf,UC=Inf,
      no_truncate=FALSE,no_censor=FALSE,verbose=FALSE,digits=2,
      LCresponse=FALSE,UCresponse=FALSE,LCdirection=TRUE,UCdirection=TRUE,
      pContaminant=NULL,pGuess=NULL,guess_window=NULL,rt_resolution=NULL
    )
  }
  TC
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
#' @param TC List of truncation/censoring arguments passed to \code{make_missing}
#'   (e.g. \code{list(LT=0.1, UC=2)}). NULL means no truncation or censoring is applied.
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


  post_functions <- NULL
  precomputed_design <- NULL
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(is(parameters, "emc")){
    if(is.null(design)) design <- get_design(parameters)[[1]] # Currently not supported for multiple designs
    if(is.null(data)) data <- get_data(parameters)
    parameters <- do.call(rbind, credint(parameters, probs = 0.5, selection = "alpha", by_subject = TRUE))
  }

  # Resolve censoring and truncation settings from the data, then the design.
  # Run this after the emc block so fitted settings are available.
  # Otherwise check_missing() would silently fall back to no-censoring defaults.
  TC <- check_missing(TC,design=design,data=data)

  # Make sure parameters are in the right format, either matrix or vector.
  # predict.emc() can provide the already-mapped parameter matrix; when the
  # parameters are named, avoid rebuilding the data-aware sampled-parameter
  # mapping just to recover names.
  mapped_parameters <- optionals$mapped_parameters
  sampled_p_names <- optionals$sampled_p_names
  if (is.null(sampled_p_names)) {
    sampled_p_names <- if (is.null(dim(parameters))) names(parameters) else colnames(parameters)
  }
  if (is.null(sampled_p_names)) sampled_p_names <- names(sampled_pars(design))
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
    # Keep the stochastic function input in the same subject order that
    # design_model() uses after accumulator expansion.  In particular, grouped
    # staircase metadata is row-aligned and must not be captured before a
    # later subject sort changes that order.
    data <- data[order(data$subjects), , drop = FALSE]
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
  if (is.null(mapped_parameters)) {
    pars <- t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
    pars <- add_constants(pars,design$constants)
  } else {
    if (!is.matrix(mapped_parameters)) {
      stop("mapped_parameters must be a matrix")
    }
    pars <- mapped_parameters
  }
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
    if (is.null(precomputed_design)) {
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
    } else {
      # predict.emc() can reuse this parameter-independent scaffold across
      # posterior draws.  It is valid only for the conditional-on-data path;
      # trialwise/unconditional simulation still rebuilds its design because
      # feedback covariates change after each simulated trial.
      data <- precomputed_design
    }
    if (is.null(mapped_parameters)) pars <- get_pars_oo(parameters, data, model())
    if (return_trialwise_parameters)
      trialwise_parameters <- pars

    pars <- model()$Ttransform(pars, data)
    if (!is.null(optionals$nobound)) attr(pars,"ok") <- rep(TRUE,nrow(pars)) else
      pars <- fix_bound(pars, model()$bound, data$lR,fix=!is.null(optionals$shrink2bound))
    pars_ok <- attr(pars, 'ok')
    if(mean(!pars_ok) > .1){
      warning("More than 10% of parameter values fall out of model bounds, see <model_name>$bounds()")
      if (!isTRUE(optionals$check_bounds)) {
        return(FALSE)
      }
    }
  }
  if (expand>1) {
      data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
                    data.frame(lapply(data,rep,times=expand)))
      # apply() drops attributes, and some rfuns (e.g. the correlated BAwL C++
      # simulator) require the "ok" bound flag, so carry it across the replication.
      pars_ok_attr <- attr(pars, "ok")
      pars <- apply(pars,2,rep,times=expand)
      if (!is.null(pars_ok_attr)) attr(pars, "ok") <- rep(pars_ok_attr, times = expand)
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
    ssd_meta$UC <- TC$UC
    attr(data, "staircase") <- ssd_meta
    attr(pars, "staircase") <- ssd_meta
  }
  c_name <- model()$c_name
  if (any(names(data)=="RACE")) {
      Rrt <- RACE_rfun(data, pars, model)
  } else if (any(names(data)=="LogicalRule") && !is.null(c_name) && grepl("LogicalRules", c_name)) {
    Rrt <- LogicalRules_rfun(data, pars, model)
  } else Rrt <- model()$rfun(data,pars)

  # One value per trial: for a race the parameter matrix has one row per
  # accumulator, so take the first accumulator's row; DDM-family models have no
  # lR column and are already one row per trial.
  first_acc <- if (is.null(data$lR)) rep(TRUE, nrow(pars)) else
    data$lR == levels(data$lR)[1]
  # unname(): these are per-trial vectors, and their pars row names would make
  # get_missing() mistake them for a subject lookup (see get_missing()).
  if (is.null(TC$pContaminant) & any(dimnames(pars)[[2]]=="pContaminant"))
    TC$pContaminant <- unname(pars[,"pContaminant"][first_acc])
  if (is.null(TC$pGuess) & any(dimnames(pars)[[2]]=="pGuess"))
    TC$pGuess <- unname(pars[,"pGuess"][first_acc])

  dropNames <- c("lR","lM","winner")
  if (!return_functions && !is.null(design$Ffunctions))
    dropNames <- c(dropNames,names(design$Ffunctions))
  # Preserve SSD columns needed by staircase simulation.
  if ("SSD"%in%dropNames) {dropNames=dropNames[!grepl("SSD",dropNames)]}
  if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
  data <- data[,!(names(data) %in% dropNames)]
  if (!is.null(ssd_meta)) {
    attr(data, "staircase") <- ssd_meta
    attr(pars, "staircase") <- ssd_meta
  }
  for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[, i]



  # Same reasoning as in design_model(): a model that cannot be binned must not
  # be binned by a TC carried on the design or handed in by the caller either,
  # or predict() would return data on a lattice the fit never used.
  if (!is.null(TC$rt_resolution) && !is.null(model) && !model_compress_ok(model))
    TC$rt_resolution <- NULL

  data <- make_missing(data,LT=TC$LT,LC=TC$LC,UC=TC$UC,UT=TC$UT,
    LCresponse = TC$LCresponse, UCresponse = TC$UCresponse,
    LCdirection = TC$LCdirection, UCdirection = TC$UCdirection,
    pContaminant=TC$pContaminant,pGuess=TC$pGuess,guess_window=TC$guess_window,
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
  lR_levels <- levels(data$lR)
  n_acc <- length(lR_levels)
  n_trials <- dim(data)[1]/n_acc
  Rrt <- data.frame(
    R = rep(NA_real_, n_trials),
    rt = rep(NA_real_, n_trials)
  )
  RACE <- data[data$lR==lR_levels[1],"RACE"]
  ok <- as.numeric(data$lR) <= as.numeric(as.character(data$RACE))
  for (i in levels(RACE)) {
    pick <- data$RACE==i
    rows <- which(RACE == i)
    data_in <- data[pick & ok,]
    data_in$lR <- factor(data$lR[pick & ok])
    tmp <- pars[pick & ok,]
    attr(tmp, "ok") <- rep(T, nrow(tmp))
    if (!is.null(attr(pars, "staircase"))) attr(tmp, "staircase") <- attr(pars, "staircase")
    Rrti <- model()$rfun(data_in,tmp)
    Rrti <- .apply_timed_guess_winner(Rrti, lR_levels)
    Rrti$R <- as.numeric(Rrti$R)
    for (nm in setdiff(names(Rrti), names(Rrt))) {
      Rrt[[nm]] <- NA
    }
    for (nm in setdiff(names(Rrt), names(Rrti))) {
      Rrti[[nm]] <- NA
    }
    Rrt[rows, names(Rrt)] <- Rrti[, names(Rrt), drop = FALSE]
  }
  Rrt$R <- factor(Rrt$R, labels = lR_levels, levels = 1:n_acc)
  return(Rrt)
}

apply_logical_rules <- function(LogicalRule, A_t, nA_t, B_t, nB_t) {
  A_yes <- A_t < nA_t
  B_yes <- B_t < nB_t

  # Equal finite finishing times are genuine ties.  Two Inf values mean that
  # neither accumulator finished, not that they tied at a finite time; in
  # particular, GNG must keep those trials as withheld responses.
  tieA <- is.finite(A_t) & is.finite(nA_t) & (A_t == nA_t)
  tieB <- is.finite(B_t) & is.finite(nB_t) & (B_t == nB_t)
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
  # OR_DETECTION_GNG is the full four-horse OR task whose "no" outcome is a
  # withheld response (rt = Inf, R = NA) rather than an overt "no".
  is_gng <- LogicalRule == "OR_DETECTION_GNG"

  if (any(is_or)) {
    R[is_or] <- ifelse(A_yes[is_or] | B_yes[is_or], "yes", "no")
    RT_yes_or <- pmin(ifelse(A_yes, A_t, Inf), ifelse(B_yes, B_t, Inf))
    RT_no_or <- pmax(nA_t, nB_t)
    RT[is_or] <- ifelse(R[is_or] == "yes", RT_yes_or[is_or], RT_no_or[is_or])
  }

  if (any(is_gng)) {
    go <- A_yes | B_yes
    RT_yes_gng <- pmin(ifelse(A_yes, A_t, Inf), ifelse(B_yes, B_t, Inf))
    R[is_gng] <- ifelse(go[is_gng], "yes", NA_character_)
    RT[is_gng] <- ifelse(go[is_gng], RT_yes_gng[is_gng], Inf)
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

# Normalized stimulus condition codes (NN/AN/NB/AB) for LogicalRules trials,
# or NULL when the data carry no stimulus column.  The compiled likelihood
# uses condition code 0 for a missing stimulus, which is the no-stimulus (NN)
# condition; keep the R simulator consistent so posterior prediction also
# works for fits whose data contain an absent/NA stimulus value.
.lr_stimulus_cond <- function(data, races) {
  stim_col <- if ("S" %in% names(data)) {
    "S"
  } else if ("stimulus" %in% names(data)) {
    "stimulus"
  } else if ("condition" %in% names(data)) {
    "condition"
  } else {
    return(NULL)
  }
  cond <- as.character(data[data$lR == races[1], stim_col])
  cond[is.na(cond)] <- "NN"
  cond[cond == "none"] <- "NN"
  cond[cond %in% c("BA", "A+B", "B+A")] <- "AB"
  cond[cond == "A"] <- "AN"
  cond[cond == "B"] <- "NB"
  cond
}

# Shared-capacity LogicalRules finishing-time sampler. One latent standard-normal
# factor per AB trial adds a shared shift to both target drift means:
# V_i = v_i + kappa + tau*z + eps_i for i in {A, B}; nontarget and
# time racers never load on the factor.  With posdrift the correlated pair is
# jointly conditioned on both target drifts being positive (rejection on
# (z, V_A, V_B)); independent racers keep their ordinary univariate
# truncation, which factorises out of the joint law.  Trials that are not AB
# (or whose capacity parameters sit at kappa = 0, tau = 0) reproduce the
# ordinary independent LBA drift draws exactly.
.lr_capacity_finish_times <- function(data, pars, races, posdrift) {
  n_trials <- nrow(data) / length(races)
  cond <- .lr_stimulus_cond(data, races)
  if (is.null(cond)) {
    stop("LogicalRules capacity simulation requires stimulus column `S` (or `stimulus`/`condition`).")
  }
  rowsA <- which(data$lR == "A")
  rowsB <- which(data$lR == "B")
  kappa <- pars[rowsA, "kappa"]
  tau <- pars[rowsA, "tau"]
  if (any(pars[rowsB, "kappa"] != kappa) || any(pars[rowsB, "tau"] != tau)) {
    stop("LogicalRules capacity requires kappa and tau shared by the A and B rows within each trial.")
  }
  if (any(!is.finite(kappa)) || any(!is.finite(tau) | tau < 0)) {
    stop("LogicalRules capacity requires finite kappa and tau >= 0.")
  }

  if (.use_cpp_rfun()) {
    return(logicalrules_capacity_finish_cpp(pars, races, cond, posdrift))
  }

  pair_active <- cond == "AB" & (kappa != 0 | tau != 0)

  trial_idx <- integer(nrow(data))
  for (r in races) trial_idx[data$lR == r] <- seq_len(n_trials)
  row_pair <- (data$lR %in% c("A", "B")) & pair_active[trial_idx]

  lower <- if (posdrift) 0 else -Inf
  drifts <- rep(NA_real_, nrow(pars))
  if (any(!row_pair)) {
    drifts[!row_pair] <- msm::rtnorm(sum(!row_pair),
                                     mean = pars[!row_pair, "v"],
                                     sd = pars[!row_pair, "sv"], lower = lower)
  }

  act <- which(pair_active)
  if (length(act)) {
    rA <- rowsA[act]
    rB <- rowsB[act]
    vA <- pars[rA, "v"]; svA <- pars[rA, "sv"]
    vB <- pars[rB, "v"]; svB <- pars[rB, "sv"]
    kap <- kappa[act]; tu <- tau[act]
    VA <- rep(NA_real_, length(act))
    VB <- rep(NA_real_, length(act))
    todo <- seq_along(act)
    for (iter in seq_len(10000L)) {
      z <- rnorm(length(todo))
      shift <- kap[todo] + tu[todo] * z
      dA <- rnorm(length(todo), mean = vA[todo] + shift, sd = svA[todo])
      dB <- rnorm(length(todo), mean = vB[todo] + shift, sd = svB[todo])
      okd <- !posdrift | (dA > 0 & dB > 0)
      VA[todo[okd]] <- dA[okd]
      VB[todo[okd]] <- dB[okd]
      todo <- todo[!okd]
      if (!length(todo)) break
    }
    if (length(todo)) {
      stop("LogicalRules capacity jointly positive drift rejection exceeded ",
           10000L, " sweeps; check that the target drift means are not far below zero.")
    }
    drifts[rA] <- VA
    drifts[rB] <- VB
  }

  # LBA finishing time per row (k = 0): t = t0 + (b - A*U)/V for V > 0,
  # otherwise the racer never finishes.
  dt <- (pars[, "b"] - pars[, "A"] * runif(nrow(pars))) / drifts
  dt[!is.finite(dt) | dt < 0] <- Inf
  ft <- pars[, "t0"] + dt
  out <- matrix(NA_real_, nrow = n_trials, ncol = length(races),
                dimnames = list(NULL, races))
  for (r in races) out[, r] <- ft[data$lR == r]
  out
}

LogicalRules_rfun <- function(data, pars, model) {
  n_trials <- dim(data)[1] / length(levels(data$lR))
  races <- levels(data$lR)
  if (!all(c("A", "B") %in% races)) stop("LogicalRules requires accumulator roles A and B.")

  logical_rule <- as.character(data$LogicalRule[data$lR == races[1]])
  allowed_rules <- c("OR", "AND", "XOR", "ID", "OR_DETECTION_ANALYTIC", "OR_DETECTION_GNG")
  if (!all(logical_rule %in% allowed_rules)) {
    stop("LogicalRule must be one of: OR, AND, XOR, ID, OR_DETECTION_ANALYTIC, OR_DETECTION_GNG")
  }

  # Simulate finishing times for each accumulator role.  The capacity variant
  # draws the target drifts jointly (shared factor on AB trials) and computes
  # the LBA finishing times directly; ordinary models delegate to the model's
  # per-role rfun.
  has_capacity <- all(c("kappa", "tau") %in% colnames(pars))
  if (has_capacity) {
    posdrift <- !grepl("IO", model()$c_name)
    Rrti <- .lr_capacity_finish_times(data, pars, races, posdrift)
  } else {
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
  }

  # Only the analytic detector activates detectors by stimulus condition.
  # OR_DETECTION_GNG is the full four-horse OR task (a withheld response in
  # place of the overt "no") and is simulated through the legacy subrace path.
  is_detection <- logical_rule == "OR_DETECTION_ANALYTIC"
  if (any(is_detection)) {
    cond <- .lr_stimulus_cond(data, races)
    if (is.null(cond)) {
      stop("LogicalRules detection rules require stimulus column `S` (or `stimulus`/`condition`).")
    }
    if (!all(cond[is_detection] %in% c("NN", "AN", "NB", "AB"))) {
      stop("LogicalRules detection rules require stimulus NN/AN/NB/AB (or A/B/AB).")
    }
  }

  out <- data.frame(R = rep(NA_character_, n_trials), rt = rep(Inf, n_trials))

  legacy <- !is_detection
  if (any(legacy)) {
    required_racers <- c("A", "n_A", "B", "n_B")
    if (!all(required_racers %in% races)) stop("LogicalRules OR/AND/XOR/ID requires A, n_A, B, n_B.")
    out_legacy <- apply_logical_rules(
      LogicalRule = logical_rule[legacy],
      A_t = Rrti[legacy, "A"],
      nA_t = Rrti[legacy, "n_A"],
      B_t = Rrti[legacy, "B"],
      nB_t = Rrti[legacy, "n_B"]
    )
    if ("time" %in% races) {
      tT <- Rrti[legacy, "time"]
      emit_T <- is.finite(tT) & (tT < out_legacy$rt)
      if (any(emit_T)) {
        rule_T <- logical_rule[legacy][emit_T]
        is_id_T <- rule_T == "ID"
        is_gng_T <- rule_T == "OR_DETECTION_GNG"
        guess_R <- character(sum(emit_T))
        if (any(is_id_T))
          guess_R[is_id_T] <- sample(c("NN", "AN", "NB", "AB"), sum(is_id_T), replace = TRUE)
        # A GNG guess is always an overt (go) "yes"; the withheld outcome cannot
        # be produced by the guess racer.
        if (any(is_gng_T)) guess_R[is_gng_T] <- "yes"
        other_T <- !is_id_T & !is_gng_T
        if (any(other_T))
          guess_R[other_T] <- sample(c("yes", "no"), sum(other_T), replace = TRUE)
        out_legacy$rt[emit_T] <- tT[emit_T]
        out_legacy$R[emit_T] <- guess_R
      }
    }
    out[legacy, ] <- out_legacy
  }

  if (any(is_detection)) {
    tA <- Rrti[is_detection, "A"]
    tB <- Rrti[is_detection, "B"]
    cond_d <- cond[is_detection]
    tT <- if ("time" %in% races) Rrti[is_detection, "time"] else rep(Inf, sum(is_detection))

    go_t <- rep(Inf, sum(is_detection))
    go_t[cond_d == "AN"] <- tA[cond_d == "AN"]
    go_t[cond_d == "NB"] <- tB[cond_d == "NB"]
    go_t[cond_d == "AB"] <- pmin(tA[cond_d == "AB"], tB[cond_d == "AB"])

    emit <- is.finite(go_t) & (go_t < tT)
    det_idx <- which(is_detection)
    out$rt[det_idx[emit]] <- go_t[emit]
    out$R[det_idx[emit]] <- "yes"

    if ("time" %in% races) {
      emit_T <- is.finite(tT) & (tT < go_t)
      out$rt[det_idx[emit_T]] <- tT[emit_T]
      out$R[det_idx[emit_T]] <- "yes"
    }

    if (any(!emit)) {
      out$rt[det_idx[!emit]] <- Inf
      out$R[det_idx[!emit]] <- NA_character_
    }
  }

  has_id <- any(logical_rule == "ID")
  has_binary <- any(logical_rule %in% c("OR", "AND", "XOR", "OR_DETECTION_ANALYTIC", "OR_DETECTION_GNG"))
  if (has_id && has_binary) {
    out$R <- factor(out$R, levels = c("no", "yes", "NN", "AN", "NB", "AB"))
  } else if (has_id) {
    out$R <- factor(out$R, levels = c("NN", "AN", "NB", "AB"))
  } else {
    out$R <- factor(out$R, levels = c("no", "yes"))
  }

  out
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
#' @param variance_proportion A non-negative double. Optional. If ``covariances``
#'   are not specified, this controls the default relative spread of the
#'   subject-level effects. Identity-transformed parameters use a standard
#'   deviation of ``variance_proportion * max(abs(mean), 1)``; exponential
#'   parameters use a log-normal spread with this coefficient of variation;
#'   probit parameters use this value as the latent-scale standard deviation.
#'   The covariances are 0.
#' @param covariances A covariance matrix. Optional. Specify the intended covariance matrix.
#' @param max_tries An integer. Maximum number of redraw rounds used to replace
#'   subject-level parameter draws that fall outside the model bounds.
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

make_random_effects <- function(design, group_means, n_subj = NULL,
                                variance_proportion = .2, covariances = NULL,
                                max_tries = 1000L){
  if(is.null(n_subj)){
    n_subj <- length(design$Ffactors$subjects)
    if (n_subj < 1L)
      stop("Could not determine number of subjects from design; specify n_subj explicitly")
    subnames <- as.character(design$Ffactors$subjects)
  } else{
    n_subj <- as.integer(n_subj)[1L]
    if (is.na(n_subj) || n_subj < 1L)
      stop("n_subj must be a positive integer")
    subnames <- as.character(1:n_subj)
  }
  sampled_names <- names(sampled_pars(design))
  if(length(group_means) != length(sampled_names))
    stop("You must specify as many means as parameters in your design")
  if (is.null(names(group_means))) {
    names(group_means) <- sampled_names
  } else {
    if (!all(sampled_names %in% names(group_means)))
      stop("names of group_means do not match sampled_pars(design)")
    group_means <- group_means[sampled_names]
  }
  if (any(!is.finite(group_means)))
    stop("group_means must be finite")
  max_tries <- as.integer(max_tries)[1L]
  if (is.na(max_tries) || max_tries < 1L)
    stop("max_tries must be a positive integer")

  if(is.null(covariances)) {
    if (length(variance_proportion) != 1L || !is.finite(variance_proportion) ||
        variance_proportion < 0)
      stop("variance_proportion must be a single non-negative finite number")
    model <- design$model()
    base_names <- get_p_types(sampled_names)
    if (is.null(model$transform) || is.null(model$transform$func)) {
      transform <- rep("identity", length(sampled_names))
    } else {
      transform <- model$transform$func[base_names]
      if (length(transform) != length(sampled_names) || anyNA(transform))
        stop("Could not determine parameter transforms for the design")
    }

    # Draws are made on EMC2's sampled scale.  Use a relative natural-scale
    # spread rather than treating a log parameter's numeric value as a scale.
    # For exp transforms, log1p(cv^2) gives a log-normal coefficient of
    # variation equal to variance_proportion exactly.
    vars <- vapply(seq_along(group_means), function(i) {
      if (transform[[i]] == "exp") {
        log1p(variance_proportion^2)
      } else if (transform[[i]] == "pnorm") {
        variance_proportion^2
      } else {
        (variance_proportion * max(abs(group_means[[i]]), 1))^2
      }
    }, numeric(1))
    covariances <- diag(vars)
  } else {
    covariances <- as.matrix(covariances)
    if (!all(dim(covariances) == length(group_means)))
      stop("covariances must be a square matrix with one row/column per parameter")
    if (!is.null(rownames(covariances)) && all(sampled_names %in% rownames(covariances)) &&
        !is.null(colnames(covariances)) && all(sampled_names %in% colnames(covariances))) {
      covariances <- covariances[sampled_names, sampled_names, drop = FALSE]
    }
  }

  # Build the same parameter map used by make_data(), so the guard checks
  # natural, cell-level parameters rather than individual regression
  # coefficients (which may legitimately be negative).
  bound_checker <- NULL
  model <- design$model()
  if (!is.null(model$bound) && !is.null(design$Ffactors)) {
    guard_design <- design
    guard_design$Ffactors$subjects <- subnames
    guard_dadm <- tryCatch({
      guard_data <- minimal_design(
        guard_design, drop_subjects = FALSE, n_trials = 1,
        add_acc = FALSE, drop_R = FALSE, drop_R_levels = FALSE,
        do_functions = FALSE, verbose = FALSE
      )
      if (!is.data.frame(guard_data))
        stop("joint designs are not supported")
      guard_data <- add_accumulators(
        guard_data, guard_design$matchfun, simulate = TRUE,
        type = model$type, Fcovariates = guard_design$Fcovariates,
        fixed_accumulator_roles = guard_design$fixed_accumulator_roles
      )
      design_model(
        guard_data, guard_design, guard_design$model, add_acc = FALSE,
        compress = FALSE, verbose = FALSE, rt_check = FALSE
      )
    }, error = function(e) {
      stop("Could not construct the parameter-bounds guard: ",
           conditionMessage(e), call. = FALSE)
    })

    guard_dadm_by_subj <- lapply(subnames, function(s) {
      row_idx <- which(guard_dadm$subjects == s)
      dadm_s <- .oo_subset_dadm(guard_dadm, row_idx)
      attr(dadm_s, "designs") <- .oo_expanded_designs(guard_dadm, row_idx, expand = FALSE)
      dadm_s$subjects <- factor(dadm_s$subjects, levels = s)
      dadm_s
    })
    names(guard_dadm_by_subj) <- subnames

    bound_checker <- function(draws) {
      sub_ids <- rownames(draws)
      vapply(seq_len(nrow(draws)), function(i) {
        s <- sub_ids[i]
        dadm_s <- guard_dadm_by_subj[[s]]
        if (is.null(dadm_s)) return(TRUE)
        p_mat <- draws[i, , drop = FALSE]
        mapped <- get_pars_oo(p_mat, dadm_s, guard_design$model)
        mapped <- model$Ttransform(mapped, dadm_s)
        # Use the same inclusive bound convention as make_data() itself.
        mapped <- fix_bound(mapped, model$bound, dadm_s$lR, fix = FALSE)
        ok <- attr(mapped, "ok")
        if (is.null(ok)) TRUE else all(ok)
      }, logical(1))
    }
  }

  random_effects <- matrix(
    NA_real_, nrow = n_subj, ncol = length(group_means),
    dimnames = list(subnames, sampled_names)
  )
  remaining <- seq_len(n_subj)
  for (attempt in seq_len(max_tries)) {
    if (!length(remaining)) break
    draws <- mvtnorm::rmvnorm(
      n = length(remaining), mean = group_means, sigma = covariances
    )
    draws <- as.matrix(draws)
    colnames(draws) <- sampled_names
    rownames(draws) <- subnames[remaining]
    valid <- if (is.null(bound_checker)) rep(TRUE, nrow(draws)) else
      bound_checker(draws)
    if (any(valid))
      random_effects[remaining[valid], ] <- draws[valid, , drop = FALSE]
    remaining <- remaining[!valid]
  }
  if (length(remaining)) {
    stop("Could not draw in-bound random effects for subject(s) ",
         paste(subnames[remaining], collapse = ", "), " after ",
         max_tries, " attempts")
  }
  return(random_effects)
}
