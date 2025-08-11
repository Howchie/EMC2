make_missing <- function(data,LT=0,UT=Inf,LC=0,UC=Inf,
                         LCresponse=TRUE,UCresponse=TRUE,LCdirection=TRUE,UCdirection=TRUE)
{
  std <- c("names", "row.names", "class", "dim")
  censor <- function(data,L=0,U=Inf,Ld=TRUE,Ud=TRUE,Lr=TRUE,Ur=TRUE)
  {
    if (Ld) Ld <- -Inf else Ld <- NA
    if (Ud) Ud <- Inf else Ud <- NA
    snams <- levels(data$subjects)
    if (length(L)==1) L <- rep(L,nrow(data))
    if (length(U)==1) U <- rep(U,nrow(data))
    pick <- data$rt < L
    pick[is.na(pick)] <- FALSE
    data$rt[pick] <- Ld
    if (!Lr) data$R[pick] <- NA
    pick <- data$rt > U
    pick[is.na(pick)] <- FALSE
    data$rt[pick] <- Ud
    if (!Ur) data$R[pick] <- NA
    
    data
  }

  pick <- is.infinite(data$rt) | (data$rt>LT & data$rt<UT)
  pick[is.na(pick)] <- TRUE
  out <- censor(data[pick,],L=LC[pick],U=UC[pick],Lr=LCresponse,Ur=UCresponse,Ld=LCdirection,Ud=UCdirection)

  out
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
                      functions = NULL, LT=NULL,LC=NULL,UT=NULL,UC=NULL,...)
{
  # #' @param LT lower truncation bound below which data are removed (scalar or subject named vector)
  # #' @param UT upper truncation bound above which data are removed (scalar or subject named vector)
  # #' @param LC lower censoring bound (scalar or subject named vector)
  # #' @param UC upper censoring bound (scalar or subject named vector)
  # #' @param LCresponse Boolean, default TRUE, if false set LC response to NA
  # #' @param UCresponse Boolean, default TRUE, if false set UC response to NA
  # #' @param LCdirection Boolean, default TRUE, set LC rt to 0, else to NA
  # #' @param UCdirection Boolean, default TRUE, set LC rt to Inf, else to NA
  # #' @param force_direction Boolean, take direction from argument not data (default FALSE)
  # #' @param force_response Boolean, take response from argument not data (default FALSE)
  # #' @param rtContaminantNA Boolean, TRUE sets contaminant trial rt to NA, if FALSE
  # #' (the default) direction is taken from data or LCdirection or UCdirection (NB
  # #' if both of these are false an error occurs as then contamination is not identifiable).
  # #' @param return_Ffunctions if false covariates are not returned

  if (!is.null(staircase)){
    staircase <- check_staircase(staircase)
  }
  
  # #' @param Fcovariates either a data frame of covariate values with the same
  # #' number of rows as the data or a list of functions specifying covariates for
  # #' each trial. Must have names specified in the design Fcovariates argument.
  check_bounds <- FALSE
  
  if(is.null(LT)){if("LT"%in%colnames(data)) LT=data$LT else{LT <- attr(data,"LT")}}; if (is.null(LT)) LT <- 0
  if(is.null(UT)){if("UT"%in%colnames(data)) UT=data$UT else{UT <- attr(data,"UT")}}; if (is.null(UT)) UT <- Inf
  if(is.null(LC)){if("LC"%in%colnames(data)) LC=data$LC else{LC <- attr(data,"LC")}}; if (is.null(LC)) LC <- 0
  if(is.null(UC)){if("UC"%in%colnames(data)) UC=data$UC else{UC <- attr(data,"UC")}}; if (is.null(UC)) UC <- Inf
  # if(!is.null(LCresponse)){LCresponse<-TRUE}
  # if(!is.null(UCresponse)){UCresponse<-TRUE}
  # if(!is.null(LCdirection)){LCdirection<-TRUE}
  # if(!is.null(UCdirection)){UCdirection<-TRUE}
  # if(!is.null(force_direction)){force_direction<-TRUE}
  # if(!is.null(force_response)){force_response<-TRUE}
  LCresponse<-FALSE
  UCresponse<-FALSE
  LCdirection<-FALSE
  UCdirection<-TRUE
  force_direction<-TRUE
  force_response<-TRUE
  rtContaminantNA<-FALSE
  return_Ffunctions <- FALSE
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
  # ZH added transform for single vector pars to be 1-row matrix, preventing crash on single trial simulation
  if(is.null(dim(parameters))){
    par_names <- names(parameters)
    nsubs = ifelse(is.null(design$Ffactors$subjects),1,length(unique(design$Ffactors$subjects)))
    parameters <- matrix(parameters, nrow = nsubs, 
                         ncol=length(parameters), dimnames = list(as.character(seq(nsubs)), par_names),byrow=TRUE)
    if(is.null(colnames(parameters))) colnames(parameters) <- sampled_p_names
  } else{
    if(!(nrow(parameters) == length(unique(design$Ffactors$subjects)))){
      stop("input parameter matrix must have number of rows equal to number of subjects specified in design")
    }
    if(is.null(colnames(parameters))) colnames(parameters) <- sampled_p_names
    if(is.null(rownames(parameters))) rownames(parameters) <-1:length(unique(design$Ffactors$subjects))
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
    
    data <- minimal_design(design_in, covariates = list(...)$covariates,
                             drop_subjects = F, n_trials = n_trials, add_acc=F,
                           drop_R = F,UC=UC,UT=UT,LC=LC,LT=LT)
  } else {
    if(length(LT)==1) {data$LT = rep(LT,nrow(data))}
    else{data$LT=LT}
    if(length(UT)==1) {data$UT = rep(UT,nrow(data))}
    else{data$UT=UT}
    if(length(LC)==1) {data$LC = rep(LC,nrow(data))}
    else{data$LC=LC}
    if(length(UC)==1) {data$UC = rep(UC,nrow(data))}
    else{data$UC=UC}

    if (!force_direction) {
      ok <- data$rt==-Inf; ok[is.na(ok)] <- FALSE
      LCdirection <- any(ok)
      ok <- data$rt==Inf; ok[is.na(ok)] <- FALSE
      UCdirection=any(ok)
    }
    if (!force_response) {
      if (!any(is.infinite(data$rt)) & any(is.na(data$R))) {
        LCresponse <- UCresponse <- FALSE
      } else {
        ok <- data$rt==-Inf
        bad <- is.na(ok)
        LCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
        ok <- data$rt==Inf
        bad <- is.na(ok)
        UCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
      }
    }
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
  data <- design_model(
      add_accumulators(data,design$matchfun,simulate=TRUE,type=model()$type,Fcovariates=design$Fcovariates,fixed_accumulator_roles = design$fixed_accumulator_roles),
      design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
      rt_check=FALSE)
  pars <- t(apply(parameters, 1, do_pre_transform, model()))
  pars <- map_p(add_constants(pars,design$constants),data, model())
  if(!is.null(model()$trend) && attr(model()$trend, "pretransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(data, model()$trend, pars)
  }
  pars <- do_transform(pars, model())
  if(!is.null(model()$trend) && attr(model()$trend, "posttransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(data, model()$trend, pars)
  }
  pars <- model()$Ttransform(pars, data)
  pars <- add_bound(pars, model()$bound, data$lR)
  pars_ok <- attr(pars, 'ok')
  if(mean(!pars_ok) > .1){
    warning("More than 10% of parameter values fall out of model bounds, see <model_name>$bounds()")
    return(FALSE)
  }
  if ( any(dimnames(pars)[[2]]=="pContaminant") && any(pars[,"pContaminant"]>0) )
    pc <- pars[data$lR==levels(data$lR)[1],"pContaminant"] else pc <- NULL
  if (expand>1) {
    data$rep <- rep(1:expand, each = nrow(data))
    data <- data[ , c("rep", setdiff(names(data), "rep")) ] 
    pars <- apply(pars,2,rep,times=expand)
  }
  if (!is.null(staircase)) {
    attr(data, "staircase") <- staircase
  }

  if (any(names(data)=="RACE")) {
    Rrt <- RACE_rfun(data, pars, model)
  } else if (any(names(data)=="LogicalRule")) {
    if (grepl("substitution",model()$c_name)){
      Rrt <- LogicalRules_substitution_rfun(data, pars, model)
    } else {Rrt <- LogicalRules_rfun(data, pars, model)}
  } else {Rrt <- model()$rfun(data,pars)}
  dropNames <- c("lR","lM","lSmagnitude")
  
  if (!return_Ffunctions && !is.null(design$Ffunctions))
    dropNames <- c(dropNames,names(design$Ffunctions))
  if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
  std <- c("names", "row.names", "class", "dim")
  data <- data[, !(names(data) %in% dropNames)]
  for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[,i]
  data <- make_missing(data[,names(data)!="winner"],LT,UT,LC,UC,
    LCresponse,UCresponse,LCdirection,UCdirection)
  if ( !is.null(pc) ) {
    if (!any(is.infinite(data$rt)) & any(is.na(data$R)))
      stop("Cannot have contamination and censoring with no direction and response")
    contam <- rbinom(length(pc),1,pc)==1
    data[contam,"R"] <- NA
    if ( any(LC!=0) | any(is.finite(UC)) ) { # censoring
      if ( (LCdirection & UCdirection) &  !rtContaminantNA)
        stop("Cannot have contamination with a mixture of censor directions")
      if (rtContaminantNA & ((any(is.finite(LC)) & !LCresponse & !LCdirection) |
                              (any(is.finite(UC)) & !UCresponse & !UCdirection)))
        stop("Cannot have contamination and censoring with no direction and response")
      if (rtContaminantNA | (!LCdirection & !UCdirection)) data[contam,"rt"] <- NA else
        if (LCdirection) data[contam,"rt"] <- -Inf  else data[contam,"rt"] <- Inf
    } else data[contam,"rt"] <- NA
  }
  attr(data,"p_vector") <- parameters;
  data
}

RACE_rfun <- function(data, pars, model){
  Rrt <- matrix(ncol=2,nrow=dim(data)[1]/length(levels(data$lR)),
         dimnames=list(NULL,c("R","rt")))
  RACE <- data[data$lR==levels(data$lR)[1],"RACE"]
  ok <- as.numeric(data$lR) <= as.numeric(as.character(data$RACE))
  for (i in levels(RACE)) {
    pick <- data$RACE==i
    data_in <- data[pick,]
    data_in$lR <- factor(data$lR[pick & ok])
    tmp <- pars[pick,, drop=FALSE]
    attr(tmp, "ok") <- rep(T, ifelse(is.null(dim(tmp)),1,nrow(tmp)))
    Rrti <- model()$rfun(data_in,tmp)
    Rrti$R <- as.numeric(Rrti$R)
    Rrt[RACE==i,] <- as.matrix(Rrti)
  }
  Rrt <- data.frame(Rrt)
  Rrt$R <- factor(Rrt$R, labels = levels(data$lR), levels = 1:length(levels(data$lR)))
  return(Rrt)
}

LogicalRules_rfun <- function(data, pars, model){
  Rrti <- matrix(ncol=length(levels(data$lR)),nrow=dim(data)[1]/length(levels(data$lR)),
                 dimnames=list(NULL,levels(data$lR)))
  RACE <- levels(data$lR) # should be 4 unless using a non-standard implementation
  df = data.frame(LogicalRule = data$LogicalRule[data$lR==levels(data$lR)[1]],
                  RuleFollow = rbinom(dim(data)[1]/length(levels(data$lR)),1,pars[data$lR==levels(data$lR)[1],"p"])==1,
                  ChannelA = rbinom(dim(data)[1]/length(levels(data$lR)),1,1-pars[data$lR==levels(data$lR)[1],"q"])==1)
  for (i in RACE) {
    pick <- data$lR==i
    data_in <- data[pick,]
    data_in$lR <- factor(data$lR[pick])
    tmp <- pars[pick,, drop=FALSE]
    attr(tmp, "ok") <- rep(T, ifelse(is.null(dim(tmp)),1,nrow(tmp)))
    Rrti[,i] <- model()$rfun(data_in,tmp)$rt
  }
  Rrt = Rrti %>%
    as.data.frame() %>%
    dplyr::mutate(R = dplyr::case_when(df$RuleFollow & df$LogicalRule=="OR" & Rrti[,"A"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ "yes", # target finishes before at least one absent
                                 df$RuleFollow & df$LogicalRule=="OR" & Rrti[,"B"]<Rrti[,"A"] & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ "yes", # target finishes before at least one absent
                                 df$RuleFollow & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"] ~ "no",
                                 df$RuleFollow & df$LogicalRule=="AND" & Rrti[,"n_A"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"n_A"] | Rrti[,"B"]>Rrti[,"n_A"]) ~ "no", # absent finishes before at least one target
                                 df$RuleFollow & df$LogicalRule=="AND" & Rrti[,"n_B"]<Rrti[,"n_A"] & (Rrti[,"A"]>Rrti[,"n_B"] | Rrti[,"B"]>Rrti[,"n_B"]) ~ "no", # absent finishes before at least one target
                                 df$RuleFollow & df$LogicalRule=="AND" & Rrti[,"A"]<Rrti[,"n_A"] & Rrti[,"A"]<Rrti[,"n_B"] & Rrti[,"B"]<Rrti[,"n_A"] & Rrti[,"B"]<Rrti[,"n_B"]  ~ "yes",
                                 !df$RuleFollow & df$ChannelA & Rrti[,"A"]<Rrti[,"n_A"] ~ "yes",
                                 !df$RuleFollow & df$ChannelA & Rrti[,"n_A"]<Rrti[,"A"] ~ "no",
                                 !df$RuleFollow & !df$ChannelA & Rrti[,"B"]<Rrti[,"n_B"] ~ "yes",
                                 !df$RuleFollow & !df$ChannelA & Rrti[,"n_B"]<Rrti[,"B"] ~ "no"),
           rt = dplyr::case_when(df$RuleFollow & df$LogicalRule=="OR" & Rrti[,"A"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ Rrti[,"A"], # target finishes before at least one absent
                           df$RuleFollow & df$LogicalRule=="OR" & Rrti[,"B"]<Rrti[,"A"] & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ Rrti[,"B"], # target finishes before at least one absent
                           df$RuleFollow & df$LogicalRule=="OR" & (Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"])  ~ pmax(Rrti[,"n_A"],Rrti[,"n_B"]),
                           df$RuleFollow & df$LogicalRule=="AND" & Rrti[,"n_A"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"n_A"] | Rrti[,"B"]>Rrti[,"n_A"]) ~ Rrti[,"n_A"], # absent finishes before at least one target
                           df$RuleFollow & df$RuleFollow & df$LogicalRule=="AND" & Rrti[,"n_B"]<Rrti[,"n_A"] & (Rrti[,"A"]>Rrti[,"n_B"] | Rrti[,"B"]>Rrti[,"n_B"]) ~ Rrti[,"n_B"], # absent finishes before at least one target
                           df$RuleFollow & df$LogicalRule=="AND" & (Rrti[,"A"]<Rrti[,"n_A"] & Rrti[,"A"]<Rrti[,"n_B"] & Rrti[,"B"]<Rrti[,"n_A"] & Rrti[,"B"]<Rrti[,"n_B"])  ~ pmax(Rrti[,"A"],Rrti[,"B"]),
                           !df$RuleFollow & df$ChannelA & Rrti[,"A"]<Rrti[,"n_A"] ~ Rrti[,"A"],
                           !df$RuleFollow & df$ChannelA & Rrti[,"n_A"]<Rrti[,"A"] ~ Rrti[,"n_A"],
                           !df$RuleFollow & !df$ChannelA & Rrti[,"B"]<Rrti[,"n_B"] ~ Rrti[,"B"],
                           !df$RuleFollow & !df$ChannelA & Rrti[,"n_B"]<Rrti[,"B"] ~ Rrti[,"n_B"])) %>%
    dplyr::select(R,rt)
  if ("guess"%in%RACE) {
    Rrt = Rrt %>% 
      dplyr::mutate(R = dplyr::case_when(df$LogicalRule=="OR" & Rrti[,"guess"]<Rrti[,"A"] & Rrti[,"guess"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"guess"] | Rrti[,"n_B"]>Rrti[,"guess"]) ~ "yes", # guess finishes before either target and at least one absent,
                                          df$LogicalRule=="AND" & Rrti[,"guess"]<Rrti[,"n_A"] & Rrti[,"guess"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"guess"] | Rrti[,"B"]>Rrti[,"guess"]) ~ "no", TRUE~R), # guess finishes before either absent and at least one target),
                    rt = dplyr::case_when(df$LogicalRule=="OR" & Rrti[,"guess"]<Rrti[,"A"] & Rrti[,"guess"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"guess"] | Rrti[,"n_B"]>Rrti[,"guess"]) ~ Rrti[,"guess"], # guess finishes before either target and at least one absent
                                           df$LogicalRule=="AND" & Rrti[,"guess"]<Rrti[,"n_A"] & Rrti[,"guess"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"guess"] | Rrti[,"B"]>Rrti[,"guess"]) ~ Rrti[,"guess"],TRUE~rt))# guess finishes before either absent and at least one target
  }
  Rrt$R <- factor(Rrt$R,levels=c("no","yes"))
  return(Rrt)
}

LogicalRules_negdrift_rfun <- function(data, pars, model){
  Rrti <- matrix(ncol=length(levels(data$lR)),nrow=dim(data)[1]/length(levels(data$lR)),
                 dimnames=list(NULL,levels(data$lR)))
  RACE <- levels(data$lR) # should be 4 unless using a non-standard implementation
  for (i in RACE) {
    pick <- data$lR==i
    data_in <- data[pick,]
    data_in$lR <- factor(data$lR[pick])
    tmp <- pars[pick,, drop=FALSE]
    attr(tmp, "ok") <- rep(T, ifelse(is.null(dim(tmp)),1,nrow(tmp)))
    if (i=="A" | i=="B"){
      Rrti[,i] <- model()$rfun(data_in,tmp,FALSE)$rt
    } else {
      Rrti[,i] <- model()$rfun(data_in,tmp,TRUE)$rt
    }
  }
  df = data.frame(LogicalRule = data$LogicalRule[data$lR==levels(data$lR)[1]],
                  AFail=is.infinite(Rrti[,"A"]),
                  BFail=is.infinite(Rrti[,"B"]))
  Rrt = Rrti %>%
    as.data.frame() %>%
    dplyr::mutate(R = dplyr::case_when(!df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"A"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ "yes", # target finishes before at least one absent
                                       !df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"B"]<Rrti[,"A"] & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ "yes", # target finishes before at least one absent
                                       !df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"] ~ "no",
                                       !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_A"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"n_A"] | Rrti[,"B"]>Rrti[,"n_A"]) ~ "no", # absent finishes before at least one target
                                       !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_B"]<Rrti[,"n_A"] & (Rrti[,"A"]>Rrti[,"n_B"] | Rrti[,"B"]>Rrti[,"n_B"]) ~ "no", # absent finishes before at least one target
                                       !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"A"]<Rrti[,"n_A"] & Rrti[,"A"]<Rrti[,"n_B"] & Rrti[,"B"]<Rrti[,"n_A"] & Rrti[,"B"]<Rrti[,"n_B"]  ~ "yes",
                                       df$BFail & !df$AFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ "yes",
                                       df$AFail & !df$BFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ "yes",
                                       df$BFail & !df$AFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] ~ "no",
                                       df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"] ~ "no",
                                       df$AFail & df$BFail ~ "no",
                                       (df$BFail | df$AFail) & df$LogicalRule=="AND" ~ "no"
                                       ),
                  rt = dplyr::case_when(!df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"A"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ Rrti[,"A"], # target finishes before at least one absent
                                        !df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"B"]<Rrti[,"A"] & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ Rrti[,"B"], # target finishes before at least one absent
                                        !df$AFail & !df$BFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"])  ~ pmax(Rrti[,"n_A"],Rrti[,"n_B"]),
                                        !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_A"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"n_A"] | Rrti[,"B"]>Rrti[,"n_A"]) ~ Rrti[,"n_A"], # absent finishes before at least one target
                                        !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_B"]<Rrti[,"n_A"] & (Rrti[,"A"]>Rrti[,"n_B"] | Rrti[,"B"]>Rrti[,"n_B"]) ~ Rrti[,"n_B"], # absent finishes before at least one target
                                        !df$AFail & !df$BFail & df$LogicalRule=="AND" & (Rrti[,"A"]<Rrti[,"n_A"] & Rrti[,"A"]<Rrti[,"n_B"] & Rrti[,"B"]<Rrti[,"n_A"] & Rrti[,"B"]<Rrti[,"n_B"])  ~ pmax(Rrti[,"A"],Rrti[,"B"]),
                                        df$BFail & !df$AFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ Rrti[,"A"],
                                        df$AFail & !df$BFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ Rrti[,"B"],
                                        df$BFail & !df$AFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] ~ pmax(Rrti[,"n_A"],Rrti[,"n_B"]),
                                        df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"] ~ pmax(Rrti[,"n_A"],Rrti[,"n_B"]),
                                        df$AFail & df$BFail & df$LogicalRule=="OR" ~ pmax(Rrti[,"n_A"],Rrti[,"n_B"]),
                                        df$AFail & df$BFail & df$LogicalRule=="AND" ~ pmin(Rrti[,"n_A"],Rrti[,"n_B"]),
                                        df$BFail & !df$AFail & df$LogicalRule=="AND" ~ pmin(Rrti[,"n_A"],Rrti[,"n_B"]), 
                                        df$AFail & !df$BFail & df$LogicalRule=="AND" ~ pmin(Rrti[,"n_A"],Rrti[,"n_B"])
                                        )) %>%
    dplyr::select(R,rt)
  if ("guess"%in%RACE) {
    Rrt = Rrt %>% 
      dplyr::mutate(R = dplyr::case_when(df$LogicalRule=="OR" & Rrti[,"guess"]<Rrti[,"A"] & Rrti[,"guess"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"guess"] | Rrti[,"n_B"]>Rrti[,"guess"]) ~ "yes", # guess finishes before either target and at least one absent,
                                         df$LogicalRule=="AND" & Rrti[,"guess"]<Rrti[,"n_A"] & Rrti[,"guess"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"guess"] | Rrti[,"B"]>Rrti[,"guess"]) ~ "no", TRUE~R), # guess finishes before either absent and at least one target),
                    rt = dplyr::case_when(df$LogicalRule=="OR" & Rrti[,"guess"]<Rrti[,"A"] & Rrti[,"guess"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"guess"] | Rrti[,"n_B"]>Rrti[,"guess"]) ~ Rrti[,"guess"], # guess finishes before either target and at least one absent
                                          df$LogicalRule=="AND" & Rrti[,"guess"]<Rrti[,"n_A"] & Rrti[,"guess"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"guess"] | Rrti[,"B"]>Rrti[,"guess"]) ~ Rrti[,"guess"],TRUE~rt))# guess finishes before either absent and at least one target
  }
  Rrt$R <- factor(Rrt$R,levels=c("no","yes"))
  return(Rrt)
}

## Important - due to the correlated drift rates in this function I sample the actual drifts inside here, and the corresponding LBA rfun just uses them directly
LogicalRules_substitution_rfun <- function(data, pars, model){
  Rrti <- matrix(ncol=length(levels(data$lR)),nrow=dim(data)[1]/length(levels(data$lR)),
                 dimnames=list(NULL,levels(data$lR)))
  RACE <- levels(data$lR) # should be 4 unless using a non-standard implementation
  accs = levels(data$lR); 
  
  # 
  lower=numeric(length=length(accs));names(lower)=accs
  if(grepl("negdrift",model()$c_name)){
    lower["A"]=-Inf;lower["B"]=-Inf
  }
  sampled_v=numeric(length(pars))
  for (i in unique(data$trials)) {
    adj_vars = pars[data$trials==i,"adj_sv"]^2; adj_vs=pars[data$trials==i,"adj_v"]
    covariances=matrix(0,nrow=length(accs),ncol=length(accs),dimnames=list(accs,accs))
    diag(covariances)=adj_vars
    tauA=pars[data$trials==i&data$lR=="A","tau"];tauB=pars[data$trials==i&data$lR=="B","tau"]
    muvA=pars[data$trials==i&data$lR=="A","v"];muvB=pars[data$trials==i&data$lR=="B","v"]
    if(tauA==tauB){tau2=tauA^2}else{stop("tau parameter mismatch")}
    covariances["A","B"] = (tau2)*muvA*muvB;covariances["B","A"] = (tau2)*muvA*muvB
    sampled_v[data$trials==i]=tmvtnorm::rtmvnorm(1,mean=adj_vs,sigma=covariances,lower=lower)
  }
  pars=cbind(pars,"sampled_v"=sampled_v)
  for (i in RACE) {
    pick <- data$lR==i
    data_in <- data[pick,]
    data_in$lR <- factor(data$lR[pick])
    tmp <- pars[pick,, drop=FALSE]
    attr(tmp, "ok") <- rep(T, ifelse(is.null(dim(tmp)),1,nrow(tmp)))
    Rrti[,i] <- model()$rfun(data_in,tmp)$rt
  }
  df = data.frame(LogicalRule = data$LogicalRule[data$lR==levels(data$lR)[1]],
                  AFail=is.infinite(Rrti[,"A"]),
                  BFail=is.infinite(Rrti[,"B"]))
  Rrt = Rrti %>%
    as.data.frame() %>%
    dplyr::mutate(R = dplyr::case_when(!df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"A"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ "yes", # target finishes before at least one absent
                                       !df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"B"]<Rrti[,"A"] & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ "yes", # target finishes before at least one absent
                                       !df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"] ~ "no",
                                       !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_A"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"n_A"] | Rrti[,"B"]>Rrti[,"n_A"]) ~ "no", # absent finishes before at least one target
                                       !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_B"]<Rrti[,"n_A"] & (Rrti[,"A"]>Rrti[,"n_B"] | Rrti[,"B"]>Rrti[,"n_B"]) ~ "no", # absent finishes before at least one target
                                       !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"A"]<Rrti[,"n_A"] & Rrti[,"A"]<Rrti[,"n_B"] & Rrti[,"B"]<Rrti[,"n_A"] & Rrti[,"B"]<Rrti[,"n_B"]  ~ "yes",
                                       df$BFail & !df$AFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B_flip"]>Rrti[,"A"]) ~ "yes",
                                       df$AFail & !df$BFail & df$LogicalRule=="OR" & (Rrti[,"n_A_flip"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ "yes",
                                       df$BFail & !df$AFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B_flip"]<Rrti[,"A"] ~ "no",
                                       df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"n_A_flip"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"] ~ "no",
                                       df$AFail & df$BFail ~ "no",
                                       (df$BFail | df$AFail) & df$LogicalRule=="AND" ~ "no"
    ),
    rt = dplyr::case_when(!df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"A"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B"]>Rrti[,"A"]) ~ Rrti[,"A"], # target finishes before at least one absent
                          !df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"B"]<Rrti[,"A"] & (Rrti[,"n_A"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ Rrti[,"B"], # target finishes before at least one absent
                          !df$AFail & !df$BFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B"]<Rrti[,"A"] & Rrti[,"n_A"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"])  ~ pmax(Rrti[,"n_A"],Rrti[,"n_B"]),
                          !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_A"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"n_A"] | Rrti[,"B"]>Rrti[,"n_A"]) ~ Rrti[,"n_A"], # absent finishes before at least one target
                          !df$AFail & !df$BFail & df$LogicalRule=="AND" & Rrti[,"n_B"]<Rrti[,"n_A"] & (Rrti[,"A"]>Rrti[,"n_B"] | Rrti[,"B"]>Rrti[,"n_B"]) ~ Rrti[,"n_B"], # absent finishes before at least one target
                          !df$AFail & !df$BFail & df$LogicalRule=="AND" & (Rrti[,"A"]<Rrti[,"n_A"] & Rrti[,"A"]<Rrti[,"n_B"] & Rrti[,"B"]<Rrti[,"n_A"] & Rrti[,"B"]<Rrti[,"n_B"])  ~ pmax(Rrti[,"A"],Rrti[,"B"]),
                          df$BFail & !df$AFail & df$LogicalRule=="OR" & (Rrti[,"n_A"]>Rrti[,"A"] | Rrti[,"n_B_flip"]>Rrti[,"A"]) ~ Rrti[,"A"],
                          df$AFail & !df$BFail & df$LogicalRule=="OR" & (Rrti[,"n_A_flip"]>Rrti[,"B"] | Rrti[,"n_B"]>Rrti[,"B"]) ~ Rrti[,"B"],
                          df$BFail & !df$AFail & df$LogicalRule=="OR" & Rrti[,"n_A"]<Rrti[,"A"] & Rrti[,"n_B_flip"]<Rrti[,"A"] ~ pmax(Rrti[,"n_A"],Rrti[,"n_B_flip"]),
                          df$AFail & !df$BFail & df$LogicalRule=="OR" & Rrti[,"n_A_flip"]<Rrti[,"B"] & Rrti[,"n_B"]<Rrti[,"B"] ~ pmax(Rrti[,"n_A_flip"],Rrti[,"n_B"]),
                          df$AFail & df$BFail & df$LogicalRule=="OR" ~ pmax(Rrti[,"n_A_flip"],Rrti[,"n_B_flip"]),
                          df$AFail & df$BFail & df$LogicalRule=="AND" ~ pmin(Rrti[,"n_A_flip"],Rrti[,"n_B_flip"]),
                          df$BFail & !df$AFail & df$LogicalRule=="AND" ~ pmin(Rrti[,"n_A"],Rrti[,"n_B_flip"]),
                          df$AFail & !df$BFail & df$LogicalRule=="AND" ~ pmin(Rrti[,"n_A_flip"],Rrti[,"n_B"])
    )) %>%
    dplyr::select(R,rt)
  if ("guess"%in%RACE) {
    Rrt = Rrt %>%
      dplyr::mutate(R = dplyr::case_when(df$LogicalRule=="OR" & Rrti[,"guess"]<Rrti[,"A"] & Rrti[,"guess"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"guess"] | Rrti[,"n_B"]>Rrti[,"guess"]) ~ "yes", # guess finishes before either target and at least one absent,
                                         df$LogicalRule=="AND" & Rrti[,"guess"]<Rrti[,"n_A"] & Rrti[,"guess"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"guess"] | Rrti[,"B"]>Rrti[,"guess"]) ~ "no", TRUE~R), # guess finishes before either absent and at least one target),
                    rt = dplyr::case_when(df$LogicalRule=="OR" & Rrti[,"guess"]<Rrti[,"A"] & Rrti[,"guess"]<Rrti[,"B"] & (Rrti[,"n_A"]>Rrti[,"guess"] | Rrti[,"n_B"]>Rrti[,"guess"]) ~ Rrti[,"guess"], # guess finishes before either target and at least one absent
                                          df$LogicalRule=="AND" & Rrti[,"guess"]<Rrti[,"n_A"] & Rrti[,"guess"]<Rrti[,"n_B"] & (Rrti[,"A"]>Rrti[,"guess"] | Rrti[,"B"]>Rrti[,"guess"]) ~ Rrti[,"guess"],TRUE~rt))# guess finishes before either absent and at least one target
  }
  Rrt$R <- factor(Rrt$R,levels=c("no","yes"))
  return(Rrt)
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

make_random_effects <- function(design, group_means, n_subj = NULL, variance_proportion = .2, covariances = NULL,covariance_prop=FALSE){
  ## ZH edits -- covariances and variance prop now influence parameters on their NATURAL SCALE, bounded correctly
  # variance prop can be a scalar, or a vector length group_means
  
  if(is.null(n_subj)){
    n_subj <- length(design$Ffactors$subjects)
    subnames <- design$Ffactors$subjects
  } else{
    subnames <- as.character(1:n_subj)
  }
  if(length(group_means) != length(sampled_pars(design))) stop("You must specify as many means as parameters in your design")
  if(!length(variance_proportion)==1 && !length(variance_proportion)==length(group_means)) stop ("variance proportion must be a scalar or exactly one value per group mean")
  if (is.null(dim(group_means))) { # Check if pars is a vector
    original_names <- names(group_means)
    group_means <- matrix(group_means, nrow = 1, ncol=length(group_means), dimnames = list(NULL, original_names),byrow=TRUE)
  }
  # re-order pars
  reordered_means = matrix(NA,nrow=nrow(group_means),ncol=ncol(group_means), dimnames = list(NULL, names(sampled_pars(design))))
  for (p in names(sampled_pars(design))) {
    reordered_means[,p] = group_means[,p]
  }
  if(is.null(covariances)&length(variance_proportion )==1) {
    variance_proportion =rep(variance_proportion ,ncol(reordered_means))
  }
  model=design$model()
  reordered_means <- t(apply(reordered_means, 1, do_pre_transform, model))
  if(!is.null(model$trend) && attr(model$trend, "pretransform")){
    # This runs the trend and afterwards removes the trend parameters
    group_means <- prep_trend(design, model$trend, reordered_means)
  }
  reordered_means <- do_transform(reordered_means, model)
  if(!is.null(model$trend) && attr(model$trend, "posttransform")){
    # This runs the trend and afterwards removes the trend parameters
    reordered_means <- prep_trend(design, model$trend, reordered_means)
  }
  # From here ZH has edited functionality
  # If variance_prop and no covariances, multiply natural-scale mean * variance prop
  # perform rmvnorm on natural scale
  # transform (and bound check post-hoc)
  if(is.null(covariances)) { # ZH modified so that variance is transformed to natural scale (e.g. 0.2*natural_scale_mu) then back-converted for transformed mvtnorm. I found the original code (a) was broken for group_means of zero (zero variance) but also estimates were more varied than expected due to the conversions
    covariances = matrix(0,nrow=ncol(reordered_means),ncol = ncol(reordered_means))
    diag(covariances)=as.numeric(reordered_means)*variance_proportion
  } 
  else if (covariance_prop) { # if covariance matrix is specified, is the diagonal to be treated like variance_proportion?
    diag(covariances)=diag(covariances)*as.numeric(reordered_means)
  }
  # If covariances specified, assume it is on natural scale
  # diag(covariances) = natural-scale mean * diag(covariances)
  # perform rmvnorm on natural scale
  # transform (and bound check post-hoc)
  # Use truncated mnvorm to respect bounds
  random_effects <- tmvtnorm::rtmvnorm(n_subj,mean=as.numeric(reordered_means),sigma=covariances,
                                       lower=as.numeric(model$bound$minmax[1, get_p_types(colnames(reordered_means))]),
                                       upper=as.numeric(model$bound$minmax[2, get_p_types(colnames(reordered_means))]),
                                                                algorithm="rejection")
  colnames(random_effects) <- colnames(reordered_means)
  rownames(random_effects) <- subnames
  random_effects <- do_reverse_transform(random_effects,model)
  return(random_effects)
}


