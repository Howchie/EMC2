#### RACE RDMSWTN ----
devtools::load_all()
library(tictoc)
#Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
#source("forceSerial_TEST.R")         # your breakpoints now trigger
#RNGkind("L'Ecuyer-CMRG")
#set.seed(123)
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR) |
  (d$lR=="pm" & as.numeric(d$S)>2)
designRDMSWTN <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=RDMSWTN,constants=c(v_RACE3=0,v=log(2),"v_RACE3:lMTRUE"=log(2),s=log(1),sv=log(0),A=log(.5),t0=log(0.2)),
  formula=list(v~RACE*lM,B~1,t0~1,A~1,s~1,sv~1),
)
designRDM <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=RDM,constants=c(v_RACE3=0,v=log(2),"v_RACE3:lMTRUE"=log(2),s=log(1),A=log(.5),t0=log(0.2)),
  formula=list(v~RACE*lM,B~1,t0~1,A~1,s~1),
)

p_vector <- sampled_pars(designRDM,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(2), log(1))

# Make square data so can remove pm in RACE = 2
template <- make_data(p_vector,designRDM,n_trials=1000)
attr(template,"UC")=Inf
template <- template[!(template$RACE==2 & (template$S %in% c("leftpm","rightpm"))),]
dat <- make_data(p_vector,designRDM,data=template)
Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R) | (d$R=="pm" & as.numeric(d$S)>2)
tapply(Cfun(dat),dat[,c("S","RACE")],mean)
# dadm <- EMC2:::design_model(dat,designRDMSWTN,compress=FALSE)


# Check likelihood
dadmRDMSWTN <- EMC2:::design_model(dat,designRDMSWTN)
dadmRDM <- EMC2:::design_model(dat,designRDM)
pars <- EMC2:::get_pars_matrix(p_vector, dadmRDM, model = attr(dadmRDM, "model")())


lfun <- function(i, x, p_vector, pname, dadm, use_c) {
  p_vector[pname] <- x[i]
  if (use_c) {
    p_matrix <- matrix(p_vector,nrow=1)
    colnames(p_matrix) <- names(p_vector)
    model <- attr(dadm, "model")()
    p_types=names(model$p_types)
    designs <- list()
    for (p in p_types) {
      designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
    }
    constants <- attr(dadm,"constants")
    if (is.null(constants)) constants <- NA
      EMC2:::calc_ll(p_matrix, dadm, constants,designs,model$c_name,
                     model$bound,model$transform,model$pre_transform,p_types,log(1e-10),model$trend)
  } else {
      EMC2:::calc_ll_R(p_vector, attr(dadm, "model")(), dadm)
  }
}


profile_plot_test <- function (data, design, p_vector, range = 0.5, layout = NA, p_min = NULL,
                               p_max = NULL, use_par = NULL, n_point = 100, n_cores = 1,use_c = FALSE,
                               round = 3, true_args = list(), ...)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  dots <- list(...)
  if (!identical(names(p_min), names(p_max)))
    stop("p_min and p_max should be specified for the same parameters")
  if (!is.null(names(p_min)) & length(p_min) == length(use_par))
    names(p_min) <- use_par
  if (!is.null(names(p_max)) & length(p_max) == length(use_par))
    names(p_max) <- use_par
  if (is.null(use_par))
    use_par <- names(p_vector)
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(use_par),
                              nplots = 1))
  }
  else {
    par(mfrow = layout)
  }
  if (is.null(dots$dadm)) {
    dadm <- EMC2:::design_model(data, design, verbose = FALSE)
  }
  else {
    dadm <- dots$dadm
  }
  out <- data.frame(true = rep(NA, length(use_par)), max = rep(NA,
                                                               length(use_par)), miss = rep(NA, length(use_par)))
  rownames(out) <- use_par
  for (p in 1:length(p_vector)) {
    cur_name <- names(p_vector)[p]
    if (cur_name %in% use_par) {
      cur_par <- p_vector[p]
      pmax_cur <- cur_par + range/2
      pmin_cur <- cur_par - range/2
      if (!is.null(p_min)) {
        if (!is.na(p_min[cur_name])) {
          pmin_cur <- p_min[cur_name]
        }
      }
      if (!is.null(p_max)) {
        if (!is.na(p_max[cur_name])) {
          pmax_cur <- p_max[cur_name]
        }
      }
      x <- seq(pmin_cur, pmax_cur, length.out = n_point)
      x <- c(x, cur_par)
      x <- unique(sort(x))
      ll <- unlist(mclapply(1:length(x), lfun, dadm = dadm, use_c = use_c,
                            x = x, p_vector = p_vector, pname = cur_name,
                            mc.cores = n_cores))
      do.call(plot, c(list(x, ll), EMC2:::fix_dots_plot(EMC2:::add_defaults(dots,
                                                                            type = "l", xlab = cur_name, ylab = "LL"))))
      do.call(abline, c(list(v = cur_par), EMC2:::fix_dots_plot(EMC2:::add_defaults(true_args,
                                                                                    lty = 2))))
      out[cur_name, ] <- c(p_vector[cur_name], x[which.max(ll)],
                           p_vector[cur_name] - x[which.max(ll)])
    }
  }
  return(round(out, 3))
}

library(parallel)
tic()
profile_plot_test(dat,designRDM,p_vector,n_cores=1,layout=c(2,3)) # good
profile_plot_test(dat,designRDM,p_vector,n_cores=1,layout=c(2,3),use_c=TRUE) # ?
profile_plot_test(dat,designRDMSWTN,p_vector,n_cores=1,layout=c(2,3)) # good
profile_plot_test(dat,designRDMSWTN,p_vector,n_cores=1,layout=c(2,3),use_c=TRUE) # ?
toc()

emc <- make_emc(dat,designRDMSWTN,type="single")
emc <- fit(emc,cores_for_chains = 2,fileName = 'samples.RData')
#recovery(emc,p_vector)
plot_pars(emc)
