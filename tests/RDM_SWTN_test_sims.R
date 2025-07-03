### Test implementation of RDMSWTN
devtools::load_all(reset = TRUE)
#Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
#source("forceSerial_TEST.R")         # your breakpoints now trigger
#Sys.unsetenv("PAR_DEBUG")     # remove when finished
options(error=recover)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
nsubs=10; ntrials=500
designRDM <- design(
  factors=list(subjects=1:nsubs,S=rep("yes","no",nsubs/2)),Rlevels=c("yes","no"),
  matchfun=matchfun,
  model=RDM,constants=c(s=log(1),A=log(2),B=log(2),t0=log(0.05)),
  formula=list(v~lM,B~1,t0~1,A~1,s~1),
)
designRDMSWTN <- design(
  factors=list(subjects=1:nsubs,S=rep("yes","no",nsubs/2)),Rlevels=c("yes","no"),
  matchfun=matchfun,
  model=RDM_SWTN,constants=c(s=log(1),A=log(2),B=log(2),t0=log(0.05), sv=log(0)),
  formula=list(v~0+lM,B~1,t0~1,A~1,s~1,sv~1),
)
p_vector_RDM <- c("v_lMFALSE"=log(2),"v_lMTRUE"=log(4))
p_vector_RDMSWTN <- c("v_lMFALSE"=log(2),"v_lMTRUE"=log(4))

reRDM = make_random_effects(design = designRDM, group_means = p_vector_RDM,n_subj=nsubs, variance_proportion = 0.05)
reRDMSWTN = make_random_effects(design = designRDMSWTN, group_means = p_vector_RDMSWTN,n_subj=nsubs, variance_proportion = 0.05)
# Make square data so can remove pm in RACE = 2
dataRDM <- make_data(reRDM,designRDM,n_trials=ntrials)
dataRDMSWTN <- make_data(reRDMSWTN,designRDMSWTN,n_trials=ntrials)

# Check likelihood
dadmRDM <- EMC2:::design_model(dataRDM,designRDM)
parsRDM <- EMC2:::get_pars_matrix(reRDM, dadmRDM, model = attr(dadmRDM, "model")())
dadmRDMSWTN <- EMC2:::design_model(dataRDMSWTN,designRDMSWTN)
parsRDMSWTN <- EMC2:::get_pars_matrix(reRDMSWTN, dadmRDMSWTN, model = attr(dadmRDMSWTN, "model")())

par(mfrow=c(2,2))
plot(dadmRDM$rt[dadmRDM$lM==TRUE],designRDM$model()$dfun(dadmRDM$rt[dadmRDM$lM==TRUE],parsRDM[dadmRDM$lM==TRUE,]), col='blue', xlim=c(0,5))
plot(dadmRDM$rt[dadmRDM$lM==TRUE],designRDM$model()$pfun(dadmRDM$rt[dadmRDM$lM==TRUE],parsRDM[dadmRDM$lM==TRUE,]), col='blue', xlim=c(0,5))
plot(dadmRDMSWTN$rt[dadmRDMSWTN$lM==TRUE],designRDMSWTN$model()$dfun(dadmRDMSWTN$rt[dadmRDMSWTN$lM==TRUE],parsRDMSWTN[dadmRDMSWTN$lM==TRUE,]), col='red', xlim=c(0,5))
plot(dadmRDMSWTN$rt[dadmRDMSWTN$lM==TRUE],designRDMSWTN$model()$pfun(dadmRDMSWTN$rt[dadmRDMSWTN$lM==TRUE],parsRDMSWTN[dadmRDMSWTN$lM==TRUE,]), col='red', xlim=c(0,5))

# Fit data from standard RDM using full RDMSWTN, and the standard RDM, use standard mu=0,sd=1 priors for now
# designRDM <- design(data=dataRDM,
#   model=RDM,constants=c(s=log(1)),
#   matchfun=matchfun,
#   formula=list(v~lM,B~1,t0~1,A~1,s~1),
# )
# designRDMSWTN <- design(data=dataRDM,
#                         model=RDM_SWTN,constants=c(s=log(1),sv=log(0)),
#                         matchfun=matchfun,
#                         formula=list(v~lM,B~1,t0~1,A~1,s~1,sv~1),
# )



fitRDM = make_emc(dataRDM,designRDM,rt_resolution = 0.02)
fitRDM = fit(fitRDM)
plot_pars(fitRDM,true_pars = p_vector_RDM)
#plot(dadmRDM$rt,designRDMSWTN$model()$dfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))
#plot(dadmRDM$rt,designRDMSWTN$model()$pfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))

#p_vector_RDMSWTN <- p_vector_RDMSWTN <- c("v"=log(15),"v_lMTRUE"=log(10),"A"=log(8),"B"=log(2),"t0"=log(0.15))
#prior_Model = prior(designRDMSWTN,mu_mean=p_vector_RDMSWTN,mu_sd=0.1)
fitRDMSWTN = make_emc(dataRDM,designRDMSWTN,rt_resolution = 0.02)
fitRDMSWTN = fit(fitRDMSWTN,fileName = 'samples.RData',verboseProgress = TRUE)
plot_pars(fitRDMSWTN,true_pars = p_vector_RDMSWTN)
#plot(dadmRDM$rt,designRDMSWTN$model()$dfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))
#plot(dadmRDM$rt,designRDMSWTN$model()$pfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))

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


profile_plot_test <- function (data, design, p_vector, range = 2, layout = NA, p_min = NULL,
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
profile_plot_test(dataRDM,designRDMSWTN,p_vector_RDM,n_cores=1,layout=c(2,3)) # good
profile_plot_test(dataRDM,designRDM,p_vector_RDM,n_cores=1,layout=c(2,3)) # good
profile_plot_test(dataRDM,designRDMSWTN,p_vector_RDM,n_cores=1,layout=c(2,3),use_c=TRUE) # ?
profile_plot_test(dataRDM,designRDM,p_vector_RDM,n_cores=1,layout=c(2,3),use_c=TRUE) # ?
toc()


tmpRDMSWTN=matrix(NA,nrow=2,ncol=6)
tmpRDMSWTN[1,] = c(2,3,2,0.2,1,1); tmpRDMSWTN[2,] = c(2,3,2,0.2,1,0)
colnames(tmpRDMSWTN) = names(designRDMSWTN$model()$p_types)
tmpRDM=matrix(NA,nrow=2,ncol=5)
tmpRDM[1,] = c(2,3,2,0.2,1); tmpRDM[2,] = c(2,3,2,0.2,1)
colnames(tmpRDM) = names(designRDM$model()$p_types)
drdmswtn_c(c(1.8,0.76),tmpRDMSWTN,c(TRUE,TRUE),1e-10,c(TRUE,TRUE))
dRDM_SWTN(c(1.8,0.76),tmpRDMSWTN)
drdm_c(c(1.8,0.76),tmpRDM,c(TRUE,TRUE),1e-10,c(TRUE,TRUE))
prdmswtn_c(c(1.8,0.76),tmpRDMSWTN,c(TRUE,TRUE),1e-10,c(TRUE,TRUE))
prdm_c(c(1.8,0.76),tmpRDM,c(TRUE,TRUE),1e-10,c(TRUE,TRUE))
tmp=drdmswtn_c(dataRDM$rt,parsRDMSWTN,rep(TRUE,nrow(dataRDM)),1e-10,rep(TRUE,nrow(dataRDM)))
tmp2=drdm_c(dadmRDM$rt,parsRDM,rep(TRUE,nrow(dadmRDM)),1e-10,rep(TRUE,nrow(dadmRDM)))
