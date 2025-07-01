### Test implementation of RDMSWTN
devtools::load_all(reset = TRUE)
Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
source("forceSerial_TEST.R")         # your breakpoints now trigger
#Sys.unsetenv("PAR_DEBUG")     # remove when finished
options(error=recover)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR)
nsubs=10; ntrials=50
designRDM <- design(
  factors=list(subjects=1:nsubs,S=rep("yes","no",nsubs/2)),Rlevels=c("yes","no"),
  matchfun=matchfun,
  model=RDM,constants=c(s=log(1)),
  formula=list(v~lM,B~1,t0~1,A~1,s~1),
)
designRDMSWTN <- design(
  factors=list(subjects=1:nsubs,S=rep("yes","no",nsubs/2)),Rlevels=c("yes","no"),
  matchfun=matchfun,
  model=RDM_SWTN,constants=c(s=log(1)),
  formula=list(v~lM,B~1,t0~1,A~1,s~1,sv~1),
)
p_vector_RDM <- c("v"=log(2),"v_lMTRUE"=log(3),"A"=log(2),"B"=log(3),"t0"=log(0.15))
p_vector_RDMSWTN <- c("v"=log(2),"v_lMTRUE"=log(3),"sv"=log(1),"B"=log(3),"t0"=log(.2), "A" = log(2))

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
designRDM <- design(data=dataRDM,
  model=RDM,constants=c(s=log(1)),
  matchfun=matchfun,
  formula=list(v~lM,B~1,t0~1,A~1,s~1),
)
designRDMSWTN <- design(data=dataRDM,
                        model=RDM_SWTN,constants=c(s=log(1),sv=log(0)),
                        matchfun=matchfun,
                        formula=list(v~lM,B~1,t0~1,A~1,s~1,sv~1),
)



fitRDM = make_emc(dataRDM,designRDM,rt_resolution = 0.02)
fitRDM = fit(fitRDM)
plot_pars(fitRDM,true_pars = p_vector_RDM)
#plot(dadmRDM$rt,designRDMSWTN$model()$dfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))
#plot(dadmRDM$rt,designRDMSWTN$model()$pfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))

p_vector_RDMSWTN <- p_vector_RDMSWTN <- c("v"=log(15),"v_lMTRUE"=log(10),"A"=log(8),"B"=log(2),"t0"=log(0.15))
prior_Model = prior(designRDMSWTN,mu_mean=p_vector_RDMSWTN,mu_sd=0.1)
fitRDMSWTN = make_emc(dataRDM,designRDMSWTN,rt_resolution = 0.02)
fitRDMSWTN = fit(fitRDMSWTN,fileName = 'samples.RData',verboseProgress = TRUE)
plot_pars(fitRDMSWTN,true_pars = p_vector_RDMSWTN)
#plot(dadmRDM$rt,designRDMSWTN$model()$dfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))
#plot(dadmRDM$rt,designRDMSWTN$model()$pfun(dadmRDM$rt,p_vec), col='red', xlim=c(0,5))

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
tmp=drdmswtn_c(dadmRDMSWTN$rt,parsRDMSWTN,rep(TRUE,nrow(dadmRDMSWTN)),1e-10,rep(TRUE,nrow(dadmRDMSWTN)))
tmp2=prdmswtn_c(dadmRDMSWTN$rt,parsRDMSWTN,rep(TRUE,nrow(dadmRDMSWTN)),1e-10,rep(TRUE,nrow(dadmRDMSWTN)))
