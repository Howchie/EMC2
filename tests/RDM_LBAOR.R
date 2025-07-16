## Script to simulate from, and test recovery of, RDMSWTN
library(EMC2)
library(sft)
library(dplyr)
devtools::load_all(reset = TRUE)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
source("tests/test_likelihood_plotfuns.R")
## Compile sims manually -- no salience effect for now
ntrials = 1000
v_ST = 2; v_DT = 2; v_NT = 2; v_F=0.5; sv=1; b=2; A=0.5; t0=0.15
stim = c("NT","STA","STB","DT")
vs = matrix(NA, nrow=4, ncol=4,dimnames=list(NULL,stim)); vs[,1]=c(v_F,v_NT,v_F,v_NT);vs[,2]=c(v_ST,v_F,v_F,v_NT) 
vs[,3] = c(v_F,v_NT,v_ST,v_F); vs[,4] = c(v_DT,v_F,v_DT,v_F) 
default_pars = c("v"=NA,"sv"=1,"b"=2,"A"=0.5,"t0"=.15); data = NULL
for (i in stim) {
  pars = matrix(default_pars,nrow=4,ncol=5,byrow = TRUE,dimnames=list(NULL,names(default_pars)))
  pars[,"v"] = vs[1:4,i]
  tmp1 = as.data.frame(rtdists::rlba_norm(ntrials,pars[1,"A"],pars[1,"b"],pars[1,"t0"],mean_v=pars[1,"v"],sd_v=pars[1,"sv"]))
  tmp2 = as.data.frame(rtdists::rlba_norm(ntrials,pars[3,"A"],pars[1,"b"],pars[3,"t0"],mean_v=pars[3,"v"],sd_v=pars[3,"sv"]))
  tmp3 = as.data.frame(rtdists::rlba_norm(ntrials,pars[2,"A"],pars[1,"b"],pars[2,"t0"],mean_v=pars[2,"v"],sd_v=pars[2,"sv"]))
  tmp4 = as.data.frame(rtdists::rlba_norm(ntrials,pars[4,"A"],pars[4,"b"],pars[4,"t0"],mean_v=pars[4,"v"],sd_v=pars[4,"sv"]))
  tmpdata = data.frame(rt=rep(NA,nrow(tmp1)),R=rep(NA,nrow(tmp1)))
  tmpdata = tmpdata %>%
    mutate(R = case_when (tmp1$rt<tmp2$rt & (tmp3$rt>tmp1$rt | tmp4$rt>tmp1$rt) ~ "yes", # target finishes before at least one absent
                          tmp2$rt<tmp1$rt & (tmp3$rt>tmp2$rt | tmp4$rt>tmp2$rt) ~ "yes", # target finishes before at least one absent
                          (tmp3$rt<tmp1$rt & tmp4$rt<tmp1$rt & tmp3$rt<tmp2$rt & tmp4$rt<tmp2$rt)  ~ "no"),
           rt = case_when (tmp1$rt<tmp2$rt & (tmp3$rt>tmp1$rt | tmp4$rt>tmp1$rt) ~ tmp1$rt, # target finishes before at least one absent
                           tmp2$rt<tmp1$rt & (tmp3$rt>tmp2$rt | tmp4$rt>tmp2$rt) ~ tmp2$rt, # target finishes before at least one absent
                           (tmp3$rt<tmp1$rt & tmp4$rt<tmp1$rt & tmp3$rt<tmp2$rt & tmp4$rt<tmp2$rt)  ~ pmax(tmp3$rt,tmp4$rt)),
    S=rep(i,nrow(tmp1))
  )
  # if (i=="NT") {
  #   tmpdata$Correct = as.numeric(tmpdata$R=="no") 
  # } else {
  #   tmpdata$Correct= as.numeric(tmpdata$R=="yes")
  # }
  data = rbind(data,tmpdata)
}
data$subjects=factor(1);data$R = factor(data$R,levels=c("yes","no")); data$S=factor(data$S, levels= c("NT","STA","STB","DT"))
# Sys.setenv(PAR_DEBUG = "1")   # turn serial mode on
# if (identical(Sys.getenv("PAR_DEBUG"), "1")) {      # opt-in via env-var
#   env <- asNamespace("parallel")
#   unlockBinding("parLapply", env)                   # temporarily unlock
#   assign("parLapply", function(cl, X, FUN, ...) {
#     lapply(X, FUN, ...)                             # <- runs in master R
#   }, envir = env)
#   lockBinding("parLapply", env)
# }
## Set simulation parameters:
nsubs=1; ntrials=200
## For maximum testing purposes we'll simulate from a model with varying n_acc
# matchfun <- function(d) dplyr::case_when((as.numeric(d$Channel1)>0 | as.numeric(d$Channel2)>0) & d$R=="yes" ~ TRUE,
#                                          (as.numeric(d$Channel1)==0 | as.numeric(d$Channel2)==0) & d$R=="no" ~ TRUE,
#                                          TRUE ~ FALSE)
matchfun <- function(d) dplyr::case_when(d$S =="NT" & d$lR=="n_A" ~ TRUE,
                                         d$S =="NT" & d$lR=="n_B" ~ TRUE,
                                         d$S =="STA" & d$lR=="A" ~ TRUE,
                                         d$S =="STA" & d$lR=="n_B" ~ TRUE,
                                         d$S =="STB" & d$lR=="n_A" ~ TRUE,
                                         d$S =="STB" & d$lR=="B" ~ TRUE,
                                         d$S =="DT" & d$lR=="A" ~ TRUE,
                                         d$S =="DT" & d$lR=="B" ~ TRUE,
                                         TRUE ~ FALSE)
designLBA <- design(
  fixed_accumulator_roles = factor(c("A","n_A","B","n_B"),levels=c("A","n_A","B","n_B")),
  matchfun=matchfun,
  model=LBAOR,constants=c(sv=log(1)),
  formula=list(v~0+lM,B~1,t0~1,A~1,sv~1),
  data=data
)
dadm <- EMC2:::design_model(data, designLBA, verbose = FALSE);model=designLBA$model()
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(v_F,v_ST,log(default_pars["b"]-default_pars["A"]),log(default_pars["t0"]),log(default_pars["A"]))
pars <- get_pars_matrix(p_vector, dadm, designLBA$model());min_ll = log(1e-10)
mapped_pars(designLBA,p_vector)


profile_plot_test(data,designLBA,p_vector,n_cores=1,layout=c(2,3)) # good
profile_plot_test(data,designLBA,p_vector,n_cores=1,layout=c(2,3),use_c=TRUE) # good
