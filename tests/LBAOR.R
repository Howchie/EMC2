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
v_ST = 2; v_DT = 2; v_NT = 2; v_F=0.5; sv=1; B=2; A=0.5; t0=0.15
nsubs=1; ntrials=200
matchfun <- function(d) dplyr::case_when(d$S =="NT" & d$lR=="n_A" ~ TRUE,
                                         d$S =="NT" & d$lR=="n_B" ~ TRUE,
                                         d$S =="STA" & d$lR=="A" ~ TRUE,
                                         d$S =="STA" & d$lR=="n_B" ~ TRUE,
                                         d$S =="STB" & d$lR=="n_A" ~ TRUE,
                                         d$S =="STB" & d$lR=="B" ~ TRUE,
                                         d$S =="DT" & d$lR=="A" ~ TRUE,
                                         d$S =="DT" & d$lR=="B" ~ TRUE,
                                         TRUE ~ FALSE)
stimLMTRUE <- function(d) factor(dplyr::case_when(d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & d$S=="DT"~ "MatchDT",
                                                  d$lM =="TRUE" & (d$lR=="A"|d$lR=="B") & substr(d$S,1,2)=="ST"~ "MatchST",
                                                  d$lM =="TRUE" & (d$lR=="n_A"|d$lR=="n_B")~ "MatchNT",
                                         TRUE ~ "Mismatch"),levels=c("Mismatch","MatchNT","MatchST","MatchDT"))
designLBA <- design(
  fixed_accumulator_roles = factor(c("A","n_A","B","n_B"),levels=c("A","n_A","B","n_B")),
  matchfun=matchfun,
  model=LogicalRules,constants=c(sv=log(1),p=qnorm(1),q=qnorm(0.5)),
  formula=list(v~0+Stim,B~1,t0~1,A~1,sv~1,p~1,q~1),functions=list(Stim=stimLMTRUE),
  factors=list(S=c("DT","STA","STB","NT"),LogicalRule=rep("OR",nsubs),subjects=rep(1:nsubs,each=4)),Rlevels = c("no","yes")
)
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(v_F,v_NT,v_ST,v_DT,log(B),log(t0),log(A))
data = make_data(p_vector,designLBA,n_trials=200)
dadm <- EMC2:::design_model(data, designLBA, verbose = FALSE,add_acc = TRUE)
pars <- EMC2:::get_pars_matrix(p_vector, dadm, designLBA$model());min_ll = log(1e-10)
mapped_pars(designLBA,p_vector)



profile_plot_test(data,designLBA,p_vector,n_cores=1,layout=c(2,3)) # good
profile_plot_test(data,designLBA,p_vector,n_cores=1,layout=c(2,3),use_c=TRUE) # good
