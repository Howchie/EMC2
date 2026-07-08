rm(list = ls())
library(EMC2)
library(dplyr)
target_omission <- 0.10
target_response <- 1 - target_omission
drift <- 3
threshold <- 1.1
t0 <- 0.2

designRDMKILL1 <- design(
    factors = list(
        S = "Target", subjects = 1
    ),
    Rlevels = c("Go"),
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mK ~ 1),
    constants = c(s = log(1), sv=log(0), A=log(0)),
    model = RDMSWTN(erlang_type = "local_kill", erlang_shape=1), UC = 3
)

designRDMKILL2 <- design(
  factors = list(
    S = "Target", subjects = 1
  ),
  Rlevels = c("Go"),
  formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mK ~ 1),
  constants = c(s = log(1), sv=log(0), A=log(0)),
  model = RDMSWTN(erlang_type = "local_kill", erlang_shape=2), UC = 3
)

designRDM <- design(
  factors = list(
    S = "Target", subjects = 1 
  ),
  Rlevels = c("Go"),
  formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, pContaminant~1),
  constants = c(s = log(1), A=log(0)),
  model = RDM, UC = 3
)

designRDM_noMix <- design(
  factors = list(
    S = "Target", subjects = 1 
  ),
  Rlevels = c("Go"),
  formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
  constants = c(s = log(1), A=log(0)),
  model = RDM, UC = 3
)

designBAwL <- design(
  factors = list(
    S = "Target", subjects = 1
  ),
  Rlevels = c("Go"),
  formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1, k~1, mK~1),
  constants = c(sv = log(1), k=log(0)),
  model = BAwL(erlang_type = "local_kill", erlang_shape=1), UC = 3
)

designLBA <- design(
  factors = list(
    S = "Target", subjects = 1
  ),
  Rlevels = c("Go"),
  formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, sv ~ 1, pContaminant~1),
  constants = c(sv = log(1)),
  model = LBA, UC = 3
)

p_vec <- sampled_pars(designRDM)
p_vec["v"] = log(c(drift))
p_vec["B"] <- log(threshold)
p_vec["t0"] <- log(t0)
p_vec["pContaminant"] <- qnorm(target_omission)

dat <- make_data(p_vec, design = designRDM, n_trials = 20000)
dat_noPC <- make_data(p_vec[1:3], design = designRDM_noMix, n_trials = 20000)
p_vec <- sampled_pars(designLBA)
p_vec["v"] = drift
p_vec["B"] <- log(threshold)
p_vec["t0"] <- log(t0)
p_vec["A"] <- log(0.4)
p_vec["pContaminant"] <- qnorm(target_omission)
datLBA <- make_data(p_vec, design = designLBA, n_trials = 20000)
mean(is.finite(dat$rt))
mean(dat$rt[is.finite(dat$rt)])
mean(is.finite(datLBA$rt))
mean(datLBA$rt[is.finite(datLBA$rt)])

emc1 = make_emc(datLBA, designRDMKILL1, type="single")
emc1 = fit(emc1)
pred1 = predict(emc1)

emc2 = make_emc(datLBA, designRDMKILL2, type="single")
emc2 = fit(emc2)
pred2 = predict(emc2)

emc3 = make_emc(datLBA, designRDM, type="single")
emc3 = fit(emc3)
pred3 = predict(emc3)

emc4 = make_emc(datLBA, designBAwL, type="single")
emc4 = fit(emc4)
pred4 = predict(emc4)

emc5 = make_emc(datLBA, designLBA, type="single")
emc5 = fit(emc5)
pred5 = predict(emc5)

comp = compare(list(RDMKILL1=emc1,RDMKILL2=emc2,RDMpC=emc3,LBAKILL1=emc4,LBA=emc5),BayesFactor = TRUE)

## Get parameter estimates
p1 = credint(emc1)
p2 = credint(emc2)
p3 = credint(emc3)
p4 = credint(emc4)
p5 = credint(emc5)
tilted_drift <- drift - log(target_response) / threshold
kill_rate <- (tilted_drift^2 - drift^2) / 2

# --- True Theoretical kill rate calculation ---
# The true probability of a response when an exponential kill process (rate lambda)
# races against a Wald accumulation process (with t0 and uniform start-point variability A).
get_response_prob <- function(lambda, v, B, A, t0) {
  if (lambda == 0) return(1)
  eta <- sqrt(v^2 + 2 * lambda) - v
  # Survival during non-decision time
  prob_t0 <- exp(-lambda * t0)
  
  # Survival during evidence accumulation
  if (A <= 1e-6) {
    prob_eam <- exp(-B * eta)
  } else {
    prob_eam <- (exp(-B * eta) - exp(-(B + A) * eta)) / (A * eta)
  }
  return(prob_t0 * prob_eam)
}

get_response_prob_erlang2 <- function(lambda, v, B, t0) {                                                               
  if (lambda == 0) return(1)                                                                                            
  
  # Base Laplace transform terms                                                                                        
  nu <- sqrt(v^2 + 2 * lambda)                                                                                          
  eta <- nu - v                                                                                                         
  
  # M(lambda) is the Laplace transform of the Wald                                                                      
  M_lambda <- exp(-B * eta)                                                                                             
  
  # minus_M_prime is the exact analytical derivative: -M'(lambda)                                                       
  minus_M_prime <- (B / nu) * exp(-B * eta)                                                                             
  
  # Combine with t0 using the Erlang-2 survivor expectation                                                             
  prob_t0_base <- exp(-lambda * t0)                                                                                     
  prob_response <- prob_t0_base * ((1 + lambda * t0) * M_lambda + lambda * minus_M_prime)                               
  
  return(prob_response)                                                                                                 
}

# Find the exact rate (lambda) that produces the target response probability
true_lambda <- uniroot(function(lambda) {
  get_response_prob(lambda, v = drift, B = threshold, A = 0, t0 = t0) - target_response
}, interval = c(0, 10))$root

# Find the exact rate (lambda) that produces the 90% target response for Erlang-2                                       
true_lambda_e2 <- uniroot(function(lambda) {                                                                            
  get_response_prob_erlang2(lambda, v = drift, B = threshold, t0 = t0) - target_response                                
}, interval = c(0, 10))$root 

expected_mK <- 1 / true_lambda
cat("Naïve kill_rate (ignoring t0 and A):", kill_rate, "\n")
cat("True theoretical kill_rate (lambda):", true_lambda, "\n")
cat("Actual estimated kill estimate:", 1/exp(p1$`1`["mK","50%"]), "\n")
cat("Actual estimated kill estimate:", 2/(exp(p2$`1`["mK","50%"])), "\n")


# Plot cuts
plot(density(dat_noPC$rt[is.finite(dat_noPC$rt)]),col='blue')
lines(density(pred3$rt[is.finite(pred3$rt)]),col='red',lty='dashed')
lines(density(pred1$rt[is.finite(pred1$rt)]),col='black',lty='dashed')
lines(density(pred2$rt[is.finite(pred2$rt)]),col='green',lty='dashed')
lines(density(pred4$rt[is.finite(pred4$rt)]),col='purple',lty='dashed')
