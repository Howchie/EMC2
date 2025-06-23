# Ensure EMC2 package is loaded (or relevant functions are sourced)
#remotes::install_github("https://github.com/Howchie/EMC2/",ref="feat/censored-truncated-race-likelihood")
library(EMC2) # If running as part of package tests or after installation
# Alternatively, if developing, source the necessary files:
# source("R/likelihood.R") # Contains log_likelihood_race, log_likelihood_race_cens_trunc
# source("R/utils.R")      # For any utilities they might call
# source("R/models.R") # For LBA() if not loaded
# Rcpp::sourceCpp("src/particle_ll.cpp") # To get the C++ test wrapper
set.seed(123)
# --- 1. Define Model (Actual LBA from EMC2) ---
# Get the LBA model definition
# This model list will contain the R versions of dfun and pfun (EMC2::dlba, EMC2::plba)
lba_model_R <- EMC2::LBA() 
# For C++ tests, we'll use "LBA_test" string, which test_c_loglik_cens_trunc_wrapper_R
# should map to lba_dfun_adapter and lba_pfun_adapter using dlba_c and plba_c.

# --- 2. Define Test Parameters ---
min_ll_test <- log(1e-10)
n_acc <- 2 # Number of accumulators

data = data.frame(rt=c(1,1.5),R=factor(c("left","right")),S=factor(c("left","right")),subjects=factor(c(1,1)))
matchfun <- function(d)as.character(d$S) == tolower(as.character(d$lR))
design <- design(data=data,model=LBA,matchfun = matchfun,
                 formula=list(v~1,sv~1,B~1,A~1,t0~1))
p_vector = sampled_pars(design)
p_vector["v"] = 2.5
p_vector["B"] = log(1)
p_vector["A"] = log(1)
p_vector["sv"] = log(1)
p_vector["t0"] = log(0.1)
sim_data = make_data(p_vector,design=design, n_trials=1)
sim_data$rt[2]=Inf
attr(sim_data,"UC")=2.5
fit_data = make_emc(sim_data,design,type="single")
# Hardcoded censoring/truncation values (matching those in the likelihood functions for now)
LT_test <- 0.1 
LC_test <- 0.2
UC_test <- 2.5 
UT_test <- 3

# --- 3. Create Dummy 'dadm' (Unique Trial Conditions) ---
# Each row is a unique trial type.
# Let's create a few scenarios:
#   1. Observed RT, no censoring/truncation effectively active
#   2. Observed RT, but would be affected by LT/UT if they were tighter
#   3. Lower censored RT (fast response, R known)
#   4. Upper censored RT (slow response, R known)
#   5. Interval censored (NA RT, R known)
#   6. Lower censored RT (fast response, R unknown)
#   7. NA RT (interval censored, R unknown)

dadm_test_list <- list(
  # Scenario 1: Standard observed RT
  data.frame(rt = 0.5, R = factor(1, levels=1:n_acc), lR = factor(1, levels=1:n_acc)),
  # Scenario 2: Another observed RT
  data.frame(rt = 1.5, R = factor(2, levels=1:n_acc), lR = factor(1, levels=1:n_acc)),
  # Scenario 3: Lower censored, Response 1 known (e.g. rt set to -Inf by experiment)
  data.frame(rt = -Inf, R = factor(1, levels=1:n_acc), lR = factor(1, levels=1:n_acc)),
  # Scenario 4: Upper censored, Response 2 known
  data.frame(rt = Inf, R = factor(2, levels=1:n_acc), lR = factor(1, levels=1:n_acc)),
  # Scenario 5: Interval censored (NA RT), Response 1 known
  data.frame(rt = NA_real_, R = factor(1, levels=1:n_acc), lR = factor(1, levels=1:n_acc)),
  # Scenario 6: Lower censored, Response unknown
  data.frame(rt = -Inf, R = factor(NA, levels=1:n_acc), lR = factor(1, levels=1:n_acc)),
  # Scenario 7: Interval censored (NA RT), Response unknown
  data.frame(rt = NA_real_, R = factor(NA, levels=1:n_acc), lR = factor(1, levels=1:n_acc))
)
dadm_test <- fit_data[[1]]$data$`1`
n_unique_trials <- nrow(dadm_test)/n_acc

# Define 'expand' attribute (each unique trial happens once for this test)
attr(dadm_test, "expand") <- 1:n_unique_trials
# Or, if you want to test expansion, make some unique trials repeat:
# attr(dadm_test, "expand") <- c(1,1,2,3,3,3,4,5,6,7) # if dadm_test had fewer rows initially

# --- 4. Create Dummy 'pars' Matrix ---
# nrow(pars) = n_unique_trials * n_acc
# Columns should match what the LBA model expects on its *estimation scale*.
# LBA()$p_types gives: v, A, b, t0, sv. sv is often fixed.
# Let's assume p_names are v, A, b, t0 (transformed values)
# For LBA, 'b' is typically B-A (on natural scale), so log(B-A) if B-A > 0.
# 't0' is often log(t0_natural). 'A' might be log(A_natural). 'v' is often natural scale.

# Get p_names from a sample LBA design to be sure about transformations
# This is a bit circular if we don't have a full design object, but for LBA they are standard.
# Typical estimation scale: v (natural), A (log), b (log of B-A or similar), t0 (log)
p_names_lba <- c("v","sv","A", "b", "t0", "B") 
n_params_per_acc <- length(p_names_lba)

pars_test <- matrix(0.0, nrow = n_unique_trials * n_acc, ncol = n_params_per_acc)
colnames(pars_test) <- p_names_lba

# Fill with some plausible (but not necessarily realistic) parameter values
# These are on the ESTIMATION scale.
# Example: For LBA, v=drift, A=upper bound for start point noise (log scale), 
# b=threshold (log scale, e.g. log(B-A)), t0=non-decision time (log scale)
base_pars_acc1 <- c(v=2.5, sv=1,  A=0.5, b=1.0, t0=0.2, B=0.5+0.5)
base_pars_acc2 <- c(v=2.0, sv=1,  A=0.4, b=0.8, t0=0.15, B=0.4+0.8)

for (i in 1:n_unique_trials) {
  row_idx_acc1 <- (i - 1) * n_acc + 1
  row_idx_acc2 <- (i - 1) * n_acc + 2
  
  # Slightly vary pars per condition for more diverse testing if desired
  pars_test[row_idx_acc1, ] <- base_pars_acc1 + c( (i-1)*0.1, 0, 0, 0, 0, 0) 
  pars_test[row_idx_acc2, ] <- base_pars_acc2 - c( (i-1)*0.05, 0, 0, 0, 0, 0)
}

# attr(pars, "ok"): all parameters are initially ok
attr(pars_test, "ok") <- rep(TRUE, nrow(pars_test))
# Example: make one parameter set for one accumulator in one trial invalid
# attr(pars_test, "ok")[3] <- FALSE 

# --- 5. Define your original log_likelihood_race_missing ---
# For comparison, you'll need this function definition available.
# [Copy your original log_likelihood_race_missing function here]
log_likelihood_race_missing <- function(pars,dadm,model,min_ll=log(1e-10))
  # Race model summed log likelihood
{
  if (any(names(dadm)=="RACE")){ # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA
  }
  if (is.null(attr(pars,"ok"))) {
    
    f <- function(t,p,dfun,pfun) {
      # Called by integrate to get race density for vector of times t given
      # matrix of parameters where first row is the winner.
      
      out <- dfun(t,
                  matrix(rep(p[1,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]])))
      if (dim(p)[1]>1) for (i in 2:dim(p)[1])
        out <- out*(1-pfun(t,
                           matrix(rep(p[i,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]]))))
      out
    }
    
    pr_pt <- function(LT,UT,ps,dadm)
      # p(untrucated response)/p(truncated response), > 1, multiplicative truncation correction
    {
      pr <- try(integrate(f,lower=0,upper=Inf,p=ps,
                          dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
      if (pr$value==0) return(0)
      pt <- try(integrate(f,lower=LT,upper=UT,p=ps,
                          dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
      out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
      if (is.infinite(out)) return(NA)
      out
    }
    
    pLU <- function(LT,LC,UC,UT,ps,dadm)
      # Probability from LT-LC + UC-UT
    {
      pL <- try(integrate(f,lower=LT,upper=LC,p=ps,
                          dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if (inherits(pL,"try-error") || suppressWarnings(is.nan(pL$value))) return(NA)
      pU <- try(integrate(f,lower=UC,upper=UT,p=ps,
                          dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if (inherits(pU,"try-error") || suppressWarnings(is.nan(pU$value))) return(NA)
      pmax(0,pmin(pL$value,1))+pmax(0,pmin(pU$value,1))
    }
    
    ok <- !logical(dim(pars)[1])
  } else{ ok <- attr(pars,"ok") }
  
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(model$dfun(rt=dadm$rt[dadm$winner],
                                     pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <- log(1-model$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- min_ll
  
  # Calculate truncation?
  LT <- attr(dadm,"LT")
  UT <- attr(dadm,"UT")
  dotrunc <- (!is.null(LT) | !is.null(UT))
  if (is.null(LT)) LT <- 0
  if (is.null(UT)) UT <- Inf
  
  # Calculate censoring
  LC <- attr(dadm,"LC")
  UC <- attr(dadm,"UC")
  
  # Response known
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      tmp <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,][order(!winner[,i]),],
                           dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if ( !inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))) )
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      tmp <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,][order(!winner[,i]),],
                           dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # No direction
  nort <- is.na(dadm$rt) & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      tmp <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,][order(!winner[,i]),],
                           dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
        p <- tmp$value
        tmp <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,][order(!winner[,i]),],
                             dfun=model$dfun,pfun=model$pfun),silent=TRUE)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
          p <- p + tmp$value else p <- 0
      } else p <- 0
      lds[tofix][i] <- log(pmax(0,pmin(p,1)))
    }
  }
  
  # Response unknown
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixfast <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      pc <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,],
                          dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
        if (!is.na(p)) {
          if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
          if (!is.na(cf)) p <- p*cf
        }
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                              dfun=model$dfun,pfun=model$pfun),silent=TRUE)
          if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
            p <- NA; break
          }
          if (pc$value != 0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
            if (!is.na(cf)) p <- p + pc$value*cf
        }
        lp <- log(p)
        if (!is.nan(lp) & !is.na(lp)) lds[tofixfast][i] <- lp else lds[tofixfast][i] <- -Inf
    }
  } else tofixfast <- NA
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixslow <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      pc <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,],
                          dfun=model$dfun,pfun=model$pfun),silent=TRUE)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
        if (!is.na(p)) {
          if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
          if (!is.na(cf)) p <- p*cf
        }
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                              dfun=model$dfun,pfun=model$pfun),silent=TRUE)
          if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
            p <- NA; break
          }
          if (pc$value != 0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
            if (!is.na(cf)) p <- p + pc$value*cf
        }
        lp <- log(p)
        if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
    }
  } else tofixslow <- NA
  # no direction
  nort <- is.na(dadm$rt) & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      pc <- pLU(LT,LC,UC,UT,mpars[,i,],dadm)
      if (is.na(pc)) p <- NA else {
        if (pc!=0 & !(LT==0 & UT==Inf)) cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- pc*cf else p <- NA
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- pLU(LT,LC,UC,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm)
          if (is.na(pc)) {
            p <- NA; break
          }
          if (pc!=0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
            if (is.na(cf)) {
              p <- NA; break
            } else p <- p +  pc*cf
        }
      }
      lp <- log(pmax(0,pmin(p,1)))
      if (!is.nan(lp) & !is.na(lp)) lds[tofix][i] <- lp
    }
  }
  
  
  # Truncation where not censored or censored and response known
  ok <- is.finite(lds[attr(dadm,"unique_nort") & dadm$winner])
  alreadyfixed <- is.na(dadm$R[attr(dadm,"unique_nort") & dadm$winner])
  ok <- ok & !alreadyfixed
  if ( dotrunc & any(ok) ) {
    tpars <- pars[attr(dadm,"unique_nort"),]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
                   dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nort")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(ok))
    for (i in 1:length(ok)) if (ok[i]) {
      cf[i] <- pr_pt(LT,UT,tpars[,i,][order(!winner[,i]),],dadm)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nort")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf)
    if (any(fix)) lds[fix] <- lds[fix] + cf[fix]
    badfix <- dadm$winner & (is.na(cf) | is.nan(cf) | is.infinite(cf))
    if (!all(is.na(tofixfast))) badfix <- badfix & !tofixfast
    if (!all(is.na(tofixslow))) badfix <- badfix & !tofixslow
    if (any(badfix)) lds[badfix] <- min_ll
  }
  
  lds[is.na(lds) | !ok] <- min_ll
  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
    ll[is.na(ll)] <- min_ll
    return(sum(pmax(min_ll,ll[attr(dadm,"expand")])))
  } else return(sum(pmax(min_ll,lds[attr(dadm,"expand")])))
}


# --- 6. Perform Calls and Comparisons ---
cat("Running Likelihood Tests:

")

results <- list()

for (i in 1:n_unique_trials) {
  cat("--- Testing Unique Trial Condition:", i, "---
")
  indices <- ((i - 1) * n_acc + 1):(i * n_acc)
  current_dadm <- dadm_test[indices, , drop = FALSE]
  attr(current_dadm, "expand") <- 1 # This trial expands to itself once
  
  
  current_pars <- pars_test[indices, , drop = FALSE]
  attr(current_pars, "ok") <- attr(pars_test, "ok")[par_indices]
  
  cat("dadm:
"); print(current_dadm)
  cat("pars:
"); print(current_pars)
  
  # A. Standard EMC2::log_likelihood_race (only for observed, non-censored/truncated)

  # B. Your original log_likelihood_race_missing
  ll_original_missing <- NA
  if (exists("log_likelihood_race_missing")) {
    # Ensure LT, LC, UC, UT are set as attributes on current_dadm if your function expects them
    attr(current_dadm, "LT") <- LT_test
    attr(current_dadm, "LC") <- LC_test
    attr(current_dadm, "UC") <- UC_test
    attr(current_dadm, "UT") <- UT_test
    # Your original function might have specific needs for 'model' object
    ll_original_missing <- log_likelihood_race_missing(
      pars = current_pars,
      dadm = current_dadm,
      model = lba_model_R, # Or whatever model structure your original func expects
      min_ll = min_ll_test
    )
  }
  cat("Original Censored LL (To be filled by user):", ll_original_missing, "
")
  
  # C. New R version: log_likelihood_race_cens_trunc
  # Make sure the actual function definition is available
  ll_old_R <- EMC2:::log_likelihood_race(
    pars = current_pars,
    dadm = current_dadm,
    model = EMC2::LBA(), 
    min_ll = min_ll_test
  )
  ll_new_R <- EMC2:::log_likelihood_race_cens_trunc(
    pars = current_pars,
    dadm = current_dadm,
    model = lba_model_R, 
    min_ll = min_ll_test
  )
  cat("New R Censored LL:", ll_new_R, "
")
  
  # D. New C++ version (via Rcpp wrapper)
  ll_new_Cpp <- NA
  # if (exists("test_c_loglik_cens_trunc_wrapper_R")) {
  #   # The C++ wrapper takes the full 'dadm' and 'pars' for *all unique trials*
  #   # and then iterates. So, we should call it once with full dadm_test, pars_test.
  #   # However, for per-trial comparison, it's easier if it could also take single trial data.
  #   # The current C++ wrapper returns a vector of ll_unique. Let's test it that way.
  #   # This call is to the C++ code, which has its own hardcoded LT/LC/UC/UT for now.
  cpp_output_vector <- EMC2:::test_c_loglik_cens_trunc_wrapper_R(
    pars = current_pars, # Full pars matrix for all unique trials
    dadm = current_dadm,  # Full dadm for all unique trials
    model_type_str = "LBA_test", 
    min_ll = min_ll_test,
    ok_params = attr(current_pars, "ok"),
    n_acc = n_acc
  )
  ll_new_Cpp <- cpp_output_vector # Get the result for the j-th unique trial
  # }
  cat("New C++ Censored LL (via wrapper):", ll_new_Cpp, "

")
  
  results[[paste0("trial_", i)]] <- list(
    ll_standard_race = ll_old_R,
    ll_original_missing = ll_original_missing,
    ll_new_R = ll_new_R,
    ll_new_Cpp = ll_new_Cpp
  )
}

# --- Comparison for Scenario 1 (Observed RT, no effective censoring/truncation) ---
# For this, we'd call log_likelihood_race_cens_trunc with LT=0, UT=Inf, LC=0, UC=Inf
# (or modify the hardcoded values temporarily in that function for this test run)
# and compare its output for dadm_test[1,] with EMC2::log_likelihood_race.

cat("--- Testing Scenario 1 (Observed RT) against standard log_likelihood_race ---
")
dadm_scen1 <- dadm_test[1, , drop = FALSE]
attr(dadm_scen1, "expand") <- 1
pars_scen1 <- pars_test[1:n_acc, , drop = FALSE]
attr(pars_scen1, "ok") <- attr(pars_test, "ok")[1:n_acc]

# To make this comparison, log_likelihood_race_cens_trunc needs non-restrictive bounds
# Easiest is to temporarily modify its hardcoded LT/LC/UC/UT for this specific call/test block
# For example:
# assignInNamespace("LT", 0, ns="EMC2", envir=environment(EMC2:::log_likelihood_race_cens_trunc))
# ... and so on for LC, UC, UT, then reset them. This is hacky.
# A better way: pass LT/LC/UC/UT as arguments to log_likelihood_race_cens_trunc.
# Since that's a future step, this comparison is harder with current hardcoding.

# Call standard log_likelihood_race
# It needs dadm where each row is an accumulator line, and pars matching that.
# It also needs dadm$winner (logical vector, TRUE for the winning accumulator row in pars)
dadm_for_std_race <- data.frame(
  rt = rep(dadm_scen1$rt, n_acc),
  R = rep(dadm_scen1$R, n_acc), # Factor
  lR = rep(dadm_scen1$lR, n_acc) # Factor
  # Any other columns needed by model$dfun/pfun if they use dadm directly
)
attr(dadm_for_std_race, "expand") <- 1:n_acc # Each accumulator line expands to itself
winner_std_race <- rep(FALSE, n_acc)
if(!is.na(dadm_scen1$R)) winner_std_race[as.integer(dadm_scen1$R)] <- TRUE
dadm_for_std_race$winner <- winner_std_race


ll_std_race_val <- EMC2:::log_likelihood_race(
  pars = pars_scen1, 
  dadm = dadm_for_std_race, 
  model = lba_model_R, 
  min_ll = min_ll_test
)
cat("Standard log_likelihood_race for Scen1:", ll_std_race_val, "")
cat("New R (with non-restrictive bounds - manual change needed) for Scen1:", 
    results$trial_1$ll_new_R, "(ensure bounds were 0,Inf,0,Inf)
")

print(results)
