calc_ll_R <- function(p_vector, model, dadm){
  if(!is.null(model$transform)){
    pars <- get_pars_matrix_oo(p_vector, dadm, model)
  } else{
    pars <- p_vector
  }
  ll <- model$log_likelihood(pars, dadm, model)
  return(ll)
}

pr_pt <- function(LT,UT,ps,dadm,model)
  # Probability of a response between LT and UT
{
  if (is.infinite(UT))
    out <- prod(1-model$pfun(rep(LT,nrow(ps)),ps)) else
      out <- c(prod(1-model$pfun(rep(LT,nrow(ps)),ps))-prod(1-model$pfun(rep(UT,nrow(ps)),ps)))
    if (!(is.nan(out) | is.na(out) | out==0)) return(1/out)
    idx <- 1:nrow(ps)
    pr <- 0
    for (i in 1:nrow(ps)) {
      pri <- my_integrate(f,lower=LT,upper=UT,p=ps[c(i,idx[-i]),],dfun=model$dfun,pfun=model$pfun)
      if (inherits(pri, "try-error") || suppressWarnings(is.nan(pri$value))) return(NA)
      pr <- pr + pri$value
    }

    1/pmax(0,pmin(pr,1))
    if (is.finite(pr)) return(pr) else return(NA)
}


pLU <- function(LT,LC,UC,UT,ps,dadm,model)
  # Probability from LT-LC + UC-UT
{
  pL <- my_integrate(f,lower=LT,upper=LC,p=ps,
                     dfun=model$dfun,pfun=model$pfun)
  if (inherits(pL,"try-error") || suppressWarnings(is.nan(pL$value))) return(NA)
  pU <- my_integrate(f,lower=UC,upper=UT,p=ps,
                     dfun=model$dfun,pfun=model$pfun)
  if (inherits(pU,"try-error") || suppressWarnings(is.nan(pU$value))) return(NA)
  pmax(0,pmin(pL$value,1))+pmax(0,pmin(pU$value,1))
}

log_surv_race <- function(t, ps, model, is_defective = FALSE) {
  # log S(t) for a race model, S(t) = prod_k (1 - F_k(t))
  if (is.infinite(t)) {
    if (is_defective) {
      return(sum(pnorm(0, ps[, "v"], ps[, "sv"], log.p = TRUE)))
    }
    return(-Inf)
  }
  p <- model$pfun(rep(t, nrow(ps)), ps)
  if (any(is.na(p)) || any(!is.finite(p))) return(NA_real_)
  p <- pmax(0, pmin(p, 1))
  if (any(p >= 1)) return(-Inf)
  sum(log1p(-p))
}

log_diff_exp_R <- function(log_hi, log_lo) {
  # Stable log(exp(log_hi) - exp(log_lo)) for log_hi >= log_lo.
  if (is.na(log_hi) || is.na(log_lo)) return(NA_real_)
  if (log_hi == -Inf) return(-Inf)
  if (log_lo == -Inf) return(log_hi)
  if (log_lo > log_hi) return(NA_real_)
  if (log_lo == log_hi) return(-Inf)
  log_hi + log1p(-exp(log_lo - log_hi))
}

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

my_integrate <- function(...,upper=Inf,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=upper,rel.tol=1e-5),silent=TRUE)
  if (!is(out,"try-error") && upper==Inf && out$subdivisions==1)
    out <- try(integrate(...,upper=big,rel.tol=1e-5),silent=TRUE)
  out
}

# doesn't have GNG branch?
log_likelihood_race_missing <- function(pars,dadm,model,min_ll=log(1e-10))
  # Race model summed log likelihood for models allowing missing values
{

  if (any(names(dadm)=="RACE")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

    is_gng <- !is.null(dadm$lR) && ("nogo" %in% levels(dadm$lR))
    nogo_idx <- if (is_gng) match("nogo", levels(dadm$lR)) else NA_integer_
    # Defective distribution: total response probability < 1 (e.g. LBAIO with
    # posdrift=FALSE). The intrinsic-omission mass at T=Inf must be added to the
    # no-response likelihood. Currently only LBAIO triggers this path.
    is_defective <- isTRUE(!is.null(model$c_name) && grepl("IO", model$c_name, fixed=TRUE))

    lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
    lds[dadm$winner] <- -Inf # Cases where winner is due only to contamination
    fix <- dadm$winner & is.finite(dadm$rt) & ok
    lds[fix] <- log(model$dfun(rt=dadm$rt[fix],pars=pars[fix,,drop=FALSE]))
    n_acc <- length(levels(dadm$R))
    fix <- !dadm$winner & is.finite(dadm$rt) & ok
    if (n_acc>1) lds[fix] <-
      log(1-model$pfun(rt=dadm$rt[fix],pars=pars[fix,,drop=FALSE]))
    lds[is.na(lds) | !ok] <- -Inf

    # Calculate truncation?
    LT <- dadm$LT
    UT <- dadm$UT
    dotrunc <- !(all(LT ==0) & all(is.infinite(UT)))

    # Calculate censoring
    LC <- dadm$LC
    UC <- dadm$UC

    # Response known
    # Fast
    nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R) & ok
    if ( any(nort) ) {
      LCs <- LC[nort & dadm$lR==levels(dadm$lR)[1]]; LTs <- LT[nort & dadm$lR==levels(dadm$lR)[1]]
      mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                     dimnames = list(NULL,NULL,colnames(pars)))
      winner <- matrix(dadm$winner[nort],nrow=n_acc)
      tofix <- dadm$winner & nort
      for (i in 1:dim(mpars)[2]) {
        if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
          pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=LTs[i],upper=LCs[i],p=pi,
                            dfun=model$dfun,pfun=model$pfun)
        if ( !inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))) )
          lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
      }
    }
    # Slow
    nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R) & ok
    if ( any(nort) ) {
      UCs <- UC[nort & dadm$lR==levels(dadm$lR)[1]]; UTs <- UT[nort & dadm$lR==levels(dadm$lR)[1]]
      mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                     dimnames = list(NULL,NULL,colnames(pars)))
      winner <- matrix(dadm$winner[nort],nrow=n_acc)
      tofix <- dadm$winner & nort
      for (i in 1:dim(mpars)[2]) if (UCs[i] < UTs[i]) {
        if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
          pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=UCs[i],upper=UTs[i],p=pi,
                            dfun=model$dfun,pfun=model$pfun)
        if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
          lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
      }
    }
    # No direction
    nort <- is.na(dadm$rt) & !is.na(dadm$R) & ok
    if ( any(nort) ) {
      LCs <- LC[nort & dadm$lR==levels(dadm$lR)[1]]; LTs <- LT[nort & dadm$lR==levels(dadm$lR)[1]]
      UCs <- UC[nort & dadm$lR==levels(dadm$lR)[1]]; UTs <- UT[nort & dadm$lR==levels(dadm$lR)[1]]
      mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                     dimnames = list(NULL,NULL,colnames(pars)))
      winner <- matrix(dadm$winner[nort],nrow=n_acc)
      tofix <- dadm$winner & nort
      for (i in 1:dim(mpars)[2]) {
        if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
          pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=LTs[i],upper=LCs[i],p=pi,
                            dfun=model$dfun,pfun=model$pfun)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
          p <- tmp$value
          if (dim(mpars)[[1]]==1) pi <- mpars[,i,,drop=FALSE] else
            pi <- mpars[,i,][order(!winner[,i]),]
          tmp <- my_integrate(f,lower=UCs[i],upper=UTs[i],p=pi,
                              dfun=model$dfun,pfun=model$pfun)
          if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
            p <- p + tmp$value else p <- 0
        } else p <- 0
        lds[tofix][i] <- log(pmax(0,pmin(p,1)))
      }
    }

    # Response unknown.
    # Fast
    nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R) & ok
    if ( any(nort) ) {
      LCs <- LC[nort & dadm$lR==levels(dadm$lR)[1]]; LTs <- LT[nort & dadm$lR==levels(dadm$lR)[1]]
      UCs <- UC[nort & dadm$lR==levels(dadm$lR)[1]]; UTs <- UT[nort & dadm$lR==levels(dadm$lR)[1]]
      mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                     dimnames = list(NULL,NULL,colnames(pars)))
      winner <- matrix(dadm$winner[nort],nrow=n_acc)
      tofixfast <- dadm$winner & nort
      for (i in 1:dim(mpars)[2]) {
        if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
          pi <- mpars[,i,]
        if (is_gng && !is.na(nogo_idx)) {
          # rt==-Inf codes an observed (fast) response: the nogo accumulator
          # cannot have won (nogo winning means inhibition, not a response).
          p <- 0
          for (j in 1:n_acc) {
            if (j == nogo_idx) next
            pj_pi <- mpars[,i,][c(j, c(1:n_acc)[-j]),]
            pc <- my_integrate(f, lower=LTs[i], upper=LCs[i], p=pj_pi,
                               dfun=model$dfun, pfun=model$pfun)
            if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
              p <- NA; break
            } else p <- pmax(0, pmin(pc$value, 1))
            # p_j <- pmax(0, pmin(pc$value, 1))
            # cf <- if (p_j != 0 && !(LTs[i]==0 & UTs[i]==Inf))
            #   pr_pt(LTs[i], UTs[i], mpars[,i,], dadm, model) else 1
            # if (is.na(cf)) { p <- NA; break }
            # p <- p + p_j * cf
          }
        } else {
          logP <- log_diff_exp_R(log_surv_race(LTs[i], mpars[,i,], model, is_defective),
                                 log_surv_race(LCs[i], mpars[,i,], model, is_defective))
          if (is.na(logP)) {
            p <- NA
          } else {
            if (!(LTs[i]==0 & UTs[i]==Inf))
              cf <- pr_pt(LTs[i],UTs[i],mpars[,i,],dadm,model) else cf <- 1
            if (!is.na(cf)) p <- exp(logP + log(cf)) else p <- NA
          }
        }
        lp <- log(p)
        if (!is.nan(lp) & !is.na(lp)) lds[tofixfast][i] <- lp else lds[tofixfast][i] <- -Inf
      }
    } else tofixfast <- NA
    # Slow
    nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R) & ok
    if ( any(nort) ) {
      LCs <- LC[nort & dadm$lR==levels(dadm$lR)[1]]; LTs <- LT[nort & dadm$lR==levels(dadm$lR)[1]]
      UCs <- UC[nort & dadm$lR==levels(dadm$lR)[1]]; UTs <- UT[nort & dadm$lR==levels(dadm$lR)[1]]
      mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                     dimnames = list(NULL,NULL,colnames(pars)))
      winner <- matrix(dadm$winner[nort],nrow=n_acc)
      tofixslow <- dadm$winner & nort
      for (i in 1:dim(mpars)[2]) {
        if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
          pi <- mpars[,i,]
        if (is_gng && !is.na(nogo_idx)) {
          pc <- my_integrate(f,lower=LTs[i],upper=UCs[i],p=pi[c(nogo_idx,c(1:n_acc)[-nogo_idx]),],
                             dfun=model$dfun,pfun=model$pfun)
          if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
            p <- NA
          } else p <- pmax(0,pmin(pc$value,1))
          # else {
          #   p <- pmax(0,pmin(pc$value,1))
          #   psurv <- prod(1 - model$pfun(rt=rep(UCs[i],n_acc),pars=pi))
          #   if (!is.na(psurv) & !is.nan(psurv)) p <- p + pmax(0,pmin(psurv,1)) else p <- NA
          # }
          # if (!is.na(p)) {
          #   if (p != 0 && !(LTs[i]==0 & UTs[i]==Inf))
          #     cf <- pr_pt(LTs[i],UTs[i],mpars[,i,],dadm,model) else cf <- 1
          #   if (!is.na(cf)) p <- p*cf
          # }
        } else {
          logP <- log_diff_exp_R(log_surv_race(UCs[i], mpars[,i,], model, is_defective),
                                 log_surv_race(UTs[i], mpars[,i,], model, is_defective))
          if (is.na(logP)) {
            p <- NA
          } else {
            if (!(LTs[i]==0 & UTs[i]==Inf))
              cf <- pr_pt(LTs[i],UTs[i],mpars[,i,],dadm,model) else cf <- 1
            if (!is.na(cf)) p <- exp(logP + log(cf)) else p <- NA
          }
          # Defective distributions (e.g. LBAIO, posdrift=FALSE): add probability
          # mass for intrinsic omissions (all accumulators have negative drift, T=Inf).
          # pr_pt already includes this mass in its denominator via 1-pfun(LT) for
          # each accumulator; we must add it to the numerator here.
          if (is_defective && !is.na(p)) {
            p_IO <- prod(pnorm(0, pi[, "v"], pi[, "sv"]))
            if (is.finite(p_IO) && p_IO > 0) {
              cf_d <- if (!(LTs[i]==0 & UTs[i]==Inf))
                pr_pt(LTs[i], UTs[i], mpars[,i,], dadm, model) else 1
              if (!is.na(cf_d)) p <- p + p_IO * cf_d
            }
          }
        }
        lp <- log(p)
        if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
      }
    } else tofixslow <- NA
    # no direction
    nort <- is.na(dadm$rt) & is.na(dadm$R) & ok
    # nort <- nort & (pars[,"pContaminant"] == 0) # Otherwise not identifiable
    if ( any(nort) ) {
      LCs <- LC[nort & dadm$lR==levels(dadm$lR)[1]]; LTs <- LT[nort & dadm$lR==levels(dadm$lR)[1]]
      UCs <- UC[nort & dadm$lR==levels(dadm$lR)[1]]; UTs <- UT[nort & dadm$lR==levels(dadm$lR)[1]]
      mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                     dimnames = list(NULL,NULL,colnames(pars)))
      winner <- matrix(dadm$winner[nort],nrow=n_acc)
      tofix <- dadm$winner & nort
      for (i in 1:dim(mpars)[2]) {
        if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
          pi <- mpars[,i,]
        pc <- pLU(LTs[i],LCs[i],UCs[i],UTs[i],pi,dadm,model)
        if (is.na(pc)) p <- NA else {
          if (pc!=0 & !(LTs[i]==0 & UTs[i]==Inf))
            cf <- pr_pt(LTs[i],UTs[i],pi,dadm,model) else cf <- 1
            if (!is.na(cf)) p <- pc*cf else p <- NA
            if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
              pc <- pLU(LTs[i],LCs[i],UCs[i],UTs[i],mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model)
              if (is.na(pc)) {
                p <- NA; break
              }
              if (pc!=0 & !(LTs[i]==0 & UTs[i]==Inf))
                cf <- pr_pt(LTs[i],UTs[i],mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model) else cf <- 1
                if (is.na(cf)) {
                  p <- NA; break
                } else p <- p +  pc*cf
            }
        }
        lp <- log(pmax(0,pmin(p,1)))
        if (!is.nan(lp) & !is.na(lp)) lds[tofix][i] <- lp
      }
    }


    # Truncation
    if ( dotrunc ) {
      okT <- attr(dadm,"unique_nortR") & dadm$winner & ok
      if (any(okT)) {
        LTs <- LT[okT]; UTs <- UT[okT]
        tpars <- pars[attr(dadm,"unique_nortR") & ok,,drop=FALSE]
        tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
                       dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
        winner <- matrix(dadm$winner[attr(dadm,"unique_nortR") & ok],nrow=n_acc)[,,drop=FALSE]
        cf <- rep(NA,length(LTs))
        for (i in 1:length(LTs)) {
          if (dim(tpars)[[1]]==1) pi <- t(as.matrix(tpars[,i,])) else
            pi <- tpars[,i,][order(!winner[,i]),]
          cf[i] <- pr_pt(LTs[i],UTs[i],pi,dadm,model)
        }
        cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nortR")]
        fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf) & !is.na(dadm$R) & ok
        if (any(fix)) {
          lds[fix] <- lds[fix] + cf[fix]
        }
      }
    }


    if (n_acc>1) {
      ll <- lds[dadm$winner]
      if (n_acc==2) {
        ll <- ll + lds[!dadm$winner]
      } else {
        # Optimization: Replace explicit matrix sum apply(..., 2, sum) with colSums
        ll <- ll + colSums(matrix(lds[!dadm$winner],nrow=n_acc-1))
      }
    } else ll <- lds
    ll[is.na(ll) | is.nan(ll)] <- -Inf

    #  llR <<- ll

    # Non-process (contaminant) miss.
    ispContaminant <- pars[dadm$winner & ok,"pContaminant"]>0
    if ( any(ispContaminant) ) {
      use <- dadm$winner & ok
      p <- exp(ll[ispContaminant])
      pc <- pars[use,"pContaminant"]
      isMiss <- c(dadm$rt==Inf)[use]
      p[isMiss] <- pc[isMiss] + (1-pc[isMiss])*p[isMiss]
      p[!isMiss] <- (1-pc[!isMiss])*p[!isMiss]
      ll[ispContaminant] <- log(p)
    }
    ll <- pmax(min_ll,ll)
    return(sum(ll[attr(dadm,"expand")]))
}


log_likelihood_ddm <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood, with protection against numerical issues
{
  ok <- attr(pars,"ok")
  like <- numeric(dim(dadm)[1])
  if (any(ok))
    like[ok] <- model$dfun(dadm$rt[ok],dadm$R[ok],
                                        pars[ok,,drop=FALSE])
  like[ok][is.na(like[ok])] <- 0

  # Truncation: renormalise each retained trial's density by the probability
  # mass falling inside the truncation window,
  #   P(LT < RT < UT) = sum over both responses of ( F_R(UT) - F_R(LT) ),
  # using the defective DDM cdf (model$pfun). Without this the DDM/WDM
  # likelihood is the un-truncated density and SBC is miscalibrated.
  LT <- dadm$LT; UT <- dadm$UT
  if (!is.null(LT) && !is.null(UT) && !(all(LT==0) & all(is.infinite(UT)))) {
    Rlevs <- levels(dadm$R)
    Rlo <- factor(rep(Rlevs[1],nrow(pars)),levels=Rlevs)
    Rhi <- factor(rep(Rlevs[2],nrow(pars)),levels=Rlevs)
    Fhi <- rep(1,nrow(pars))            # F(UT)=1 when UT=Inf (mass sums to 1)
    sel <- ok & is.finite(UT)
    if (any(sel)) Fhi[sel] <-
      model$pfun(UT[sel],Rlo[sel],pars[sel,,drop=FALSE]) +
      model$pfun(UT[sel],Rhi[sel],pars[sel,,drop=FALSE])
    Flo <- numeric(nrow(pars))          # F(LT)=0 when LT=0
    sel <- ok & (LT>0)
    if (any(sel)) Flo[sel] <-
      model$pfun(LT[sel],Rlo[sel],pars[sel,,drop=FALSE]) +
      model$pfun(LT[sel],Rhi[sel],pars[sel,,drop=FALSE])
    Pin <- Fhi - Flo
    fix <- ok & like>0 & is.finite(Pin) & Pin>0
    like[fix] <- like[fix]/Pin[fix]
  }
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}

log_likelihood_ddmgng <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood for go/nogo model
{
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok"))) {
    isna <- is.na(dadm$rt)|is.infinite(dadm$rt) # allow Inf like other censored models
    ok <- attr(pars,"ok") & !isna
    like[ok] <- model$dfun(dadm$rt[ok],dadm$R[ok],pars[ok,,drop=FALSE])
    ok <- attr(pars,"ok") & isna
    like[ok] <- # dont terminate on go boundary before timeout
      pmax(0,pmin(1,(1-model$pfun(dadm$TIMEOUT[ok],dadm$Rgo[ok],pars[ok,,drop=FALSE]))))

  }
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}



#### sdt choice likelihoods ----

log_likelihood_sdt <- function(pars,dadm, model,lb=-Inf, min_ll=log(1e-10))
  # probability of ordered discrete choices based on integrals of a continuous
  # distribution between thresholds, with fixed lower bound for first response
  # lb. Upper bound for last response is a fixed value in threshold vector
{
  first <- dadm$lR==levels(dadm$lR)[1]
  last <- dadm$lR==levels(dadm$lR)[length(levels(dadm$lR))]
  pars[last,"threshold"] <- Inf
  # upper threshold
  ut <- pars[dadm$winner,"threshold"]
  # lower threshold fixed at lb for first response
  pars[first &  dadm$winner,"threshold"] <- lb
  # otherwise threshold of response before one made
  notfirst <- !first &  dadm$winner
  pars[notfirst,"threshold"] <- pars[which(notfirst)-1,"threshold"]
  lt <- pars[dadm$winner,"threshold"]
  # log probability
  ll <- numeric(sum(dadm$winner))
  if (!is.null(attr(pars,"ok"))) { # Bad parameter region
    ok <- attr(pars,"ok")
    okw <- ok[dadm$winner]
    ll[ok] <- log(model$pfun(lt=lt[okw],ut=ut[okw],pars=pars[dadm$winner & ok,,drop=FALSE]))
  } else ll <- log(model$pfun(lt=lt,ut=ut,pars=pars[dadm$winner,,drop=FALSE]))
  ll <- ll[attr(dadm,"expand")]
  ll[is.na(ll)] <- 0
  sum(pmax(min_ll,ll))
}

# Two options:
# 1) component = NULL in which case we do all likelihoods in one block
# 2) component = integer, in which case we are blocking the ll and only want that one
log_likelihood_joint <- function(proposals, dadms, model_list, component = NULL){
  parPreFixs <- unique(gsub("[|].*", "", colnames(proposals)))
  i <- 0
  k <- 0
  total_ll <- 0
  for(dadm in dadms){
    i <- i + 1
    if(is.null(component) || component == i){
      if(is.data.frame(dadm)){
        k <- k + 1
        # Sometimes designs are the same across models (typically in fMRI)
        # Instead of storing multiple, we just store a pointer to the first original design
        if(is.numeric(attr(dadm, "designs"))){
          ref_idx <- attr(dadm, "designs")
          attr(dadm, "designs") <- attr(dadms[[ref_idx]], "designs")
        }
        parPrefix <- parPreFixs[k]
        # Optimization: Use startsWith to match parameter prefixes
        columns_to_use <- startsWith(colnames(proposals), paste0(parPrefix, "|"))
        currentPars <- proposals[,columns_to_use, drop = F]
        colnames(currentPars) <- gsub(".*[|]", "", colnames(currentPars))
        total_ll <- total_ll +  calc_ll_manager(currentPars, dadm, model_list[[i]])
      }
    }
  }
  return(total_ll)
}
