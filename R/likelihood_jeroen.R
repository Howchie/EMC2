#### Missing ----

my_integrate <- function(...,upper=Inf,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=upper),silent=TRUE)
  if (!is(out,"try-error") && upper==Inf && out$subdivisions==1)
    out <- try(integrate(...,upper=big),silent=TRUE)
  out
}


# pars <- cbind(a=c(1,1),v=c(1,1),t0=c(.3,.3),z=c(.5,.5),d=c(0,0),sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
# LC <- .4
# UC=2
#
# EMC2:::pDDM(c(LC,LC),"lower",pars)
# EMC2:::pDDM(UC,"lower",pars)
#
# EMC2:::pDDM(LC,"upper",pars)
# EMC2:::pDDM(UC,"upper",pars)


log_likelihood_ddm_missing <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood, with protection against numerical issues
{

  pr_pt <- function(LT,UT,R,p,model)
  # p(untruncated response)/p(truncated response), > 1, multiplicative truncation correction
  {
    pr <- model$pfun(rep(Inf,length(R)),R,p)
    out <- rep(1,length(R))
    if (!any(is.na(pr))) {
      ok <- pr>0
      if (any(ok)) {
        pt <- model$pfun(UT[ok],R[ok],p[ok,,drop=FALSE]) -
            model$pfun(LT[ok],R[ok],p[ok,,drop=FALSE])
        if (!any(is.na(pt))) out[ok] <- pr[ok]/pt
      }
    }
    out[is.na(out) | is.nan(out) | !is.finite(out) | out < 1 ] <- 1
    out
  }

  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok"))) {
    rt <- dadm$rt[attr(pars,"ok")]
    R <- dadm$R[attr(pars,"ok")]
    p <- pars[attr(pars,"ok"),,drop=FALSE]

    # Calculate truncation?
    LT <- attr(dadm,"LT")
    UT <- attr(dadm,"UT")
    dotrunc <- (!is.null(LT) | !is.null(UT))
    if (is.null(LT)) LT <- 0
    if (is.null(UT)) UT <- Inf

    # Calculate censoring
    LC <- attr(dadm,"LC")
    UC <- attr(dadm,"UC")

    likeok <- rep(NA,sum(attr(pars,"ok")))

    # Response known
    # Fast
    nort <- rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,model$pfun(rep(LC,sum(nort)),R[nort],p[nort,,drop=FALSE]))
    }
    # Slow
    nort <- rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,model$pfun(rep(Inf,sum(nort)),R[nort],p[nort,,drop=FALSE])-
                      model$pfun(rep(UC,sum(nort)),R[nort],p[nort,,drop=FALSE]))
    }
    # No direction
    nort <- is.na(rt) & !is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,model$pfun(rep(LC,sum(nort)),R[nort],p[nort,,drop=FALSE]) +
                      (model$pfun(rep(Inf,sum(nort)),R[nort],p[nort,,drop=FALSE])-
                         model$pfun(rep(UC,sum(nort)),R[nort],p[nort,,drop=FALSE])))
    }

    # Response unknown.
    # Fast
    nort <- rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,model$pfun(rep(LC,sum(nort)),"lower",p[nort,,drop=FALSE]) +
                      model$pfun(rep(LC,sum(nort)),"upper",p[nort,,drop=FALSE]))
    }
    # Slow
    nort <- rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(R)
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,(model$pfun(rep(Inf,sum(nort)),"lower",p[nort,,drop=FALSE])-
                       model$pfun(rep(UC,sum(nort)),"lower",p[nort,,drop=FALSE])) +
                      (model$pfun(rep(Inf,sum(nort)),"upper",p[nort,,drop=FALSE])-
                       model$pfun(rep(UC,sum(nort)),"upper",p[nort,,drop=FALSE])))
    }
    # no direction
    nort <- is.na(rt) & is.na(R)
    likeok[nort] <- 0
    nort <- nort & (p[,"pContaminant"] == 0) # Otherwise not identifiable
    if ( any(nort) ) {
      likeok[nort] <- pmax(0,
        model$pfun(rep(LC,sum(nort)),"lower",p[nort,,drop=FALSE]) +
        model$pfun(rep(LC,sum(nort)),"upper",p[nort,,drop=FALSE]) +
        (model$pfun(rep(Inf,sum(nort)),"lower",p[nort,,drop=FALSE])-
         model$pfun(rep(UC,sum(nort)),"lower",p[nort,,drop=FALSE])) +
        (model$pfun(rep(Inf,sum(nort)),"upper",p[nort,,drop=FALSE])-
         model$pfun(rep(UC,sum(nort)),"upper",p[nort,,drop=FALSE])))
    }

    # Truncation where not censored or censored and response known
    ok <- is.na(likeok) & !is.na(R)
    mult <- rep(1,length(likeok))
    if ( dotrunc & any(ok) )
      mult[ok] <- pr_pt(rep(LT,sum(ok)),rep(UT,sum(ok)),R[ok],p[ok,,drop=FALSE],model)

    # Usual non-missing update x truncation ratio
    ok <- is.na(likeok)
    likeok[ok] <- mult[ok]*model$dfun(rt[ok],R[ok],p[ok,,drop=FALSE])

    # Non-process (contaminant) miss.
    ispContaminant <- p[,"pContaminant"]>0
    if ( any(ispContaminant) ) {
      pc <- p[,"pContaminant"]
      isMiss <- is.na(R)
      likeok[isMiss] <- pc[isMiss] + (1-pc[isMiss])*likeok[isMiss]
      likeok[!isMiss] <- (1-pc[!isMiss])*likeok[!isMiss]
    }

    like[attr(pars,"ok")] <- likeok
  }

  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
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


log_likelihood_race_missing <- function(pars,dadm,model,min_ll=log(1e-10))
  # Race model summed log likelihood for models allowing missing values
{

  if (any(names(dadm)=="RACE")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(model$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <-
    log(1-model$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
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
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
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
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
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
  nort <- is.na(dadm$rt) & !is.na(dadm$R)
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
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
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
      pc <- my_integrate(f,lower=LTs[i],upper=LCs[i],p=pi,
                         dfun=model$dfun,pfun=model$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LTs[i]==0 & UTs[i]==Inf))
          cf <- pr_pt(LTs[i],UTs[i],mpars[,i,],dadm,model) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=LTs[i],upper=LCs[i],p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=model$dfun,pfun=model$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LTs[i]==0 & UTs[i]==Inf))
          cf <- pr_pt(LTs[i],UTs[i],mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model) else cf <- 1
          if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixfast][i] <- lp else lds[tofixfast][i] <- -Inf
    }
  } else tofixfast <- NA
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    LCs <- LC[nort & dadm$lR==levels(dadm$lR)[1]]; LTs <- LT[nort & dadm$lR==levels(dadm$lR)[1]]
    UCs <- UC[nort & dadm$lR==levels(dadm$lR)[1]]; UTs <- UT[nort & dadm$lR==levels(dadm$lR)[1]]
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixslow <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) if (UCs[i] < UTs[i]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- my_integrate(f,lower=UCs[i],upper=UTs[i],p=pi,
                         dfun=model$dfun,pfun=model$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LTs[i]==0 & UTs[i]==Inf))
          cf <- pr_pt(LTs[i],UTs[i],mpars[,i,],dadm,model) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=UCs[i],upper=UTs[i],p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=model$dfun,pfun=model$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LTs[i]==0 & UTs[i]==Inf))
          cf <- pr_pt(LTs[i],UTs[i],mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model) else cf <- 1
          if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
    }
  } else tofixslow <- NA
  # no direction
  nort <- is.na(dadm$rt) & is.na(dadm$R)
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
    okT <- attr(dadm,"unique_nortR") & dadm$winner
    LTs <- LT[okT]; UTs <- UT[okT]
    tpars <- pars[attr(dadm,"unique_nortR"),,drop=FALSE]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
                   dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nortR")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(LTs))
    for (i in 1:length(LTs)) {
      if (dim(tpars)[[1]]==1) pi <- t(as.matrix(tpars[,i,])) else
        pi <- tpars[,i,][order(!winner[,i]),]
      cf[i] <- pr_pt(LTs[i],UTs[i],pi,dadm,model)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nortR")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf) & !is.na(dadm$R)
    if (any(fix)) {
      lds[fix] <- lds[fix] + cf[fix]
    }
  }


  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
  } else ll <- lds
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  # Non-process (contaminant) miss.
  ispContaminant <- pars[,"pContaminant"]>0
  if ( any(ispContaminant) ) {
    p <- exp(ll)
    pc <- pars[dadm$winner,"pContaminant"]
    isMiss <- dadm$rt[dadm$winner]==Inf
    p[isMiss] <- pc[isMiss] + (1-pc[isMiss])*p[isMiss]
    p[!isMiss] <- (1-pc[!isMiss])*p[!isMiss]
    ll <- log(p)
  }


  return(sum(pmax(min_ll,ll[attr(dadm,"expand")])))
}


# NOT CHECKED !!!!!

log_likelihood_race_missing_LBAU <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood for an LBA allowing negative rates and
  # missing values
{

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

  pr_pt <- function(LT,UT,ps,dadm,model)
    # p(untrucated response)/p(truncated response), > 1, multiplicative truncation correction
  {
    pr <- my_integrate(f,lower=0,upper=Inf,p=ps,
                       dfun=model$dfun,pfun=model$pfun)
    if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
    if (pr$value==0) return(0)
    pt <- my_integrate(f,lower=LT,upper=UT,p=ps,
                       dfun=model$dfun,pfun=model$pfun)
    if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
    out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
    if (is.infinite(out)) return(NA)
    out
  }

  # pr_pt <- function(LT,UT,ps,dadm,model)
  #   # Sum over responses of the probability of a response with rt between LT and UT
  # {
  #   for ()
  #   pr <- my_integrate(f,lower=0,upper=Inf,p=ps,
  #                      dfun=model$dfun,pfun=model$pfun)
  #   if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
  #   if (pr$value==0) return(0)
  #   pt <- my_integrate(f,lower=LT,upper=UT,p=ps,
  #                      dfun=model$dfun,pfun=model$pfun)
  #   if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
  #   out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
  #   if (is.infinite(out)) return(NA)
  #   out
  # }


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


  pars <- get_pars_matrix(p_vector,dadm)

  if (any(names(dadm)=="RACE")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(model$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <-
    log(1-model$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- -Inf

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
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
                          dfun=model$dfun,pfun=model$pfun)
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
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
                          dfun=model$dfun,pfun=model$pfun)
      # if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
      #   lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
      if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
        lds[tofix][i] <- pmax(0,pmin(tmp$value,1))
      lds[tofix][i] <- log(lds[tofix][i] + (1-lds[tofix][i])*prod(pnorm(0,mpars[,i,"v"],mpars[,i,"sv"])))
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
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
                          dfun=model$dfun,pfun=model$pfun)
      if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
        p <- tmp$value
        if (dim(mpars)[[1]]==1) pi <- mpars[,i,,drop=FALSE] else
          pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
                            dfun=model$dfun,pfun=model$pfun)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
          p <- p + tmp$value else p <- 0
      } else p <- 0
      lds[tofix][i] <- log(pmax(0,pmin(p,1)))
    }
  }

  # Response unknown.
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixfast <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- my_integrate(f,lower=LT,upper=LC,p=pi,
                         dfun=model$dfun,pfun=model$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm,model) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=LT,upper=LC,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=model$dfun,pfun=model$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model) else cf <- 1
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
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- my_integrate(f,lower=UC,upper=UT,p=pi,
                         dfun=model$dfun,pfun=model$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm,model) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=UC,upper=UT,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
                           dfun=model$dfun,pfun=model$pfun)
        if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model) else cf <- 1
          if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p + (1-p)*prod(pnorm(0,mpars[,i,"v"],mpars[,i,"sv"])))
      if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
    }
  } else tofixslow <- NA
  # no direction
  nort <- is.na(dadm$rt) & is.na(dadm$R)
  nort <- nort & (pars[,"pContaminant"] == 0) # Otherwise not identifiable
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
                   dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
        pi <- mpars[,i,]
      pc <- pLU(LT,LC,UC,UT,pi,dadm,model)
      if (is.na(pc)) p <- NA else {
        if (pc!=0 & !(LT==0 & UT==Inf)) cf <- pr_pt(LT,UT,pi,dadm,model) else cf <- 1
        if (!is.na(cf)) p <- pc*cf else p <- NA
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- pLU(LT,LC,UC,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model)
          if (is.na(pc)) {
            p <- NA; break
          }
          if (pc!=0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm,model) else cf <- 1
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
    tpars <- pars[attr(dadm,"unique_nort"),,drop=FALSE]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
                   dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nort")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(ok))
    for (i in 1:length(ok)) if (ok[i]) {
      if (dim(tpars)[[1]]==1) pi <- t(as.matrix(tpars[,i,])) else
        pi <- tpars[,i,][order(!winner[,i]),]
      cf[i] <- pr_pt(LT,UT,pi,dadm,model)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nort")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf)
    if (any(fix)) lds[fix] <- lds[fix] + cf[fix]
    badfix <- dadm$winner & (is.na(cf) | is.nan(cf) | is.infinite(cf))
    if (!all(is.na(tofixfast))) badfix <- badfix & !tofixfast
    if (!all(is.na(tofixslow))) badfix <- badfix & !tofixslow
    if (any(badfix)) lds[badfix] <- -Inf
  }

  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
  } else ll <- lds
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  # Non-process (contaminant) miss.
  ispContaminant <- pars[,"pContaminant"]>0
  if ( any(ispContaminant) ) {
    p <- exp(ll)
    pc <- pars[dadm$winner,"pContaminant"]
    isMiss <- is.na(dadm$R[dadm$winner])
    p[isMiss] <- pc[isMiss] + (1-pc[isMiss])*p[isMiss]
    p[!isMiss] <- (1-pc[!isMiss])*p[!isMiss]
    ll <- log(p)
  }



  return(sum(pmax(min_ll,ll[attr(dadm,"expand")])))
}
