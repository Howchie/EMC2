### Single accumulator model

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with:
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review. Comments and changes added by Andrew Heathcote. Trish's
# code is for k = threshold, a = half width of uniform threshold variability,
# l = rate of accumulation. Note that Wald mean = k/l and shape = k^2.

# Andrew Heathcote: Following functions use a different parameterization in terms of v=l (rate),
# uniform start point variability from 0-A (A>=0), threshold b (>0) and hence
# B=b-A (>=0) as a threshold gap. Hence k = b-A/2 = B + A/2 and a=A

# Guy Hawkins: Addition of fS_E, fs_T and modified rWaldRace, dWaldRace. Distributed with:
#   Hawkins & Heathcote. Racing against the clock: Evidence-based vs timed-based decisions

rWald <- function(n,v,B,A)
  # random function for single acumulator
{

  rwaldt <- function(n,k,l,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss

    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }

    flag <- l>tiny
    x <- rep(NA,times=n)

    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- k^2

    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]

    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)

    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }

  # Act as if negative v never terminates, cluge to do single accumulator
  # case by passing negative v
  if (length(v)!=n) v <- rep(v,length.out=n)
  if (length(B)!=n) B <- rep(B,length.out=n)
  if (length(A)!=n) A <- rep(A,length.out=n)

  # Kluge to return -Inf for negative rates, so can implment one accumulator case
  out <- numeric(n)
  ok <- v>0
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok])
  out[!ok] <- Inf
  out
}


dWald <- function(t,v,B,A)
  # density for single accumulator
{

  digt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
    # pdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns digt.0 if a<1e-10

    digt.0 <- function(t,k=1,l=1) {
      # pdf of inverse gaussian at t with no k variability
      # much faster than statmod's dinvgauss funciton

      lambda <- k^2
      l0 <- l==0
      e <- numeric(length(t))
      if ( any(!l0) ) {
        mu <- k[!l0]/l[!l0]
        e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
      }
      if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
      x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
      x[t<=0] <- 0
      x
    }

    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))

    tpos <- t<=0

    atiny <- a<=tiny & !tpos
    a[atiny] <- 0

    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0

    x <- numeric(length(t))

    # No threshold variability
    if ( any(atiny) )
      x[atiny] <- digt.0(t=t[atiny],k=k[atiny],l=l[atiny])

    # Threshold variability
    if ( any(!atiny) ) {

      if ( any(notltiny) ) { # rate non-zero

        sqr.t <- sqrt(t[notltiny])

        term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])

        term.2a <- log(.5)+log(l[notltiny])
        term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
        term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.2d <- term.2b+term.2c
        term.2 <- exp(term.2a)*term.2d

        term.3 <- term.1+term.2
        term.4 <- log(term.3)-log(2)-log(a[notltiny])
        x[notltiny] <- exp(term.4)
      }

      if ( any(ltiny) ) {  # rate zero
        log.t <- log(t[ltiny])
        term.1 <- -.5*(log(2)+log(pi)+log.t)
        term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
        term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
        term.4 <- (exp(-term.2)-exp(-term.3))
        term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
        x[ltiny] <- exp(term.5)
      }

    }

    x[x<0 | is.nan(x) ] <- 0
    x
  }

  out <- numeric(length(t))
  ok <- v>0
  out[ok] <- digt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2)
  out[!ok] <- 0
  out
}



pWald <- function(t,v,B,A)
  # cumulative density for single accumulator
{
  pigt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
    # cdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns pigt.0 if a<=0

    pigt.0 <- function(t,k=1,l=1) {
      # cdf of inverse gaussian at t with no k variability
      # much faster than statmod's pinvgauss funciton

      mu <- k/l
      lambda <- k^2

      e <- exp(log(2*lambda) - log(mu))
      add <- sqrt(lambda/t) * (1 + t/mu)
      sub <- sqrt(lambda/t) * (1 - t/mu)

      p.1 <- 1 - pnorm(add)
      p.2 <- 1 - pnorm(sub)
      x <- exp(e + log(p.1)) + p.2

      x[t<0] <- 0
      x
    }

    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))

    tpos <- t<=0

    atiny <- a<=tiny & !tpos
    a[atiny] <- 0

    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0

    x <- numeric(length(t))

    # No threshold variability
    if ( any(atiny) )
      x[atiny] <- pigt.0(t[atiny],k[atiny],l[atiny])

    # Threshold variability
    if ( any(!atiny) ) {

      if ( any(notltiny) ) { # rate non-zero

        log.t <- log(t[notltiny])
        sqr.t <- sqrt(t[notltiny])

        term.1a <- .5*log.t-.5*log(2*pi)
        term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1 <- exp(term.1a)*(term.1b-term.1c)

        term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) +
                         log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) +
                         log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])

        term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
        term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
        term.4 <- term.4c*term.4a + term.4d*term.4b

        x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
      }

      if ( any(ltiny) ) {  # rate zero
        sqr.t <- sqrt(t[ltiny])
        log.t <- log(t[ltiny])
        term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
        term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
        term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])

        term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6 <- 1 + exp(term.6b) - exp(term.6a)

        x[ltiny] <- term.5 + term.6
      }

    }

    x[x<0 | is.nan(x) ] <- 0
    x
  }

  out <- numeric(length(t))
  ok <- v>0
  out[ok] <- pigt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2)
  out[!ok] <- 0
  out

}



# density and survivor for evidence accumulators
fS_E <- function(t, x, drift, threshold, n.acc) {
  f <- S <- matrix(nrow=n.acc, ncol=length(t))
  for(i in 1:n.acc) {
    f[i,] <- dWald(t, v=drift[i], B=threshold[i], A=0)
    S[i,] <- 1-pWald(t, v=drift[i], B=threshold[i], A=0)
  }
  list(f=f, S=S)
}

# density and survivor for timing accumulator
fS_T <- function(t, x) {
  f <- dWald(t, v=x["v_T"], B=x["b_T"], A=0)
  S <- 1-pWald(t, v=x["v_T"], B=x["b_T"], A=0)
  list(f=f, S=S)
}


### random generation for N racing Walds for N choice with additional timing accumulator
rWaldRace <- function(n, x, drift, threshold) {
  # random function for Wald race with timer
  n.acc <- ifelse(is.null(dim(drift)), length(drift), dim(drift)[1])
  # first, sample responses from evidence accumulation process
  ttf_E <- matrix(x["t0_E"] + rWald(n*n.acc, v=drift, B=threshold, A=0), nrow=n.acc)
  # and store responses and RT
  resp <- apply(ttf_E, 2, which.min)
  out <- data.frame(R=resp, RT=ttf_E[cbind(resp,1:n)],timer.response=FALSE)

  # second, sample timing RTs then resample responses if timer terminated first (to 1 of N alternatives)
  ttf_T <- x["t0_T"] + rWald(n, v=x["v_T"], B=x["b_T"], A=0)
  timer.replace <- ttf_T < out$RT

  # resample response only if there were >0 timer responses
  if(any(timer.replace)) {
    out$timer.response[timer.replace] <- TRUE
    out$RT[timer.replace] <- ttf_T[timer.replace]
    # then determine which response was given
    # unbiased guessing probablity of 1/N
    p.guess <- rep(1/n.acc, n.acc)
    out$R[timer.replace] <- sample(1:n.acc, size=sum(timer.replace), replace=TRUE, prob=p.guess)
  }
  out
}



### density for N racing Walds for N choice with additional timing accumulator
dWaldRace <- function(t, x, drift, threshold) {
  # Generates defective PDF for responses on node=1, t (decison time) is a vector of times
  # get density for drift[1] (observed response) vs drift[2] (non-response) vs timer

  # race equation for evidence process
  RaceLike <- function(E, n.acc) {
    if(n.acc==2) {
      ES <- E$S[n.acc,]
    } else {
      ES <- apply(E$S[-1,], 2 ,prod)
    }
    # density of node 1 response (E$f) multiplied by survivor of remaining
    #   evidence accumulator nodes (ES)
    E$f[1,] * ES
  }

  # race equation for evidence process vs time process
  TimerRaceLike <- function(E, T, p.guess, n.acc) {
    # two ways to generate a response:
    # 1. density of node 1 response (E$f) multiplied by survivor of remaining
    #   evidence accumulator nodes (ES) and timing node (T$S)
    # 2. density of timing node (T$f) and probability of response (p.guess)
    #   multiplied by survivor of all evidence accumulator nodes (E$S)
    (RaceLike(E, n.acc) * T$S) +
    (T$f * p.guess * apply(E$S, 2, prod))
  }

  # helper function for simpler code below
  GetLike <- function(t_E, t_T, x, drift, threshold, p.guess, n.acc) {
    E <- fS_E(t=t_E, x, drift, threshold, n.acc)
    T <- fS_T(t=t_T, x)
    TimerRaceLike(E, T, p.guess, n.acc)
  }

  n.acc <- ifelse(is.null(dim(drift)), length(drift), dim(drift)[1])
  # assume unbiased guessing probablity of 1/N
  p.guess <- 1/n.acc
  out <- numeric(length(t))

  ### timing process starts earlier than evidence process ###
  if(x["t0_T"] < x["t0_E"]) {
    dt <- t - x["t0_T"]
    # assume evidence process starts at d=t0_E-t0_T after the timing process
    d <- x["t0_E"] - x["t0_T"]
    u <- dt < d
    if(any(u)) {  # evidence process hasn't started so must be a timer response with a guess
      out[u] <- fS_T(dt[u], x)$f * p.guess
    }
    if(any(!u)) { # could be a response from the timer or evidence process
      out[!u] <- GetLike(t_E=dt[!u]-d, t_T=dt[!u], x, drift, threshold, p.guess, n.acc)
    }

  ### evidence process starts earlier than timing process ###
  } else if(x["t0_E"] < x["t0_T"]) {
    dt <- t - x["t0_E"]
    # assume timing process starts at d=t0_T-t0_E after the evidence process
    d <- x["t0_T"] - x["t0_E"]
    u <- dt < d
    if(any(u)) {  # timing process hasn't started so must be an evidence response
      E <- fS_E(t=dt[u], x, drift, threshold, n.acc)
      out[u] <- RaceLike(E, n.acc)
    }
    if(any(!u)) { # could be a response from the timer or evidence process
      out[!u] <- GetLike(t_E=dt[!u], t_T=dt[!u]-d, x, drift, threshold, p.guess, n.acc)
    }

    ### evidence process and timing process start at the same time ###
  } else if(x["t0_E"] == x["t0_T"]) {
    dt <- t - x["t0_E"]
    out <- GetLike(t_E=dt, t_T=dt, x, drift, threshold, p.guess, n.acc)
  }
  out
}

