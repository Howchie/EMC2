# Isolates the cause of the fixed alpha recovery offset: rt_resolution.
# make_emc floors rt to 1/60 by default and the likelihood then evaluates the
# density at the floored time instead of integrating it over the bin, and alpha
# absorbs the mismatch.  Env: ATRUE (generating alpha), RTRES (rt_resolution,
# empty for NULL).  With RTRES empty recovery is unbiased at every alpha; at
# 1/60 it is biased high by +0.012 to +0.023, worst near alpha 1.5.
.libPaths(c(Sys.getenv("RLFLIB"), .libPaths())); suppressMessages(library(EMC2))
B <- 1.5; A <- 0.5; T0 <- 0.2
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
des <- design(factors=list(subjects=1,S=c("left","right")), Rlevels=c("left","right"),
  matchfun=matchfun, model=RLF(), formula=list(v~lM,B~1,A~1,t0~1,s~1,alpha~1),
  constants=c(s=log(1)))
p0 <- sampled_pars(des, doMap=FALSE)
p0["B"]<-log(B); p0["A"]<-log(A); p0["t0"]<-log(T0)
p0["v"]<-log(1); p0["v_lMTRUE"]<-log(2)
m <- RLF()
options(emc2.rlf_simd_batch=FALSE, emc2.rlf_horizon_split=TRUE)
tg <- c(seq(1e-3, 8, length.out=4000), seq(8.01, 22, length.out=1500))

cdf_at <- function(a, v, nx, rich) {
  o <- options(emc2.rlf_nx=nx, emc2.rlf_richardson=rich,
               emc2.rlf_richardson_ratio=1.5, emc2.rlf_dt=4e-3)
  on.exit(options(o))
  pars <- cbind(v=rep(v, length(tg)), B=B, A=A, t0=0, s=1, alpha=a)
  cummax(EMC2:::pRLF(tg, pars))
}
draw <- function(n, C) {
  u <- runif(n); out <- rep(Inf, n); ok <- u <= C[length(C)]
  out[ok] <- approx(C, tg, xout=u[ok], ties="ordered", rule=2)$y
  out
}
fit_peak <- function(dat, a, nx, rich, step=0.01) {
  o <- options(emc2.rlf_nx=nx, emc2.rlf_richardson=rich,
               emc2.rlf_richardson_ratio=1.5, emc2.rlf_dt=8e-3)
  on.exit(options(o))
  emc <- suppressMessages(make_emc(dat, des, type="single", rt_resolution=RTRES))
  dd <- EMC2:::.cache_ll_data_attrs(emc[[1]]$data[[1]])
  dg <- EMC2:::.oo_expanded_designs(dd)
  cons <- attr(dd,"constants"); if (is.null(cons)) cons <- NA
  g <- seq(max(1.005, a-0.25), min(1.985, a+0.25), by=step)
  P <- t(sapply(g, function(x){q<-p0; q["alpha"]<-qnorm(x-1); q}))
  colnames(P) <- names(p0)
  ll <- EMC2:::calc_ll_oo(P, dd, constants=cons, designs=dg, type="RLF",
    bounds=m$bound, transforms=m$transform, pretransforms=m$pre_transform,
    p_types=names(m$p_types), min_ll=log(1e-10), trend=m$trend, marginalise=NULL)
  i <- which.max(ll); v <- g[i]
  if (i > 1 && i < length(ll)) {
    y <- ll[(i-1):(i+1)]; d <- (y[1]-y[3])/(2*(y[1]-2*y[2]+y[3]))
    if (is.finite(d) && abs(d) <= 1) v <- v + d*step
  }
  v
}
mkdat <- function(seed, n, cdfs) {
  set.seed(seed)
  S <- sample(c("left","right"), n, TRUE)
  fin <- cbind(ifelse(S=="left", draw(n, cdfs[[2]]), draw(n, cdfs[[1]])),
               ifelse(S=="right", draw(n, cdfs[[2]]), draw(n, cdfs[[1]])))
  w <- max.col(-fin, ties.method="first")
  d <- data.frame(subjects=factor(1), S=factor(S, levels=c("left","right")),
                  R=factor(c("left","right")[w], levels=c("left","right")),
                  rt=fin[cbind(1:n, w)] + T0)
  d[is.finite(d$rt), ]
}

A_TRUE <- if (nzchar(Sys.getenv("ATRUE"))) as.numeric(Sys.getenv("ATRUE")) else 1.5; N <- 8000; SEEDS <- 1:12
RTRES <- if (nzchar(Sys.getenv("RTRES"))) as.numeric(Sys.getenv("RTRES")) else NULL
cfg <- list(
  list("gen 192 pair -> fit 192 pair", 192L, TRUE, 192L, TRUE))
cat(sprintf("=== matched vs mismatched generator/estimator, alpha %.1f, n %d ===\n",
            A_TRUE, N))
cat(sprintf("%-30s %9s %9s %8s\n", "configuration", "mean", "bias", "se"))
for (k in cfg) {
  cdfs <- lapply(c(1, 2), function(v) cdf_at(A_TRUE, v, k[[2]], k[[3]]))
  pk <- sapply(SEEDS, function(s) fit_peak(mkdat(7000+s, N, cdfs), A_TRUE,
                                           k[[4]], k[[5]]))
  cat(sprintf("%-30s %9.4f %+9.4f %8.4f\n", k[[1]], mean(pk),
              mean(pk)-A_TRUE, sd(pk)/sqrt(length(pk)))); flush.console()
}
cat("\nmatcheddone\n")
