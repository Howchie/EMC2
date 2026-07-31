# Removing the simulator entirely (data drawn by inversion from the solver's own
# cdf) leaves the alpha = 1.5 bias intact: +0.033 self-generated vs +0.035
# simulated.  So it is not a simulator/likelihood mismatch, and the solver is
# converged and invariant.  That leaves a finite-sample property of the
# estimator.  Discriminate the two remaining stories by scaling n: an O(1/n)
# estimator bias falls with sample size, a model or likelihood error does not.
.libPaths(c(Sys.getenv("RLFLIB"), .libPaths()))  # RLFLIB: library to test
suppressMessages(library(EMC2))
source("refcdf.R")
B <- 1.5; A <- 0.5; T0 <- 0.2; b0 <- B + A; z0 <- A
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
des <- design(factors=list(subjects=1,S=c("left","right")), Rlevels=c("left","right"),
  matchfun=matchfun, model=RLF(), formula=list(v~lM,B~1,A~1,t0~1,s~1,alpha~1),
  constants=c(s=log(1)))
p0 <- sampled_pars(des, doMap=FALSE)
p0["B"]<-log(B); p0["A"]<-log(A); p0["t0"]<-log(T0)
p0["v"]<-log(1); p0["v_lMTRUE"]<-log(2)
m <- RLF()
options(emc2.rlf_dt=8e-3, emc2.rlf_simd_batch=FALSE, emc2.rlf_horizon_split=TRUE,
        emc2.rlf_nx=192L, emc2.rlf_richardson=TRUE, emc2.rlf_richardson_ratio=1.5)
tg <- c(seq(1e-3, 8, length.out=4000), seq(8.01, 22, length.out=1500))
draw <- function(n, C) {
  u <- runif(n); out <- rep(Inf, n); ok <- u <= C[length(C)]
  out[ok] <- approx(C, tg, xout=u[ok], ties="ordered", rule=2)$y
  out
}
one <- function(a, seed, cdfs, n, step=0.01) {
  set.seed(seed)
  S <- sample(c("left","right"), n, TRUE)
  fin <- cbind(ifelse(S=="left", draw(n, cdfs[[2]]), draw(n, cdfs[[1]])),
               ifelse(S=="right", draw(n, cdfs[[2]]), draw(n, cdfs[[1]])))
  w <- max.col(-fin, ties.method="first")
  dat <- data.frame(subjects=factor(1), S=factor(S, levels=c("left","right")),
                    R=factor(c("left","right")[w], levels=c("left","right")),
                    rt=fin[cbind(1:n, w)] + T0)
  dat <- dat[is.finite(dat$rt), ]
  emc <- suppressMessages(make_emc(dat, des, type="single"))
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
cat("=== does the bias fall like 1/n? (self-generated, no simulator) ===\n")
cat(sprintf("%-7s %8s %10s %8s %8s %9s\n", "alpha", "n", "mean", "bias", "sd", "se"))
for (a in c(1.5, 1.1)) {
  cdfs <- lapply(c(1, 2), function(v) ref_cdf(a, v, tg, B, A))
  for (n in c(2000, 8000, 32000)) {
    pk <- sapply(1:12, function(s) one(a, 5000 + s, cdfs, n))
    cat(sprintf("%-7.1f %8d %10.4f %+8.4f %8.4f %9.4f\n", a, n, mean(pk),
        mean(pk) - a, sd(pk), sd(pk)/sqrt(length(pk)))); flush.console()
  }
}
cat("\nnscaledone\n")
