# The +0.035 recovery bias at alpha = 1.5 survives a 5x finer simulator step and
# is invariant to every discretisation knob, so it is neither the Euler
# simulator's step nor the PDE's grid.  Remaining possibilities: the simulator
# and the solver describe slightly different processes, or the profile estimator
# is itself biased at n = 2000.  Cut between them by drawing data from the
# solver's OWN cdf by inversion -- no simulator anywhere -- and recovering.
.libPaths(c(Sys.getenv("RLFLIB"), .libPaths()))  # RLFLIB: library to test
suppressMessages(library(EMC2))
pde <- EMC2:::rlf_fht_pdf_cdf_vec
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
source("refcdf.R")
cdf_for <- function(a) lapply(c(1, 2), function(v) ref_cdf(a, v, tg, B, A))
draw <- function(n, C) {
  u <- runif(n); out <- rep(Inf, n); ok <- u <= C[length(C)]
  out[ok] <- approx(C, tg, xout=u[ok], ties="ordered", rule=2)$y
  out
}
one <- function(a, seed, cdfs) {
  set.seed(seed); n <- 2000
  S <- sample(c("left","right"), n, TRUE)
  fin <- cbind(ifelse(S=="left", draw(n, cdfs[[2]]), draw(n, cdfs[[1]])),
               ifelse(S=="right", draw(n, cdfs[[2]]), draw(n, cdfs[[1]])))
  w <- max.col(-fin, ties.method="first")
  rt <- fin[cbind(1:n, w)] + T0
  dat <- data.frame(subjects=factor(1),
                    S=factor(S, levels=c("left","right")),
                    R=factor(c("left","right")[w], levels=c("left","right")),
                    rt=rt)
  dat <- dat[is.finite(dat$rt), ]
  emc <- suppressMessages(make_emc(dat, des, type="single"))
  dd <- EMC2:::.cache_ll_data_attrs(emc[[1]]$data[[1]])
  dg <- EMC2:::.oo_expanded_designs(dd)
  cons <- attr(dd,"constants"); if (is.null(cons)) cons <- NA
  g <- seq(max(1.01, a-0.30), min(1.98, a+0.30), by=0.02)
  P <- t(sapply(g, function(x){q<-p0; q["alpha"]<-qnorm(x-1); q}))
  colnames(P) <- names(p0)
  ll <- EMC2:::calc_ll_oo(P, dd, constants=cons, designs=dg, type="RLF",
    bounds=m$bound, transforms=m$transform, pretransforms=m$pre_transform,
    p_types=names(m$p_types), min_ll=log(1e-10), trend=m$trend, marginalise=NULL)
  i <- which.max(ll); v <- g[i]
  if (i > 1 && i < length(ll)) {
    y <- ll[(i-1):(i+1)]; d <- (y[1]-y[3])/(2*(y[1]-2*y[2]+y[3]))
    if (is.finite(d) && abs(d) <= 1) v <- v + d*0.02
  }
  v
}
cat("=== recovery from the solver's own cdf (inversion, no simulator) ===\n")
for (a in c(1.1, 1.5, 1.7)) {
  cdfs <- cdf_for(a)
  pk <- sapply(1:8, function(s) one(a, 2000+s, cdfs))
  cat(sprintf("alpha %.1f: %s\n  mean %.3f  bias %+0.3f  sd %.3f  se %.3f\n",
      a, paste(sprintf("%.3f", pk), collapse=" "), mean(pk), mean(pk)-a,
      sd(pk), sd(pk)/sqrt(8))); flush.console()
}
cat("\nselfdone\n")
