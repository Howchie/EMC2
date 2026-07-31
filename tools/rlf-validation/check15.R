# The alpha = 1.5 data set's converged peak sits at 1.581, +0.081 from truth,
# larger and opposite in sign to the alpha = 1.1 offset.  Three candidates:
#   (a) one unlucky seed  -> scatter across seeds should cover it
#   (b) simulator bias    -> shrinking emc2.rlf_sim_dt should move it
#   (c) the solver        -> the peak should be invariant to knobs that only
#                            change the discretisation (domain width, dt), and
#                            all three reference rungs agreeing does not prove
#                            that, since they share one scheme.
.libPaths(c(Sys.getenv("RLFLIB"), .libPaths()))  # RLFLIB: library to test
suppressMessages(library(EMC2))
matchfun <- function(d) as.numeric(d$S) == as.numeric(d$lR)
des <- design(factors=list(subjects=1,S=c("left","right")), Rlevels=c("left","right"),
  matchfun=matchfun, model=RLF(), formula=list(v~lM,B~1,A~1,t0~1,s~1,alpha~1),
  constants=c(s=log(1)))
p0 <- sampled_pars(des, doMap=FALSE)
p0["B"]<-log(1.5); p0["A"]<-log(0.5); p0["t0"]<-log(0.2)
p0["v"]<-log(1); p0["v_lMTRUE"]<-log(2)-log(1)
m <- RLF()
options(emc2.rlf_dt=8e-3, emc2.rlf_simd_batch=FALSE, emc2.rlf_horizon_split=TRUE,
        emc2.rlf_nx=192L, emc2.rlf_richardson=TRUE, emc2.rlf_richardson_ratio=1.5)

gen <- function(atrue, seed, simdt=1e-3) {
  p <- p0; p["alpha"] <- qnorm(atrue-1)
  old <- getOption("emc2.rlf_sim_dt"); options(emc2.rlf_sim_dt=simdt)
  set.seed(seed); dat <- make_data(p, des, n_trials=2000)
  options(emc2.rlf_sim_dt=old)
  dat
}
prep <- function(dat, atrue) {
  emc <- suppressMessages(make_emc(dat, des, type="single"))
  dd <- EMC2:::.cache_ll_data_attrs(emc[[1]]$data[[1]])
  dg <- EMC2:::.oo_expanded_designs(dd)
  cons <- attr(dd,"constants"); if (is.null(cons)) cons <- NA
  g <- seq(max(1.01, atrue-0.30), min(1.98, atrue+0.30), by=0.02)
  p <- p0; p["alpha"] <- qnorm(atrue-1)
  P <- t(sapply(g, function(a){q<-p; q["alpha"]<-qnorm(a-1); q}))
  colnames(P) <- names(p0)
  list(dd=dd, dg=dg, cons=cons, P=P, g=g)
}
peak <- function(D) {
  ll <- EMC2:::calc_ll_oo(D$P, D$dd, constants=D$cons, designs=D$dg, type="RLF",
    bounds=m$bound, transforms=m$transform, pretransforms=m$pre_transform,
    p_types=names(m$p_types), min_ll=log(1e-10), trend=m$trend, marginalise=NULL)
  i <- which.max(ll); v <- D$g[i]
  if (i > 1 && i < length(ll)) {
    y <- ll[(i-1):(i+1)]; d <- (y[1]-y[3])/(2*(y[1]-2*y[2]+y[3]))
    if (is.finite(d) && abs(d) <= 1) v <- v + d*(D$g[2]-D$g[1])
  }
  v
}

cat("=== (a) scatter of the converged peak across seeds (nx 192+288) ===\n")
for (atrue in c(1.1, 1.5, 1.7)) {
  pk <- numeric(0)
  for (sd in 1:8) pk <- c(pk, peak(prep(gen(atrue, 1000+sd), atrue)))
  cat(sprintf("alpha %.1f: %s\n  mean %.3f  bias %+0.3f  sd %.3f  se(mean) %.3f\n",
      atrue, paste(sprintf("%.3f", pk), collapse=" "), mean(pk),
      mean(pk)-atrue, sd(pk), sd(pk)/sqrt(length(pk)))); flush.console()
}

cat("\n=== (b) simulator step: same seeds, emc2.rlf_sim_dt 1e-3 vs 2e-4 ===\n")
for (atrue in c(1.1, 1.5)) {
  a3 <- a4 <- numeric(0)
  for (sd in 1:4) {
    a3 <- c(a3, peak(prep(gen(atrue, 1000+sd, 1e-3), atrue)))
    a4 <- c(a4, peak(prep(gen(atrue, 1000+sd, 2e-4), atrue)))
  }
  cat(sprintf("alpha %.1f  dt 1e-3 mean %.3f | dt 2e-4 mean %.3f | shift %+0.3f\n",
      atrue, mean(a3), mean(a4), mean(a4)-mean(a3))); flush.console()
}

cat("\n=== (c) solver invariance on the original alpha=1.5 seed ===\n")
dat <- gen(1.5, 42+150); D <- prep(dat, 1.5)
probe <- function(lbl, ...) {
  o <- options(...); on.exit(options(o))
  cat(sprintf("  %-34s peak %.3f\n", lbl, peak(D))); flush.console()
}
probe("192+288, dt 8e-3, width x1 (base)")
probe("  dt 4e-3", emc2.rlf_dt=4e-3)
probe("  dt 2e-3", emc2.rlf_dt=2e-3)
probe("  domain width x1.5", emc2.rlf_width_scale=1.5)
probe("  domain width x2.5", emc2.rlf_width_scale=2.5)
probe("  horizon split off", emc2.rlf_horizon_split=FALSE)
probe("  nx 288+432", emc2.rlf_nx=288L)
probe("  nx 384+576", emc2.rlf_nx=384L)
probe("  nx 384 raw (no extrapolation)", emc2.rlf_nx=384L, emc2.rlf_richardson=FALSE)
cat("\ncheckdone\n")
