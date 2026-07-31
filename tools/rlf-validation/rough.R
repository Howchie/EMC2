# Where does the roughness of ll(alpha) come from?  sched.R found 1.2-2.9 nats
# of RMS residual about a local quadratic under an alpha-dependent nx, against
# 0.05-0.33 at fixed nx -- but one of those rows had nx constant at 48 across
# the whole window, so quantised nx cannot be the only source.  Candidates:
#   (a) nx steps as the schedule crosses a rounding boundary,
#   (b) coarse grids being intrinsically jagged in alpha,
#   (c) horizon bucketing: t_crit = 2 b^alpha moves with alpha, so individual
#       trials change bucket mid-profile and their domain jumps.
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
options(emc2.rlf_dt=8e-3, emc2.rlf_simd_batch=FALSE)

mk <- function(atrue, agrid) {
  p <- p0; p["alpha"] <- qnorm(atrue-1)
  set.seed(42 + round(atrue*100)); dat <- make_data(p, des, n_trials=2000)
  emc <- suppressMessages(make_emc(dat, des, type="single"))
  dd <- EMC2:::.cache_ll_data_attrs(emc[[1]]$data[[1]])
  dg <- EMC2:::.oo_expanded_designs(dd)
  cons <- attr(dd,"constants"); if (is.null(cons)) cons <- NA
  P <- t(sapply(agrid, function(a){q<-p; q["alpha"]<-qnorm(a-1); q}))
  colnames(P) <- names(p0)
  list(dd=dd, dg=dg, cons=cons, P=P)
}
run <- function(D) EMC2:::calc_ll_oo(D$P, D$dd, constants=D$cons, designs=D$dg,
  type="RLF", bounds=m$bound, transforms=m$transform, pretransforms=m$pre_transform,
  p_types=names(m$p_types), min_ll=log(1e-10), trend=m$trend, marginalise=NULL)

cfg <- list(
  list("raw   nx= 48        ",  48L, FALSE, 0,    TRUE),
  list("raw   nx= 64        ",  64L, FALSE, 0,    TRUE),
  list("raw   nx= 96        ",  96L, FALSE, 0,    TRUE),
  list("raw   nx=128        ", 128L, FALSE, 0,    TRUE),
  list("pair  nx= 48        ",  48L, TRUE,  0,    TRUE),
  list("pair  nx=128        ", 128L, TRUE,  0,    TRUE),
  list("pair  nx=128 nosplit", 128L, TRUE,  0,    FALSE),
  list("raw   nx=128 nosplit", 128L, FALSE, 0,    FALSE),
  list("sched nx=128 k=1.84 ", 128L, TRUE,  1.84, TRUE),
  list("sched nx=128 k=1.84 nosplit", 128L, TRUE, 1.84, FALSE))

for (a in c(1.3, 1.7)) {
  g <- seq(a-0.05, a+0.05, by=0.005)
  D <- mk(a, g)
  cat(sprintf("\n== alpha %.1f, profile step 0.005, quadratic detrended ==\n", a))
  cat(sprintf("%-28s %9s %9s %9s\n", "config", "rms", "max", "d(ll) step"))
  for (k in cfg) {
    options(emc2.rlf_nx=k[[2]], emc2.rlf_richardson=k[[3]],
            emc2.rlf_richardson_ratio=1.5, emc2.rlf_nx_alpha=k[[4]],
            emc2.rlf_horizon_split=k[[5]], emc2.rlf_nx_floor=48L)
    ll <- run(D)
    r <- resid(lm(ll ~ poly(g, 2)))
    cat(sprintf("%-28s %9.4f %9.4f %9.4f\n", k[[1]], sqrt(mean(r^2)),
                max(abs(r)), max(abs(diff(r))))); flush.console()
  }
}
cat("\nroughdone\n")
