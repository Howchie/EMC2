# Run from the repository root. Research validation, not a production test.
Rcpp::sourceCpp("WorkingTests/rdmswtn_analytic_series.cpp")
set.seed(20260921)
n <- 2000L
p <- cbind(t=runif(n,.1,3),mu=runif(n,.05,3),b=0,A=runif(n,.01,.8),
           s=1,sv=runif(n,.05,1.5),pos=1)
p[,"b"] <- runif(n,.1,2)+p[,"A"]
for (pos in c(1,0)) {
  p[,"pos"] <- pos
  if (!pos) p[,"mu"] <- p[,"mu"]-1.5
  ref <- quadrature_cdf(p,80)
  stopifnot(max(abs(backward_cdf(p)-ref))<1e-10)
  for (degree in c(16,24,32,48)) {
    ans <- series_cdf(p,degree)
    cat("pos",pos,"degree",degree,"max abs",max(abs(ans-ref)),
        "max relative",max(abs(ans-ref)/pmax(ref,1e-12)),"\n")
  }
  if(pos) {
    worst <- order(abs(ans-ref),decreasing=TRUE)[1:3]
    print(cbind(p[worst,],series=ans[worst],quad=ref[worst]))
    worst_rows <- p[worst,,drop=FALSE]
    cat("Subdivided series max abs",max(abs(series_cdf(p,32,0)-ref)),"\n")
    cat("Backward series max abs",max(abs(backward_cdf(p)-ref)),"\n")
  }
}
p[,"pos"] <- 1; p[,"mu"] <- abs(p[,"mu"])
# Independent order of integration: average analytic fixed-drift SPV CDF.
ref_drift <- function(row,density=FALSE) {
  norm <- pnorm(row["mu"]/row["sv"])
  kernel <- if(density) EMC2:::dwald else EMC2:::pwald
  integrate(function(z) vapply(z,function(zz) {
    drift <- row["mu"]+row["sv"]*zz
    kernel(row["t"], drift, row["b"], row["A"], s=row["s"],
                 posdrift=FALSE)*dnorm(zz)/norm
  },numeric(1)),lower=max(-row["mu"]/row["sv"],-12),upper=12,
  rel.tol=1e-10,abs.tol=1e-12,subdivisions=300)$value
}
ii <- seq(1,n,length.out=200)
independent <- apply(p[ii,,drop=FALSE],1,ref_drift)
cat("Independent drift integration max abs",max(abs(series_cdf(p[ii,],32)-independent)),"\n")
cat("Independent backward max abs",max(abs(backward_cdf(p[ii,])-independent)),"\n")
stopifnot(max(abs(backward_cdf(p[ii,])-independent))<1e-8)
dt <- 1e-5
plus <- minus <- p[ii,,drop=FALSE]
plus[,1] <- plus[,1]+dt; minus[,1] <- minus[,1]-dt
derivative <- (backward_cdf(plus)-backward_cdf(minus))/(2*dt)
pdf <- apply(p[ii,,drop=FALSE],1,function(r)
  EMC2:::drdmswtn(r[1],r[2],r[3],r[4],s=r[5],sv=r[6]))
cat("Time derivative vs existing analytic density max abs",max(abs(derivative-pdf)),"\n")
pdf_ref <- apply(p[ii,,drop=FALSE],1,ref_drift,density=TRUE)
cat("Time derivative vs independent density max abs",max(abs(derivative-pdf_ref)),"\n")
cat("Existing analytic density vs independent max abs",max(abs(pdf-pdf_ref)),"\n")
stopifnot(max(abs(derivative-pdf_ref))<1e-7)
print(cbind(worst_rows,reference=apply(worst_rows,1,ref_drift),
            series=series_cdf(worst_rows,32),split=series_cdf(worst_rows,32,0),
            quad=quadrature_cdf(worst_rows,80)),digits=15)
# Fair compiled-loop timing; include a one-threshold baseline.
timing <- function(f) median(replicate(5,system.time(for(j in 1:20) f())["elapsed"]))
baseline <- p; baseline[,"A"] <- 0
times <- c(fixed=timing(function() quadrature_cdf(baseline,1)),
           original_fixed=timing(function() quadrature_cdf(baseline,1,TRUE)),
           original_gl20=timing(function() quadrature_cdf(p,20,TRUE)),
           gl20=timing(function() quadrature_cdf(p,20)),
           series32=timing(function() series_cdf(p,32)),
           split32=timing(function() series_cdf(p,32,0)),
           backward24=timing(function() backward_cdf(p)))
print(times); print(times/times["series32"])
cat("Original algebra max difference vs backward",
    max(abs(quadrature_cdf(p,20,TRUE)-backward_cdf(p))),"\n")
# Stress dimensions: these intentionally expose when a single expansion fails.
stress <- as.matrix(expand.grid(t=c(.005,.05,.5,5),mu=c(.05,1,5),
                               B=c(.05,.5,3),A=c(.05,.5,2),s=c(.5,1),
                               sv=c(.1,1,3),pos=1))
stress[,"B"] <- stress[,"B"]+stress[,"A"]
colnames(stress)[3] <- "b"
ref <- quadrature_cdf(stress,80); ans <- series_cdf(stress,32)
cat("Stress cases",nrow(stress),"nonfinite reference",sum(!is.finite(ref)),
    "single-series failures abs>1e-8",sum(abs(ans-ref)>1e-8,na.rm=TRUE),"\n")
back <- backward_cdf(stress)
cat("Stress backward nonfinite",sum(!is.finite(back)),"difference>1e-8",
    sum(abs(back-ref)>1e-8,na.rm=TRUE),"max finite difference",
    max(abs(back-ref),na.rm=TRUE),"\n")
ind_stress <- apply(stress,1,ref_drift)
cat("Stress backward vs independent max abs",max(abs(back-ind_stress)),
    "failures>1e-8",sum(abs(back-ind_stress)>1e-8),"\n")
stopifnot(all(is.finite(back)),max(abs(back-ind_stress))<1e-8)
