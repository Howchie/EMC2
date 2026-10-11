pkgload::load_all(".", quiet=TRUE, compile=FALSE)
emc <- EMC2:::restore_custom_kernel_pointers(EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")), quiet=TRUE)
Tn <- dim(emc[[1]]$samples$alpha)[3]; all_idx <- seq_len(Tn); p <- emc[[1]]$n_pars
chain <- emc[[2]]
alpha <- matrix(chain$samples$alpha[,,chain$samples$idx],p)
ll <- as.numeric(chain$samples$subj_ll[,chain$samples$idx])
mu <- chain$samples$theta_mu[,chain$samples$idx]
V <- as.matrix(chain$samples$theta_var[,,chain$samples$idx])
chart <- EMC2:::.emc_collapsed_chart(EMC2:::.emc_shear_v_spec(chain),p)
phi <- EMC2:::.emc_collapsed_to_phi(mu,V,chart)
fixed_idx <- EMC2:::.emc_collapsed_fixed_idx(chain)
fixed_value <- alpha[fixed_idx,,drop=FALSE]
cat("fixed parameters:",paste(chain$par_names[fixed_idx],collapse=", "),"\n")

fit_ridge <- function(A,y,lambda) {
  ok <- is.finite(y) & stats::complete.cases(A); A<-A[ok,,drop=FALSE]; y<-y[ok]
  p<-ncol(A); iu<-which(upper.tri(diag(p),diag=TRUE),arr.ind=TRUE)
  if(nrow(A)<2L*(1L+p+nrow(iu))) return(NULL)
  abar<-colMeans(A); e<-eigen(stats::cov(A),symmetric=TRUE)
  if(any(!is.finite(e$values)) || any(e$values<=1e-12*max(e$values))) return(NULL)
  W<-t(e$vectors)/sqrt(e$values); Z<-sweep(A,2,abar)%*%t(W)
  Q<-Z[,iu[,1]]*Z[,iu[,2],drop=FALSE]; X<-cbind(1,Z,Q)
  scl<-sqrt(colMeans(X^2)); scl[1]<-1; Xs<-sweep(X,2,scl,"/")
  XtX<-crossprod(Xs); penalty<-diag(ncol(Xs)); penalty[1,1]<-0
  b<-tryCatch(solve(XtX+lambda*nrow(Xs)*penalty,crossprod(Xs,y)),error=function(e)NULL)
  if(is.null(b)) return(NULL); b<-drop(b)/scl
  cq<-b[-seq_len(p+1L)]; B<-matrix(0,p,p); B[iu]<- -cq; B<-B+t(B); diag(B)<- -2*cq[iu[,1]==iu[,2]]
  eb<-eigen(B,symmetric=TRUE); Bp<-eb$vectors%*%(pmax(eb$values,0)*t(eb$vectors))
  g<-qr.coef(qr(cbind(1,Z)),y+0.5*rowSums((Z%*%Bp)*Z))[-1L]
  tol<-max(abs(eb$values),.Machine$double.xmin)*1e-10; keep<-which(eb$values>tol)
  if(length(keep)) g<-drop(eb$vectors[,keep,drop=FALSE]%*%crossprod(eb$vectors[,keep,drop=FALSE],g)) else g[]<-0
  H<-t(W)%*%Bp%*%W; H<-(H+t(H))/2; lin<-drop(H%*%abar+t(W)%*%g)
  if(any(!is.finite(H))||any(!is.finite(lin))) return(NULL)
  list(prec=H,lin=lin)
}
fit_ridge_all <- function(lambda) lapply(seq_len(chain$n_subjects),function(s) {
  A<-do.call(rbind,lapply(emc,function(ch)t(matrix(ch$samples$alpha[,s,all_idx],p))))
  y<-unlist(lapply(emc,function(ch)ch$samples$subj_ll[s,all_idx]),use.names=FALSE)
  fit_ridge(A,y,lambda)
})
base_fits <- EMC2:::.emc_collapsed_fits(emc,all_idx)
fit_sets <- list(empirical=base_fits,ridge01=fit_ridge_all(0.1))
base_plan <- chain$collapsed$plan; C<-base_plan$C
ns<-asNamespace("EMC2"); binding<-".emc_wpool_subject_ll"; old<-get(binding,envir=ns,inherits=FALSE); capture<-new.env(parent=emptyenv())
unlockBinding(binding,ns)
assign(binding,function(pool,rows,ctx,part=NULL,source) { z<-old(pool,rows,ctx,part,source); if(identical(source,"collapsed")){capture$rows<-rows;capture$ll<-as.numeric(z$ll)}; z },envir=ns)
lockBinding(binding,ns)
on.exit({unlockBinding(binding,ns);assign(binding,old,envir=ns);lockBinding(binding,ns)},add=TRUE)
results<-list()
for(nm in names(fit_sets)) {
  chain2<-chain; chain2$collapsed$fits<-fit_sets[[nm]]
  ctx<-EMC2:::.emc_collapsed_context(chain2,fixed_idx,fixed_value)
  eps<-EMC2:::.emc_collapsed_tune_eps(phi,C,base_plan$eps,ctx)
  plan<-base_plan; plan$C<-C; plan$eps<-eps
  cur<-EMC2:::.emc_shear_state(mu,V,ctx$fits,chain$prior,ctx$spec)
  x0<-list(cur=cur,alpha=alpha,ll=ll,rounds=0L,wpool=NULL,wpool_part=NULL,moved=FALSE)
  set.seed(44127)
  acc<-jump<-logcorr<-numeric(20)
  started<-proc.time()[["elapsed"]]
  for(i in seq_len(20)) {
    capture$rows<-NULL; capture$ll<-NULL
    xx<-EMC2:::.emc_collapsed_round(x0,plan,chain2,ctx)
    acc[i]<-xx$acc
    jump[i]<-sum(xx$jump^2)
    if(!is.null(capture$rows)) {
      anew<-matrix(vapply(seq_len(chain$n_subjects),function(s)as.numeric(capture$rows[[s]][1L,chain$par_names]),numeric(p)),p)
      logcorr[i]<-sum(capture$ll-EMC2:::.emc_collapsed_log_lhat(anew,ctx$fits))-
        sum(ll-EMC2:::.emc_collapsed_log_lhat(alpha,ctx$fits))
    } else logcorr[i]<-NA_real_
  }
  sec<-proc.time()[["elapsed"]]-started
  results[[nm]]<-data.frame(model=nm,eps=eps,seconds=sec,mean_acceptance=mean(acc),median_log_correction=median(logcorr,na.rm=TRUE),
    accepted_movement_proxy_per_second=mean(acc*jump)/sec,mean_exact_calls_per_second=xx$rounds*20/sec,
    log_correction_q10=quantile(logcorr,.1,na.rm=TRUE),log_correction_q90=quantile(logcorr,.9,na.rm=TRUE))
  cat(nm,":",paste(capture.output(print(results[[nm]],row.names=FALSE)),collapse=" "),"\n")
}
out<-do.call(rbind,results);write.csv(out,"/tmp/collapsed_ridge_screen_summary.csv",row.names=FALSE)
