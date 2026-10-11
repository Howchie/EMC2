pkgload::load_all(".", quiet=TRUE, compile=FALSE)
chains <- EMC2:::restore_custom_kernel_pointers(EMC2:::restore_duplicates(readRDS("/tmp/lba_sampling_audit/snapshot.rds")), quiet=TRUE)
stopifnot(length(chains)==3L)
p <- chains[[1]]$n_pars
Tn <- dim(chains[[1]]$samples$alpha)[3]
fold_metrics <- list(); k <- 0L
fit_ridge <- function(A,y,lambda) {
  ok <- is.finite(y) & stats::complete.cases(A); A<-A[ok,,drop=FALSE]; y<-y[ok]
  p<-ncol(A); iu<-which(upper.tri(diag(p),diag=TRUE),arr.ind=TRUE)
  abar<-colMeans(A); e<-eigen(stats::cov(A),symmetric=TRUE)
  if(any(!is.finite(e$values)) || any(e$values<=1e-12*max(e$values))) return(NULL)
  W<-t(e$vectors)/sqrt(e$values); Z<-sweep(A,2,abar)%*%t(W)
  Q<-Z[,iu[,1]]*Z[,iu[,2],drop=FALSE]
  X<-cbind(1,Z,Q)
  scl<-sqrt(colMeans(X^2)); scl[1]<-1; Xs<-sweep(X,2,scl,"/")
  XtX<-crossprod(Xs); penalty<-diag(ncol(Xs)); penalty[1,1]<-0
  b<-tryCatch(solve(XtX + lambda*nrow(Xs)*penalty, crossprod(Xs,y)),error=function(e)NULL)
  if(is.null(b)) return(NULL); b<-drop(b)/scl
  cq<-b[-seq_len(p+1L)]; B<-matrix(0,p,p); B[iu]<- -cq; B<-B+t(B); diag(B)<- -2*cq[iu[,1]==iu[,2]]
  eb<-eigen(B,symmetric=TRUE); Bp<-eb$vectors%*%(pmax(eb$values,0)*t(eb$vectors))
  g<-qr.coef(qr(cbind(1,Z)), y + 0.5*rowSums((Z%*%Bp)*Z))[-1L]
  tol<-max(abs(eb$values),.Machine$double.xmin)*1e-10; keep<-which(eb$values>tol)
  if(length(keep)) g<-drop(eb$vectors[,keep,drop=FALSE]%*%crossprod(eb$vectors[,keep,drop=FALSE],g)) else g[]<-0
  H<-t(W)%*%Bp%*%W; H<-(H+t(H))/2; lin<-drop(H%*%abar+t(W)%*%g)
  list(prec=H,lin=lin)
}
predict_fit <- function(fit,A,train_y) {
  if(is.null(fit)) return(NULL)
  base <- -0.5*rowSums((A%*%fit$prec)*A) + drop(A%*%fit$lin)
  # Fit only the additive constant on training data.
  # Caller passes training predictions separately through an attribute.
  base
}
for (held in 1:3) {
  tr <- setdiff(1:3,held)
  for (s in seq_len(chains[[1]]$n_subjects)) {
    Atrain <- do.call(rbind,lapply(tr,function(j) t(matrix(chains[[j]]$samples$alpha[,s,seq_len(Tn)],nrow=p))))
    ytrain <- unlist(lapply(tr,function(j) chains[[j]]$samples$subj_ll[s,seq_len(Tn)]),use.names=FALSE)
    Atest <- t(matrix(chains[[held]]$samples$alpha[,s,seq_len(Tn)],nrow=p))
    ytest <- as.numeric(chains[[held]]$samples$subj_ll[s,seq_len(Tn)])
    fits <- list(full=EMC2:::.emc_collapsed_fit_subject(Atrain,ytrain))
    for(lam in c(0.001,0.01,0.1,1)) fits[[paste0("ridge",lam)]]<-fit_ridge(Atrain,ytrain,lam)
    for(nm in names(fits)) {
      f<-fits[[nm]]
      if(is.null(f)) next
      f<-EMC2:::.emc_collapsed_safe_fit(f)
      rawtr<- -0.5*rowSums((Atrain%*%f$prec)*Atrain)+drop(Atrain%*%f$lin)
      rawte<- -0.5*rowSums((Atest%*%f$prec)*Atest)+drop(Atest%*%f$lin)
      c0<-mean(ytrain-rawtr)
      err<-ytest-(rawte+c0)
      k<-k+1L
      fold_metrics[[k]]<-data.frame(held=held,subject=s,model=nm,
        train_rmse=sqrt(mean((ytrain-(rawtr+c0))^2)),
        test_rmse=sqrt(mean(err^2)),test_mae=mean(abs(err)),
        test_sd_error=sd(err),test_r2=1-sum(err^2)/sum((ytest-mean(ytrain))^2))
    }
  }
  cat("finished heldout chain",held,"\n")
}
d<-do.call(rbind,fold_metrics)
# Summarize per subject first so long/high-variance subjects do not dominate.
by_subject<-aggregate(cbind(train_rmse,test_rmse,test_mae,test_sd_error,test_r2)~model+subject,d,FUN=mean)
summary<-do.call(rbind,lapply(split(by_subject,by_subject$model),function(z)data.frame(model=z$model[1],
  subjects=nrow(z),median_test_rmse=median(z$test_rmse),q90_test_rmse=quantile(z$test_rmse,.9),
  median_test_mae=median(z$test_mae),median_test_r2=median(z$test_r2),
  median_train_rmse=median(z$train_rmse),test_wins_vs_full=NA_integer_)))
base<-by_subject[by_subject$model=="full",c("subject","test_rmse")]; names(base)[2]<-"base_rmse"
by_subject<-merge(by_subject,base,by="subject",suffixes=c("","_base"))
win<-aggregate(I(test_rmse<base_rmse)~model,by_subject,FUN=sum); names(win)[2]<-"win_subjects"
summary<-merge(summary,win,by="model",all.x=TRUE)
summary<-summary[order(summary$median_test_rmse),]
write.csv(summary,"/tmp/collapsed_surrogate_cv_summary.csv",row.names=FALSE)
write.csv(by_subject,"/tmp/collapsed_surrogate_cv_subject.csv",row.names=FALSE)
print(summary,row.names=FALSE)
