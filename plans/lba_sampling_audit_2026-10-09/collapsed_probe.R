library(EMC2)
x<-EMC2:::restore_duplicates(readRDS('/tmp/lba_sampling_audit/snapshot.rds'));x<-EMC2:::restore_custom_kernel_pointers(x,quiet=TRUE)
set.seed(sum(utf8ToInt('lba-collapsed-audit-2026-10-09')))
records<-list(); contributions<-list()
for(j in seq_along(x)){
ch<-x[[j]];p<-ch$n_pars;ns<-ch$n_subjects;idx<-ch$samples$idx;ctx<-EMC2:::.emc_collapsed_context(ch)
a0<-ch$samples$alpha[,,idx];ll0<-ch$samples$subj_ll[,idx];mu0<-ch$samples$theta_mu[,idx];V0<-ch$samples$theta_var[,,idx]
cur<-EMC2:::.emc_shear_state(mu0,V0,ctx$fits,ch$prior,ctx$spec)
phi0<-EMC2:::.emc_collapsed_to_phi(mu0,V0,ctx$chart);lp0<-EMC2:::.emc_collapsed_surrogate(phi0,ctx)$lp
r<-sapply(seq_len(ns),function(s)drop(cur$frames[[s]]$R %*% (a0[,s]-cur$frames[[s]]$m)))
props<-list();new_states<-list();trs<-list();scales<-rep(c(1,.3,.1,.03),each=4)
for(k in seq_along(scales)){
tr<-EMC2:::collapsed_trajectories_cpp(phi0,ch$collapsed$plan$C*scales[k],ch$collapsed$plan$eps,3L,1.5,ctx$cpp)
st<-EMC2:::.emc_collapsed_from_phi(tr$phi,ctx$chart)
ne<-EMC2:::.emc_shear_state(st$mu,st$V,ctx$fits,ch$prior,ctx$spec)
props[[k]]<-sapply(seq_len(ns),function(s)ne$frames[[s]]$m+drop(backsolve(ne$frames[[s]]$R,r[,s])))
new_states[[k]]<-ne;trs[[k]]<-tr
}
saveRDS(list(chain=j,props=props,scales=scales,a0=a0),paste0('/tmp/lba_sampling_audit/proposals_chain',j,'.rds'))
ll<-sapply(seq_len(ns),function(s){A<-do.call(rbind,lapply(props,function(a)a[,s]));colnames(A)<-ch$par_names;EMC2:::calc_ll_manager(A,dadm=ch$data[[s]],model=ch$model)})
for(k in seq_along(scales)){
a<-props[[k]];ne<-new_states[[k]];tr<-trs[[k]]
lhat_delta<-vapply(seq_len(ns),function(s){fit<-ctx$fits[[s]];if(is.null(fit))return(0);d<-a[,s]-a0[,s];sum(d*(fit$lin-.5*drop(fit$prec %*% (a[,s]+a0[,s]))))},numeric(1))
delta<-ll[k,]-ll0-lhat_delta
la<-sum(delta)
# Independent full target/Jacobian/proposal correction.
chart_delta<-EMC2:::.emc_collapsed_from_phi(tr$phi,ctx$chart)$log_jac-EMC2:::.emc_collapsed_from_phi(phi0,ctx$chart)$log_jac
full<-sum(ll[k,]-ll0)+(ne$lp-cur$lp)+(EMC2:::.emc_shear_log_hier(a,ne)-EMC2:::.emc_shear_log_hier(a0,cur))+(cur$logdet-ne$logdet)+chart_delta-(tr$lp-lp0)
implemented<-sum(ll[k,]-EMC2:::.emc_collapsed_log_lhat(a,ctx$fits))-sum(ll0-EMC2:::.emc_collapsed_log_lhat(a0,ctx$fits))
records[[length(records)+1]]<-data.frame(chain=j,proposal=k,scale=scales[k],inner_acc=mean(tr$acc),surrogate_delta=tr$lp-lp0,true_LL_delta=sum(ll[k,]-ll0),approx_LL_delta=sum(lhat_delta),log_accept=la,accept=exp(min(0,la)),full_error=full-la,implemented_error=implemented-la,max_subject_delta=max(abs(delta)),worst_subject=dimnames(ch$samples$alpha)[[2]][which.min(delta)])
contributions[[length(contributions)+1]]<-data.frame(chain=j,proposal=k,scale=scales[k],subject=dimnames(ch$samples$alpha)[[2]],true_LL_delta=ll[k,]-ll0,approx_LL_delta=lhat_delta,correction=delta)
}
cat('chain',j,'finished\n');flush.console()
}
d<-do.call(rbind,records);con<-do.call(rbind,contributions);write.csv(d,'/tmp/lba_sampling_audit/collapsed_probes.csv',row.names=FALSE);write.csv(con,'/tmp/lba_sampling_audit/collapsed_subject_corrections.csv',row.names=FALSE)
print(d,row.names=FALSE,digits=4);cat('SCALE SUMMARY\n');print(aggregate(cbind(inner_acc,accept,log_accept)~scale,d,mean));cat('MAX RATIO ERRORS',max(abs(d$full_error)),max(abs(d$implemented_error)),'\n');cat('WORST SUBJECTS full scale\n');q<-aggregate(correction~subject,con[con$scale==1,],mean);print(head(q[order(q$correction),],15),row.names=FALSE)
