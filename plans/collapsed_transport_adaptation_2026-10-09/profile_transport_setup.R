pkgload::load_all('.',quiet=TRUE,compile=FALSE)
x<-EMC2:::restore_custom_kernel_pointers(EMC2:::restore_duplicates(readRDS('/tmp/lba_sampling_audit/calibrated_snapshot.rds')),quiet=TRUE)
ch<-x[[1]]; i<-ch$samples$idx; ch$collapsed$plan$rounds<-1
pars<-list(tmu=ch$samples$theta_mu[,i],tvar=ch$samples$theta_var[,,i],tvinv=solve(ch$samples$theta_var[,,i]))
proposals<-rbind(ch$samples$alpha[,,i],ch$samples$subj_ll[,i])
ctx<-EMC2:::.emc_collapsed_context(ch)
phi<-EMC2:::.emc_collapsed_to_phi(pars$tmu,pars$tvar,ctx$chart)
print(c(p=ch$n_pars,n=ch$n_subjects,dim=ctx$chart$dim,eps=ch$collapsed$plan$eps,scale=ch$collapsed$plan$scale))
bench<-function(f,n=20L){t<-proc.time();for(k in seq_len(n))f();(proc.time()-t)/n}
print(bench(function() EMC2:::.emc_collapsed_context(ch)))
print(bench(function() EMC2:::.emc_shear_state(pars$tmu,pars$tvar,ctx$fits,ch$prior,ctx$spec)))
print(bench(function() EMC2:::collapsed_trajectories_cpp(phi,ch$collapsed$plan$C*ch$collapsed$plan$scale,ch$collapsed$plan$eps,3L,1.5,ctx$cpp)))
print(bench(function() EMC2:::.emc_collapsed_move(ch,pars,proposals,'sample')))
