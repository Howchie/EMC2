pkgload::load_all('/data/work/EMC2_dev_oo',quiet=TRUE,compile=FALSE)
x<-EMC2:::restore_custom_kernel_pointers(readRDS('/tmp/lba_sampling_audit/calibrated_snapshot.rds'),quiet=TRUE)
options(emc2.subject_hmc=FALSE)
# Exercise calibration through allocation with likelihood workers, including a
# previously disabled move and the persisted Gibbs stream.
x<-x[1];class(x)<-'emc';x[[1]]$collapsed$plan$calibrate<-TRUE;x[[1]]$collapsed$plan$rounds<-0
before<-x[[1]]$samples
seed<-x[[1]]$rng$gibbs
out<-EMC2:::.emc_alloc_rebuild(x,'sample',n_cores=2L,B=100L)
stopifnot(identical(before,out[[1]]$samples),!identical(seed,out[[1]]$rng$gibbs),!out[[1]]$collapsed$plan$calibrate)
print(out[[1]]$alloc$last[c('collapsed_rounds','collapsed_scale','collapsed_cost','collapsed_calibration')])
# One retained transition uses the calibrated transport and preserves its scale.
ch<-out[[1]];i<-ch$samples$idx
ch$collapsed$plan$rounds<-1
pars<-list(tmu=ch$samples$theta_mu[,i],tvar=ch$samples$theta_var[,,i],tvinv=solve(ch$samples$theta_var[,,i]))
res<-EMC2:::.emc_collapsed_move(ch,pars,rbind(ch$samples$alpha[,,i],ch$samples$subj_ll[,i]),'sample')
stopifnot(identical(res$sampler$collapsed$plan$scale,ch$collapsed$plan$scale))
fresh<-vapply(seq_len(ch$n_subjects),function(s)EMC2:::calc_ll_manager(matrix(res$proposals[seq_len(ch$n_pars),s],1,dimnames=list(NULL,ch$par_names)),ch$data[[s]],ch$model),0)
cat('post-transition max likelihood error',max(abs(fresh-res$proposals[ch$n_pars+1,])),'\n')
stopifnot(max(abs(fresh-res$proposals[ch$n_pars+1,]))<1e-8)
