pkgload::load_all('/data/work/EMC2_dev_oo',quiet=TRUE,compile=FALSE)
x<-EMC2:::restore_custom_kernel_pointers(EMC2:::restore_duplicates(readRDS('/tmp/lba_sampling_audit/snapshot.rds')),quiet=TRUE)
results<-list()
for(j in seq_along(x)){
 ch<-x[[j]];before<-ch$samples
 ch<-EMC2:::.emc_collapsed_set_plan(ch,ch$collapsed$plan$C,'sample')
 ch<-EMC2:::.emc_collapsed_calibrate(ch,n_cores=1)
 stopifnot(identical(before,ch$samples))
 cat('CHAIN',j,'scale',ch$collapsed$plan$scale,'validation acceptance',ch$collapsed$plan$calibration$acceptance,'\n')
 print(ch$collapsed$plan$calibration);flush.console()
 ctx<-EMC2:::.emc_collapsed_context(ch);ii<-ch$samples$idx;p<-ch$n_pars
 cur<-EMC2:::.emc_shear_state(ch$samples$theta_mu[,ii],ch$samples$theta_var[,,ii],ctx$fits,ch$prior,ctx$spec)
 base<-list(cur=cur,alpha=matrix(ch$samples$alpha[,,ii],p),ll=ch$samples$subj_ll[,ii],rounds=0L,wpool=NULL,wpool_part=NULL)
 set.seed(sum(utf8ToInt('collapsed-transport-heldout'))+j)
 acc<-numeric(20);score<-numeric(20)
 t0<-proc.time()[['elapsed']]
 for(k in seq_along(acc)){o<-EMC2:::.emc_collapsed_round(base,ch$collapsed$plan,ch,ctx);acc[k]<-o$acc;score[k]<-o$acc*sum(forwardsolve(ch$collapsed$plan$C,o$phi_jump)^2)}
 results[[j]]<-list(chain=j,scale=ch$collapsed$plan$scale,calibration=ch$collapsed$plan$calibration,acceptance=mean(acc),esjd=mean(score),cost=(proc.time()[['elapsed']]-t0)/length(acc))
 print(results[[j]][c('chain','scale','acceptance','esjd','cost')]);flush.console()
 x[[j]]<-ch
}
saveRDS(results,'/tmp/lba_sampling_audit/adaptation_validation.rds')
saveRDS(x,'/tmp/lba_sampling_audit/calibrated_snapshot.rds')
