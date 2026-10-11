library(EMC2)
e <- new.env(); load('/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_lba7r_sepA_sepL_sepV_t0T_mt_v2.RData', e)
x <- e$emc
saveRDS(x,'/tmp/lba_sampling_audit/snapshot.rds')
pn <- x[[1]]$par_names; ns <- x[[1]]$n_subjects; n <- x[[1]]$samples$idx
cat('snapshot draws',n,'transitions',x[[1]]$sample_transitions,'\n'); print(x[[1]]$production_options)
old <- function(z) EMC2:::gelman_diag_robust(coda::mcmc.list(lapply(seq_len(ncol(z)),function(i)coda::mcmc(matrix(z[,i],ncol=1)))),omit_mpsrf=TRUE)[1]
new <- EMC2:::rhat_new
rows <- list()
for(p in seq_along(pn)) for(s in 0:ns){
 z<-sapply(x,function(ch) if(s==0) ch$samples$theta_mu[p,seq_len(n)] else ch$samples$alpha[p,s,seq_len(n)])
 rows[[length(rows)+1]]<-data.frame(parameter=pn[p],subject=if(s==0)'mu' else as.character(x[[1]]$subjects[s]),old=old(z),new=new(z),move_min=min(colMeans(z[-1,,drop=FALSE]!=z[-n,,drop=FALSE])),acf_max=max(apply(z,2,function(y)acf(y,lag.max=1,plot=FALSE)$acf[2])))
}
d<-do.call(rbind,rows);write.csv(d,'/tmp/lba_sampling_audit/diagnostics.csv',row.names=FALSE)
for(group in c('mu','alpha')){ q<-d[if(group=='mu')d$subject=='mu' else d$subject!='mu',];cat('\nTOP',group,'\n');print(head(q[order(-q$old),],12),row.names=FALSE)}
top<-head(d[d$subject!='mu',][order(-d$old[d$subject!='mu']),],4)
for(i in seq_len(nrow(top))){p<-match(top$parameter[i],pn);s<-match(top$subject[i],as.character(x[[1]]$subjects));z<-sapply(x,function(ch)ch$samples$alpha[p,s,seq_len(n)]);cat('\nSUBJECT',top$subject[i],pn[p],'\n');print(rbind(mean=colMeans(z),sd=apply(z,2,sd),min=apply(z,2,min),max=apply(z,2,max),first100=colMeans(head(z,100)),last100=colMeans(tail(z,100)))); for(k in c(200,400,600,n)){zz<-head(z,k);cat('prefix',k,'old',old(zz),'new',new(zz),'\n')}; for(k in c(100,200,400)){zz<-tail(z,k);cat('tail',k,'old',old(zz),'new',new(zz),'\n')}}
cat('\nMU CHAIN MOMENTS\n');for(p in seq_along(pn)){z<-sapply(x,function(ch)ch$samples$theta_mu[p,seq_len(n)]);cat(pn[p],'mean',round(colMeans(z),3),'sd',round(apply(z,2,sd),3),'first-last',round(colMeans(tail(z,100))-colMeans(head(z,100)),3),'\n')}
for(j in seq_along(x)){ch<-x[[j]];cat('\nMOVES CHAIN',j,'\n');print(ch$alloc$last[c('rounds','collapsed_rounds','collapsed_g','collapsed_a','collapsed_cost','hmc_rates')]);str(ch$subject_hmc,max.level=2,list.len=5);cat('rng unique',length(unique(vapply(ch$rng$subjects,paste,collapse=',',FUN.VALUE=''))),'\n');cat('LL chain sums range ',range(colSums(ch$samples$subj_ll[,seq_len(n)])),'\n');cat('LL first-last ',mean(colSums(ch$samples$subj_ll[,tail(seq_len(n),100)]))-mean(colSums(ch$samples$subj_ll[,head(seq_len(n),100)])),'\n')}
pdf('/tmp/lba_sampling_audit/traces.pdf',width=11,height=8.5)
par(mfrow=c(3,2),mar=c(3,4,3,1))
for(i in seq_len(nrow(top))){p<-match(top$parameter[i],pn);s<-match(top$subject[i],as.character(x[[1]]$subjects));z<-sapply(x,function(ch)ch$samples$alpha[p,s,seq_len(n)]);matplot(z,type='l',lty=1,col=c('#0072B2','#D55E00','#009E73'),main=paste('Subject',top$subject[i],pn[p]),xlab='Retained draw',ylab='Untransformed coefficient')}
for(p in match(c('A','B'),pn)){z<-sapply(x,function(ch)ch$samples$theta_mu[p,seq_len(n)]);matplot(z,type='l',lty=1,col=c('#0072B2','#D55E00','#009E73'),main=paste('Group mean',pn[p]),xlab='Retained draw',ylab='Untransformed coefficient')}
dev.off()
