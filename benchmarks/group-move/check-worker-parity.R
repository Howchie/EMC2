# Investigative whole-fit diagnostic, independent of the analytic mixing suite.
# Run with an isolated installed EMC2 on R_LIBS to exercise spawned workers.
suppressPackageStartupMessages(library(EMC2))
stopifnot(EMC2:::.emc_wpool_backend(list())=="spawn")
options(emc2.group_move_method="am",emc2.group_move_proposals=1L,
        emc2.group_move_warmup=20L,emc2.group_move_settle=5L)
data(forstmann)
ids <- unique(forstmann$subjects)[1:2]
dat <- droplevels(do.call(rbind,lapply(ids,function(id) head(forstmann[forstmann$subjects==id,],30))))
des <- design(data=dat,model=LNR,formula=list(m~1,s~1),constants=c(t0=log(.2)))
set.seed(8931)
x <- make_emc(dat,des,type="diagonal-gamma",compress=FALSE,n_chains=2L)
for(i in seq_along(x)) x[[i]] <- EMC2:::run_stages(x[[i]],stage="preburn",iter=5,
             verbose=FALSE,verboseProgress=FALSE,particle_factor=10,n_cores=1)
seed <- .Random.seed
rows <- list()
for(method in c("am","disabled")) {
  initial <- x
  if(method=="disabled") for(i in seq_along(initial)) {
    initial[[i]]$group_move_config <- list(method="legacy",enabled=FALSE,K=1L,scale=TRUE)
    initial[[i]]["group_move"] <- list(NULL)
  }
  run <- function(cores) fit(initial,cores_for_chains=1L,cores_per_chain=cores,
    stop_criteria=list(preburn=list(iter=10L),burn=list(iter=15L),adapt=list(iter=15L,min_unique=0L),sample=list(iter=20L)),
    verbose=FALSE,particle_factor=10,step_size=10L,max_tries=5L,trim=FALSE)
  assign(".Random.seed",seed,globalenv()); a <- EMC2:::restore_duplicates(run(1L))
  assign(".Random.seed",seed,globalenv()); b <- EMC2:::restore_duplicates(run(2L))
  for(i in seq_along(a)) {
    delta <- colSums(abs(a[[i]]$samples$theta_mu-b[[i]]$samples$theta_mu))
    first <- which(delta>1e-12)[1]
    rows[[length(rows)+1L]] <- data.frame(method=method,chain=i,first_difference=first,
      stage=if(is.na(first)) NA_character_ else a[[i]]$samples$stage[first],max_abs_difference=max(delta),
      identical_draws=identical(a[[i]]$samples,b[[i]]$samples),
      identical_rng=identical(a[[i]]$rng,b[[i]]$rng))
  }
}
print(do.call(rbind,rows))
write.csv(do.call(rbind,rows),"benchmarks/group-move/results/worker-parity.csv",row.names=FALSE)
stopifnot(all(vapply(rows,function(z) z$identical_draws && z$identical_rng,logical(1))))
