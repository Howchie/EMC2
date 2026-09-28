# Execute the actual package transition on an analytically solvable hierarchy.
suppressPackageStartupMessages(pkgload::load_all(quiet=TRUE))
if (!requireNamespace("posterior", quietly=TRUE)) stop("Install posterior for benchmark diagnostics")
args <- commandArgs(TRUE)
suite <- if(length(args)) args[1] else "development"
stopifnot(suite %in% c("development","holdout","smoke"))
preburn <- if(suite=="smoke") 25L else 500L
outdir <- file.path("benchmarks/group-move/results",suite)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
fingerprint <- tools::md5sum(c(list.files("R",full.names=TRUE),
                              "benchmarks/group-move/run.R"))
ll <- function(proposals, dadm, model, ...) {
  z <- sweep(proposals,2,dadm$y)
  -.5*rowSums((z %*% dadm$prec)*z)
}
run_case <- function(p, obs_max, angle, seed, method, warm=6000L, settle=500L, draws=4000L, chains=4L) {
  set.seed(seed)
  n <- 8L
  Q <- if(p==2) matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2) else qr.Q(qr(matrix(rnorm(p*p),p)))
  Obs <- Q %*% diag(exp(seq(log(.01),log(obs_max),length.out=p)),p) %*% t(Q)
  Sigma <- diag(.01,p); Si <- diag(100,p); Oi <- solve(Obs)
  reference <- solve(diag(.01,p)+n*solve(Sigma+Obs))
  eig <- eigen(reference,symmetric=TRUE)
  Vmu <- solve(diag(.01,p)+n*Si); Lmu <- t(chol(Vmu))
  Va <- solve(Si+Oi); La <- t(chol(Va))
  proposal_count <- if (method == "legacy") 4L else 1L
  options(emc2.group_move_method=if(method %in% c("disabled","oracle")) "am" else method,
          emc2.group_move_ram_batch=as.integer(Sys.getenv("RAM_BATCH","16")),
          emc2.group_move=method!="disabled",emc2.group_move_proposals=proposal_count,
          emc2.group_move_scale=FALSE,emc2.group_move_warmup=if(method=="ram") NULL else warm,
          emc2.group_move_settle=if(method=="ram") NULL else settle)
  arr <- array(NA_real_,c(draws,chains,p))
  history <- array(NA_real_,c(p,chains,preburn+warm+settle+draws))
  states <- vector("list",chains)
  started <- proc.time()[["elapsed"]]
  for(ch in seq_len(chains)) {
    set.seed(seed+ch)
    sampler <- list(type="standard",n_pars=p,n_subjects=n,par_names=paste0("p",seq_len(p)),
      nuisance=rep(FALSE,p),marginalised_idx=rep(FALSE,p),is_blocked=rep(TRUE,p),par_group=rep(1,p),
      prior=list(theta_mu_mean=rep(0,p),theta_mu_invar=diag(.01,p)),
      data=rep(list(list(y=rep(0,p),prec=Oi)),n))
    mu <- as.numeric(eig$vectors %*% (rep(c(-2,2),length.out=p)*sqrt(eig$values)*if(ch%%2) 1 else -1))
    if(method=="oracle") {
      sampler$group_move_config <- EMC2:::.emc_group_move_config()
      spec <- EMC2:::.emc_group_move_spec(sampler, scale=FALSE)
      gm <- EMC2:::.emc_group_move_initialize(
        sampler, spec, list(tvar=Sigma), sampler$group_move_config)
      oracle_covariance <- solve(diag(.01,p) + n*Oi)
      oracle_factor <- t(chol(oracle_covariance)) * (2.38/sqrt(p))
      gm$baseline <- oracle_covariance
      gm$shape_factor <- oracle_factor
      gm$factor <- oracle_factor
      gm$status <- "ready"
      sampler$group_move <- gm
    }
    states[[ch]] <- list(sampler=sampler,alpha=matrix(mu,p,n),rng=.Random.seed)
  }
  # Lockstep iterations reproduce the package's cross-chain legacy adaptation.
  # New adapters still own independent state; no other chain enters their update.
  for(t in seq_len(dim(history)[3])) {
    for(ch in seq_len(chains)) {
      state <- states[[ch]]
      assign(".Random.seed",state$rng,globalenv())
      sampler <- state$sampler; alpha <- state$alpha
      # These are the original sampler's exact group-mean and subject Gibbs
      # updates. The "disabled" arm below skips only the added group proposal.
      mu <- as.numeric(Vmu %*% (Si %*% rowSums(alpha)) + Lmu %*% rnorm(p))
      alpha <- matrix(as.numeric(Va %*% Si %*% mu),p,n) + La %*% matrix(rnorm(p*n),p,n)
      pars <- list(tmu=mu,tvar=Sigma,tvinv=Si)
      old_ll <- -.5*colSums(alpha*(Oi %*% alpha))
      sampler$rng$gibbs <- .Random.seed
      stage <- if(method=="oracle") "sample" else
        if(t<=preburn) "preburn" else
          if(t<=preburn+warm+settle) "adapt" else "sample"
      step <- EMC2:::.emc_group_move(sampler,pars,rbind(alpha,old_ll),stage,
                NULL,NULL,NULL)
      sampler <- step$sampler; mu <- step$pars$tmu; alpha <- step$proposals[seq_len(p),,drop=FALSE]
      states[[ch]] <- list(sampler=sampler,alpha=alpha,rng=.Random.seed)
      history[,ch,t] <- mu
      if(t>preburn+warm+settle)
        arr[t-preburn-warm-settle,ch,] <-
          as.numeric(t(eig$vectors) %*% mu)/sqrt(eig$values)
    }
    if(method=="legacy" && t<=preburn+warm && t%%100L==0L) {
      ix <- max(1,t-249L):t
      emc <- lapply(seq_len(chains),function(ch) {
        s <- states[[ch]]$sampler
        s$samples <- list(theta_mu=matrix(history[,ch,ix],p),
                         theta_var=array(rep(Sigma,length(ix)),c(p,p,length(ix))))
        s
      })
      spec <- EMC2:::.emc_group_move_spec(emc[[1]],FALSE)
      V <- EMC2:::.emc_group_move_covariance(emc,spec,seq_along(ix))
      if(!is.null(V)) for(ch in seq_len(chains)) {
        states[[ch]]$sampler <- EMC2:::.emc_group_move_install(
          states[[ch]]$sampler,V,spec,proposal_count)
      }
    }
  }
  factors <- lapply(states,function(s) if(method %in% c("am","ram","oracle"))
    s$sampler$group_move$factor else s$sampler$group_move_cov)
  elapsed <- proc.time()[["elapsed"]]-started
  metrics <- lapply(seq_len(p),function(j) {
    z <- arr[,,j]
    data.frame(suite=suite,p=p,obs_max=obs_max,angle=angle,seed=seed,method=method,
      direction=j,mean=mean(z),variance=var(as.numeric(z)),rhat=posterior::rhat(z),
      bulk_ess=posterior::ess_bulk(z),tail_ess=posterior::ess_tail(z),
      elapsed=elapsed,likelihood_evals=if(method=="disabled") 0 else
        (preburn+warm+settle+draws)*chains*n*proposal_count)
  })
  diagnostics <- lapply(states,function(s) {
    gm <- s$sampler$group_move
    if(is.null(gm)) return(NULL)
    gm[c("baseline","n_adapt","status","factor","attempts_warmup",
         "accepts_warmup","sum_accept_prob_warmup")]
  })
  list(metrics=do.call(rbind,metrics),factors=factors,diagnostics=diagnostics,
       draws=arr,reference=reference)
}
cases <- if(suite=="holdout") list(c(2,3,.31),c(8,100,.31)) else
  list(c(2,.01,.7853981634),c(2,1,.7853981634),c(2,100,.7853981634))
seeds <- if(suite=="holdout") c(44001,54001,64001) else
  c(14001,24001,34001)
if(suite=="smoke") {cases <- cases[3]; seeds <- seeds[1]}
methods <- c("disabled","legacy","am","ram","oracle")
if(nzchar(Sys.getenv("GM_METHODS"))) methods <- strsplit(Sys.getenv("GM_METHODS"),",")[[1]]
warmups <- if(suite=="smoke") 500L else 6000L
settling <- if(suite=="smoke") 50L else 500L
production <- if(suite=="smoke") 500L else 4000L
n_chains <- if(suite=="smoke") 2L else 4L
rows <- list()
for(cs in cases) for(seed in seeds) for(warm in warmups) for(method in methods) {
  name <- sprintf("p%d-v%g-seed%d-%s-warm%d",cs[1],cs[2],seed,method,warm)
  file <- file.path(outdir,paste0(name,".rds"))
  result <- if(file.exists(file)) readRDS(file) else NULL
  if(is.null(result) || !identical(result$code_fingerprint,fingerprint)) {
    result <- testthat::with_mocked_bindings(
      run_case(as.integer(cs[1]),cs[2],cs[3],seed,method,warm=warm,
               settle=settling,draws=production,chains=n_chains),
      calc_ll_manager=ll,.package="EMC2")
    result$code_fingerprint <- fingerprint
    saveRDS(result,file)
  }
  rows[[length(rows)+1L]] <- result$metrics
  write.csv(do.call(rbind,rows),file.path(outdir,"metrics.csv"),row.names=FALSE)
  cat(name, "rhat",round(max(result$metrics$rhat),3),"min ESS",round(min(result$metrics$bulk_ess)),
      "variance",paste(round(result$metrics$variance,2),collapse=","),"seconds",round(result$metrics$elapsed[1],1),"\n")
  flush.console()
}
writeLines(c(capture.output(sessionInfo()),capture.output(system("git rev-parse HEAD",intern=TRUE))),file.path(outdir,"session.txt"))
