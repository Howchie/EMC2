# compare_versions.R
library(callr)
rm(list=ls())
dev_lib <- Sys.getenv("EMC2_DEV_LIB", "/tmp/r_libs/dev_version")
new_lib <- Sys.getenv("EMC2_NEW_LIB", "/tmp/r_libs/new_version")
upstream_lib <- Sys.getenv("EMC2_UPSTREAM_LIB", "/tmp/r_libs/upstream_version")
optimized_lib <- Sys.getenv("EMC2_OPTIMIZED_LIB", "/tmp/r_libs/optimized_version")
hybrid_lib <- Sys.getenv("EMC2_HYBRID_LIB", "/tmp/r_libs/hybrid_version")

benchmark_script_dev <- function() {
  library(EMC2)
  
  ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
  matchfun=function(d)d$S==d$lR
  
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
  dat$subjects <- droplevels(dat$subjects)
  
  design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                       formula=list(m~lM,s~1,t0~1),
                       contrasts=list(m=list(lM=ADmat)))
  LNR_s1 <- make_emc(dat, design_LNR, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector1 <- sampled_pars(LNR_s1)
  p_vector1[['s']] = 1; p_vector1[['m']] = 1; p_vector1[['m_lMd']] = .5; p_vector1[['t0']] = .2
  
  design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                       formula=list(v~lM,s~1,t0~1,A~1,B~1),
                       contrasts=list(v=list(lM=ADmat)))
  RDM_s <- make_emc(dat, design_RDM, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector2 <- sampled_pars(RDM_s)
  p_vector2[['s']] = 1; p_vector2[['v']] = 1; p_vector2[['v_lMd']] = .5; p_vector2[['t0']] = .3; p_vector2[['A']] = .5; p_vector2[['B']] = 0
  
  design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                       formula=list(v~lM,sv~1,t0~1,A~1,B~1),
                       contrasts=list(v=list(lM=ADmat)))
  LBA_s <- make_emc(dat, design_LBA, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector21 <- sampled_pars(LBA_s)
  p_vector21[['s']] = 1; p_vector21[['v']] = 1; p_vector21[['v_lMd']] = .5; p_vector21[['t0']] = .3; p_vector21[['A']] = .5; p_vector21[['B']] = 0
  
  design_WDM <- design(data = dat, model=DDM,
                       formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                       constants=c(s=1))
  WDM_s <- make_emc(dat, design_WDM, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector3 <- sampled_pars(WDM_s)
  p_vector3[1:length(p_vector3)] <- c(1, 1, .3, .5, 0, 0)
  
  design_DDM <- design(data = dat,model=DDM,
                       formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                       constants=c(s=1))
  DDM_s <- make_emc(dat, design_DDM, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector4 <- sampled_pars(DDM_s)
  p_vector4[1:length(p_vector4)] <- c(1, 1, .3, .5, .1, .1)
  
  designs_and_p_vectors <- list('LNR'=list(p_vector=p_vector1, emc=LNR_s1),
                                'RDM'=list(p_vector=p_vector2, emc=RDM_s),
                                'LBA'=list(p_vector=p_vector21, emc=LBA_s),
                                'WDM'=list(p_vector=p_vector3, emc=WDM_s),
                                'DDM'=list(p_vector=p_vector4, emc=DDM_s))

  results <- list()
  out_old <- list()
  for(i in 1:length(designs_and_p_vectors)) {
    nm <- names(designs_and_p_vectors)[i]
    emc <- designs_and_p_vectors[[nm]][['emc']]
    p_vector <- designs_and_p_vectors[[nm]][['p_vector']]
    model <- emc[[1]]$model()
    p_types <- names(model$p_types)
    dadm <- emc[[1]]$data[[1]]
  
    designs <- list()
    for(p in p_types){
      designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
    }
    constants <- attr(dadm, "constants")
    if(is.null(constants)) constants <- NA
  
    set.seed(1)
    p_mat <- matrix(rnorm(100*length(p_vector)), ncol=length(p_vector))
    colnames(p_mat) <- names(p_vector)
    if(nm == 'WDM') {
      p_mat[,'sv'] = log(0)
      p_mat[,'SZ'] = qnorm(0)
    }
    if(nm == 'RDM_noA') {
      constants = c('A'=log(0))
      c_name='RDM-A0'
      p_mat <- p_mat[,!colnames(p_mat) == c('A')]
    } else {
      c_name = model$c_name
    }
    if(nm %in% c('RDM', 'RDM_noA')) p_mat[,'s'] <- log(1)
    
    # Check if trend exists
    trend <- if (!is.null(model$trend)) model$trend else list()
    
    print(paste("c_name is:", c_name))
    out_old[[nm]] <- EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = c_name,
                                    model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                                    trend)
    bm <- microbenchmark::microbenchmark(
      f1 = EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = c_name,
                          model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                          trend),
      times=300, control=list(warmup=20)
    )
    
    results[[nm]] <- bm
  }
  return(list(results=results,lls=out_old))
}

benchmark_script_new <- function() {
  library(EMC2)
  
  ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
  matchfun=function(d)d$S==d$lR
  
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
  dat$subjects <- droplevels(dat$subjects)
  
  design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                       formula=list(m~lM,s~1,t0~1),
                       contrasts=list(m=list(lM=ADmat)))
  LNR_s <- make_emc(dat, design_LNR, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector1 <- sampled_pars(LNR_s)
  p_vector1[['s']] = 1; p_vector1[['m']] = 1; p_vector1[['m_lMd']] = .5; p_vector1[['t0']] = .2
  
  design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                       formula=list(v~lM,s~1,t0~1,A~1,B~1),
                       contrasts=list(v=list(lM=ADmat)))
  RDM_s <- make_emc(dat, design_RDM, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector2 <- sampled_pars(RDM_s)
  p_vector2[['s']] = 1; p_vector2[['v']] = 1; p_vector2[['v_lMd']] = .5; p_vector2[['t0']] = .3; p_vector2[['A']] = .5; p_vector2[['B']] = 0
  
  design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                       formula=list(v~lM,sv~1,t0~1,A~1,B~1),
                       contrasts=list(v=list(lM=ADmat)))
  LBA_s <- make_emc(dat, design_LBA, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector21 <- sampled_pars(LBA_s)
  p_vector21[['s']] = 1; p_vector21[['v']] = 1; p_vector21[['v_lMd']] = .5; p_vector21[['t0']] = .3; p_vector21[['A']] = .5; p_vector21[['B']] = 0
  
  design_WDM <- design(data = dat, model=DDM,
                       formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                       constants=c(s=1))
  WDM_s <- make_emc(dat, design_WDM, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector3 <- sampled_pars(WDM_s)
  p_vector3[1:length(p_vector3)] <- c(1, 1, .3, .5, 0, 0)
  
  design_DDM <- design(data = dat,model=DDM,
                       formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                       constants=c(s=1))
  DDM_s <- make_emc(dat, design_DDM, rt_resolution = 1/60, n_chains = 2, compress=FALSE)
  p_vector4 <- sampled_pars(DDM_s)
  p_vector4[1:length(p_vector4)] <- c(1, 1, .3, .5, .1, .1)
  
  designs_and_p_vectors <- list('LNR'=list(p_vector=p_vector1, emc=LNR_s),
                                'RDM'=list(p_vector=p_vector2, emc=RDM_s),
                                'LBA'=list(p_vector=p_vector21, emc=LBA_s),
                                'WDM'=list(p_vector=p_vector3, emc=WDM_s),
                                'DDM'=list(p_vector=p_vector4, emc=DDM_s))

  results <- list()
  out_old <- list()
  out_new <- list()
  for(i in 1:length(designs_and_p_vectors)) {
    nm <- names(designs_and_p_vectors)[i]
    emc <- designs_and_p_vectors[[nm]][['emc']]
    p_vector <- designs_and_p_vectors[[nm]][['p_vector']]
    model <- emc[[1]]$model()
    p_types <- names(model$p_types)
    dadm <- emc[[1]]$data[[1]]
  
    designs <- list()
    for(p in p_types){
      designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
    }
    constants <- attr(dadm, "constants")
    if(is.null(constants)) constants <- NA
  
    set.seed(1)
    p_mat <- matrix(rnorm(100*length(p_vector)), ncol=length(p_vector))
    colnames(p_mat) <- names(p_vector)
    if(nm == 'WDM') {
      p_mat[,'sv'] = log(0)
      p_mat[,'SZ'] = qnorm(0)
    }
    if(nm == 'RDM_noA') {
      constants = c('A'=log(0))
      c_name='RDM-A0'
      p_mat <- p_mat[,!colnames(p_mat) == c('A')]
    } else {
      c_name = model$c_name
    }
    if(nm %in% c('RDM', 'RDM_noA')) p_mat[,'s'] <- log(1)
    
    # Check if trend exists
    trend <- if (!is.null(model$trend)) model$trend else list()
    
    print(paste("c_name is:", c_name))
    
    # Check that outputs are identical before benchmarking
    out_old[[nm]] <- EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = c_name,
                        model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                        model$trend)
    out_new[[nm]] <- EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = c_name,
                        model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                        model$trend)
    print(paste("Outputs identical:", isTRUE(all.equal(out_old[[nm]], out_new[[nm]]))))
    
    bm <- microbenchmark::microbenchmark(
      f1_old_calc_ll=EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = c_name,
                        model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                        model$trend),
      f2_calc_ll_oo=EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = c_name,
                        model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                        model$trend),
      times=300, control=list(warmup=20)
    )
    
    results[[nm]] <- bm
  }
  return(list(results=results,lls=out_new))
}

cat("=========================================\n")
cat("Running DEV Branch (Original calc_ll)\n")
cat("=========================================\n")
dev_results <- callr::r(
  func = benchmark_script_dev,
  show = TRUE,
  libpath = c(dev_lib, .libPaths())
)

for (nm in names(dev_results$results)) {
  cat("\n---", nm, "---\n")
  print(dev_results[[nm]])
}

cat("\n=========================================\n")
cat("Running DEV-OO Branch (New calc_ll and calc_ll_oo)\n")
cat("=========================================\n")
new_results <- callr::r(
  func = benchmark_script_new,
  libpath = c(new_lib, .libPaths())
)

for (nm in names(new_results$results)) {
  cat("\n---", nm, "---\n")
  print(new_results$results[[nm]])
}

cat("\n=========================================\n")
cat("Running UPSTREAM OO Branch (oo_refactor_simd)\n")
cat("=========================================\n")
upstream_results <- callr::r(
  func = benchmark_script_new,
  show = TRUE,
  libpath = c(upstream_lib, .libPaths())
)

for (nm in names(upstream_results$results)) {
  cat("\n---", nm, "---\n")
  print(upstream_results$results[[nm]])
}

cat("\n=========================================\n")
cat("Running OPTIMIZED Branch (Newest build)\n")
cat("=========================================\n")
optimized_results <- callr::r(
  func = benchmark_script_new,
  show = TRUE,
  libpath = c(optimized_lib, .libPaths())
)

for (nm in names(optimized_results$results)) {
  cat("\n---", nm, "---\n")
  print(optimized_results$results[[nm]])
}

cat("\n=========================================\n")
cat("Running HYBRID Branch (Newest build)\n")
cat("=========================================\n")
hybrid_results <- callr::r(
  func = benchmark_script_new,
  show = TRUE,
  libpath = c(hybrid_lib, .libPaths())
)

for (nm in names(hybrid_results$results)) {
  cat("\n---", nm, "---\n")
  print(hybrid_results$results[[nm]])
}

res = matrix(NA,ncol=5,nrow=6,dimnames=list(c("old_dev","oo_calc_ll","upstream_calc_ll_oo","zach_calc_ll_oo","optimized_calc_ll_oo","hybrid_calc_ll_oo"),c("LBA","RDM","LNR","WDM","DDM")))
for (m in c("LBA","RDM","LNR","WDM","DDM")) {
    res[1,m] = median(dev_results$results[[m]]$time[dev_results$results[[m]]$expr=="f1"]) * 1e-6
    res[2,m] = median(upstream_results$results[[m]]$time[upstream_results$results[[m]]$expr=="f1_old_calc_ll"]) * 1e-6
    res[3,m] = median(upstream_results$results[[m]]$time[upstream_results$results[[m]]$expr=="f2_calc_ll_oo"])* 1e-6
    res[4,m] = median(new_results$results[[m]]$time[new_results$results[[m]]$expr=="f2_calc_ll_oo"]) * 1e-6  
    res[5,m] = median(optimized_results$results[[m]]$time[optimized_results$results[[m]]$expr=="f2_calc_ll_oo"]) * 1e-6
    res[6,m] = median(hybrid_results$results[[m]]$time[hybrid_results$results[[m]]$expr=="f2_calc_ll_oo"]) * 1e-6
}
print(res)

## Compare numerical accuracy
res2 = matrix(NA,ncol=5,nrow=4,dimnames=list(c("upstream vs dev","zach vs dev", "optimized vs dev", "hybrid vs dev"),
                                             c("LBA","RDM","LNR","WDM","DDM")))
for (m in names(optimized_results$results)) {
  res2[1,m] = max(abs(upstream_results$lls[[m]]-dev_results$lls[[m]]))
  res2[2,m] = max(abs(new_results$lls[[m]]-dev_results$lls[[m]]))
  res2[3,m] = max(abs(optimized_results$lls[[m]]-dev_results$lls[[m]]))
  res2[4,m] = max(abs(hybrid_results$lls[[m]]-dev_results$lls[[m]]))
}
