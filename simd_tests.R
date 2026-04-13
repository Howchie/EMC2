#remotes::install_github("Howchie/EMC2", ref="dev-oo")

rm(list=ls())
## messing around with SIMD
library(EMC2)


# Simplest design, no trend -----------------------------------------------
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR

# Drop most subjects for speed
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)

design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                     formula=list(m~lM,s~1,t0~1),
                     # transform=list(func=c('m'='identity', 's'='identity')),
                     contrasts=list(m=list(lM=ADmat)))
LNR_s1 <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, compress=FALSE)

p_vector1 <- sampled_pars(LNR_s1)
p_vector1[['s']] = 1
p_vector1[['m']] = 1
p_vector1[['m_lMd']] = .5
p_vector1[['t0']] = .2


# RDM no trend ----------------------------------------------------------
design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                     formula=list(v~lM,s~1,t0~1,A~1,B~1),
                     # transform=list(func=c('v'='identity', 's'='identity')),
                     contrasts=list(v=list(lM=ADmat)))
RDM_s <- make_emc(dat, design_RDM, rt_resolution = 0.05, n_chains = 2, compress=FALSE)

p_vector2 <- sampled_pars(RDM_s)
p_vector2[['s']] = 1
p_vector2[['v']] = 1
p_vector2[['v_lMd']] = .5
p_vector2[['t0']] = .3
p_vector2[['A']] = .5
p_vector2[['B']] = 0

# LBA no trend ----------------------------------------------------------
design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                     formula=list(v~lM,sv~1,t0~1,A~1,B~1),
                     # transform=list(func=c('v'='identity', 's'='identity')),
                     contrasts=list(v=list(lM=ADmat)))
LBA_s <- make_emc(dat, design_LBA, rt_resolution = 0.05, n_chains = 2, compress=FALSE)

p_vector21 <- sampled_pars(LBA_s)
p_vector21[['s']] = 1
p_vector21[['v']] = 1
p_vector21[['v_lMd']] = .5
p_vector21[['t0']] = .3
p_vector21[['A']] = .5
p_vector21[['B']] = 0

# WDM no trend ----------------------------------------------------------
design_WDM <- design(data = dat, model=DDM,
                     formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                     # transform=list(func=c('a'='identity', 's'='identity',
                     #                       'Z'='identity', 'sv'='identity', 'SZ'='identity')),
                     constants=c(s=1))
WDM_s <- make_emc(dat, design_WDM, rt_resolution = 0.05, n_chains = 2)
p_vector3 <- sampled_pars(WDM_s)
p_vector3[1:length(p_vector3)] <- c(1, 1, .3, .5, 0, 0)

# DDM ---
design_DDM <- design(data = dat,model=DDM,
                     formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                     # transform=list(func=c('a'='identity', 's'='identity', 'Z'='identity', 'sv'='identity', 'SZ'='identity')),
                     constants=c(s=1))
DDM_s <- make_emc(dat, design_DDM, rt_resolution = 0.05, n_chains = 2)
p_vector4 <- sampled_pars(DDM_s)
p_vector4[1:length(p_vector4)] <- c(1, 1, .3, .5, .1, .1)


# Benchmarks --------------------------------------------------------------
designs_and_p_vectors <- list('LNR (2 acc)'=list(p_vector=p_vector1, emc=LNR_s1),
                              'RDM (2 acc)'=list(p_vector=p_vector2, emc=RDM_s),
                              'RDM_noA (2 acc)'=list(p_vector=p_vector2, emc=RDM_s),
                              'LBA (2 acc)'=list(p_vector=p_vector21, emc=LBA_s),
                              'WDM'=list(p_vector=p_vector3, emc=WDM_s),
                              'DDM'=list(p_vector=p_vector4, emc=DDM_s)
)
for(i in 1:length(designs_and_p_vectors)) {#nm in names(designs_and_p_vectors)) {
  nm <- names(designs_and_p_vectors)[i]
  print(nm)
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

  type = model$c_name
  bound=model$bound
  transform=model$transform
  pre_transform=model$pre_transform
  trend=model$trend

  p_mat <- t(as.matrix(p_vector))
  colnames(p_mat) <- names(p_vector)

  pars_old <- EMC2:::get_pars_c_wrapper(p_matrix = p_mat, data = dadm, constants = constants, designs = designs,
                                        bounds = bound, transforms = transform, pretransforms = pre_transform,
                                        p_types = p_types, trend = trend
  )
  pars_new <- EMC2:::get_pars_c_wrapper_oo(particle_matrix = p_mat, data = dadm, constants = constants, designs = designs,
                                           bounds = bound, transforms = transform, pretransforms = pre_transform,
                                           trend = trend
  )
  # print('All parameters equal?')
  # print(all(pars_old == pars_new))
  if(!all(round(pars_old,8) == round(pars_new,8))) stop('Discrepancy between old and new')


  # Parameter mapping speed? ------------------------------------------------
  print('Benchmark for a single get_pars_matrix()')
  bm1 <- microbenchmark::microbenchmark(
    f1=EMC2:::get_pars_c_wrapper(p_matrix = p_mat, data = dadm, constants = constants, designs = designs, bounds = model$bound, transforms = model$transform,
                                 pretransforms = model$pre_transform,p_types = p_types,trend = trend),
    f2=EMC2:::get_pars_c_wrapper_oo(particle_matrix = p_mat, data = dadm, constants = constants, designs = designs,
                                    bounds = bound, transforms = transform, pretransforms = pre_transform,
                                    trend = trend),
    times=2000, control=list(warmup=50)
  )
  print(bm1)



  # Likelihood evaluation speed? 50 particles  ----------------------------------
  set.seed(1)
  p_mat <- matrix(rnorm(100*length(p_vector)), ncol=length(p_vector))
  colnames(p_mat) <- names(p_vector)
  if(nm == 'WDM') {
    p_mat[,'sv'] = log(0)
    p_mat[,'SZ'] = qnorm(0)
  }
  if(nm == 'RDM_noA (2 acc)') {
    constants = c('A'=log(0))
    c_name='RDM-A0'
    p_mat <- p_mat[,!colnames(p_mat) == c('A')]
  } else {
    c_name = model$c_name
  }
  if(nm %in% c('RDM_noA (2 acc)', 'RDM (2 acc)')) p_mat[,'s'] <- log(1)
  lls_old <- EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                            model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                            model$trend)
  lls_new <- EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = c_name,
                               model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                               model$trend)
  # lls_new2 <- EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
  #                              model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
  #                              model$trend, use_pt_pipeline = FALSE)

  print(all(round(lls_old,6) == round(lls_new,6)))
  # if(!all(round(lls_old,8) == round(lls_new,8))) stop('Discrepancy between old and new')
  print(sum(lls_new))

  print('Benchmark for calc_ll with 50 particles')
  bm2 <- microbenchmark::microbenchmark(
    f1=EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                      model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                      model$trend),
    f2=EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                         model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                         model$trend),
    times=300, control=list(warmup=20)
  )
  print(bm2)
  print('')
}


# for(i in 1:1000) {
#   EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
#                     model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
#                     model$trend)
# }

