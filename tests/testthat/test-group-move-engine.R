local_rng_guard()
engine_real_ll_manager <- EMC2:::calc_ll_manager

engine_fixture <- function(p = 2L, n = 4L, type = "diagonal-gamma") {
  s <- list(type = type, n_pars = p, n_subjects = n,
            par_names = paste0("p", seq_len(p)), nuisance = rep(FALSE, p),
            marginalised_idx = rep(FALSE, p), model = NULL, marginalise = NULL)
  s <- if (type == "standard") EMC2:::add_info_standard(s, prior = NULL,
                par_groups = rep(1, p)) else EMC2:::add_info_diag_gamma(s, prior = NULL)
  s$data <- lapply(seq_len(n), function(i) list(y = rep(0, p), prec = diag(p)))
  alpha <- matrix(.2, p, n, dimnames = list(s$par_names, NULL))
  pars <- list(tmu = rep(.1, p), tvar = diag(.3, p), tvinv = diag(1/.3, p),
               alpha = alpha, a_half = rep(1, p))
  list(sampler = s, pars = pars, proposals = rbind(alpha, rep(-.5*p*.2^2, n)))
}
engine_ll <- function(props, data, model, marginalise = NULL) {
  r <- sweep(props, 2, data$y)
  -.5*rowSums((r %*% data$prec)*r)
}
engine_step <- function(x, stage = "adapt") {
  out <- EMC2:::.emc_group_move(x$sampler, x$pars, x$proposals, stage, NULL, NULL, NULL)
  out[c("sampler", "pars", "proposals")]
}
engine_options <- function(method = "am", warmup = 20, settle = 5) {
  withr::local_options(list(emc2.group_move_method = method, emc2.group_move = TRUE,
    emc2.group_move_proposals = 1L, emc2.group_move_scale = TRUE,
    emc2.group_move_warmup = warmup, emc2.group_move_settle = settle),
    .local_envir = getOption("testthat_topenv", parent.frame()))
}

test_that("the group move is off by default and every method defaults to one proposal", {
  withr::local_options(list(emc2.group_move = NULL, emc2.group_move_method = "legacy",
                            emc2.group_move_proposals = NULL))
  expect_false(EMC2:::.emc_group_move_options()$enabled)
  expect_identical(EMC2:::.emc_group_move_options()$K, 1L)
  withr::local_options(list(emc2.group_move_method = "am",
                            emc2.group_move_proposals = NULL))
  expect_identical(EMC2:::.emc_group_move_options()$K, 1L)
  expect_identical(EMC2:::.emc_group_move_config()$K, 1L)
})

test_that("checkpoint round trips preserve every chain's dynamic fields", {
  x <- engine_fixture()$sampler
  x$samples <- list(idx=1); x$rng <- list(gibbs=1:7)
  x$group_move <- list(factor=diag(2), n_adapt=13)
  x$chains_mu <- list(c(1,2)); attr(x, "prop_var") <- 2
  y <- x; y$group_move$factor <- diag(3,2); y$group_move$n_adapt <- 99
  y$rng$gibbs <- 8:14; y$chains_mu <- list(c(3,4)); attr(y,"prop_var") <- 4
  emc <- list(x,y)
  expect_identical(EMC2:::restore_duplicates(emc), emc)
  stripped <- EMC2:::strip_duplicates(emc)
  expect_null(stripped[[2]]$data)
  path <- tempfile(); on.exit(unlink(path))
  saveRDS(stripped,path)
  expect_identical(EMC2:::restore_duplicates(readRDS(path)), emc)
  expect_identical(EMC2:::restore_duplicates(EMC2:::strip_duplicates(stripped)), emc)
  old <- list(x, list(samples=y$samples,rng=y$rng))
  expect_warning(migrated <- EMC2:::restore_duplicates(old), "legacy checkpoint")
  expect_null(migrated[[2]]$group_move)
  expect_match(migrated[[2]]$group_move_migration,"unavailable")
})

test_that("group-move signatures ignore reconstructed prior design closures", {
  engine_options(warmup=20, settle=5)
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll, .package="EMC2")
  new_model_closure <- function() local({
    value <- list(marker = 1L)
    function() value
  })
  design_a <- list(list(model = new_model_closure()))
  design_b <- list(list(model = new_model_closure()))

  x <- engine_fixture()
  attr(x$sampler$prior, "design") <- design_a
  spec <- EMC2:::.emc_group_move_spec(x$sampler)
  sig_a <- EMC2:::.emc_group_move_signature(x$sampler, spec)
  expect_null(attr(sig_a$prior, "design"))
  expect_identical(attr(sig_a$prior, "type"), "diagonal-gamma")
  expect_identical(attr(x$sampler$prior, "design"), design_a)

  attr(x$sampler$prior, "design") <- design_b
  expect_identical(EMC2:::.emc_group_move_signature(x$sampler, spec), sig_a)

  # Emulate a saved v1 checkpoint whose signature still contains the original
  # design closure, then restore a fresh closure as run_emc() does at a stage
  # boundary. The state should resume and rewrite its normalized signature.
  old_signature <- list(spec=spec, prior=x$sampler$prior, gd=x$sampler$gd,
                        marginal=x$sampler$marginalised_idx,
                        nuisance=x$sampler$nuisance)
  first <- engine_step(x)
  first$sampler$group_move$signature <- old_signature
  attr(first$sampler$prior, "design") <- list(list(model = new_model_closure()))
  resumed <- engine_step(first)
  expect_null(attr(resumed$sampler$group_move$signature$prior, "design"))
  expect_identical(attr(resumed$sampler$prior, "design")[[1]]$model(),
                   list(marker = 1L))

  resumed$sampler$prior$theta_mu_mean[1] <-
    resumed$sampler$prior$theta_mu_mean[1] + 0.1
  expect_error(engine_step(resumed), "Group move schema changed")
})

test_that("regular AM freezes and retains chain-local state through autosave", {
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll, .package="EMC2")
  engine_options("am", warmup=50, settle=5)
  set.seed(731); x <- engine_fixture()
  x$sampler$rng$gibbs <- .Random.seed
  a <- b <- x
  for (i in 1:55) {
    a <- engine_step(a); b <- engine_step(b)
    chain2 <- b$sampler; chain2$group_move$n_adapt <- 300L
    restored <- EMC2:::restore_duplicates(EMC2:::strip_duplicates(list(chain2,b$sampler)))
    b$sampler <- restored[[2]]
  }
  expect_identical(a,b)
  expect_identical(a$sampler$group_move$status,"ready")
  factor <- a$sampler$group_move$factor
  withr::local_options(list(emc2.group_move_method="legacy",emc2.group_move_scale=FALSE,
                            emc2.group_move_proposals=4L))
  z <- engine_step(a,"sample")
  expect_identical(z$sampler$group_move$factor,factor)
  expect_equal(z$sampler$group_move$n_adapt,50)
  expect_equal(z$sampler$group_move$attempts_sample,1)
  expect_error(engine_step(z,"adapt"),"retune")
})

test_that("new engine validates warmup, dimensions, capability and configuration", {
  engine_options()
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll,.package="EMC2")
  x <- engine_fixture()
  expect_error(engine_step(x,"sample"),"warmup")
  z <- engine_step(x)
  expect_error(engine_step(z,"sample"),"incomplete")
  z$sampler$prior$theta_mu_mean[1] <- 9
  expect_error(engine_step(z),"schema changed")
  x$sampler$type <- "factor"
  expect_identical(engine_step(x)$sampler$group_move$status,"unsupported")
  withr::local_options(list(emc2.group_move_proposals=4L))
  expect_error(EMC2:::.emc_group_move_config(),"require")
  old_checkpoint <- engine_fixture()
  old_checkpoint$sampler$group_move_config <- list(method="removed",enabled=TRUE)
  expect_error(engine_step(old_checkpoint),"method is no longer supported")
})

test_that("transport inverse and joint density cancellation are correct", {
  set.seed(618)
  for (type in c("standard","diagonal-gamma")) {
    for (p in c(1L,2L,4L)) {
      x <- engine_fixture(p=p,type=type)
      sp <- EMC2:::.emc_group_move_spec(x$sampler)
      cur <- x$pars; cur$subj_mu <- matrix(cur$tmu,p,4)
      delta <- rnorm(2*p,sd=.2)
      new <- EMC2:::.emc_group_move_apply(delta,sp,cur)
      inv <- EMC2:::.emc_group_move_apply(-delta,sp,new)
      expect_equal(inv,cur[names(inv)],tolerance=1e-12)
      prior <- x$sampler$prior
      P <- EMC2:::.emc_group_move_mu_invar(sp,prior)
      lp <- function(z) EMC2:::.emc_group_move_log_prior(z$tmu,z$tvar,sp,prior,P,cur$a_half)
      independent_prior <- function(z) {
        out <- mvtnorm::dmvnorm(z$tmu,prior$theta_mu_mean,solve(P),log=TRUE)
        if(type=="standard" && p>1) {
          out <- out + log(EMC2:::robust_diwish(z$tvar,prior$v+p-1,diag(2*prior$v/cur$a_half,p)))
        } else {
          tau <- 1/diag(z$tvar)
          shape <- if(type=="standard") prior$v/2 else prior$shape
          rate <- if(type=="standard") prior$v/cur$a_half else prior$rate
          out <- out + sum(dgamma(tau,shape,rate,log=TRUE)+2*log(tau))
        }
        out
      }
      expect_equal(lp(new)-lp(cur),independent_prior(new)-independent_prior(cur),tolerance=1e-10)
      normal <- function(z) sum(vapply(1:4,function(i) mvtnorm::dmvnorm(z$alpha[,i],z$subj_mu[,i],z$tvar,log=TRUE),0))
      reduced <- lp(new)-lp(cur)+sum(delta[p+seq_len(p)]*sp$jac_weight)
      covariance_jacobian <- if (type == "standard" && p > 1) p+1 else 2
      full <- independent_prior(new)-independent_prior(cur)+normal(new)-normal(cur)+
        sum(delta[p+seq_len(p)]*(covariance_jacobian+4))
      expect_equal(reduced,full,tolerance=1e-10)
    }
  }
})

test_that("checked evaluation turns likelihood failures into recorded rejections", {
  EMC2:::.emc_reject_reset("group_move")
  assign("announced", character(), envir=EMC2:::.emc_profile_state)
  props <- matrix(0,1,2)
  local_mocked_bindings(calc_ll_manager=function(...) -Inf,.package="EMC2")
  expect_identical(EMC2:::.emc_group_move_ll_checked(props,NULL,NULL),-Inf)
  local_mocked_bindings(calc_ll_manager=function(...) stop("non-finite likelihood"),.package="EMC2")
  expect_silent(expect_identical(EMC2:::.emc_group_move_ll_candidate(props,NULL,NULL),-Inf))
  expect_equal(EMC2:::.emc_reject_counts("group_move")[["numerical"]],1L)
  local_mocked_bindings(calc_ll_manager=function(...) stop("broken evaluator"),.package="EMC2")
  expect_warning(expect_identical(EMC2:::.emc_group_move_ll_candidate(props,NULL,NULL),-Inf),
                 "group_move update failed")
  expect_equal(EMC2:::.emc_reject_counts("group_move")[["unknown"]],1L)
  local_mocked_bindings(calc_ll_manager=function(...) NaN,.package="EMC2")
  expect_silent(expect_identical(EMC2:::.emc_group_move_ll_candidate(props,NULL,NULL),-Inf))
})

test_that("worker requests use checked evaluator without consuming RNG", {
  set.seed(8); before <- .Random.seed
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll,
                       .emc_group_move_ll_one=function(...) stop("legacy"),.package="EMC2")
  props <- list(matrix(c(1,2),1),matrix(c(3,4),1))
  ctx <- list(data=list(list(y=c(0,0),prec=diag(2)),list(y=c(1,1),prec=diag(2))),model=NULL)
  msg <- list(kind="group_move_ll",checked=TRUE,subs=c(2L,1L),props=props)
  out <- EMC2:::.emc_wpool_compute(msg,ctx)
  expect_equal(as.numeric(out$ll),c(-.5,-12.5))
  expect_identical(.Random.seed,before)
})

test_that("support rejection leaves parameter and likelihood caches untouched", {
  engine_options(warmup=3,settle=1)
  local_mocked_bindings(.emc_group_move_ll_checked=function(...) -Inf,.package="EMC2")
  set.seed(43); x <- engine_fixture(); z <- engine_step(x)
  expect_identical(z$pars,x$pars)
  expect_identical(z$proposals,x$proposals)
  expect_equal(z$sampler$group_move$accepts_warmup,0)
  expect_equal(z$sampler$group_move$n_adapt,1)
  expect_equal(z$sampler$group_move$sum_accept_prob_warmup,0)
})

test_that("AM ignores preburn and burn coordinates and applies scale feedback immediately", {
  engine_options(warmup=40, settle=2)
  local_mocked_bindings(.emc_group_move_ll_checked=function(...) -Inf,.package="EMC2")
  set.seed(44)
  x <- engine_fixture()
  x$pars$tmu[] <- 30
  preburn <- engine_step(x, "preburn")
  burn <- engine_step(preburn, "burn")
  expect_equal(burn$sampler$group_move$n_adapt, 0L)
  expect_equal(burn$sampler$group_move$am_n, 0L)

  burn$pars$tmu[] <- .1
  adapted <- engine_step(burn, "adapt")
  gm <- adapted$sampler$group_move
  expect_equal(gm$n_adapt, 1L)
  expect_equal(gm$am_n, 1L)
  expect_equal(gm$am_mean, c(.1, .1, rep(.5*log(.3), 2)), tolerance=1e-12)
  expect_equal(gm$factor, gm$shape_factor * exp(gm$am_log_scale), tolerance=1e-12)
  expect_false(identical(gm$factor, gm$shape_factor))
})

test_that("a failed group likelihood rejects the move without aborting the fit", {
  engine_options(warmup=3, settle=1)
  EMC2:::.emc_reject_reset("group_move")
  local_mocked_bindings(calc_ll_manager=function(...) stop("non-finite likelihood"),
                        .package="EMC2")
  set.seed(45)
  x <- engine_fixture()
  moved <- expect_silent(engine_step(x))
  expect_identical(moved$pars, x$pars)
  expect_identical(moved$proposals, x$proposals)
  expect_equal(moved$sampler$group_move$accepts_warmup, 0L)
  expect_gt(EMC2:::.emc_reject_counts("group_move")[["numerical"]], 0L)
})

test_that("new engine composes with fit stages and saves distinct chain state", {
  skip_on_os("windows")
  engine_options(warmup=30,settle=5)
  data("forstmann",package="EMC2")
  ids <- unique(forstmann$subjects)[1:3]
  dat <- droplevels(forstmann[forstmann$subjects %in% ids,])
  des <- design(data=dat,model=LNR,formula=list(m~1,s~1),constants=c(t0=log(.2)))
  local_mocked_bindings(calc_ll_manager=function(proposals,...) {
    matrix(-.5*rowSums(proposals^2)/10,ncol=1)
  }, .package="EMC2")
  set.seed(9371)
  emc <- suppressMessages(make_emc(dat,des,type="diagonal-gamma",compress=FALSE,n_chains=2L))
  before_fit <- .Random.seed
  path <- tempfile(fileext=".RData"); on.exit(unlink(path))
  fitted <- fit(emc,cores_for_chains=1L,stop_criteria=list(preburn=list(iter=10L),
    burn=list(iter=15L),adapt=list(iter=25L,min_unique=0L),sample=list(iter=30L)),
    verbose=FALSE,particle_factor=5,step_size=10L,max_tries=5L,trim=FALSE,fileName=path)
  fitted <- EMC2:::restore_duplicates(fitted)
  expect_true(all(vapply(fitted,function(s) s$group_move$status=="sampling",TRUE)))
  expect_true(all(vapply(fitted,function(s) s$group_move$n_adapt==30,TRUE)))
  expect_false(identical(fitted[[1]]$group_move$factor,fitted[[2]]$group_move$factor))
  expect_true(file.exists(path))
  assign(".Random.seed", before_fit, globalenv())
  unsaved <- fit(emc,cores_for_chains=1L,stop_criteria=list(preburn=list(iter=10L),
    burn=list(iter=15L),adapt=list(iter=25L,min_unique=0L),sample=list(iter=30L)),
    verbose=FALSE,particle_factor=5,step_size=10L,max_tries=5L,trim=FALSE)
  unsaved <- EMC2:::restore_duplicates(unsaved)
  for (i in seq_along(fitted)) {
    expect_identical(fitted[[i]]$samples, unsaved[[i]]$samples)
    expect_identical(fitted[[i]]$group_move, unsaved[[i]]$group_move)
    expect_identical(fitted[[i]]$rng, unsaved[[i]]$rng)
  }
})

test_that("AM preserves unknown-variance hierarchical posteriors", {
  skip_model_validation()
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll,.package="EMC2")
  run <- function(type,method,seed,iters=20000L) {
    engine_options(method,warmup=2000,settle=500)
    set.seed(seed); x <- engine_fixture(type=type)
    s <- x$sampler; p <- s$n_pars; n <- s$n_subjects
    y <- matrix(c(-.8,.3,.1,-.2,.7,1.1,-.3,.8),p)
    s$data <- lapply(seq_len(n),function(i) list(y=y[,i],prec=diag(.5,p)))
    s$samples <- list(theta_mu=matrix(0,p,1),theta_var=array(diag(p),c(p,p,1)),
                     last_theta_var_inv=diag(p),a_half=matrix(1,p,1),idx=1)
    alpha <- y
    out <- matrix(NA_real_,iters,7)
    for(t in seq_len(iters)) {
      pars <- EMC2:::gibbs_step(s,alpha,type)
      V <- solve(pars$tvinv+diag(.5,p)); L <- t(chol(V))
      alpha <- V %*% (matrix(pars$tvinv %*% pars$tmu,p,n)+.5*y)+L %*% matrix(rnorm(p*n),p,n)
      if(method!="legacy") {
        s$rng$gibbs <- .Random.seed
        ll <- vapply(1:n,function(i) engine_ll(matrix(alpha[,i],1),s$data[[i]],NULL),0)
        step <- engine_step(list(sampler=s,pars=pars,proposals=rbind(alpha,ll)),
                            if(t<=2500) "adapt" else "sample")
        s <- step$sampler; pars <- step$pars; alpha <- step$proposals[1:p,,drop=FALSE]
      }
      s$samples$theta_mu[,1] <- pars$tmu
      s$samples$theta_var[,,1] <- pars$tvar
      s$samples$last_theta_var_inv <- pars$tvinv
      if(!is.null(pars$a_half)) s$samples$a_half[,1] <- pars$a_half
      out[t,] <- c(pars$tmu,log(diag(pars$tvar)),alpha[,1],cov2cor(pars$tvar)[1,2])
    }
    out[-(1:2500),,drop=FALSE]
  }
  se <- function(z) {
    k <- floor(length(z)/50)
    sd(colMeans(matrix(z[seq_len(50*k)],k)))/sqrt(50)
  }
  for(type in c("standard","diagonal-gamma")) {
    ref <- run(type,"legacy",9721)
    for(method in "am") {
      x <- run(type,method,9722)
      for(j in seq_len(ncol(x))) {
        if(sd(ref[,j])==0) next
        # Means and second moments; a narrow but apparently precise chain fails.
        for(power in 1:2) {
          a <- x[,j]^power; b <- ref[,j]^power
          z <- abs(mean(a)-mean(b))/sqrt(se(a)^2+se(b)^2)
          expect_lt(z,6,label=paste(type,method,j,power))
        }
      }
    }
  }
})

test_that("group move diagnostics distinguish tuning and production", {
  engine_options(warmup=3,settle=1)
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll,.package="EMC2")
  x <- engine_fixture()
  for(i in 1:4) x <- engine_step(x)
  x <- engine_step(x,"sample")
  d <- group_move_diagnostics(list(x$sampler))
  expect_identical(d$status,"sampling")
  expect_identical(d$method,"am")
  expect_equal(d$n_adapt,3)
  expect_equal(d$sample_attempts,1)
  expect_true(d$sample_acceptance %in% c(0,1))
})

test_that("real worker likelihood pool preserves serial group-move trajectories", {
  skip_on_os("windows")
  engine_options(warmup=20,settle=5)
  withr::local_options(list(emc2.worker_backend="fork"))
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll,.package="EMC2")
  set.seed(981)
  a <- b <- engine_fixture()
  a$sampler$rng$gibbs <- b$sampler$rng$gibbs <- .Random.seed
  ctx <- list(data=a$sampler$data,model=NULL,marginalise=NULL)
  pool <- EMC2:::.emc_wpool_start(2L,ctx)
  skip_if(is.null(pool),"worker pool unavailable")
  on.exit(EMC2:::.emc_wpool_stop(pool),add=TRUE)
  for(i in 1:35) {
    stage <- if(i<=25) "adapt" else "sample"
    a <- engine_step(a,stage)
    out <- EMC2:::.emc_group_move(b$sampler,b$pars,b$proposals,stage,pool,ctx,
                                 list(c(1L,3L),c(2L,4L)))
    pool <- out$wpool; b <- out[c("sampler","pars","proposals")]
  }
  expect_identical(a,b)
})

test_that("failed group likelihood replies retry the same candidate on master", {
  ctx <- list(data=list(list(y=c(0,0),prec=diag(2)),list(y=c(1,1),prec=diag(2))),
              model=NULL,group_move_checked=TRUE)
  props <- list(matrix(c(1,2),1),matrix(c(3,4),1))
  local_mocked_bindings(.emc_group_move_ll_checked=engine_ll,
    .emc_wpool_request=function(...) TRUE,
    .emc_wpool_reply=function(...) list(failed="injected worker error"),
    .package="EMC2")
  set.seed(23); before <- .Random.seed
  out <- EMC2:::.emc_wpool_group_move_ll(list(n=2L,alive=TRUE),list(1L,2L),props,ctx)
  expect_equal(as.numeric(out$ll),c(-2.5,-6.5))
  expect_identical(.Random.seed,before)
})

test_that("worker likelihood errors become counted master-side rejections", {
  ctx <- list(data=list(list(y=c(0,0),prec=diag(2)),list(y=c(1,1),prec=diag(2))),
              model=NULL,group_move_checked=TRUE)
  props <- list(matrix(c(1,2),1),matrix(c(3,4),1))
  EMC2:::.emc_reject_reset("group_move")
  local_mocked_bindings(.emc_group_move_ll_checked=function(...) {
      stop("non-finite likelihood")
    }, .emc_wpool_request=function(...) TRUE,
    .emc_wpool_reply=function(...) list(failed="worker numerical error"),
    .package="EMC2")
  set.seed(24); before <- .Random.seed
  out <- EMC2:::.emc_wpool_group_move_ll(list(n=2L,alive=TRUE),list(1L,2L),props,ctx)
  expect_true(all(out$ll == -Inf))
  expect_equal(EMC2:::.emc_reject_counts("group_move")[["numerical"]],2L)
  expect_identical(.Random.seed,before)
})

test_that("installed spawned workers match serial group transitions on a real likelihood", {
  skip_on_os("windows")
  withr::local_options(list(emc2.worker_backend="spawn"))
  skip_if(EMC2:::.emc_wpool_backend(list())!="spawn","requires isolated installed package")
  engine_options(warmup=20,settle=5)
  # Explicitly use the real manager: preceding integration fixtures replace it.
  local_mocked_bindings(calc_ll_manager=engine_real_ll_manager,.package="EMC2")
  data("forstmann",package="EMC2")
  ids <- unique(forstmann$subjects)[1:2]
  dat <- droplevels(do.call(rbind,lapply(ids,function(id) head(forstmann[forstmann$subjects==id,],30))))
  des <- design(data=dat,model=LNR,formula=list(m~1,s~1),constants=c(t0=log(.2)))
  set.seed(8931)
  x <- suppressMessages(make_emc(dat,des,type="diagonal-gamma",compress=FALSE,n_chains=2L))
  s <- EMC2:::run_stages(x[[1]],stage="preburn",iter=5L,verbose=FALSE,
                         verboseProgress=FALSE,particle_factor=10,n_cores=1L)
  j <- s$samples$idx
  pars <- list(tmu=s$samples$theta_mu[,j],tvar=s$samples$theta_var[,,j],
               tvinv=s$samples$last_theta_var_inv)
  props <- rbind(s$samples$alpha[,,j],s$samples$subj_ll[,j])
  ctx <- list(data=s$data,model=EMC2:::.emc_wpool_slim_model(s$model),marginalise=s$marginalise)
  pool <- EMC2:::.emc_wpool_start(2L,ctx)
  skip_if(is.null(pool),"spawned pool unavailable")
  on.exit(EMC2:::.emc_wpool_stop(pool),add=TRUE)
  expect_identical(pool$backend,"spawn")
  a <- b <- list(sampler=s,pars=pars,proposals=props)
  for(i in 1:35) {
    stage <- if(i<=20) "adapt" else "sample"
    a <- engine_step(a,stage)
    out <- EMC2:::.emc_group_move(b$sampler,b$pars,b$proposals,stage,pool,ctx,list(1L,2L))
    pool <- out$wpool; b <- out[c("sampler","pars","proposals")]
  }
  expect_identical(a,b)
})

test_that("conditional proposal regularization is deterministic and preserves RNG", {
  set.seed(96); before <- .Random.seed
  singular <- matrix(1,3,3)
  a <- EMC2:::condMVN(rep(0,3),singular,dependent.ind=1L,given.ind=2:3,X.given=c(.2,.2))
  expect_identical(.Random.seed,before)
  set.seed(197)
  b <- EMC2:::condMVN(rep(0,3),singular,dependent.ind=1L,given.ind=2:3,X.given=c(.2,.2))
  expect_identical(a,b)
  expect_true(all(eigen(a$condVar,symmetric=TRUE)$values>0))
  S <- matrix(c(2,.5,.5,1),2)
  c <- EMC2:::condMVN(c(0,0),S,dependent.ind=1L,given.ind=2L,X.given=.3)
  expect_equal(c$condMean,.15)
  expect_equal(c$condVar,matrix(1.75,1,1))
})

test_that("conditional regularization respects parameter units", {
  S <- matrix(1,3,3)
  scales <- c(1e-6,1e3,2)
  a <- EMC2:::condMVN(rep(0,3),S,1L,2:3,c(.2,.2))
  b <- EMC2:::condMVN(rep(0,3),S*outer(scales,scales),1L,2:3,c(.2,.2)*scales[2:3])
  expect_equal(b$condMean/scales[1],a$condMean,tolerance=1e-7)
  expect_equal(b$condVar/scales[1]^2,a$condVar,tolerance=1e-7)
  tiny <- EMC2:::condMVN(c(0,0),diag(1e-12,2),1L,2L,0)
  expect_equal(tiny$condVar,matrix(1e-12,1,1),tolerance=1e-20)
})
