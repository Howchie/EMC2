# Extracted from test-variant_funs.R:29

# prequel ----------------------------------------------------------------------
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)
design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                          formula=list(m~lM,s~1,t0~1),
                          contrasts=list(m=list(lM=ADmat)))
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
N <- 25
idx <- N + 1

# test -------------------------------------------------------------------------
skip_on_os("windows")
LNR_diag <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "diagonal")
LNR_diag <- init_chains(LNR_diag, cores_for_chains = 1, particles = 5)
