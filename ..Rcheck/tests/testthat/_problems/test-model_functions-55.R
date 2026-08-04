# Extracted from test-model_functions.R:55

# prequel ----------------------------------------------------------------------
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)
design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                            formula=list(m~lM + E,s~1,t0~1),
                            contrasts=list(m=list(lM=ADmat)))
p_LNR <- c(m=-1,m_lMd=2, m_Eneutral = .1, m_Eaccuracy = .1, s = .5, t0=log(.2))
design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                            formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
p_LBA <- c(v=1.5,v_lMd=1, sv_lMTRUE = log(.5), B = log(2), B_Eneutral = .5,
           B_Eaccuracy = .3, B_lRright = 0, A = log(.4), t0=log(.2))
design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                            formula=list(v~lM,s~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))
p_RDM <- c(v=.5,v_lMd=1, s_lMTRUE = log(.5), B = log(2), B_Eneutral = .5,
           B_Eaccuracy = .3, B_lRright = 0, A = log(.4), t0=log(.2))
design_DDM <- design(data = dat,model=DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))
p_DDM <- c(v_Sleft = -2, v_Sright = 2, a = log(2), a_Eneutral = .3, a_Eaccuracy = .3, t0 = log(.2),
           Z = qnorm(.5), sv = 1, SZ = qnorm(.3))
LNR_s <- make_emc(dat, design_LNR, compress = F, n_chains = 1)
LBA_s <- make_emc(dat, design_LBA, compress = F, n_chains = 1)
RDM_s <- make_emc(dat, design_RDM, compress = F, n_chains = 1)
DDM_s <- make_emc(dat, design_DDM, compress = F, n_chains = 1)

# test -------------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
expect_snapshot(init_chains(LBA_s, particles = 10, cores_per_chain = 1)[[1]]$samples)
