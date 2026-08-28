devtools::load_all(".")

pars <- matrix(c(v=1.5, sv=1, B=1, A=1, t0=0.2, k=1, ell=0.1, alpha=0.5, eta=0), nrow=1)
colnames(pars) <- c("v", "sv", "b", "A", "t0", "k", "ell", "alpha", "eta")
attr(pars, "ok") <- TRUE

rt <- 0.5
cat("PDF:\n")
print(dBAwDD(rt, pars, launch=4L, posdrift=TRUE))
cat("CDF:\n")
print(pBAwDD(rt, pars, launch=4L, posdrift=TRUE))

