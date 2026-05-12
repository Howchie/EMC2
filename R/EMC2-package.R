#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom abind abind
#' @importFrom Brobdingnag as.brob
#' @importFrom coda acfplot
#' @importFrom coda as.mcmc
#' @importFrom coda cumuplot
#' @importFrom coda effectiveSize
#' @importFrom coda gelman.diag
#' @importFrom coda gelman.plot
#' @importFrom coda mcmc
#' @importFrom coda mcmc.list
#' @importFrom coda traceplot
#' @importFrom colorspace diverging_hcl
#' @importFrom corrplot colorlegend
#' @importFrom corrplot corrplot
#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @importFrom graphics axis
#' @importFrom graphics hist pairs rect text
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics matplot
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot.default
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics segments
#' @importFrom graphics title
#' @importFrom grDevices adjustcolor
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom lpSolve lp.assign
#' @importFrom magic adiag
#' @importFrom MASS ginv
#' @importFrom Matrix nearPD
#' @importFrom matrixcalc is.positive.definite
#' @importFrom methods is
#' @importFrom mvtnorm dmvnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom parallel mclapply
#' @importFrom psych fa.diagram
#' @importFrom Rcpp sourceCpp
#' @importFrom stats acf
#' @importFrom stats aggregate
#' @importFrom stats approxfun
#' @importFrom stats as.dist
#' @importFrom stats as.formula
#' @importFrom stats ave
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats cov2cor
#' @importFrom stats cutree
#' @importFrom stats dbinom
#' @importFrom stats density
#' @importFrom stats density.default
#' @importFrom stats dexp
#' @importFrom stats dgamma
#' @importFrom stats dnorm
#' @importFrom stats ecdf
#' @importFrom stats hclust
#' @importFrom stats integrate
#' @importFrom stats lm
#' @importFrom stats median
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats na.pass
#' @importFrom stats optimize
#' @importFrom stats pexp
#' @importFrom stats pnorm
#' @importFrom stats predict
#' @importFrom stats qbinom
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats rchisq
#' @importFrom stats reformulate
#' @importFrom stats reshape
#' @importFrom stats residuals
#' @importFrom stats rexp
#' @importFrom stats rgamma
#' @importFrom stats rmultinom
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom stats smooth.spline
#' @importFrom stats terms
#' @importFrom stats update
#' @importFrom stats var
#' @importFrom stats varimax
#' @importFrom utils combn
#' @importFrom utils head
#' @importFrom WienR dWDM
#' @importFrom WienR pWDM
#' @importFrom WienR rWDM
#' @useDynLib EMC2, .registration = TRUE
## usethis namespace: end
NULL

.emc2_samples_lnr_path <- function(libname, pkgname) {
  path <- file.path(libname, pkgname, "extdata", "samples_LNR.rds")
  if (!file.exists(path)) {
    path <- system.file("extdata", "samples_LNR.rds", package = pkgname)
  }
  return(path)
}

.emc2_load_samples_lnr <- function(libname, pkgname) {
  readRDS(.emc2_samples_lnr_path(libname, pkgname))
}

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  delayedAssign(
    "samples_LNR",
    .emc2_load_samples_lnr(libname, pkgname),
    assign.env = ns
  )
}

.onAttach <- function(libname, pkgname) {
  pkg_env <- as.environment(paste0("package:", pkgname))
  delayedAssign(
    "samples_LNR",
    get("samples_LNR", envir = asNamespace(pkgname), inherits = FALSE),
    assign.env = pkg_env
  )
}
