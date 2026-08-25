\name{GOM}
\alias{GOM}
\alias{GOMPERTZ}
\alias{Gompertz}
\title{Gompertz Growth Process race model}
\description{
Racing Gompertz growth processes. Each accumulator follows
\deqn{dX = \alpha X\log(K/X)dt + \beta X dW}
from a physical-space start range \eqn{X(0)\sim U(1,1+A)}, with response
threshold \eqn{1+B+A}. The likelihood is evaluated after the exact
\eqn{Y=\log(X)} reduction to an Ornstein--Uhlenbeck process.
}
\usage{
GOM(boundary_collapse = c("fixed", "exponential", "linear_additive",
                          "linear_multiplicative", "weibull"))
GOMPERTZ(boundary_collapse = c("fixed", "exponential", "linear_additive",
                               "linear_multiplicative", "weibull"))
Gompertz(boundary_collapse = c("fixed", "exponential", "linear_additive",
                              "linear_multiplicative", "weibull"))
}
\arguments{
  \item{boundary_collapse}{Boundary family. Non-fixed variants estimate
  physical-space \code{Binf} and \code{tau}; the Weibull form also estimates
  \code{pw}.}
}
\details{
The transformed OU parameters are
\eqn{k=\alpha}, \eqn{s=\beta}, and
\eqn{v=\alpha\log(K)-\beta^2/2}. Physical-space start variability and
moving boundaries are retained through the coordinate change.

\tabular{llllll}{
\strong{Parameter} \tab \strong{Transform} \tab \strong{Natural scale} \tab \strong{Default} \tab \strong{Mapping} \tab \strong{Interpretation} \cr
\emph{alpha} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Gompertz growth-rate parameter. \cr
\emph{beta} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Multiplicative diffusion scale. \cr
\emph{K} \tab log \tab [0, Inf] \tab log(5) \tab  \tab Carrying capacity of the Gompertz process. \cr
\emph{B} \tab log \tab [0, Inf] \tab log(1) \tab \emph{b} = 1 + \emph{B} + \emph{A} \tab Distance from the baseline start level to the response boundary. \cr
\emph{A} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Start-point range above the baseline level 1. \cr
\emph{t0} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Non-decision time. \cr
\emph{Binf} \tab log \tab [0, Inf] \tab log(1.5) \tab  \tab Asymptotic boundary value; used for non-fixed boundary collapse. \cr
\emph{tau} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Boundary-collapse time scale; used for non-fixed collapse. \cr
\emph{pw} \tab log \tab [0, Inf] \tab log(1) \tab  \tab Weibull boundary-collapse exponent; used only for \code{boundary_collapse = "weibull"}. \cr
}

\code{Binf} and \code{tau} are included for every non-fixed boundary-collapse
form; \code{pw} is included only for the Weibull form. The generic
\code{pContaminant} and \code{pGuess} nuisance parameters are appended by the
data pipeline when requested.
}
\value{A list defining an EMC2 race model.}
\seealso{\code{ROU}}
