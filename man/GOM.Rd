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
}
\value{A list defining an EMC2 race model.}
\seealso{\code{ROU}}
