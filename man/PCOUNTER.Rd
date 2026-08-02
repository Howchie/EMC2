% Generated from R/model_LNR.R documentation.
\name{PCOUNTER}
\alias{PCOUNTER}
\title{The gamma-mixed self-exciting Poisson Counter race model}
\usage{
PCOUNTER()
}
\value{A model list compatible with \code{design()}.}
\description{
Each accumulator is a pure-birth counter with transition rate
\code{nu + gamma * n}. The trialwise input rate is gamma distributed with mean
\code{nu} and standard deviation \code{sv}. The threshold is \code{k} plus a
geometric excess with mean \code{omega}; \code{t0} is a nondecision-time shift.
}
\details{
The model integrates analytically over both the continuous input-rate
distribution and the infinite geometric threshold distribution. Setting
\code{sv = 0}, \code{gamma = 0}, or \code{omega = 0} selects the corresponding
fixed-rate, ordinary Poisson, or fixed-threshold branch. The criterion is
rounded to the nearest positive integer when evaluating the process.

For estimated parameters, \code{sv}, \code{gamma}, and \code{omega} have lower
bound \code{1e-4}; exact zero is available when they are fixed explicitly.

The model parameters are:

\tabular{llll}{
\strong{Parameter} \tab \strong{Transform} \tab \strong{Default} \tab \strong{Interpretation} \cr
\emph{nu} \tab log \tab log(10) \tab Mean external input rate. \cr
\emph{sv} \tab log \tab log(0) \tab Trialwise input-rate standard deviation. \cr
\emph{gamma} \tab log \tab log(0) \tab Self-excitation increment per event. \cr
\emph{k} \tab log \tab log(5) \tab Fixed criterion location. \cr
\emph{omega} \tab log \tab log(0) \tab Mean geometric criterion excess. \cr
\emph{t0} \tab log \tab log(0) \tab Nondecision-time shift. \cr
}
}
