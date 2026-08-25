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
\code{nu} and standard deviation \code{sv}. The sampled \code{k} parameter is a
non-negative criterion offset. It is rounded to the nearest non-negative
integer with half-integers rounded upward, then shifted by two, so the
canonical event threshold is \eqn{K = 2 + floor(k + 1/2)}. A geometric excess
with mean \code{omega} is added to
that threshold; \code{t0} is a nondecision-time shift.
}
\details{
Each accumulator is a pure-birth counter with event rate
\code{nu + gamma * n}. The trialwise input rate has mean \code{nu} and
standard deviation \code{sv}; the response threshold is
\eqn{K = 2 + floor(k + 1/2)} and has a geometric excess with mean
\code{omega}. Exact-zero \code{sv}, \code{gamma}, and \code{omega} select
simpler nested branches.

The closed form evaluates its Stirling recurrence row by row, so finite
\code{k} values are not restricted by a substantive threshold cap. The
implementation only checks that the integer index is representable on the
target machine. The proper \code{sv = 0}, \code{omega = 0},
\code{gamma > 0} member overlaps with the non-defective boundary of FRQ:
for integer \code{alpha = K}, \code{beta = nu/gamma}, \code{h = 1}, and
\code{lambda = gamma}, both have
\eqn{F(t) = I_{1-exp(-gamma*t)}(K, nu/gamma)}. FRQ is not otherwise
redundant: it supplies a finite-reservoir defect and threshold-shape
relaxation, while PCOUNTER supplies trialwise input variation,
self-excitation, and geometric threshold variation.
For estimated parameters, \code{k} has lower bound zero, while \code{sv},
\code{gamma}, and \code{omega} have lower bound \code{1e-4}; exact zero is
available when they are fixed explicitly.

The model uses the following parameter matrix:

\tabular{llllll}{
\strong{Parameter} \tab \strong{Transform} \tab \strong{Natural scale} \tab \strong{Default} \tab \strong{Mapping} \tab \strong{Interpretation} \cr
\emph{nu} \tab log \tab [0, Inf] \tab log(10) \tab  \tab Mean trialwise input/event rate. \cr
\emph{sv} \tab log \tab [0, Inf] \tab log(0) \tab  \tab SD of the trialwise input-rate draw. \cr
\emph{gamma} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Self-excitation increment per accumulated event. \cr
\emph{k} \tab log \tab [0, Inf] \tab log(3) \tab \emph{K} = 2 + floor(\emph{k} + 0.5) \tab Threshold offset determining the baseline event count. \cr
\emph{omega} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Mean geometric excess above the baseline threshold. \cr
\emph{t0} \tab log \tab [0, Inf] \tab log(0) \tab  \tab Non-decision time. \cr
\emph{pContaminant} \tab probit \tab [0, 1] \tab qnorm(0) \tab  \tab Optional omission contaminant probability. \cr
\emph{pGuess} \tab probit \tab [0, 1] \tab qnorm(0) \tab  \tab Optional uniform outlier probability. \cr
}

\code{k} is sampled continuously but converted to the integer threshold shown
in the mapping column. Exact zero values of \code{sv}, \code{gamma}, \code{omega},
and \code{k} select the corresponding nested boundary models.
}
