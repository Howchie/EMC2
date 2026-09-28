\name{group_move_diagnostics}
\alias{group_move_diagnostics}
\title{Inspect group-move adaptation and frozen sampling state}
\usage{group_move_diagnostics(emc)}
\arguments{\item{emc}{An EMC fit or list of sampler chains.}}
\value{A data frame with one row per chain.}
\description{Reports algorithm state, not posterior convergence. Use posterior diagnostics for means, scales, correlations and subject parameters separately.}
\details{The adaptive method is regularized adaptive Metropolis (AM), selected before fitting with \code{options(emc2.group_move_method = "am")}. The legacy method remains the default while AM release gates are in progress and retains its four-proposal default. AM collects adaptation moments during the discarded \code{adapt} stage, then freezes its proposal after the configured warmup and settling periods. AM uses one proposal per iteration. See \code{benchmarks/group-move/README.md} in the source repository for the validation status.}
