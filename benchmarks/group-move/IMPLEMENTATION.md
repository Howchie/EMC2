# Group-move implementation

The legacy method remains the default. To test regularized adaptive Metropolis
(AM) on a fresh fit:

```r
options(
  emc2.group_move_method = "am",
  emc2.group_move_proposals = 1L,
  emc2.group_move_warmup = 6000L,
  emc2.group_move_settle = 500L
)
```

AM uses the Gibbs-step covariance as its initial proposal shape. During the
discarded `adapt` stage it updates a regularized empirical covariance and a
scalar log scale. The empirical covariance is blended with the initial metric
using weight `n / (n + 5*d)`. Shape is refreshed every 25 adaptation
opportunities; scalar feedback is applied to the proposal factor on every
opportunity. Preburn and burn transitions use the initial metric but do not
enter AM's moments or scale adaptation. AM freezes after its configured
adaptation and settling budgets. Production requires a ready, frozen state.

The AM state and its mathematical prior signature are stored per chain and
survive checkpoint serialization. The prior's auxiliary `design` attribute is
excluded from the schema signature because its model closures and C++ pointers
are reconstructed across stage and worker boundaries.

Candidate likelihood failures are retried on the master after a worker reports
an evaluation error. The failure policy records the cause and treats the
candidate as a zero-density rejection by default; `options(emc2.failure_policy
= "strict")` stops on nonnumerical failures. Numerical failures remain ordinary
rejections under every policy. A failed group Gibbs transition remains fatal.

The transport is currently supported for standard and diagonal-gamma
hierarchies with an intercept for translated coordinates. It can include
residual-scale coordinates when `emc2.group_move_scale = TRUE`. Unsupported
representations retain their existing sampler. AM requires one candidate per
iteration and has a configurable dense-dimension limit.

See [README.md](README.md) for the remaining release gates.
