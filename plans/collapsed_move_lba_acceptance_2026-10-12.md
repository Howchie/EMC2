# Collapsed group move: acceptance on LBA, and a request for expert review

## Minimal summary

The collapsed move (R/collapsed_move.R, src/collapsed_move.cpp) proposes a
jump of the whole group state phi = (mu, log-Cholesky of V, log sds) with
every subject carried along in its conditional frame (whitened residual
r_s fixed, a_s = m_s(phi) + R_s(phi)^-1 r_s). The proposal is m HMC
trajectories on the Laplace-collapsed surrogate phat(phi) (no likelihood
calls), followed by one exact Metropolis correction with
log w = sum_s [log L_s(a_s) - log Lhat_s(a_s)] at the transported subjects,
one likelihood round. Each subject's Lhat_s is a concave quadratic fitted by
least squares to its log likelihood over the window draws (Laplace fit when
too few draws). The move is exact for any surrogate; the surrogate only sets
the acceptance. It is the one mechanism aimed at the group log-sd funnel:
under exact Gaussian likelihoods phi and r are independent, so the move
walks the log sds with the subjects re-centred for free.

Where it works: the T30 fixture (30 subjects, log-normal race, 2026-10-08):
outer acceptance ~30%, allocator holds it at its cap of 4 moves per
iteration, min ESS of log sds 3-6x better than without it, stop rule met in
every seed.

Where it fails, measured:

| fixture | subjects | rounds tried per chain | accepted | summed acceptance prob | step size eps |
|---|---|---|---|---|---|
| N-back (2026-10-11 checkpoint) | 93 | 324-360 | 0 | ~0 | at floor |
| lbaT seed 1 (2026-10-11, current build) | 80 | 88 / 124 / 170 | 0 / 0 / 0 | 2e-16 / 0 / 2e-12 | 0.010 / 0.059 / 0.077 |
| LBA screen (plans/collapsed_transport_adaptation_2026-10-09) | 80 | - | - | expected ~5e-34 | - |

The allocator measures a factor of exactly one and buys zero rounds, which
is correct. On lbaT the slowest posterior direction is the group log sd of
the start-point range A (with a little of its mean): chains traverse the
posterior range of that coordinate in excursions of hundreds of iterations,
giving ~5 effective draws per 1700 and R-hat 1.3-1.6 at the 20-try cap on
every build. Only the shear scale rounds touch that direction and their
measured factors are near one, so the move that should carry it is the one
that never accepts.

Known mechanisms of the failure (diagnosed on N-back, 2026-10-11, scratch
scripts collapsed_diag.R, wall_emul.R, curv_check.R):

1. Declared bounds: t0 >= 0.05 floors ~17 of 93 N-back subjects' cells; a
   transported subject pushed through the wall costs -800 nats (finite, so
   non-finite checks miss it). Fixed by the wall-aware transport
   (truncated-quantile map along the tightest half-space, committed
   4be15aa8). Acceptance on N-back stayed ~0.
2. Surrogate scale: the quadratic surrogates have 2-6 zero eigenvalues per
   subject after PSD projection, so phat is far wider than the posterior.
   Inner HMC proposes ~2 posterior-sd jumps that phat accepts at 93% and
   the exact correction rejects at -13000 nats. With a wall-aware transport
   the non-floored subjects alone sum to -70 nats at 1.1-sd jumps and -15
   at 0.6 sd; acceptance only at ~0.2-sd jumps. A joint accept over n
   subjects needs per-subject surrogate error below ~1/n nats where the
   posterior lives; quadratic surrogates do not reach that at n = 80-93 on
   LBA (curved A/B/t0 ridges, one-sided plateau likelihoods).
3. LBA Laplace fits are poor: far whole-group jumps with residuals fixed
   give log w differences of -30 to -80 on an earlier LBA fixture.

Tried and removed: a subject-block prototype (accepting per subject block)
recovered acceptance ~0.47 but gave worse ESS per second over a short
continuation; exact-pilot scale searches and global scale calibration
(reverted, see plans/collapsed_transport_adaptation_2026-10-09/report.md);
bounding tricks beyond the wall transport (judged not worth retrying).

Constraints from the project: the sampler stays model-independent (no
per-model code), exactness must hold, no serial work that idles a chain's
likelihood workers (inner chains are C++ in the chain process; one
likelihood round per outer correction is the budget), and one default
sampler with no options framework.

## The sampler, for context

Model: y_s | a_s ~ L_s (per-subject likelihood over p parameters on an
unconstrained scale, e.g. log A, log B, log t0); a_s ~ N(mu, V) across
subjects; priors on mu and V (V full, inverse-Wishart style, or a
diagonal/factor structure). Three chains run in parallel, each with its own
pool of likelihood workers (9 per chain here); every move that evaluates
likelihoods is a batch over subjects on that pool. One iteration of one
chain, in order (R/sampling.R):

1. **Gibbs on the group level.** (mu, V) | {a_s} by their conjugate
   conditionals.
2. **Particle (importance) step per subject.** For each subject, K = 50 sqrt(p)
   particles from a mixture of its own adapted proposal and the conditional
   prior, weighted by L_s and N(a_s | mu, V), one resampled (an
   independence Metropolis step on a_s, exact). Each subject runs it with
   probability q_s set by the allocator (default q_s = 1; it stands down
   only where HMC replaces it).
3. **Ensemble group update.** A move of (mu, V) using draws pooled across
   chains as an ensemble proposal, with the subjects fixed; exact by
   Metropolis.
4. **Shear rounds (the frame-transport moves).** Each round steps one group
   coordinate: a mean direction d (from the Gibbs-conditional covariance's
   slow eigen-directions) or one log sd. The subjects are not fixed: each is
   transported in its conditional frame, a_s = m_s(mu, V) + R_s(mu, V)^-1 r_s
   with r_s fixed, where N(m_s, (R_s' R_s)^-1) is the subject's conditional
   posterior under a quadratic surrogate Lhat_s of its likelihood (a Laplace
   fit at the conditional mode, or on curved ridges a per-subject atlas of
   conditional-mode fits regressed on the slow directions so the frame
   follows the ridge as mu moves). Declared bounds enter as half-space walls
   with a truncated-quantile map along the wall normal. Eight candidates
   (Tjelmeland multi-proposal) are scored with one likelihood round; the
   Jacobian of the transport makes the step exact. The allocator sets the
   expected number of rounds per direction per iteration (0 to 4).
5. **Collapsed move** (the subject of this document): a move of the whole
   group state with the subjects transported as in 4, proposed by HMC on the
   Laplace-collapsed group posterior and corrected with one likelihood
   round. Allocated like a round; currently 0 on LBA.
6. **Pooled parameters** (if the design declares any): a move of parameters
   shared across subjects, one likelihood round.
7. **Subject HMC round.** For each subject with an allocated rate r_s, with
   probability r_s one HMC transition on a_s | mu, V (gradient by finite
   differences of L_s on the worker, whitened in the subject's frame, step
   size from a bracketed warm-up and pilots). Under the "refresh" option a
   subject must be refreshed by step 2 or step 7 each iteration (q_s + r_s
   >= 1) and the cheaper kernel that meets the target is chosen; by default
   the particle step is a floor and HMC is only added on top.

Allocation (R/move_allocation.R): every move is priced by its measured CPU
per use and its measured effect on the lag-one autocorrelation of each
coordinate (group means along the shear directions, group log sds, each
subject's whitened coordinates), from expected squared jumps against the
pooled production variance. At each production rebuild (every 100
iterations, on an increasingly rare schedule) the plan minimises
max over coordinates of tau x cost with tau = (1 + rho) / (1 - rho), over
rounds per direction, collapsed moves, subject HMC rates and particle rates,
with pilots for subjects lacking HMC evidence. Coordinates no move can lower
no longer bound the target (changed 2026-10-11). The lag-one surrogate is an
AR(1) model and understates the integrated time on slow directions.

Stages: pre-burn, burn (mode-start from a batched Hessian at the Laplace
mode, random-walk steps, stop on mean R-hat), adapt (builds the subjects'
adapted proposals), sample (production, stop rule on R-hat < 1.1 and
minimum ESS over up to 20 tries of 100 iterations).

What the pieces are for: 1 and 3 move the group given the subjects; 2 and 7
move the subjects given the group; 4 and 5 are the only moves that change a
group coordinate and the subjects together, which is what a hierarchical
funnel needs. On lbaT, 4 moves the log sd of A only a little per round
(factor near one) and 5 never accepts.

## Prompt for an expert review

You are reviewing a hierarchical MCMC sampler for cognitive models (EMC2,
R with C++ kernels). Per subject s = 1..n (n = 30-100), a likelihood
L_s(a_s) over p = 8-20 parameters (race models such as the LBA: the
likelihood is smooth but strongly ridged in start-point range A, threshold
B and non-decision time t0, with one-sided plateaus and declared bounds
such as t0 >= 0.05). Group level: a_s ~ N(mu, V) with V a full covariance
(inverse-Wishart style prior). Within an iteration the sampler already has:
Gibbs on (mu, V) given the a_s; a particle (importance) step and/or HMC per
subject; an ensemble move on the group; and "shear" moves that step mu or a
log sd while each subject is transported in its conditional frame under a
quadratic likelihood surrogate (exact by a Metropolis correction).

The slow direction that remains on LBA data is the group log sd of a weakly
identified parameter (the start-point range A at ~150 trials per subject),
coupled to the subjects' values of that parameter: a hierarchical funnel.
A collapsed move was built for it: the whole group state (mu and V) is moved
by HMC on the Laplace-collapsed group posterior with subjects carried in
their conditional frames, then corrected exactly with one likelihood round,
log w = sum over subjects of [log L_s - log Lhat_s]. It works on an easy
fixture (acceptance ~30%, n = 30) and has acceptance ~1e-12 to ~0 on LBA
fixtures with n = 80-93: the per-subject quadratic surrogate errors add up
across subjects, and a joint accept needs per-subject error well below 1/n
nats. Per-subject block acceptance was tried and gave worse ESS per second.

Questions:

1. Is there a way to move the group log sd of A jointly with the subjects
   on this model that keeps exactness and a cost of order one likelihood
   round per move? Candidates we have not evaluated: delayed-acceptance or
   multiple-try variants over the inner trajectories; tempered or
   annealed transitions between phat and the target; surrogates beyond a
   concave quadratic per subject (e.g. per-subject emulators along the
   ridge coordinates, or using the exact 1-D conditional profile along the
   moved direction); restricting the move to the one slow direction with
   exact per-subject line fits; control variates or subsampled
   pseudo-marginal corrections that reduce the sum's variance.
2. Is a different parameterisation the answer instead: a non-centred
   treatment of the weakly identified parameter(s) only, chosen by the
   sampler from the measured posterior geometry rather than per model, in a
   way that composes with the existing centred steps for well-identified
   parameters?
3. Which of these is most likely to work at n = 80-100 with ridged,
   plateaued likelihoods, and what diagnostic on a saved checkpoint would
   show whether it can before implementing it?
4. Is there evidence in the description above that the problem is not the
   surrogate but the chart (log-Cholesky of V with the subjects' whitened
   residuals fixed), such that a different conditional frame would make the
   quadratic surrogate adequate?

Please answer with mechanisms and the conditions under which they hold,
not with a list of papers; name a paper only where its construction is the
one to use. Assume the implementer can run a diagnostic on saved draws
within minutes but a full fit costs 20-60 minutes on 27 cores.

## Pointers

- Code: R/collapsed_move.R (design header), src/collapsed_move.cpp (lane
  batched surrogate), R/shear_move.R (frames, curved atlas, wall
  transport), R/move_allocation.R (allocation, factor measurement).
- Fits: scratch `bench/P_mech/lbaT_fullc_1_emc.rds` (lbaT, 80 subjects, all
  three chains with collapsed records), `ckpt.rds` (N-back 93),
  `bench/A/nback31_fullc_1_emc.rds` (N-back 31).
- Earlier reports: plans/collapsed_transport_adaptation_2026-10-09/report.md,
  plans/lba_sampling_audit_2026-10-09/report.md,
  plans/sampler_mechanisms_review_2026-10-11.md (+ response).
