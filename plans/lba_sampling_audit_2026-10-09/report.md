# LBA checkpoint sampling investigation, 2026-10-09

Checkpoint: `/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_lba7r_sepA_sepL_sepV_t0T_mt_v2.RData`. The live fit continued during inspection. The frozen diagnostic snapshot contains 859 retained sample draws per chain and 1567 production transitions, through sample try 16. It is at `/tmp/lba_sampling_audit/snapshot.rds`. Installed EMC2 3.4.0; installed `check_progress()` body agrees with the current repository source. No sampler or live checkpoint changes were made.

## Rising Rhat and trimming

`check_progress()` captures the draw count before calling `check_gd()`. The latter may discard the first third when doing so strictly improves the largest Rhat divided by its threshold. The logged Rhats are after that operation, whereas the iteration count is before it. Thus the lower count appears on the next line. The comparison is against the untrimmed window of the same check, not against the previous try. It also optimizes the worst rule, rather than requiring every individual selection to improve.

Reconstructed relevant windows:

| Try | Reported count | Retained count | Alpha max Rhat | Mu max Rhat |
|---|---:|---:|---:|---:|
| 10 | 933 | 622 | 1.074 | 1.079 |
| 11 | 689 | 689 | 1.163 | 1.113 |
| 14 | 989 | 659 | 1.474 | 1.158 |
| 15 | 759 | 759 | 1.405 | 1.182 |
| 16 | 859 | 859 | 1.306 | 1.169 |

Try 11 adds 67 draws; the other listed follow-up blocks add 100. The higher Rhats on try 11 come from new draws and a different split of the window after try 10's accepted trim.

The later alpha maxima are participant 52091, internal subject index 14, chiefly `B_lR_target`. In chain 3, this participant moves into a different region of B, threshold contrasts, A and t0 for hundreds of draws, then returns toward the region visited by chains 1 and 2. Their B_lR_target chain means in the frozen window are -0.066, -0.062 and -0.288. The group maxima are mainly v.e_novel and v.g. This is compatible with a difficult posterior ridge or multiple regions; the traces alone do not establish separate modes.

Fresh likelihood calculations at 11 selected draws of participant 52091 in each chain match the 33 stored likelihoods exactly, including the excursion. This checks state/likelihood consistency in those selected states, not correctness of every transition.

## Collapsed move

The checkpoint uses diagonal-gamma, with shape 2 on precisions and parameter-specific rates. It does not use the script's diagHW flag. `test-collapsed-move.R` explicitly exercises diagonal-gamma prior/Jacobian/gradient/correction identities and joint invariance; it also exercises standard Huang-Wand blocks plus an unblocked diagonal coordinate. These small Gaussian unit-test fixtures are not evidence for efficient transport on this nonlinear LBA posterior. The subject-HMC benchmark runner also calls make_emc without specifying a hierarchy type, using its standard default, rather than explicitly configuring an all-diagonal Huang-Wand or diagonal-gamma hierarchy. Thus the benchmark panel does not by itself establish efficiency for this hierarchy.

Stored collapsed counters show zero accepted iterations across about 484-491 likelihood-evaluated proposals per chain. The allocator sets its rate to zero by production transition 100. Epsilon calibration targets the inner surrogate acceptance only; allocation can remove an unsuccessful move, but does not calibrate its outer transport distance.

Diagnostic probes used each chain's final saved state and frozen proposal fits, with four independently generated proposals per chain per transport scale (12 per scale). They do not update the checkpoint. The scale multiplies the mass square root C; epsilon, three inner trajectories and trajectory length 1.5 are retained. Acceptance entries below are mean computed acceptance probabilities, not empirical accepted transition counts or efficiency estimates.

| C multiplier | Inner acceptance probability | Exact acceptance probability |
|---|---:|---:|
| 1.00 | 0.926 | 0 |
| 0.30 | 0.997 | 1.14e-18 |
| 0.10 | 1.000 | 0.0757 |
| 0.03 | 1.000 | 0.700 |

At current scale, true summed likelihood changes are about -8,665 to -28,446 while the quadratic predicts changes from +16 to +168. Mean exact log acceptance is -16,395. The Gaussian transport crosses the RT support: all 12 current-scale candidates give t0 >= an observed finite RT for at least two subjects; the mean is 8.25 subjects per candidate. One proposal gives participant 52091 nondecision times around 3.36-3.60 seconds despite recorded RTs around one second. Smaller scales avoid this t0 crossing in these probes. Explicit model bounds, with their zero exceptions respected, are also violated for a mean of 10.75 subjects at full scale, 3.33 at 0.3 scale, 0.42 at 0.1 scale and none at 0.03 scale. Other nonlinear likelihood errors also remain.

An independent full joint target ratio includes the subject Gaussian priors, group mean/variance priors, transport determinant, chart Jacobian and reverse surrogate density ratio. It agrees with the implemented likelihood residual correction to max absolute error 3.8e-11 across the 48 candidates. A cancellation-resistant evaluation of quadratic differences agrees with the implemented correction to 1.1e-11. Thus these checked ratios give no evidence of a missing diagonal-gamma prior/Jacobian term or catastrophic cancellation.

The verified failure is proposal efficiency: the surrogate permits large transports that are unusable under the true nonlinear likelihood/support, and current tuning does not shrink the outer transport in response. Shrinking C demonstrates recoverable acceptance, but does not establish an optimal scale or effective sampling. A principled next change would separately calibrate transport distance using exact correction/jump evidence during discarded warmup, and ensure changed surrogates can receive fresh evidence after allocation has disabled the move. Changing the prior is not supported as the remedy by these checks.

## Artifacts

- `traces.png`, `traces.pdf`: selected chain traces (PNG includes the group parameters driving later maxima).
- `diagnostics.csv`: original and rank-normalized folded split Rhat for every coefficient, with participant IDs and internal indices; move frequencies and lag-one autocorrelation.
- `verify.txt`: fresh likelihood comparisons and reconstructed trim counts.
- `collapsed_probes.csv`: per-candidate exact/surrogate correction and independent ratio checks.
- `collapsed_subject_corrections.csv`: per-subject decomposition.
- `collapsed_bounds.csv`: mapped-parameter support checks.
- R scripts: reproducible diagnostic calculations, using the frozen temporary snapshot where applicable.
