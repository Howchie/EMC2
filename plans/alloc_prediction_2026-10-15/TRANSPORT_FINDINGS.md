# Transport-scored allocator: findings (2026-10-15)

The patch `transport_allocator.patch` holds the whole working-tree change. It is not committed.

## What it does
- Each iteration records each move family's change to a state statistic: the group mean, the log
  group sds, the subjects' mean and the subjects' log spread.
- A group coordinate is scored on transport: the Bartlett lag-sum within L iterations of each
  family's change against the total change, over the pooled window variance.
- The coordinate's rho is that of the AR(1) with the same transport. Subject moves are credited at
  the plan's mean rates.

## Results
- **The attribution works.** On PM, shear carries most of the ridge's persistent transport, and
  subject HMC and particle steps carry little. T30 shows the same pattern, and Forstmann (the null
  check) is clean.
- **Prediction does not improve** (`eval_transport/horizon.R`, scored on the saved arms):

  | horizon L | median ratio, predicted / run tau, all coordinates | same, slowest 8 | Spearman |
  |---|---|---|---|
  | 1 | 0.17–0.41 | 0.12–0.27 | 0.58–0.80 |
  | 10–25 | 0.48–0.92 | 0.26–0.48 | 0.62–0.91 |
  | 50–99 | noise | noise | ≈0 |

  The slow tail stays 2.5–5x underpredicted at every horizon a 100-draw block supports.
- **PM screening arm, seed 1** (`arm_transport_1`):
  - CPU 19.6k s; HMC regime 1 every block (40,400 transitions).
  - Ridge R-hat 1.37–1.47, min ESS 6, median ESS 53.
  - new_1 on the same seed: CPU 11.0k s, ridge R-hat 1.20–1.26.
- **Across all arms the ridge outcome is mostly noise.** Shear rounds were about the same in every
  arm (1,300–1,800), yet ridge R-hat ranged 1.20–1.74. What separates the arms is CPU: any arm on
  HMC regime 1 pays about +80%.

## Why it chose regime 1
- At L = 25 the group looked well mixed (predicted tau about 22, realized about 150), so the binding
  units were the slowest subjects' lag-one rho.
- HMC also earns a measurable share of transport (0.005–0.04 per unit rate) on some mean
  coordinates, so regime 1 won the objective by a few percent.

## Conclusion
Not an improvement over bf677e14 (new_1) on this evidence. The block-local signal does not see the
drift between chains. The allocator therefore cannot be steered to it by better prediction or by
attribution alone.
