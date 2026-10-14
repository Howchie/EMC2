I ran the offline emulation on the lbaE replay. It fails the gate we set, so I stopped before building anything: no C++, no sampler changes, no tests.

  What worked: subject transport no longer has a summed-residual problem. Each subject moves exactly along its own line, so with 80 subjects the correction step accepts about 90% of moves. The
  leftover error is about 0.02 nats per subject, and the forward and reverse maps match to 1e-13. Subject cliffs no longer matter, because the node values capture them.

  What didn't: holding the other dimensions of each subject still pins the group sd. Measured as sd of lsd(B_Espeed) in window sds, over 3 chains × 8 draws:

  ┌────────────────────────────────────────────────────────┬──────────────────┐
  │                    Each subject...                     │ sd in window sds │
  ├────────────────────────────────────────────────────────┼──────────────────┤
  │ fully held (centred Gibbs)                             │ 0.11–0.18        │
  ├────────────────────────────────────────────────────────┼──────────────────┤
  │ collapsed along one line (this move)                   │ 0.30–0.39        │
  ├────────────────────────────────────────────────────────┼──────────────────┤
  │ fully collapsed under its Gaussian frame (upper bound) │ 0.60–0.97        │
  └────────────────────────────────────────────────────────┴──────────────────┘

  So the information that pins this coordinate is spread over the other 10 dimensions of each subject, not just the cliff direction.

  Against the gate:
  - Jump acceptance at −1/−0.5/+0.5/+1 window sd: ±1 is no better than the refit ceiling (0.10/0.37/0.45/0.16). A hybrid that moves the other dimensions by the Gaussian frame and the line exactly
    gets 0.16/0.42/0.45/0.14. Fixed frames on the same draws get 0.04/0.18/0.34/0.13.
  - Distance moved per round: about 1.7× a fixed-frame shear round. Each round costs about 2× in likelihood calls, so per CPU it's roughly a wash, well short of the 3× target.
  - Line direction: the frame-based direction beats the plain parameter axis (0.34 vs 0.24).
  - Likelihood cost on lbaE: 0.30 ms per subject for 1 row, 1.06 ms for 16, 6.8 ms for 128. The fixed cost per call dominates up to about 16 rows.

  I didn't run T30 or N-back, since lbaE already fails the bar.

  What this means for next steps: a fix for this geometry has to free more than one dimension per subject, and still handle the non-Gaussian shape in the dimensions it frees. Any 1-D cliff fix,
  curved or exact, tops out at the refit ceiling. The next idea to screen would be a 2-D collapse per subject, with the plane chosen to maximise the Gaussian full-collapse width. That costs about
  64–128 rows per subject per round, so it only pays if it gets most of the way to the 0.6–1.0 bound. The emulation script already measures that, so it's a cheap test before building anything; I
  haven't run it yet.