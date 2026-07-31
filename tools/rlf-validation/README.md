# RLF solver validation scripts

Diagnostics for the RLF nonlocal Fokker--Planck likelihood. They are not part of
the test suite: each takes minutes to hours and several are studies rather than
pass/fail checks. Run them against an explicitly built library so you are not
measuring a stale install:

```sh
R CMD INSTALL --no-docs --no-byte-compile -l /path/to/lib .
cd tools/rlf-validation && RLFLIB=/path/to/lib Rscript nxlaw.R
```

`refcdf.R` is a shared helper the others `source()`, so run them from this
directory.

## What each script establishes

| script | question |
| --- | --- |
| `nxlaw.R` | how much resolution an extrapolated pair needs, as a function of alpha |
| `rough.R` | where roughness of the profile log-likelihood in alpha comes from |
| `check15.R` | is a recovery offset the data, the simulator, or the solver |
| `selfgen15.R` | recovery when data are drawn from the solver's own cdf |
| `nscale.R` | does a recovery offset fall like 1/n |
| `matched.R` | does the offset survive matching generator and estimator resolution |
| `cdfcmp.R` | simulator empirical cdf against the solver's, no estimator involved |

## Findings (2026-07-30/31)

**Resolution needed falls steeply with alpha.** Scored as mean `|d log pdf|` per
trial over the central 96% against an nx = 1024 reference, averaged over
v in {0.5, 1, 3} and b0 in {1, 2}, the nx at which an extrapolated pair holds
0.01 per trial is

    alpha  1.05 1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80+
    nx      172  146  115   99   84   70   57   48    48

log-linear in alpha at a rate of about 1.84 until it bottoms out.

**But nx must not be scheduled on alpha.** Implemented and reverted, for two
independent reasons. It makes the discretisation bias a function of alpha, which
adds a spurious term to the score for alpha: a schedule reaching nx = 96 at the
alpha = 1.3 peak recovers -0.050 where a *flat* nx = 96 recovers -0.013, the
same grid at the peak for four times the error. And quantised nx steps the
log-likelihood by ~1 nat per cell crossed: RMS residual about a local quadratic
at alpha = 1.3 is 0.054 nats at flat nx = 128 against 1.25 nats under a schedule
spanning nx 81-97. Independently, grids below about 96 are jagged in alpha on
their own account (alpha = 1.7: 0.13 nats at nx = 128, 0.14 at 96, 0.60 at 64,
1.54 at 48; extrapolation amplifies it because it differences two solves), which
is what the cheap end of any schedule is made of. The best schedule found was
dominated by simply setting a lower flat nx.

**The horizon split is still required after Richardson extrapolation.**
Extrapolation assumes an asymptotic error expansion in h, and a ballooned domain
leaves too few cells in the boundary layer for one to hold. With the split off,
on 2000-trial data (converged peak in brackets):

    alpha 1.10 [1.111]   raw 128  1.010 -> 1.330   pair 128+192  1.092 -> 1.245
    alpha 1.50 [1.581]   raw 128  1.508 -> 1.517   pair  64+96   1.619 -> 1.568

and profile roughness at alpha = 1.3 goes 0.054 -> 0.478 nats. Its cost is
self-limiting: it only engages once t_max passes t_crit, so at alpha = 1.7 with
5 s of data no bucket is occupied and it is a no-op. Loosening the threshold to
`2 (c b)^alpha` (c = 1.3, 1.6, 2.0) saves 5-20% and is erratic in both
directions; rejected.

**Open: a fixed +0.03 recovery offset at alpha >= 1.5.** Recovering alpha from
2000-trial data with a converged solver gives -0.003 +- 0.006 at alpha = 1.1 but
+0.035 +- 0.009 at 1.5 and +0.025 +- 0.013 at 1.7 — the opposite alpha
dependence from the discretisation error, so a second effect. Excluded so far:

* *the data being one unlucky seed* — it is a mean over 8 seeds;
* *the simulator's step* — refining `emc2.rlf_sim_dt` 5x moves alpha-hat by
  +0.001 at alpha = 1.5;
* *the discretisation* — the peak is invariant across `emc2.rlf_dt` 8e-3 to
  2e-3, domain width x1 to x2.5, nx 192+288 to 384+576, and split on/off;
* *the simulator's formulation* — data drawn by inversion from the solver's own
  cdf, no simulator anywhere, reproduce it (+0.033 +- 0.016 at alpha = 1.5);
* *a finite-sample estimator bias* — it does not fall with n: +0.021 at
  n = 2000, +0.026 at 8000, +0.028 at 32000, while the sd falls as 1/sqrt(n)
  (0.038, 0.016, 0.010) exactly as it should.

A fixed offset between a generator and an estimator that are both meant to be
the same model means they are not the same distribution. The remaining suspect
is a residual difference between the resolution the data were generated at and
the one they are fitted at — including the extrapolation *order*: p = 1 is
assumed, and where the true order is higher the pair family converges somewhere
slightly different from the raw family. On the alpha = 1.5 seed the raw ladder
(128, 192, 256, 384 -> 1.508, 1.525, 1.534, 1.545) extrapolates to about 1.565
while the pair family sits at 1.579-1.581. `matched.R` is the discriminator:
generate and fit at the same resolution, with and without extrapolation.

Practical consequence either way: at alpha >= 1.5 this offset is 1-3x the
per-fit discretisation error, so it, not the grid, sets the floor on alpha
recovery there.

## Two traps these scripts exist to avoid

**Never build a reference cdf from one `rlf_fht_pdf_cdf_vec` call at a long
t_max.** `rlf_lower_extent` sizes the domain from the horizon, so a 22 s grid at
nx = 768 leaves ~12 cells between the start point and the absorbing boundary —
coarser there than the nx = 128 default is at 2 s. It is the worst-resolved
configuration available, not the best. It produced an impossible result (all
eight alpha = 1.1 seeds pinned at the profile edge) before it was caught. Build
references through `pRLF` (see `refcdf.R`), which buckets each time by its own
horizon.

**Never score recovery against the generating value alone.** A data set's own
MLE is displaced from truth by sampling and by whatever the generator does; at
alpha = 1.9 every configuration, raw and paired, agrees on 1.878 while truth is
1.90. Score configurations against a converged solver *on the same data* to
isolate discretisation, and treat the converged-versus-truth gap as a separate
question. Relatedly, time a profile at the generating alpha, not at the bottom
of the profile grid, or configurations whose cost depends on alpha all report
the same number.
