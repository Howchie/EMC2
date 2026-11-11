# Benchmark protocol and gates

The paper fixture fixes data generation at seed 123 and uses sampler seeds
101, 202 and 303. EXG uses 10,000 trials, a 90th-percentile censoring bound,
1 ms resolution, and three single-subject chains. LNR uses the paper's
truncation and 90th-percentile censoring construction. Fixed-transition runs
request 50 preburn, 200 burn, 200 adapt and 400 sample iterations, with
100-iteration blocks; check the saved chain stages to determine actual counts.
Time-to-precision runs use the same setup but require minimum alpha ESS 800
across chains and at least 400 sampling iterations. They must report whether
the target was reached rather than assuming the stop condition succeeded.

Run baseline and candidate in frozen installed libraries, with the benchmark
lock, CPUs 0–26 and OMP, OpenBLAS and MKL threads set to one. Run one fresh R
process per case, seed and arm. Capture `/usr/bin/time -v`, fitted objects and
session information. Controlled profiles disable the fit host and enable the
sampler profiler; time-to-precision and fixed-transition comparisons use the
default host. `strace` runs are separate from timing runs because tracing
changes runtime. Repeat arm pairs rather than comparing only best times.

For setup-only changes, compare each stored sample field, subject and chain
RNG states, and every continuation field, interpreting process pointers and
closures separately. No draw or RNG differences are accepted. Compare
paired elapsed times and inspect system CPU, faults and block profiles before
attributing a gain. Record negative or inconclusive results as such.

For later changes that intentionally alter streams, compare each 5%, 50%,
and 95% posterior quantile using the combined Monte Carlo standard error of
the two runs; require the absolute difference to be within four combined
MCSEs for every parameter, with more draws or seeds when MCSE is unstable.
For ESS per draw, the proposed non-inferiority margin is 10%: the lower end
of an across-seed confidence interval for candidate/baseline must exceed 0.90.
Three seeds are a pilot, not enough to make a reliable non-inferiority claim;
add seeds until interval precision supports a decision. Report minimum bulk and
tail ESS per draw and per second, not just averages. These tolerances are
specified before evaluating a change that alters the transition stream.
