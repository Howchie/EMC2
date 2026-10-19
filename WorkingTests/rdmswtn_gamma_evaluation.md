# Evaluation of `RDMSWTN(drift_dist = "gamma")` Implementation

**Date:** 2026-10-19  
**Commit evaluated:** `8529ae34` (`RDMSWTN: Gamma drift law (drift_dist = "gamma")`) and worktree files (`src/rdmswtn_gamma.h`, `src/model_RDMSWTN_gamma.h`, `src/model_RDMSWTN_gamma.cpp`, `R/model_RDM.R`, `Math/RDMSWTN_Derivation.tex` §4, `tests/testthat/test-rdmswtn-gamma.R`).

---

## 1. Mathematical Derivation

The mathematical architecture in `Math/RDMSWTN_Derivation.tex` (§4, "Gamma Drift Variability", lines 283–664) and `src/rdmswtn_gamma.h` is sound and well-structured:

1. **Dimensionless reduction and cylinder moments (`src/rdmswtn_gamma.h:7-25`)**:
   - With $V \sim \operatorname{Gamma}(k, \text{rate }\tau)$ where $k = (v/s_v)^2$ and $\tau = v/s_v^2$, scaling drift to $Y = V\sqrt{t}/s \sim \operatorname{Gamma}(k, \text{rate }p)$ with $p = \tau s/\sqrt{t}$, $\delta = d/(s\sqrt{t})$, $z = p - \delta$, and $q = p/\delta = \tau s^2/d$ reduces the 5-parameter integral to $(k, p, \delta)$.
   - Completing the square yields $f(t) = (\delta/t)\mathcal{D}_k$ with
     $$\mathcal{D}_\nu = \mathbb{E}_{\operatorname{Gamma}(\nu, p)}[\phi(Y - \delta)] = \frac{p^\nu e^{-\delta^2/2}}{\sqrt{2\pi}} G_\nu(z), \qquad G_\nu(w) = \frac{J_\nu(w)}{\Gamma(\nu)} = e^{w^2/4} D_{-\nu}(w).$$
   - Integration by parts on $J_\nu(x) = \int_0^\infty y^{\nu-1} e^{-y^2/2 - x y}\,dy$ gives the division-free recurrence
     $$G_\nu(x) = (\nu + 1)G_{\nu+2}(x) + x G_{\nu+1}(x).$$
     For $x \ge 0$, every coefficient is nonnegative and $G_\nu$ is the minimal solution, so downward recursion (`GWalk::step`, lines 122–125) is unconditionally stable.
   - **Implicit invariant in `logG_asym` (lines 90–108)**: Line 107 calls `stirling(nu)` directly, which is only accurate for $\nu \ge 12$ (unlike `lgam(a)` on lines 73–78, which shifts $a < 12$ up to $\ge 12$). Because every caller of `logG_asym` has $x \ge 0$ and $\kappa = \nu + y_0^2 \ge \text{KAPPA\_MIN} = 64$, and $x \ge 0 \implies y_0 = \frac{2\nu}{x + \sqrt{x^2 + 4\nu}} \le \sqrt{\nu} \implies \kappa \le 2\nu$, the condition $\nu \ge \kappa/2 \ge 32 \ge 12$ holds automatically.

2. **Complementary positive series across $q = 1$**:
   - **NB branch ($q \le 1$, `nb_branch`, lines 270–420)**:
     - At $q = 1$ ($z = 0$), Legendre's duplication formula reduces $G_k(0) = \sqrt{\pi}/(2^{k/2}\Gamma((k+1)/2))$ and the decision time is exactly $\operatorname{InvGamma}((k+1)/2, d^2/(2s^2))$.
     - For $q < 1$, expanding $\operatorname{Gamma}(k, \tau) = \sum_{n\ge 0} w_n \operatorname{Gamma}(k+n, d/s^2)$ with negative-binomial weights $w_n = \frac{(k)_n}{n!} q^k (1-q)^n$ expresses $f$, $F$, and $S$ as positive mixtures of inverse-gamma laws with shape $\alpha_n = (k+n+1)/2$ and $A = \delta^2/2$.
     - Splitting $n$ by parity makes $\alpha_n$ step by $1$ within each parity stream. Using $P(\alpha, A) = \iota(\alpha) + P(\alpha+1, A)$ and summing by parts with cumulative same-parity weights $W_m = \sum_{n_0 \le n \le m,\, n\equiv m} w_n$ (lines 332, 406–407) evaluates $S$ in a single pass with only two `lower_ratio` calls at the top of the two parity streams.
   - **Moment branch ($q > 1$, `moment_logS` and `moment_logF`, lines 422–581)**:
     - Gamma integration by parts $\mathbb{E}_{\nu+1}[h] - \mathbb{E}_\nu[h] = p^{-1}\mathbb{E}_{\nu+1}[h']$ applied to $F_W(Y)$ cancels the Gaussian derivatives and leaves $\partial_Y F_W(Y) = 2\delta e^{2\delta Y}\Phi(-Y-\delta)$. Applying integration by parts a second time to the reflection term yields a first-order recurrence in shape whose backward solution telescopes to
       $$S = \frac{1}{p}\sum_{m=2}^\infty \bigl(1 - r^{m-1}\bigr)\mathcal{D}_{k+m}, \qquad r = 1 - \frac{2}{q} \in (-1, 1).$$
       Because $|r| < 1$ for all $q > 1$, $1 - r^{m-1} > 0$ for all $m \ge 2$ even when $r \in (-1, 0)$ ($1 < q < 2$).
     - For small $F$, the symmetric Mills identity $e^{-\delta y} F_W(y) = \sqrt{2\pi}\,\phi(y)\phi(\delta)[R(\delta-y) + R(\delta+y)]$ and $R(\delta\mp y) = \int_0^\infty e^{-\delta u - u^2/2 \pm y u}\,du$ cancel odd powers of $y$ and give the positive even-order series
       $$F = \frac{2 p^k e^{-\delta^2/2}}{\sqrt{2\pi}} \sum_{n\text{ even}} W_n I_n(\delta), \qquad W_n = \frac{J_{k+n}(z)}{\Gamma(k)}, \quad I_n(\delta) = G_{n+1}(\delta).$$
   - **Narrow branch (`narrow_branch`, lines 594–672)**:
     - Expands $\phi(y-\delta)$ and $\partial_y S_W(y) = -2\delta\phi(y-\delta)R(y+\delta)$ around $y_0 = k/p$ against the exact central moments $\mu_n$ of $\operatorname{Gamma}(k, p)$, forming the Taylor coefficients of $\phi(x_0+h)R(w_0+h)$ via the Cauchy product of Hermite coefficients $a_i$ and $(-1)^j I_j(w_0)$.
   - **Wide / tilted EMG branch (`wide_branch` and `tilt_branch`, lines 674–847)**:
     - Uses $-\partial_y S_W(y) = \int_0^\infty 2\delta e^{-2\delta u}\phi(y - \delta + u)\,du$, the density of $\delta - U + Z$ with $U \sim \operatorname{Exp}(2\delta) \perp Z \sim \mathcal{N}(0, 1)$, so $S_W(y) = P(U - Z < \delta - y)$ for all $y \in \mathbb{R}$ and $S = \mathbb{E}_X[G_Y(\delta - X)]$ with $X = U - Z$.
     - When the log-slope $\psi_1 = (k-1)/\delta - p$ of $g_Y$ at $\delta$ is steep ($|\psi_1| > 1$), `tilt_branch` factors $g_Y(\delta+h) = g_Y(\delta)e^{\psi_1 h}\varphi(h)$ and exponentially tilts $A = Z - U$ by $e^{\psi_1 A}$ so $Z \sim \mathcal{N}(\psi_1, 1)$ and $U \sim \operatorname{Exp}(2\delta + \psi_1)$.

---

## 2. Accuracy & Numerical Issues Found

### A. Critical bug in simulator `rRDMSWTN_gdrift` for $k = (v/s_v)^2 \ll 1$ (spurious omissions & mantissa loss)

In `R/model_RDM.R` (lines 1065–1083):
```r
v <- p$v[ok]; sv <- p$sv[ok]
vary <- is.finite(sv) & sv > 1e-10 * v
V <- v
V[vary] <- rgamma(sum(vary), shape = (v[vary] / sv[vary])^2, rate = v[vary] / sv[vary]^2)
dt[ok] <- rWald(sum(ok), B = pmax(p$b[ok], 0), v = V, A = rep(0, sum(ok)), s = p$s[ok]) + p$t0[ok]
```

1. **Spurious `Inf` (omissions) when `rgamma` underflows to `0.0`**:
   When $k = (v/s_v)^2 < 1$ (for example, row 4 of the test grid in `tests/testthat/test-rdmswtn-gamma.R:80-82`: $v = 0.05, s_v = 4 \implies k = 1.5625 \times 10^{-4}$), `rgamma` in R frequently underflows to `0.0` in IEEE-754 double precision:
   $$P(V = 0) \approx (\text{DBL\_TRUE\_MIN} \cdot \tau)^k \approx (5 \times 10^{-324})^{1.5625 \times 10^{-4}} \approx 88.7\%.$$
   Because `rWald` (`R/model_RDM.R:30-99`) defaults to `posdrift = TRUE`, line 70 sets `pos <- v > 0` and treats every `V == 0` draw as a defective non-finishing accumulator (`rt = Inf`), instead of routing `V == 0` to `rlevy` (the driftless Brownian hitting time $\operatorname{InvGamma}(1/2, b^2/(2s^2))$).
   - **Reproducer**:
     ```r
     lR <- factor(rep("a", 10000), levels = "a")
     pars <- cbind(v = 0.05, B = 1.0, t0 = 0.0, s = 1.0, sv = 4.0)[rep(1, 10000), ]
     sim <- EMC2:::rRDMSWTN_gdrift(lR, pars)
     mean(is.infinite(sim$rt))
     # [1] 0.8873  (88.73% spurious omissions; exact model is proper with 0% omissions)
     mean(sim$rt <= 5.0)
     # [1] 0.0740  (vs exact pRDMSWTN_gdrift_cpp(5.0, 0.05, 1.0, 0.0, 1.0, 4.0) = 0.6551)
     ```
2. **Severe mantissa loss in `rWald` / `rwald_fpt_r` for small positive $V$**:
   In `rWald` (`R/model_RDM.R:55-56`) and `rwald_fpt_r` (`src/model_rng.h:71-72`), the Michael–Schucany–Haas step evaluates:
   ```r
   x.0 <- mu.0 + mu.0^2 * y / (2 * lambda.0) -
     sqrt(4 * mu.0 * lambda.0 * y + mu.0^2 * y^2) * mu.0 / (2 * lambda.0)
   ```
   When $10^{-6} < V \ll 1$, $\mu_0 = b/V$ is large and $\mu_0^2 y / (2\lambda_0) \propto V^{-2}$, so subtracting the square root loses up to 39 of 53 mantissa bits (at $V = 1.1\times 10^{-6}$, simulated RTs are quantized to integer multiples of $2^{-13}$).
   - **Fix**: Rationalize the smaller root of the Michael–Schucany–Haas quadratic:
     $$x_0 = \frac{2 b^2}{2 b V + s^2 y + \sqrt{s^2 y (4 b V + s^2 y)}}, \qquad \text{test} = \frac{b}{b + V x_0}.$$
     This algebraically identical form has no subtraction, no division by $V$, needs no arbitrary `tiny = 1e-6` cutoff or `qnorm` call, and is smooth and accurate to full double precision for all $V \ge 0$ down to and including $V = 0$ (where $x_0 = b^2/(s^2 Z^2)$ and $\text{test} = 1$).

### B. Precision loss & failure in `moment_logS` for large $q$ under non-`-ffast-math` builds (`compatibility` profile)

In `src/rdmswtn_gamma.h` (lines 434, 463–465):
```cpp
const double r = 1.0 - 2.0 / q, lp = std::log(p);
...
const double rp = std::pow(r, mlo - 2);
const double head = (r > 0.0) ? -std::expm1((mlo - 2) * std::log(r)) : 1.0 - rp;
const double sum = T * head + (1.0 - r) * rp * H;
```
- When $q = p/\delta = \tau s^2/d$ is large, `r = 1.0 - 2.0 / q` loses precision, and `1.0 - r` on line 465 (as well as `std::log(r)` on line 464 when `mlo > 2`) suffers catastrophic cancellation.
- Under the `legacy` build profile (`-ffast-math`), GCC rewrites `1.0 - (1.0 - 2.0 / q)` as `2.0 / q` (though `std::log(r)` still loses precision when `mlo > 2`). Under the `compatibility` build profile or standard IEEE-754 compilation (e.g., `Rcpp::sourceCpp` or CRAN flags), `1.0 - r` is not simplified:
  - At $q = 10^{13}$ ($t=1, v=0.5, s_v=\sqrt{5}/10, b=10^{-12}, s=1$), `logS` loses 12 digits (`-28.51695` vs exact `-28.516703348382`).
  - At $q = 10^{16}$ ($b = 10^{-15}$), `S` has $>10\%$ relative error.
  - At $q \ge 10^{17}$ ($b = 10^{-16}$), `r == 1.0`, `sum == 0.0`, all 6 window attempts fail, and `moment_logS` returns `NaN` / `-Inf` after 1,370 terms.
- **Fix**: Define `const double omr = 2.0 / q;` (`1 - r`) and `const double lr = (r > 0.0) ? std::log1p(-omr) : NAN_D;`, using `omr` directly in `sum = T * head + omr * rp * H` and `lr` in `head = -std::expm1((mlo - 2) * lr)`.

### C. Unclamped `o.logS` / `o.logF` before `log1m_exp` in `gamma_wald`

In `nb_branch` (lines 402–413), `o.logF` and `o.logS` are clamped via `std::min(0.0, ...)` before calling `log1m_exp`.
In `gamma_wald` (lines 942–958, `MOMENT` branch):
```cpp
o.logS = moment_logS(k, p, delta, z, q, want & WANT_PDF, lDk, t1);
o.logf = ldt + lDk;
o.logF = log1m_exp(o.logS);
if ((want & WANT_CDF) && o.logF < std::log(TAIL_SMALL)) { ... }
```
If $S$ is very close to $1$ and `moment_logS` rounds 1 ULP above `0.0` (e.g., `+2e-16`), `log1m_exp(o.logS)` evaluates `std::log(-std::expm1(a))` on a negative argument and returns `NaN`. Because `o.logF` is `NaN`, `o.logF < std::log(TAIL_SMALL)` on line 958 evaluates to `false`, skipping the direct `moment_logF` call and returning `NaN`.
- **Fix**: Clamp `o.logS = std::min(0.0, ...)` and `o.logF = std::min(0.0, ...)` before `log1m_exp`.

---

## 3. Optimisation of Implementation

### A. Analytic elimination of $O(k\varepsilon)$ cancellation in $\log\mathcal{D}_\nu$ (removes `DENS_DIRECT_KMAX = 2000`)

In `gamma_wald` (lines 910–917), direct $O(1)$ evaluation of `WANT_PDF` for $q > 1$ is restricted to `k <= DENS_DIRECT_KMAX` ($2000$):
```cpp
if (want == WANT_PDF && force == 0 && q > 1.0 && k <= DENS_DIRECT_KMAX) {
  o.logf = std::log(delta / t) + k * std::log(p) - 0.5 * delta * delta - LOG_SQRT_2PI + logG(k, z);
  o.branch = MOMENT;
  return o;
}
```
Evaluating $k\log p - \frac{1}{2}\delta^2$ and $\log G_k(z) = k\log(y_0/k) + \frac{1}{2}y_0^2 + \dots$ separately loses $O(k\varepsilon)$ absolute accuracy (and the same cancellation occurs in the seed `ns * lp + logG_asym(ns, z) + dconst` of `moment_logS` at line 445). Furthermore, when $k > 2000$ and `want == WANT_PDF` falls outside `narrow_branch`'s gate (e.g., $\operatorname{sd}(Y) > 0.3$ or $|y_0 - \delta| > p/3$), `gamma_wald` falls through to `moment_logS`, summing the entire survivor window `[mlo, mhi]` even though only the PDF was requested, and when `mlo > 48` (line 485) still falls back to `k * lp + logG(k, z) + dconst`.

The cancellation between $\nu\log p - \frac{1}{2}\delta^2$ and $\log G_\nu(z)$ can be removed **analytically**:
Since the saddle point $y_0 = \operatorname{saddle}(\nu, z)$ satisfies $y_0^2 + (p - \delta)y_0 = \nu$, we have $p y_0 = \nu + y_0(\delta - y_0)$, so
$$\frac{p y_0}{\nu} = 1 + u, \qquad u := \frac{y_0(\delta - y_0)}{\nu}.$$
Using $\nu u = y_0(\delta - y_0)$, the linear part of $\nu\log(1 + u)$ combines with $\frac{1}{2}(y_0^2 - \delta^2)$ into an exact square:
$$\nu\log\left(\frac{p y_0}{\nu}\right) + \frac{1}{2}(y_0^2 - \delta^2) = -\frac{1}{2}(y_0 - \delta)^2 + \nu\bigl[\log(1 + u) - u\bigr].$$
Therefore, for $\kappa(\nu, z) \ge \text{KAPPA\_MIN}$:
$$\log\mathcal{D}_\nu = -\frac{1}{2}(y_0 - \delta)^2 + \nu\bigl[\log(1 + u) - u\bigr] - \frac{1}{2}\log(2\pi) + \frac{1}{2}\log\alpha + \log(1 + \text{acc}) - \operatorname{stirling}(\nu).$$
As $s_v \to 0$ ($k \to \infty$), $y_0 \to v\sqrt{t}/s = O(1)$, $-\frac{1}{2}(y_0 - \delta)^2 \to -\frac{(vt - d)^2}{2s^2 t}$ (the exact Wald exponent), $u = O(1/k)$, and $k[\log(1+u)-u] = O(1/k) \to 0$ with **zero cancellation for any $k > 0$**. Using this in `gamma_wald` and `moment_logS` removes `DENS_DIRECT_KMAX` and makes all $q > 1$ `WANT_PDF` evaluations $O(1)$.

### B. `narrow_branch` ignores `want == WANT_PDF` and uses an unscaled stopping threshold for `tc`

In `src/rdmswtn_gamma.h` (lines 594–672):
1. `narrow_branch` ignores `want` (`(void)want;` on line 670). When called with `want == WANT_PDF`, it still runs `walk_from(NT + 2.0, w)`, populates `I[0..NT+1]`, computes the $O(\text{NT}^2)$ Cauchy dot product `for (int i = 0; i < n; ++i) g += a[i] * bb[i];` at every step, and calls `R::pnorm` in `mills(x)`, even though `o.logf` depends only on the $O(\text{NT})$ 3-term recurrences for `a[n]` and `mu`.
2. On line 639:
   ```cpp
   if (n >= 6 && std::fabs(td) <= EPS * std::fabs(dens) && std::fabs(tc) <= EPS * 1e-3 * std::max(1.0, std::fabs(C)))
   ```
   Because `C` is subtracted from `SW` (when $x > 0$) or added to `Fs0 = mills(-x) + I[0]` (when $x \le 0$), and $\text{SW} = R(x) - R(w) \approx 2\delta/(x w)$ scales with $\delta$, comparing `tc` against `1e-3 * std::max(1.0, std::fabs(C))` is overly tight by $10^3$ when $|C| \ge 1$ and unscaled when $\delta \ll 1$. Computing `base_c = (x > 0.0) ? SW : (mills(-x) + I[0])` before the loop allows testing `std::fabs(tc) <= EPS * std::fabs((x > 0.0) ? base_c - C : base_c + C)` directly.

### C. Tightening and accelerating `far_tail` via the Chernoff–Laplace saddle bound

In a 10,000-draw random stress test across $t \in [10^{-3}, 20]$, $v \in [10^{-3}, 20]$, $s_v \in [10^{-4}, 20]$, $b \in [0.05, 20]$, $s \in [0.05, 5]$, the 10 slowest calls to `gamma_wald` (requiring 5,760–11,041 terms and up to 235 µs) had true log-tails between `-1018` and `-3510` (well below `LOG_NEGLIGIBLE = -1000`), yet `far_tail` (lines 861–892) did not trigger on any of them.
- **Cause**: `far_tail` splits the expectation at $c$ into $\Phi(-h) + P(Y \gtrless c)$, which replaces one factor by $1$ on each half-line (losing at least half of the true exponential decay rate) and chooses $c$ using Gaussian-like weights `w_lo = sdY / (1 + sdY)` that do not reflect the linear tail decay $e^{-p c}$ of the Gamma distribution. It also makes up to 4 special-function calls (`R::pgamma` and `R::pnorm` twice).
- **Exact Chernoff–Laplace saddle bound**:
  For all $h > 0$ and $w \in \mathbb{R}$, the normal tail and density satisfy
  $$2\Phi(-w) \le e^{-w^2/2} \le \exp\left(\frac{1}{2}h^2 - h w\right) \quad (w \ge 0), \qquad 1 < \exp\left(\frac{1}{2}h^2 - h w\right) \quad (w < 0),$$
  $$\phi(w) = \frac{1}{\sqrt{2\pi}}\exp\left(\frac{1}{2}h^2 - h w - \frac{1}{2}(w - h)^2\right) \le \frac{1}{\sqrt{2\pi}}\exp\left(\frac{1}{2}h^2 - h w\right).$$
  Taking expectation over $Y \sim \operatorname{Gamma}(k, p)$ at the saddle point
  $$c = \begin{cases} \dfrac{2k}{z + \sqrt{z^2 + 4k}}, & z \ge 0, \\[6pt] \dfrac{-z + \sqrt{z^2 + 4k}}{2}, & z < 0, \end{cases}$$
  which always lies strictly between $y_0 = k/p$ and $\delta$ and satisfies $p - (\delta - c) = c + z = k/c > 0$, yields the global upper bound valid for all $k, p, \delta > 0$:
  $$\log\left(\frac{t\sqrt{2\pi}}{\delta} f(t)\right) \le \Psi(k, p, \delta), \qquad \log F(t) \le \Psi(k, p, \delta) \ (y_0 < \delta), \qquad \log S(t) \le \Psi(k, p, \delta) \ (y_0 > \delta),$$
  where
  $$\Psi(k, p, \delta) = -\frac{1}{2}(c - \delta)^2 + k\bigl[\log(1 + u) - u\bigr], \qquad u = \frac{c(\delta - c)}{k}, \quad 1 + u = \frac{p c}{k} > 0.$$
  This bound costs only one `std::sqrt` and one `std::log1p` (no `pgamma` or `pnorm` calls) and is tight to within a few log-units (e.g., giving `-3504.0` where the true $\log F = -3510.1$, whereas the split bound gave `-350`).

### D. Minor parity waste in `moment_F_length`

In `moment_F_length` (lines 504–516), the convergence check runs when `(nn & 7) == 7` (odd `nn`) and returns `nn + 8`, which is always **odd**. Because `moment_logF` (lines 552–560) only sums even indices (`if (n & 1) return;`), on `attempt == 0` the top step `N` is always odd and discarded. Returning `nn + 9` ensures `N` is even.

---

## 4. Completeness & Package Wiring

| Area | Status | Details |
|---|---|---|
| **Constructor & option validation** (`R/model_RDM.R:852-866, 1016-1043`) | Complete | `RDMSWTN(drift_dist = "gamma")` rejects `erlang_type != "none"`, `erlang_shape != 1L`, `posdrift = FALSE`, and `correlated = TRUE`. Includes `pContaminant` and `pGuess` via `add_nuisance_pars()`. |
| **Column registry & race adapter** (`src/col_registry.h:580-588`, `src/race_dispatch.cpp:167-176`) | Complete | `RDMSWTN_GDRIFT` is matched before `RDMSWTN` in `resolve_race_model_adapter`, wiring scalar (`drdmswtn_gdrift_scalar`, `prdmswtn_gdrift_scalar`), raw (`drdmswtn_gdrift_raw`, `prdmswtn_gdrift_raw`), and truncation (`rdmswtn_gdrift_logS_at_t`) kernels. |
| **Particle likelihood & WAIC/LOO** (`src/particle_ll.cpp:1371-1372`) | Complete | Excluded from the Gaussian `RDMSWTN` column-index check in `calc_ll_oo`; works in both `calc_ll_oo` and `calc_ll_oo_pw`. |
| **R-level `dfun` / `pfun` recycling** (`R/model_RDM.R:1051-1063`) | Minor gap | `dRDMSWTN_gdrift` and `pRDMSWTN_gdrift` call `rep(rt, length.out = nrow(pars))`. Unlike `dRDMSWTN` (`R/model_RDM.R:1830-1840`) and `dRDMGBM` (`R/model_RDM.R:242-249`), passing a vector `rt` (`length(rt) > 1`) with a 1-row or vector `pars` truncates `rt` to length 1 instead of expanding `pars` to `length(rt)`. |
| **Simulation (`rRDMSWTN_gdrift`)** (`R/model_RDM.R:1065-1083`) | Bug for $k \ll 1$; R-only | 1. Suffers from the `V == 0` underflow omission bug and small-$V$ cancellation in `rWald` described in §2A.<br>2. Does not call `.apply_timed_guess_winner(out, levels(lR))` at the end like `rRDM` (`R/model_RDM.R:133`).<br>3. Implemented in R only (no `rrdmswtn_gdrift_cpp` in `src/model_rng.cpp` / `R/model_rng.R`). |
| **Git tracking of `Math/RDMSWTN_Derivation.tex`** | Note | Commit `8529ae34` references `Math/RDMSWTN_Derivation.tex, "Gamma Drift Variability"` in its commit message and in `src/rdmswtn_gamma.h:5`, and the section is present in the worktree (`Math/RDMSWTN_Derivation.tex:283-664`), but `Math/RDMSWTN_Derivation.tex` is ignored by `.gitignore` and was not part of the commit. |
