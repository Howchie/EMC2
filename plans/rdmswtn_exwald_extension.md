# RDMSWTN ex-Wald extension

For accumulator $i$, let $D_i$ be a Wald first-passage time with threshold $a$, drift $v>0$, and diffusion SD $s$. Add independent hesitation $H_i\sim\operatorname{Exp}(\text{mean }\tau)$, with $\tau=hF$, and set $T_i=D_i+H_i$. Write $\lambda=1/\tau$.

The Wald density and CDF are

\[
f_W(t\mid a,v,s)=\frac{a}{s\sqrt{2\pi t^3}}\exp\!\left[-\frac{(a-vt)^2}{2s^2t}\right],
\]
\[
F_W(t\mid a,v,s)=\Phi\!\left(\frac{vt-a}{s\sqrt t}\right)+e^{2av/s^2}\Phi\!\left(-\frac{vt+a}{s\sqrt t}\right).
\]

Define \(\kappa=\sqrt{v^2-2s^2\lambda}\). For real \(\kappa\), the ex-Wald density and CDF are

\[
f_X(t)=\lambda e^{-\lambda t+a(v-\kappa)/s^2}F_W(t\mid a,\kappa,s),
\qquad
F_X(t)=F_W(t\mid a,v,s)-\tau f_X(t).
\]

When \(v^2<2s^2\lambda\), evaluate the density using the real Faddeeva-function branch. For RDMSWTN variability, average these conditional densities and CDFs over its drift and threshold distributions. For independent accumulators, racer \(i\)'s winning density is

\[
g_i(t)=f_{X_i}(t)\prod_{j\ne i}\bigl[1-F_{X_j}(t)\bigr].
\]

## References

- Schwarz, W. (2001). The ex-Wald distribution as a descriptive model of response times. *Behavior Research Methods, 33*, 457–469. https://doi.org/10.3758/BF03195403
- Heathcote, A. (2004). Fitting Wald and ex-Wald distributions to response time data. *Behavior Research Methods, 36*, 678–694. https://doi.org/10.3758/BF03206550
