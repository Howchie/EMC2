# joint A + sv density, and assessment of guess/kill cases

## 1. setup

Let the Wald density be

\[
f_W(t;d,\mu,\sigma)
=
\frac{d}{\sqrt{2\pi\sigma^2 t^3}}
\exp\!\left[-\frac{(d-\mu t)^2}{2\sigma^2 t}\right],
\qquad t>0.
\]

Assume start-point variability in threshold distance

\[
d \sim U[d_L,d_H],
\qquad d_L=b-A,
\qquad d_H=b,
\]

and truncated-normal drift variability

\[
\mu \sim TN(\mu_0,sv^2,c),
\qquad
h(\mu)=\frac{\phi((\mu-\mu_0)/sv)}{sv\,Z},
\qquad
Z=\Phi\!\left(\frac{\mu_0-c}{sv}\right),
\qquad \mu\ge c.
\]

Then the joint density is

\[
f_{A,sv}(t)
=
\frac{1}{A}
\int_{d_L}^{d_H}
\int_c^{\infty}
 f_W(t;d,\mu,\sigma)
 h(\mu)
\,d\mu\,dd.
\]

Using

\[
f_W(t;d,\mu,\sigma)=\frac{d}{t}\,\varphi(d;\mu t,\sigma^2 t),
\]

where \(\varphi(\cdot;m,s^2)\) is a normal density, this becomes

\[
f_{A,sv}(t)
=
\frac{1}{A Z t}
\int_{d_L}^{d_H}
\int_c^{\infty}
 d\,\varphi(d;\mu t,\sigma^2 t)
\varphi(\mu;\mu_0,sv^2)
\,d\mu\,dd.
\]

## 2. latent bivariate-normal representation

Define the latent Gaussian pair

\[
\mu \sim N(\mu_0,sv^2),
\qquad
D_t=\mu t + \sigma\sqrt{t}\,\varepsilon,
\qquad
\varepsilon\sim N(0,1),
\]

with \(\varepsilon\) independent of \(\mu\). Then \((D_t,\mu)\) is bivariate normal with

\[
m_D = E[D_t] = \mu_0 t,
\qquad
s_D^2 = \operatorname{Var}(D_t) = sv^2 t^2 + \sigma^2 t,
\qquad
s_D=\sqrt{sv^2 t^2 + \sigma^2 t},
\]

and correlation

\[
\rho
=
\frac{\operatorname{Cov}(D_t,\mu)}{s_D\,sv}
=
\frac{sv^2 t}{s_D\,sv}
=
\frac{sv\sqrt{t}}{\sqrt{\sigma^2 + sv^2 t}}.
\]

Therefore

\[
\boxed{
f_{A,sv}(t)
=
\frac{1}{A Z t}
E\left[D_t\,\mathbf 1\{d_L\le D_t\le d_H,\ \mu\ge c\}\right].
}
\]

So the density is a first moment of a truncated bivariate normal over a rectangle.

## 3. standardised form

Standardise

\[
Y = \frac{D_t-m_D}{s_D},
\qquad
X = \frac{\mu-\mu_0}{sv},
\]

so \((Y,X)\) is standard bivariate normal with correlation \(\rho\).

Let

\[
\alpha_L = \frac{d_L-m_D}{s_D},
\qquad
\alpha_H = \frac{d_H-m_D}{s_D},
\qquad
\gamma = \frac{c-\mu_0}{sv}.
\]

Then

\[
D_t = m_D + s_D Y,
\]

and hence

\[
f_{A,sv}(t)
=
\frac{1}{A Z t}\left[m_D P_R + s_D M_R\right],
\]

where

\[
P_R = P(\alpha_L \le Y \le \alpha_H,\ X \ge \gamma),
\]

and

\[
M_R = E\left[Y\,\mathbf 1\{\alpha_L \le Y \le \alpha_H,\ X \ge \gamma\}\right].
\]

## 4. probability term

Using the standard bivariate normal CDF \(\Phi_2(\cdot,\cdot;\rho)\),

\[
\boxed{
P_R
=
\bigl[\Phi(\alpha_H)-\Phi(\alpha_L)\bigr]
-
\bigl[\Phi_2(\alpha_H,\gamma;\rho)-\Phi_2(\alpha_L,\gamma;\rho)\bigr].
}
\]

This is just the probability of the strip \(\alpha_L \le Y \le \alpha_H\) with the truncation \(X\ge\gamma\).

## 5. first-moment term

Conditioning on \(Y=y\),

\[
P(X\ge\gamma\mid Y=y)
=
\Phi\!\left(\frac{\rho y-\gamma}{\sqrt{1-\rho^2}}\right).
\]

Therefore

\[
M_R
=
\int_{\alpha_L}^{\alpha_H}
 y\phi(y)
 \Phi\!\left(\frac{\rho y-\gamma}{\sqrt{1-\rho^2}}\right)
\,dy.
\]

Integrating by parts gives the endpoint form

\[
\boxed{
M_R
=
\Biggl[
-\phi(y)
\Phi\!\left(\frac{\rho y-\gamma}{\sqrt{1-\rho^2}}\right)
+
\rho\,\phi(\gamma)
\Phi\!\left(\frac{y-\rho\gamma}{\sqrt{1-\rho^2}}\right)
\Biggr]_{y=\alpha_L}^{\alpha_H}.
}
\]

The second endpoint coefficient is \(\rho\phi(\gamma)\), not
\(\rho\phi(\gamma)/\sqrt{1-\rho^2}\): the derivative contributes
\(\rho/\sqrt{1-\rho^2}\), while the product-of-Gaussians integral contributes
the cancelling \(\sqrt{1-\rho^2}\).

So the full density is

\[
\boxed{
f_{A,sv}(t)
=
\frac{1}{A Z t}
\left[
 m_D P_R + s_D M_R
\right],
}
\]

with the explicit forms for \(P_R\) and \(M_R\) above.

This is a fully analytic density for the joint \(A>0, sv>0\), no-clock case.
No GL quadrature is required for this density.

## 6. implementation sketch

```cpp
double dens_A_sv(double t, double b, double A,
                 double mu0, double sv, double c,
                 double sigma) {

    double dL = b - A;
    double dH = b;

    double Z = Phi((mu0 - c) / sv);

    double mD = mu0 * t;
    double sD = sqrt(sv*sv*t*t + sigma*sigma*t);
    double rho = (sv * sqrt(t)) / sqrt(sigma*sigma + sv*sv*t);
    double s = sqrt(1.0 - rho*rho);

    double aL = (dL - mD) / sD;
    double aH = (dH - mD) / sD;
    double gam = (c - mu0) / sv;

    double PR = (Phi(aH) - Phi(aL))
              - (Phi2(aH, gam, rho) - Phi2(aL, gam, rho));

    auto G = [&](double y) {
        return -phi(y) * Phi((rho*y - gam) / s)
             + rho * phi(gam) * Phi((y - rho*gam) / s);
    };

    double MR = G(aH) - G(aL);

    return (mD * PR + sD * MR) / (A * Z * t);
}
```

## 7. assessment of the other cases

### (a) A = 0, sv > 0, no Erlang

Yes: your current approach is already right.

- SWTN density is analytic.
- SWTN CDF is analytic via the current transform / bivariate-normal trick.

### (b) A > 0, sv > 0, no Erlang

Split by target:

- **Density**: there is an analytic form, namely the one above. So here GL is not fundamentally necessary.
- **CDF**: your current GL strategy still looks right. The image term

\[
\exp\left(\frac{2\mu d}{\sigma^2}\right)
\Phi\left(-\frac{d+\mu t}{\sigma\sqrt t}\right)
\]

introduces a bilinear \(\mu d\) interaction, which destroys the clean Gaussian closure.

### (c) A = 0, sv > 0, Erlang kill / guess variants

Yes: your current strategy is still the right one.

For fixed \((d,\mu)\), the killed Wald pieces are analytic. But once you integrate over truncated-normal drift, the tilt introduces

\[
q(\mu)=\sqrt{\mu^2 + 2\sigma^2 r},
\]

and terms like

\[
\exp\left[\frac{a}{\sigma^2}(\mu-q(\mu))\right]
\Phi\left(\frac{q(\mu)U-a}{\sigma\sqrt U}\right),
\]

which no longer have the Gaussian-completion structure that made SWTN analytic.

So for the finite-window kill / guess CDF pieces with \(sv>0\), numerical integration over drift is still the correct call.

### (d) A > 0, sv > 0, Erlang kill / guess variants

Yes: your current approach is still the right one.

Once both start-point variability and drift variability are active and the Erlang tilt is present, the integrand contains both the bilinear \(\mu d\) structure from the Wald image term and the nonlinear square-root tilt structure from

\[
q(\mu)=\sqrt{\mu^2+2\sigma^2 r}.
\]

That is well beyond the clean Gaussian / bivariate-normal family. So GL (or similar quadrature) remains the appropriate bridge.

## 8. bottom line

The only obvious missed analytic opportunity is:

\[
\boxed{
A>0,\ sv>0,\ \text{density only}
}
\]

which can be done in closed form using a truncated bivariate-normal first moment.
This is implemented for the standard `drdmswtn` density path when
\(\lambda_g=\lambda_k=0\). The CDF and Erlang guess/kill variants still use the
existing quadrature paths.

For the rest:

- non-killed SWTN density and CDF: you already do the right thing;
- joint \(A+sv\) CDF: GL is still justified;
- killed / guess variants with \(sv>0\): GL is still justified.
