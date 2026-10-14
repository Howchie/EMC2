
I think the reorganisation is a genuine improvement. Retiring PCOUNTER in favour of the revised FRQ is mathematically defensible, and I prefer the CV parameterization for the way you normally use rates in a race model.

The October 14 document establishes something stronger than simply noticing that the distributions resemble one another: it constructs an exact gamma-mixed Poisson representation of FRQ, including the continuous-threshold case, the ordinary Poisson limit, and linear self-excitation.

I wouldn't undo the restructuring. However, there are some substantive qualifications, particularly about what the counter representation means mechanistically, the new threshold-variability parameter , and the identification claims.

### 1\. The CV parameterization is the right default

Your change from to , with

\\beta=cv^{-2},\\qquad s_v=cv\\,v,

is particularly appropriate because it makes a pure time-scale parameter when are held fixed:

F(u;v)=F(vu;1).

This is exact, not approximate.

It means that changing the processing rate doesn't inadvertently change the distribution's relative variability. If you're fitting a design where experimental conditions alter evidence-processing efficiency, that's a desirable structural restriction.

By contrast, with fixed , increasing also increases , thereby making the rate distribution relatively less variable. That is a legitimate psychological hypothesis, but shouldn't happen implicitly just because you introduced an experimental effect on .

I'd therefore treat `parameterization = "cv"` as the primary implementation, retaining `sv` as an alternative when there is a substantive reason to assume constant absolute variability.

Your exact branch also matters. It gives

T_K\\sim\\operatorname{Gamma}(K,\\text{rate}=v)

for , regardless of the chosen clock shape. That is exactly the continuous-threshold counter we were looking for previously.

The qualification to your time-scale argument is that a common scaling of all rates preserves race choice probabilities only when the accumulators share a common nondecision time. Equation (1) allows different , in which case scaling decision times can change which accumulator finishes first.

### 2\. The counter representation is exact, but the two underlying processes are not identical

Proposition 1 is the strongest part of the revision:

F(u)=I_{q(u)}(K,\\beta)
=
P\\left\[G_K\\le W\\frac{q(u)}{1-q(u)}\\right\],

where and are independent.

This means FRQ can be represented as a gamma-mixed Poisson counter operating on the clock

\\phi(u)=\\frac{q(u)}{1-q(u)}.

For integer , this is an exact Poisson-event threshold. For real , the gamma threshold provides a continuous extension.

The conditional intensity result in Proposition 2 is also correct. In fact, for , it simplifies to

\\boxed{
\\lambda_n(u)=
\\frac{v(1+n/\\beta)}
{1+(c^2v/\\beta)u}.
}

This immediately shows your two familiar cases:

* : , the linear self-exciting birth process.

* : , the gamma-mixed homogeneous Poisson process after conditioning on the observed count.

But there is an important mechanistic distinction.

Consider . In the original finite-reservoir representation, three independent exponentially delayed units register with successive total rates

3\\lambda,\\quad2\\lambda,\\quad\\lambda.

In the counter representation, the successive event rates are

\\lambda,\\quad2\\lambda,\\quad3\\lambda.

Both give exactly

F(u)=(1-e^{-\\lambda u})^3,

because their completion times sum the same three independent exponential waiting-time distributions, just in a different order.

Yet their intermediate event dynamics are quite different: one is depletion, the other reinforcement.

This isn't a defect in the revised model. It establishes a stronger conclusion: the reservoir and contagious-counter mechanisms are observationally equivalent when only completion times are observed.

I'd make that distinction explicit in Section 3. The representation is an exact equality of first-passage laws, not an assertion that the reservoir's literal registrations follow the same counting process.

It also affects Section 4's interpretation of . In the counter representation, the initial mean event rate is . In a literal integer reservoir, it is

Np\\lambda=(K+\\beta-1)p\\lambda.

These can differ appreciably. For the example above, the rates are versus . Thus should specifically be called the counter-equivalent initial rate, not the actual aggregate initial registration rate under both mechanistic stories.

### 3\. There is a mathematical error involving

Section 4 states that at ,

f(0^+)=v

for every .

That's true when , but not for the new threshold-variability extension.

You define

F_\\delta(u)=H_\\delta(F_0(u)),

and derive

H_\\delta'(z)=
\\frac{\\sinh(\\delta)/\\delta}
{1+2z(1-z)(\\cosh\\delta-1)}.

At the onset , hence

\\boxed{
f_\\delta(0^+)=v\\frac{\\sinh\\delta}{\\delta},
\\qquad K=1.
}

At , the observed initial density is approximately .

The underlying counter's mean initial rate is still , but the marginal threshold variability changes the response-time distribution at onset.

This is repairable by qualifying the claim in Section 4, but it has an interpretive consequence: changes early response behaviour even at fixed , not just the distribution's spread.

There's a related correction to Table 1. To reproduce the listed classical counter distributions exactly, the table needs to specify . Otherwise its stated CDFs are transformed by .

### 4\. I'd be more cautious about the new

The new construction in Section 5 is mathematically neat. The uniform perturbation of the log-odds scale is exactly integrable, its inverse is analytical, and it works with continuous .

However, it is not quite the same kind of mechanism as varying the number of required evidence events.

You are modifying the cumulative completion probability:

F_\\delta(u)=H_\\delta(F_0(u)).

This produces a valid generative distribution, but represents variability in a probability-scale criterion, rather than explicitly varying an event-count threshold. In particular, even for integer , the resulting first-passage distribution is no longer generally the time of the th Poisson event.

My bigger concern is identification.

For small , your transformation expands as

H_\\delta(z)=z+
\\frac{\\delta^2}{6}z(1-z)(1-2z)
+O(\\delta^4).

So the effect is second-order near zero. Small threshold variability will be intrinsically difficult to identify.

Furthermore, the transform increases early and late completion-time mass at the expense of the middle of the distribution. That's precisely the broad qualitative effect for which you're retaining rate CV.

They aren't mathematically equivalent. For example, with , rate CV changes the power-law tail exponent , whereas only changes its leading multiplicative constant. But that distinction might be difficult to identify in ordinary RT samples.

Section 7.2 provides a convincing reason to drop geometric threshold excess in favour of rate CV. It doesn't yet provide evidence that the newly retained is independently identifiable from CV.

I'd want a recovery comparison with jointly free, using the complete race likelihood. Until then, keeping as the default is sensible.

### 5\. The clock-shape identification claim is too strong

In Section 7.4, you report effectively perfect posterior correlation between clock shape and rate variability. That's good evidence for fixing in your present application.

But they're not generally structurally non-identifiable.

Using the simplified intensity above, at ,

F(u)=1-
\\left(1+\\frac{c^2v}{\\beta}u\\right)^{-\\beta/c^2}.

Here only matters. Consequently, and really are exactly confounded. At , the distribution reduces to an exponential and disappears altogether.

For , however, the dependence on in the successive event rates provides additional information. The equivalence no longer holds.

I checked this numerically: two parameter configurations with identical produce identical CDFs but materially different CDFs.

So I'd describe Section 7.4 as establishing severe practical confounding under the tested designs, rather than claiming that both parameters cannot be estimated in principle. Your decision to fix the clock remains well justified.

### 6\. FRQfade has an exact overlap with ordinary FRQ

I also like FRQfade as a process interpretation: evidence-registration opportunities diminish over time, and the model's defectiveness follows from a finite operational clock instead of a directly fitted .

But there is an exact equivalence worth documenting.

When , FRQfade gives

F_{\\mathrm{fade}}(u)=
H_\\delta\\left\[
P_K\\left(
\\frac v\\kappa(1-e^{-\\kappa u})
\\right)\\right\].

Now take ordinary defective FRQ, again with , but set its registration clock to . Its CDF is

F_{\\mathrm{FRQ}}(u)=
H_\\delta\\left\[P_K\\left(E(1-e^{-\\lambda u})\\right)\\right\].

Setting

E=\\frac v\\kappa,\\qquad
\\lambda=\\kappa,\\qquad
h=H_\\delta\\left\[P_K(v/\\kappa)\\right\]

makes these distributions exactly identical at every time.

That means FRQfade is already contained in the defective FRQ family when rate variability is absent and the appropriate registration clock is allowed.

This doesn't make FRQfade pointless. It gives a more direct explanation of the fading mechanism, and the two families need not remain equivalent when . But any claim that fading and finite evidence availability have been empirically distinguished must account for this overlap.

### 7\. Was deleting PCOUNTER justified?

For the independent-accumulator models you actually want to fit, I think yes.

FRQ now gives you a coherent continuum from the ordinary Poisson counter through gamma-mixed Poisson timing and linear self-excitation. You no longer need separate model implementations for these cases, and the rate/CV representation is easier to use with experimental designs.

One wording change matters, though. Section 8 acknowledges that the combined self-excitation plus additional rate variability model and the correlated-gamma PCOUNTER model are not retained.

Therefore, the statement in the abstract that PCOUNTER was removed as redundant is slightly too broad. Its central independent-counter submodels are redundant. The complete former model family isn't.

I would also soften the absolute assertion that no closed form exists for jointly varying rate and self-excitation. Your earlier finite-sum solution fails for real , but that does not establish the nonexistence of an analytical representation using more general special functions.

### Overall

The revised architecture is more unified than the previous FRQ/PCOUNTER split, and I'd keep it. The mathematical counter construction is sound, the continuous threshold is principled as a gamma-level extension, and CV gives the clean relationship between experimental rate manipulations and RT distribution shape that you want.

The two changes I'd prioritise before calling the document settled are correcting the interpretation of the counter-equivalent rate and qualifying the claims about . The clock identification wording and the FRQfade overlap are comparatively straightforward additions.

The unresolved modelling decision, to my mind, is whether you actually need two variability mechanisms, and . You've demonstrated that CV adds meaningful RT behaviour and can be estimated when sufficiently large. You haven't yet established the same for .

For the N-back application, I'd start with the revised FRQ at , , shared across racers unless the data justify otherwise, and . That's a relatively interpretable model, and its stronger restrictions would make any subsequent improvement from freeing much more informative.
