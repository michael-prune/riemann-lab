# Xi exchange-monotonicity probe

## Outcome

The candidate

\[
Q_{r,k+1}<Q_{r,k},\qquad
Q_{r,k}=\frac{D_{r,k}D_{r-2,k}}
{D_{r-1,k-1}D_{r-1,k+1}}
=\frac{D_{r-1,k}^2}{D_{r-1,k-1}D_{r-1,k+1}}-1
\]

survives a 130-decimal numerical stress test for xi at ranks 2--15 and
shifts through 40.  It is false for generic PF-infinity sequences: a seeded
degree-70 polynomial with only negative real zeros violates it at every tested
rank 2--8.

There is also a rigorous but non-effective asymptotic result:

- the rank-2 inequality holds for all sufficiently large `k`;
- the rank-3 inequality holds for all sufficiently large `k`.

Those tail statements follow from published all-orders asymptotics, not from
the numerical script.  They do not reach the finite/diagonal regime needed for
RH.

## Rank 2: exact reduction

Put `f(x)=log(a(x))`, where

\[
a(x)=\frac{1}{\Gamma(2x+1)}\int_0^\infty u^{2x}\Phi(u)\,du.
\]

Then

\[
Q_{2,k}=\frac{a_k^2}{a_{k-1}a_{k+1}}-1=e^{\tau_k}-1,
\quad
\tau_k=-\log\frac{a_{k-1}a_{k+1}}{a_k^2}.
\]

Thus `Q_(2,k+1)<Q_(2,k)` is exactly

\[
a_k^3a_{k+2}>a_{k-1}a_{k+1}^3.
\]

A sufficient continuous condition is `f'''(x)>0`.  Under the tilted
probability law proportional to `u^(2x) Phi(u) du`, with `L=log(u)`,

\[
f'''(x)=8\left(\kappa_{3,x}(L)-\psi_2(2x+1)\right).
\]

Both terms on the right are negative in the samples, but the gamma
polygamma term has larger magnitude, making `f'''` positive.

Cormac O'Sullivan's all-orders Lambert-`W` expansion for the xi coefficients
has an explicit smooth main term whose third derivative is
`(2+o(1))/x^2`.  The rigorous discrete conclusion does not require
differentiating the asymptotic remainder: take sufficiently many terms and
apply the four-point difference directly to obtain

\[
\log q_{k+1}-\log q_k=\frac{2+o(1)}{k^2}>0.
\]

This proves eventual rank-2 exchange monotonicity, though no effective
starting index is extracted here.

## Rank 3: inherited curvature

Let

\[
b_k=D_{2,k}=a_k^2-a_{k-1}a_{k+1},\qquad g_k=\log b_k.
\]

Then `Q_(3,k+1)<Q_(3,k)` is the same four-point inequality for `b`.
For the continuous interpolation, put

\[
h(x)=f(x-1)+f(x+1)-2f(x),\qquad
g(x)=2f(x)+\log(1-e^{h(x)}).
\]

The probe evaluates `g'''` from the first three tilted cumulants, without a
finite-difference approximation.  It is positive at every sampled point
`x=1,2,3,5,10,20,40,80`.

The same coefficient expansion gives, at integer arguments,

\[
h(k)=-\frac{2+o(1)}k,
\qquad
3g_k+g_{k+2}-g_{k-1}-3g_{k+1}
=\frac{4+o(1)}{k^2}>0.
\]

Using one more asymptotic order controls the division by
`1-exp(h(k)) ~ 2/k`, proving eventual rank-3 monotonicity.

## Numerical stress test

The committed result uses 130 decimal digits, 18 theta terms, coefficients
through `a_62`, ranks 2--15, and shifts through 40.  Every sampled exchange
coordinate is positive and strictly decreasing.  For ranks 10--15, the
smallest relative one-step decrease is between `0.0252` and `0.0265`, so the
new high-rank signs are not near numerical zero at this precision.

Selected continuous values are:

| `x` | `f'''(x)` | `(log D_2)'''(x)` |
|---:|---:|---:|
| 1 | 0.417534394991 | 0.600061098642 |
| 5 | 0.0368825817488 | 0.0656029306981 |
| 20 | 0.00302342492961 | 0.00585599972103 |
| 80 | 0.000219202583832 | 0.000434961981622 |

These are high-precision observations, not interval enclosures.

## Prior-art boundary

The checked primary sources establish nearby but weaker facts:

- [Guo, higher Turan inequalities](https://arxiv.org/abs/2012.03530)
  gives coefficient inequalities for real-rooted polynomials, but no ordering
  of consecutive Turan quotients.
- [Planat and Sole, second-level concavity](https://arxiv.org/abs/2608.19160)
  proves a double-Turan inequality for normalized xi coefficients, not this
  four-point ratio monotonicity.
- [Michalowski, explicit cubic wedge](https://arxiv.org/abs/2607.16795)
  proves positivity of `D_(r,k)` for `k >= 10^18 r^3` and supplies analytic
  derivative bounds, but those bounds do not determine the sign of `f'''`.
- [O'Sullivan, xi coefficient asymptotics](https://arxiv.org/abs/2007.13582)
  supplies the all-orders expansion used for the two eventual statements.

The PF-infinity counterexample proves that the candidate is not a formal
consequence of total positivity or the Laguerre-Polya property.  The literature
audit found no exact global theorem matching it; that is an audit result, not
a novelty claim.

## Reproduce

```bash
uv run --with mpmath \
  experiments/xi-exchange-monotonicity-20260827/probe.py \
  --output experiments/xi-exchange-monotonicity-20260827/results/probe-v1.json
```
