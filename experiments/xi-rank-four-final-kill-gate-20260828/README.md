# Xi rank-four final kill gate

## Decision

The angle passes the stated final kill gate, narrowly but for a substantive
reason.  Rank four retains the inherited `x^-2` sign, its conservative
effective threshold is only twice the rank-three threshold, and the exact
first rank-dependent asymptotic correction is quadratic rather than
factorial or exponential.

Conditional on the same unreplayed saddle and contour lemmas used in the
rank-two and rank-three artifacts, the directed 256-bit Arb calculation gives

\[
(\log D_3(x))'''>0\qquad(x\ge2{,}000{,}000),
\]

and hence

\[
Q_{4,k+1}<Q_{4,k}\qquad(k\ge2{,}000{,}001).
\]

This is still only a fixed-rank tail result.  The favorable rank scaling is a
formal asymptotic identity, not a remainder bound uniform in rank.  Therefore
the result says **do not abandon the mechanism yet**, not that it is close to
RH.

## Exact all-rank recursion

Put `F_r(x)=log D_r(x)`, with `F_0=0`, `F_1=f=log a`, and

\[
\tau_r(x)=2F_r(x)-F_r(x-1)-F_r(x+1),\qquad
J(t)=\log(1-e^{-t}).
\]

Desnanot--Jacobi gives the exact recursion, wherever the displayed
determinants are positive,

\[
F_{r+1}=2F_r-F_{r-1}+J(\tau_r).
\]

For rank four, write `g=F_2`.  Then

\[
F_3'''=2g'''-f'''+J(\tau_g)'''
       =3f'''+2J(\tau_f)'''+J(\tau_g)'''.
\]

The inherited term is positive and of order `x^-2`; both nonlinear terms are
of order `x^-3`.

## Effective rank-four budget

The preceding source-conditional certificates supply

\[
x^2f'''(x)>1.23,
\]

and, after nested Cauchy estimates,

\[
|f^{(3)}|<\frac{37}{n^2},\quad
|f^{(4)}|<\frac{1850}{n^3},\quad
|f^{(5)}|<\frac{185000}{n^4},
\]

with the next two coefficients `27,750,000` and `5,550,000,000`.
Differentiating `J(tau_f)` explicitly through order five and normalizing each
derivative by its natural power of `n` gives, for `n>=10^6`,

\[
|g'''|<\frac{76.011}{n^2},\quad
|g''''|<\frac{4332.771}{n^3},\quad
|g'''''|<\frac{631300.899}{n^4}.
\]

The second-difference identity also gives `tau_g>.98/n`.  Consequently,

\[
|J(\tau_g)'''|<\frac{2{,}606{,}140}{n^3},\qquad
2|J(\tau_f)'''|<\frac{4{,}021{,}535}{n^3}.
\]

At `n=2,000,000`, including the nearest-integer adjustment, the total
negative correction has scaled coefficient below `3.31384`, while the
inherited term is above `3.69`.  The remaining coefficient margin is greater
than `0.376`.

## Why the rank-growth signal is better than the raw constants

Suppose the base asymptotic begins

\[
f'''(x)=\frac2{x^2}+\frac{d_1}{x^3}+\cdots.
\]

The exact recursion and `J(t)=log(t)+O(t)` force, for every fixed rank,

\[
F_r'''(x)=\frac{2r}{x^2}
 +\frac{r d_1-r(r-1)}{x^3}+\cdots.
\]

Indeed the leading coefficients obey a homogeneous second-difference
recursion, while the universal next coefficients obey the same recursion with
constant forcing `-2`.  Thus the first dangerous rank term is `-r(r-1)`, not
`r!` or `C^r`.  At `x=lambda*r`, its ratio to the main term tends
`1/(2*lambda)`.

That is the serious positive result of the gate.  It identifies a plausible
normalized target such as `x*tau_r/(2r)` for a uniform wedge argument.  It
does not prove that all omitted terms remain polynomial in rank.

## Numerical rank-four diagnostic

Independent saddle-centered 60-decimal quadrature is positive at every sample
from `x=100` through `x=1,000,000`.  The scaled value `x^2 F_3'''` rises from
approximately `4.23` to `5.46`, consistent with the predicted limiting value
`6`.  These quadrature signs are not interval enclosures.

## Calibrated conclusion

- Fixed-rank technical result: real progress, conditional on unreplayed source
  lemmas.
- Structural result: the first two asymptotic orders support polynomial rank
  scaling and pass the predefined kill criterion.
- Missing result: a uniform remainder theorem in the regime `x/r` bounded
  below, followed by a bridge from a tail wedge to all shifts.
- RH significance today: still low.  No global rank-four, all-rank, or RH
  theorem follows.

Primary analytic source: Wojciech Michałowski,
[An explicit uniform cubic wedge for consecutive Toeplitz minors of the
Riemann xi coefficients](https://arxiv.org/abs/2607.16795).

## Reproduce

```bash
uv run --with python-flint \
  experiments/xi-rank-four-final-kill-gate-20260828/certify_tail.py \
  --output experiments/xi-rank-four-final-kill-gate-20260828/results/tail-v1.json

uv run --with mpmath \
  experiments/xi-rank-four-final-kill-gate-20260828/probe.py \
  --output experiments/xi-rank-four-final-kill-gate-20260828/results/wide-probe-v1.json
```
