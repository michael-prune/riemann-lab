# Xi rank-three kill test

## Decision

Rank three passes the predefined kill test, but narrowly in the strategic
sense: the xi-specific saddle mechanism lifts effectively, while generic
rank-two monotonicity does not.

Conditional on the same unreplayed analytic saddle and contour lemmas used by
the rank-two certificate, the directed calculation gives

\[
(\log D_2(x))'''>0\qquad(x\ge10^6),
\]

and hence

\[
Q_{3,k+1}<Q_{3,k}\qquad(k\ge1{,}000{,}001).
\]

This is an effective fixed-rank tail statement.  It is not a global rank-three
theorem and does not control the RH-critical regime in which rank grows with
shift.

## Exact rank-three transform

Let `f=log a`,

\[
\tau(x)=2f(x)-f(x-1)-f(x+1)>0,
\qquad J(t)=\log(1-e^{-t}).
\]

Since

\[
D_2(x)=a(x)^2(1-e^{-\tau(x)}),
\]

we have

\[
g(x):=\log D_2(x)=2f(x)+J(\tau(x))
\]

and

\[
g'''=2f'''+J'''(\tau)\tau'^3
      +3J''(\tau)\tau'\tau''+J'(\tau)\tau'''.
\]

The inherited term is order `x^-2`.  The complete nonlinear correction is
order `x^-3`, because `tau` is order `x^-1` and its first three derivatives
are respectively order `x^-2`, `x^-3`, and `x^-4`.

## Effective directed bound

The rank-two threshold certificate gives

\[
x^2f'''(x)>1.23\qquad(x\ge10^5).
\]

On nested source disks, direct Cauchy estimates give the conservative real
derivative hierarchy

\[
|f'''|<\frac{37}{n^2},\qquad
|f''''|<\frac{1850}{n^3},\qquad
|f'''''|<\frac{185000}{n^4}.
\]

The ordinary Turan curvature window gives `tau_n>1/(2n)` at integers.
Moving by at most one half with the first derivative bound yields
`tau(x)>.499/n`.

For `J(t)=log(1-exp(-t))`, elementary hyperbolic inequalities give

\[
J'(t)\le t^{-1},\qquad |J''(t)|\le t^{-2},\qquad
J'''(t)\le2t^{-3}.
\]

Therefore the nonlinear correction is at most

\[
\frac{2{,}010{,}768}{n^3}.
\]

At `n=10^6`, after the nearest-integer adjustment, the coefficient budget is

\[
2x^2f'''(x)>2.46,qquad
x^2|J(\tau(x))'''|<2.011,
\]

leaving a directed margin greater than `0.449`.

## Wide xi diagnostic

Saddle-centered 60-decimal quadrature samples `x=100` through `2,000,000`.
Every value of `g'''` is positive.  The scaled quantity `x^2 g'''` rises from
approximately `2.84` at `x=100` to `3.66` at `x=2,000,000`.  The nonlinear
correction is already less than one percent of the inherited term at `x=100`
and decays by another power of `x`.

These samples are not interval enclosures.

## Exact generic obstruction

The lift is not formal.  Set

\[
q_k=\frac{a_{k-1}a_{k+1}}{a_k^2}
\]

to the increasing rational sequence

\[
\frac{11}{100},\frac7{50},\frac6{25},\frac{23}{50},
\frac{11}{20},\frac{73}{100},\frac{22}{25},
\]

and reconstruct positive rational `a_k` from `a_0=a_1=1`.  Because
`Q_(2,k)=1/q_k-1`, every rank-two exchange ratio strictly decreases.  But

\[
Q_{3,4}=\frac{30449}{10051}
<\frac{3911}{1089}=Q_{3,5}.
\]

Thus rank two does not imply rank three, even exactly and locally.  Xi passes
because its saddle supplies a controlled derivative hierarchy, not because of
total positivity or the exchange relation alone.

## Program decision

Do not abandon the lane solely because of rank three: it has a real,
quantitative lift using the same analytic infrastructure as rank two.

Do not upgrade this to a strong RH program either.  The mechanism presently
works only at fixed rank and in a tail.  The next mandatory gate is whether a
rank-four formula preserves an `x^-2` inherited term with corrections whose
constants grow mildly enough to say something uniform in rank.  Failure there,
or factorial/exponential growth of the constants, should kill the angle.

Primary analytic source: Wojciech Michałowski,
[An explicit uniform cubic wedge for consecutive Toeplitz minors of the
Riemann xi coefficients](https://arxiv.org/abs/2607.16795).

## Reproduce

```bash
uv run --with python-flint \
  experiments/xi-rank-three-kill-test-20260827/certify_tail.py \
  --output experiments/xi-rank-three-kill-test-20260827/results/tail-v1.json

uv run --with mpmath \
  experiments/xi-rank-three-kill-test-20260827/probe.py \
  --output experiments/xi-rank-three-kill-test-20260827/results/wide-probe-v1.json
```
