# Xi uniform-rank Toda limit and stopping boundary

## Decision

This push found a genuine rank-normalized structure, but it also located the
point where the exchange-monotonicity lane stops being a tractable tail
argument and becomes the missing all-order content of RH.

The positive result is an exact formal hydrodynamic limit.  It explains why
rank four survived and why the fixed-rank constants did not grow
factorially.  The negative result is that neither positive-moment structure,
strict kernel log-concavity, nor the known second-level kernel concavity can
control the finite-rank bridge.  Controls separate each of those properties
from the desired exchange inequality or from PF-infinity.

**Program verdict:** archive this as a strong structural candidate and stop
pushing the present angle.  Reactivate it only if a new xi-specific all-order
comparison principle appears; more fixed-rank certificates, larger numerical
grids, or additional finite concavity levels will not cross the identified
boundary.

## Exact Toda curvature recursion

Let `F_r(x)=log D_r(x)`, `F_0=0`, `F_1=f`, and

\[
\tau_r(x)=2F_r(x)-F_r(x-1)-F_r(x+1),\qquad
J(t)=\log(1-e^{-t}).
\]

Desnanot--Jacobi gives

\[
F_{r+1}=2F_r-F_{r-1}+J(\tau_r)
\]

and therefore the exact discrete Toda recursion

\[
\tau_{r+1}=2\tau_r-\tau_{r-1}+LJ(\tau_r),
\qquad Lh(x)=2h(x)-h(x-1)-h(x+1).
\]

Exchange monotonicity for rank `r+1` is equivalent to `tau_r` decreasing in
the shift while the relevant minors remain positive.

## Exact formal rank growth

The rational formal-series engine propagates the Gamma-2 saddle core
`f'''=2/x^2` through ten inverse powers.  Twelve ranks determine the
coefficient polynomials and twelve further ranks are exact held-out checks.

Writing

\[
-\tau_r'(x)=\sum_{m\ge2}\frac{c_m(r)}{x^m},
\]

the first terms are

\[
c_2=2r,
\quad c_3=-r(r-1),
\quad c_4=\frac{r(5r^2-3r+2)}4,
\]

and every checked `c_m` has degree exactly `m-1` in rank.  This degree law
also follows formally by induction: `L` raises inverse order by two, the
regular part of `J(tau)` preserves the rank-degree filtration, and solving
the second difference in rank adds two degrees.  Thus the expansion organizes
in powers of `r/x`, rather than producing a factorial rank ladder.

This is exact formal algebra, not a bound on the infinite remainder.

## Hydrodynamic continuum theorem

Take `r,x -> infinity` with `lambda=x/r` fixed and suppose the base curvature
has amplitude `A`, `tau_1(x)~A/x`.  The Toda recursion has the continuum
equation

\[
T_A'(\lambda)=
-\frac{A}{\lambda^2+1/(e^{T_A(\lambda)}-1)},
\qquad T_A(\lambda)\sim\frac A\lambda.
\]

This equation is explicitly solvable parametrically.  Put

\[
y=1-e^{-T},\qquad
u_A(T)=y\,{}_2F_1\!\left(1-\frac1A,1+\frac1A;2;y\right),
\]

then

\[
\lambda(T)=A\frac{u_A'(T)}{u_A(T)}.
\]

Indeed the inverse Riccati equation linearizes to

\[
u''+\frac{u}{A^2(e^T-1)}=0,
\]

which becomes the displayed hypergeometric equation after `y=1-e^-T`.
For `A>=1`, both hypergeometric parameters are nonnegative, so `u_A>0` and
`u_A'>0`.  The Riccati equation then makes `lambda(T)` strictly decrease
from infinity to zero.  Consequently `T_A(lambda)>0` and
`T_A'(lambda)<0` for every positive `lambda`.

The solid-minor exchange range needs only `lambda>=1`.  The boundary case is

\[
T_1(\lambda)=\log(1+1/\lambda).
\]

Comparison with this solution gives the useful continuum margin

\[
-T_A'(\lambda)\ge\frac{A}{\lambda^2+\lambda}
\ge\frac{A}{2\lambda^2}
\qquad(A\ge1,\ \lambda\ge1).
\]

This continuum theorem is an exact internal result.  It is conditional as a
description of xi because convergence of the finite Toda lattice has not been
proved uniformly in rank.

## Numerical controls

### Gamma core

For `a_k=1/Gamma(2k+1)`, 350-decimal condensation through rank 120 is strictly
exchange-monotone on the solid range.  At rank 120, relative errors against
the continuum curve are approximately `1.8e-6` at `lambda=1`, `2.1e-7` at
`lambda=2`, and `6.5e-8` at `lambda=3`.

### Xi diagonal

Independent 120-decimal theta quadrature through rank 20 is strictly
exchange-monotone for shifts from `r-1` through `3r+1`.  It also satisfies the
stronger finite comparison

\[
\tau_r^\xi(k)>\log(1+r/k)
\]

and its one-step decrease exceeds that of the explicit Gamma-1 model on every
tested point.  At rank 20 the minimum extra decrease is greater than
`0.00225`.

These are numerical diagnostics, not interval enclosures.

## Why the bridge fails with the available structure

### Positive moments are insufficient

The exact rational control

\[
a_k=\frac{8\,4^k+9\,5^k+3\,9^k+11^k+3\,14^k}{(2k)!}
\]

has all minors involved in the test positive, but its rank-four exchange
ratio increases by about `5.8%` from shift 3 to shift 4.  Thus xi's positive
moment representation and even-factorial normalization do not imply the
candidate.

### Finite kernel concavity is insufficient

The xi kernel is now source-backed as strictly log-concave, and a separate
source proves a second-level concavity inequality.  But the classical smooth
kernel `exp(-u^4)` is strictly log-concave; after `s(t)=exp(-t^2)`, its first
Laguerre expression is exactly `2 exp(-2t^2)`, whose logarithm is also
strictly concave.  Nevertheless its Fourier transform has nonreal zeros, so
its coefficient sequence is not PF-infinity.

Therefore first- and second-level kernel concavity cannot yield the all-rank
exchange principle.  Large finite success for such kernels cannot remove that
exact external obstruction.

### The tail comparison is nonperturbative at the diagonal

Michałowski's source-backed positive wedge uses the natural determinant
parameter `r^3/k`.  Its q-Pascal/Banach-algebra perturbation ceases to be small
long before `k` is comparable to `r`.  The Toda continuum explains the
leading diagonal shape, but controls show that this shape does not stabilize
finite determinant signs by itself.

The missing statement would have to preserve an all-order xi-specific
comparison from the cubic positive tail down to `k~r`.  That is the hard
content, not a remaining scalar constant optimization.

## Claim boundary and reactivation condition

- The formal degree law and hypergeometric continuum solution are exact
  internal mathematics.
- Gamma and xi convergence, the Gamma-1 barrier, and the rank-100 quartic scan
  are non-interval numerical evidence.
- The rational positive-moment reversal is exact.
- No uniform finite-rank xi theorem, global exchange theorem, PF-infinity
  theorem, or RH implication has been proved.
- Reactivate only with a new theorem controlling the full theta-specific Toda
  evolution or an equivalent all-order positive representation.  Do not
  reactivate for rank five, a larger finite scan, or another finite kernel
  concavity level.

Primary sources:

- Wojciech Michałowski,
  [uniform cubic wedge](https://arxiv.org/abs/2607.16795).
- Cormac O'Sullivan,
  [xi coefficient asymptotics](https://arxiv.org/abs/2007.13582).
- Avi Gershon,
  [log-concavity of the xi kernel](https://doi.org/10.20944/preprints202604.0159.v2)
  (preprint).
- Michel Planat and Patrick Solé,
  [second-level concavity of the xi kernel](https://arxiv.org/abs/2608.19160).

## Reproduce

```bash
uv run --with sympy \
  experiments/xi-uniform-rank-wedge-20260828/formal_series.py \
  --output experiments/xi-uniform-rank-wedge-20260828/results/formal-toda-v1.json

uv run --with mpmath --with scipy \
  experiments/xi-uniform-rank-wedge-20260828/gamma_core_probe.py \
  --output experiments/xi-uniform-rank-wedge-20260828/results/gamma-core-v1.json

uv run --with mpmath --with scipy \
  experiments/xi-uniform-rank-wedge-20260828/diagonal_probe.py \
  --output experiments/xi-uniform-rank-wedge-20260828/results/diagonal-v1.json

uv run --with mpmath \
  experiments/xi-uniform-rank-wedge-20260828/controls.py \
  --output experiments/xi-uniform-rank-wedge-20260828/results/controls-v1.json
```
