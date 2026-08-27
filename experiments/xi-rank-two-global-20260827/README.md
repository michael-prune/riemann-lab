# Global rank-two xi exchange attempt

## Follow-up

The theta-transfer obligation identified below was reduced and its scalar
closure certified in
[`../xi-rank-two-certification-20260827/README.md`](../xi-rank-two-certification-20260827/README.md).
That follow-up proves a tail statement only conditional on the analytic bounds
stated in the cited source: the arXiv v1 source archive did not contain the
advertised certificate modules, so those bounds could not be replayed.  The
compact range below `10^9` remains open.

## Outcome

The xi inequality remains alive, but the obvious general theorem behind it does
not exist.

For

\[
a(x)=\frac{1}{\Gamma(2x+1)}\int_0^\infty u^{2x}\Phi(u)\,du,
\qquad f(x)=\log a(x),
\]

the continuous sufficient condition is

\[
f'''(x)=8\left(\kappa_{3,x}(\log u)-\psi_2(2x+1)\right)>0.
\]

The right side compares the third log-cumulant of the tilted xi-kernel law
with that of a Gamma law.  A dense xi scan found no crossing.  However, a
simple decreasing log-concave density violates the comparison, so Borell,
Berwald, or ordinary kernel log-concavity cannot prove the xi statement by
themselves.

The dominant theta saddle has a large explicit sign margin.  This reduces a
possible effective tail proof to one narrow missing estimate: transfer the
dominant-summand saddle bound to the full theta saddle with directed bounds on
the derivatives of the tiny theta remainder.

No global rank-two theorem, rank-three lift, or RH result is proved here.

## 1. A sharp no-go for the general log-concavity route

Consider the density on `[0,infinity)`

\[
\varphi(u)=e^{-V(u)},\qquad
V(u)=\begin{cases}2u,&0\le u\le1,\\3u-1,&u\ge1.\end{cases}
\]

The potential is convex, so `varphi` is log-concave; it is also positive,
continuous, and strictly decreasing.  At `x=0`, 80-decimal quadrature gives

\[
\kappa_3(\log U)=-2.4273387546492343288\ldots,
\qquad
\psi_2(1)=-2.4041138063191885708\ldots,
\]

and hence

\[
8(\kappa_3-\psi_2)=-0.1857995866403660641\ldots<0.
\]

Thus the desired third-derivative sign is false even for decreasing
log-concave kernels.  The kink can be smoothed while preserving convexity and
the strict negative margin, so smoothness alone does not repair the route.

This explains the limit of the Berwald--Borell theorem.  Dividing the Mellin
transform of a log-concave density by the Gamma function proves the normalized
Mellin transform is log-concave, i.e. `f''<=0`.  It does not force the next sign
`f'''>=0`.  The recent log-concave Bernstein theorem supplies stronger
representations for individual Turan differences, but not this ordering of
successive normalized curvatures.

Primary source: Klartag and Lehec,
[Poisson processes and a log-concave Bernstein theorem](https://www.math.tau.ac.il/~klartagb/papers/log_concave_bernstein.pdf).

## 2. Xi-specific numerical test

The probe uses 5,000-point Gauss--Legendre quadrature on `[0,5]`, twelve theta
terms, and 2,501 parameter values: a uniform mesh on `[0,100]` followed by a
geometric mesh through `10^4`.

Every sampled value of

\[
\kappa_{3,x}(\log u)-\psi_2(2x+1)
\]

was positive.  The largest ratio of the two negative quantities,
`kappa_3/psi_2`, occurs at `x=0` and is approximately `0.87365`.  It then
decreases steadily, reaching approximately `0.14167` at `x=10^4`.  Therefore
the xi distribution moves farther from a crossing in relative terms.

This is finite double-precision evidence, not interval certification.

## 3. Dominant-saddle sign budget

For the first theta summand, write

\[
K(u)=u\left(B-\frac92-\frac6{B-3}\right),\qquad
B=2\pi e^{4u},
\]

so the real saddle is `k=K(u)`.  If `S_0` is the saddle action, exact implicit
differentiation gives

\[
S_0'''(k)=-2\left(\frac{K''}{uK'^3}+\frac1{u^2K'^2}\right)<0.
\]

The dimensionless amount consumed from the positive Gamma contribution is

\[
C(u)=-k^2S_0'''(k).
\]

It decreases numerically from about `0.7493` at `k=16` to `0.1142` at
`k=10^9`.  More usefully, elementary bounds on `K`, `K'`, and `K''` give, for
`u>=0.98`,

\[
C(u)\le
U(u)=\frac{80000(3202u^2+2402u+197)}{(800u+197)^3}.
\]

The derivative numerator is

\[
-160000(1280800u^2+1290806u-197),
\]

so `U` is decreasing there and `U(0.98)<0.477`.

By contrast, the Gamma term obeys

\[
-8\psi_2(2k+1)>\frac{8}{(2k+1)^2},
\]

whose coefficient after multiplication by `k^2` is essentially `2`.
At `k=10^9`, the dominant saddle uses less than `0.477` of that budget.

Michałowski's certified zero-free factorization on the five-percent complex
disk bounds the curvature-logarithm and relative-error Cauchy contributions by
approximately `0.0012/k^2` and `0.000002/k^2`, respectively.  More than
`1.52/k^2` remains for the passage from the first theta summand to the full
theta saddle.

The paper already gives a theta remainder below `1.1e-14` on the relevant
tube.  The unfinished proof obligation is to propagate that remainder through
the first two derivatives of the full saddle map with directed rounding.  The
margin suggests this is realistic, but it has not been certified here.

Primary source: Michałowski,
[An explicit uniform cubic wedge for consecutive Toeplitz minors of the Riemann xi coefficients](https://arxiv.org/abs/2607.16795).

## 4. What this changes

The global rank-two problem has become a compact-plus-tail certification
problem:

1. finish the full-theta saddle perturbation to obtain an effective large-`x`
   theorem;
2. certify the remaining compact interval for the exact cumulant expression;
3. only then revisit the rank-three transform.

The earlier idea of proving a universal higher-order Berwald theorem should be
discarded.  The counterexample shows that xi's special theta structure is
essential.

## Reproduce

```bash
uv run --with mpmath,numpy,scipy \
  experiments/xi-rank-two-global-20260827/probe.py \
  --output experiments/xi-rank-two-global-20260827/results/probe-v1.json
```
