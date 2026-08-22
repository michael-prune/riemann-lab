# Rank-Two Transport Geometry

## Outcome

The most literal rank-two chamber involution does not work.  After fixing the
product of the two moment variables, some fibers have negative total mass.
Therefore no measure-preserving pairing confined to one such fiber can prove
the Turan determinant positive.

There is, however, a genuine positive geometry behind the classical rank-two
proof: a beta simplex followed by continuous Cauchy-Binet composition over an
ordered chamber.  This is a continuous Lindstrom-Gessel-Viennot mechanism.
The open extension target is sign-regularity of a theta-specific cumulative
kernel, together with a non-circular transfer back to every xi Toeplitz minor.

Nothing in this packet proves a new Turan inequality or RH.

## 1. Exact signed chamber

For `k >= 1`, write

```text
b_k = integral_0^infinity u^(2k) Phi(u) du,
a_k = b_k/(2k)!,
D_(2,k) = a_k^2 - a_(k-1)a_(k+1),
c_k = (2k)(2k-1)/((2k+1)(2k+2)).
```

Symmetrizing the product integrals gives the exact identity

```text
D_(2,k) = 1/(2(2k)!^2) double-integral
  Phi(u)Phi(v)u^(2k-2)v^(2k-2)
  [2u^2v^2 - c_k(u^4+v^4)] du dv.
```

Set `x = uv`, `y = log(u/v)`.  Since `du dv = (1/2) dx dy`, this becomes

```text
D_(2,k) = 1/(2(2k)!^2) integral_0^infinity x^(2k)
  integral_R W_x(y)[1-c_k cosh(2y)] dy dx,

W_x(y) = Phi(sqrt(x)e^(y/2)) Phi(sqrt(x)e^(-y/2)).
```

The chamber near `y=0` is positive and sufficiently unequal pairs are
negative.  More strongly, the inner integral itself has both signs as `x`
varies.  The committed 60-digit probe includes, for example:

| k | x | fiber sign |
| ---: | ---: | --- |
| 1 | 0.001 | negative |
| 1 | 0.01 | positive |
| 2 | 0.01 | negative |
| 2 | 0.03 | positive |
| 3 | 0.01 | negative |
| 3 | 0.03 | positive |
| 5 | 0.03 | negative |
| 5 | 0.1 | positive |

This kills every proposed involution that preserves `x=uv`.  Any valid
transport must mix product scales or use a different carrier.

## 2. The rank-two positive carrier already present in the classical proof

Csordas, Norfolk, and Varga introduced

```text
K(t) = integral_(sqrt(t))^infinity u Phi(u) du
```

and proved `log K(t)` strictly concave.  In the square variable this is exactly

```text
s(t) = Phi(sqrt(t)),
K(t) = (1/2) integral_t^infinity s(v) dv,
K'(t) = -s(t)/2.
```

Thus `K` is the first Karlin cumulant of the square-variable theta kernel.
Define

```text
lambda_x = 1/Gamma(x+1) integral_0^infinity t^x K(t) dt.
```

Integration by parts gives

```text
lambda_x = 2/Gamma(x+2) integral_0^infinity u^(2x+3) Phi(u) du,
lambda_(m-3/2) = 2b_m/Gamma(m+1/2),
a_m = sqrt(pi) lambda_(m-3/2)/(2 * 4^m * m!).
```

The beta identity expresses the normalized power kernel through a positive
simplex.  Continuous Cauchy-Binet then integrates over ordered variables; its
two determinants antisymmetrize crossing configurations, while the ordered
residual has one sign.  In geometric language:

```text
positive beta simplex
  -> ordered two-particle chamber
  -> crossing cancellation by determinant antisymmetry
  -> positive rank-two boundary minor.
```

This reconstructs the known 1986 proof as a positive geometry.  It is not a new
proof of a new statement.  It also explains why a pointwise swap in the
original `(u,v)` plane was the wrong calibration target: the successful
crossing cancellation occurs only after composition through an auxiliary
simplex.

## 3. Higher-rank experimental target

For ordered `x_1 < ... < x_r` and `y_1 < ... < y_r`, test the additive kernel

```text
mathcal K_r(x,y) = [K(x_i+y_j)]_(i,j=1,...,r)
```

against the reverse sign pattern

```text
epsilon_r det(mathcal K_r) > 0,
epsilon_r = (-1)^(r(r-1)/2).
```

On the six-point grid

```text
0.001, 0.003, 0.01, 0.03, 0.1, 0.3
```

all 923 minors through rank six had the predicted sign in the committed
high-precision calculation.  This consists of 922 minors through rank five and
the single full rank-six determinant.  It is a finite, non-certified probe—not
a theorem and not evidence for all inputs.  An independent 90-digit pass
matched the first 30 displayed digits of every fiber value and every
rank-minimum oriented minor from the committed 60-digit pass.

The literature makes the boundary unusually informative:

- The original additive kernel `Phi(x+y)` is proved reverse sign-regular
  through rank three, has strong numerical support at rank four, and is known
  to fail at rank five.
- Karlin cumulants cannot reduce an already established order of
  sign-regularity.  A recent Nuttall argument proves unbounded cumulative order
  for one Conrey-Ghosh analogue, while treating extension to the Riemann kernel
  as a further problem.
- The August 2026 second-level-concavity result for `s(t)=Phi(sqrt(t))` proves
  positivity inside a sheared quadratic cone.  That is compatible with a
  rank-two/rank-three cumulative geometry, but does not prove the all-rank
  statement here.

## 4. Rank-five falsification

A coalescing-point follow-up overturns the all-rank square-tail conjecture.
For additive kernels, the decisive local objects are

```text
W_p(t) = det[K^(i+j-2)(t)]_(i,j=1,...,p).
```

The required signs for ranks four, five, and six are respectively `+`, `+`,
and `-`.  At `t=0.001`, the derivative calculation gives the correct rank-four
and rank-six signs but the wrong rank-five sign:

```text
W_4(0.001) =  1.8914868092e5
W_5(0.001) = -1.1851429664e13   (required positive)
W_6(0.001) = -4.4601867506e24
```

The failure was reproduced without numerical differentiation.  Set

```text
x_i = i h,
y_j = 0.001 + j h,
M_(i,j) = K(x_i+y_j),
```

for `i,j = 0,...,4`.  For each

```text
h = 0.0001, 0.00005, 0.00002, 0.00001,
```

the actual rank-five determinant `det M` is negative, although rank five
requires a positive determinant.  After division by the positive coalescence
factor `h^20`, the values converge toward the negative derivative Wronskian.
For example, at `h=0.0001`,

```text
det M = -1.0544596288e-67,
det M / h^20 = -1.0544596288e13.
```

The calculation used 220 decimal digits and 16 theta summands.

Therefore the square-tail kernel fails the numerical `RR5` test.  It also
cannot be `RR6`, because `RR6` requires the correct compound-kernel signs at
every rank from one through six.  The isolated rank-six determinant having the
correct sign does not repair the rank-five failure.

This is a high-confidence numerical counterexample, not yet an interval proof.
The full record is in
[`results/square-tail-wronskian-v1.json`](results/square-tail-wronskian-v1.json).

The finite rank-six grid from the preceding section was therefore misleading:
widely separated points can have the desired signs even when coalescing points
expose a wrong Wronskian sign.  The square-tail kernel cannot be the desired
all-rank positive carrier.

## 5. What remains at rank four

For an additive kernel, sign-regularity through rank `r` can be reduced to the
signs of the two-way Wronskians

```text
W_p(t) = det[K^(i+j-2)(t)]_(i,j=1,...,p),  p=1,...,r.
```

The rank-four Wronskian had the correct sign at every sampled point from
`t=0.000001` through `t=0.5`, using

```text
K'(t) = -Phi(sqrt(t))/2
```

but this does not prove its sign on all `t >= 0`.  The next worthwhile question
is now narrower: certify or falsify `RR4`, and then determine whether a different
number of cumulant integrations can raise the available order as in Nuttall's
program.  The one-fold square-tail kernel should no longer be treated as an
all-rank candidate.

## Reproduce

```bash
.venv/bin/python experiments/positive-grassmannian/rank_two_transport.py \
  --output experiments/positive-grassmannian/results/rank-two-transport-v1.json

.venv/bin/python experiments/positive-grassmannian/wronskian_counterexample.py \
  --output experiments/positive-grassmannian/results/square-tail-wronskian-v1.json
```

## Primary sources

- G. Csordas, T. Norfolk, and R. Varga,
  [The Riemann Hypothesis and the Turan inequalities](https://www.math.kent.edu/~varga/pub/paper_157.pdf).
- J. Nuttall,
  [Wronskians, Cumulants, and the Riemann Hypothesis](https://publish.uwo.ca/~jnuttall/det3CAEXP_2.pdf).
- J. Nuttall,
  [Cumulants, the Riemann Hypothesis, and Similar Problems](https://publish.uwo.ca/~jnuttall/cumanal_a.pdf).
- P. Planat and P. Sole,
  [Second-Level Concavity of the Riemann Xi Kernel](https://arxiv.org/abs/2608.19160).
