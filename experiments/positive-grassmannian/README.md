# Xi Positive Geometry / Grassmannian Attempt

## Research question

Can the completed Riemann xi function be represented by a positive geometry or
positive network whose boundary measurements make the Riemann Hypothesis (RH)
positivity conditions manifest?

The exact starting point is already known: RH is equivalent to total
nonnegativity of a Toeplitz matrix built from positive Taylor coefficients of
the completed xi function. The research target here is the missing geometric
reason for those minors to be positive.

## The object map

| Object | Singularities and zeros | Role in this attempt |
| --- | --- | --- |
| `zeta(s)` | Simple pole at `s = 1`; trivial zeros at negative even integers; nontrivial zeros in the critical strip | Arithmetic source, but not the clean geometry object |
| `Lambda(s) = pi^(-s/2) Gamma(s/2) zeta(s)` | Meromorphic with poles at `0` and `1`; satisfies `Lambda(s) = Lambda(1-s)` | Makes the functional symmetry explicit |
| `xi(s) = s(s-1)Lambda(s)/2` | Entire; its zeros are exactly the nontrivial zeta zeros; `xi(s) = xi(1-s)` | Removes the pole pair and the trivial-zero/gamma-pole bookkeeping |
| `Xi(x) = xi(1/2 + ix/2)/8` | Even, real entire; RH says all its zeros are real | Fourier transform of the positive theta kernel `Phi` |
| `G(z) = xi(1/2 + sqrt(z)/2)/8 = sum a_k z^k` | Entire of order `1/2`; `a_k > 0`; a xi zero `rho` maps to `z = 4(rho-1/2)^2` | RH says every zero of `G` is on the negative real axis |
| `d log G = G'(z)/G(z) dz` | Simple poles at mapped xi zeros, with positive integer residues equal to multiplicities | Canonical-form-like meromorphic object; under RH all poles lie on one negative ray |
| `T(a) = [a_(j-i)]`, with `a_n = 0` for `n < 0` | Its minors are finite algebraic functions of the xi coefficients | RH is equivalent to total nonnegativity of `T(a)` |

This identifies the correct use of poles. The zeta pole at `1` is removed before
the positivity problem is posed. The nontrivial zeros then reappear as poles of
the logarithmic differential `d log G`, where a canonical-form comparison is
mathematically natural.

## Exact positive-Grassmannian bridge

For

```text
Phi(u) = sum_{n>=1}
  (2 pi^2 n^4 e^(9u) - 3 pi n^2 e^(5u)) e^(-pi n^2 e^(4u)),
```

the coefficients are

```text
a_k = 1/(2k)! integral_0^infinity u^(2k) Phi(u) du > 0.
```

The Aissen-Schoenberg-Whitney-Edrei characterization gives

```text
RH
<=> (a_k) is PF_infinity
<=> every minor of T(a) is nonnegative.
```

For fixed rank `r`, width `n`, and offset `q`, the matrix

```text
B(q)_(i,j) = a_(q+j-i),  0 <= i < r, 0 <= j < n
```

represents a point of `Gr(r,n)` whenever it has full rank. Its maximal minors
are its Plucker coordinates. Thus the infinite RH criterion is a compatible
family of finite nonnegative-Grassmannian conditions.

The consecutive chamber minors

```text
D_(r,k) = det[a_(k+j-i)]_(i,j=0,...,r-1)
```

satisfy the Desnanot-Jacobi / Plucker exchange relation

```text
D_(r,k) D_(r-2,k)
= D_(r-1,k)^2 - D_(r-1,k-1) D_(r-1,k+1).
```

This supplies local positive coordinates whenever all terms are positive:

```text
X_(r,k) = D_(r,k)D_(r-2,k)
          / (D_(r-1,k-1)D_(r-1,k+1)).
```

## Hydrotope-inspired carrier

Each coefficient has an exact simplex-volume representation:

```text
u^(2k)/(2k)! = Vol{0 <= t_1 <= ... <= t_(2k) <= u}.
```

Therefore `a_k` is a positive theta-weighted volume of ordered simplices. In a
permutation term of `D_(r,k)`, the simplex dimensions are

```text
2(k + sigma(i) - i),  i = 0,...,r-1.
```

Their total is independent of the permutation:

```text
sum_i 2(k + sigma(i) - i) = 2rk.
```

This is the proposed **theta-Hydrotope carrier**: permutation-labelled products
of ordered simplices all live over the same `2rk`-dimensional fiber count. The
determinant supplies orientations. The hoped-for geometry is a gluing,
sign-reversing involution, or nonintersecting-path model that cancels opposite
chambers and leaves a positive residual region.

The carrier observation is elementary and its novelty is unknown. It does not
establish the required gluing or positivity.

## Non-circular breakthrough target

Construct a planar directed network or a sliced positive domain with all local
weights explicitly positive functions of the theta data, such that its boundary
measurement matrix is `T(a)`. By the Lindstrom-Gessel-Viennot mechanism, its
minors would become sums of positive nonintersecting-path weights. If the same
construction works uniformly for every rank and shift, it would prove the
Toeplitz total-positivity criterion and hence RH.

The [rank-two transport calibration](rank-two-transport.md) shows that a pairing
confined to fixed `uv` fibers cannot work: some complete fibers have negative
mass. It also reconstructs the classical Turan proof as a beta-simplex followed
by continuous Cauchy-Binet cancellation over an ordered chamber. The higher-rank
target is now the square-tail cumulative kernel and an exact, non-circular
composition back to the xi Toeplitz minors.

## Baseline experiment

`scan.py` computes the theta moments, consecutive Toeplitz minors, finite
Plucker slices, Plucker-exchange coordinates, cancellation margins, and two
controls.

The committed baseline used 100 decimal digits, 14 theta summands, coefficients
`a_0` through `a_20`, ranks `1` through `6`, and shifts `0` through `12`.

- All 78 consecutive minors were numerically positive.
- All 756 maximal minors across six `4 x 9` Toeplitz slices were numerically
  positive.
- All 55 exchange coordinates were numerically positive.
- The weakest determinant-to-absolute-permutation-sum margin was about
  `5.78e-14`, at rank `6`, shift `12`; high precision is essential.
- The negative-real-root control had no negative consecutive minors in the
  scanned range.
- The positive-coefficient non-PF control was rejected immediately, with
  `D_(2,1) = -9`. This confirms that the scanner does not confuse positive
  coefficients with total positivity.
- An independent 80-digit, 12-summand pass produced the same signs and the same
  first 30 rendered digits for every scanned consecutive minor.

These are reproducibility checks, not new evidence for RH. Low finite orders
are already covered much more strongly by known zero verification and sector
theorems. The calculation is not interval-certified.

## Run

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --index-url https://pypi.org/simple \
  -r experiments/positive-grassmannian/requirements.txt
.venv/bin/python experiments/positive-grassmannian/scan.py \
  --output experiments/positive-grassmannian/results/baseline-v1.json
```

## Kill conditions

Stop or reclassify a proposed construction if any of these occurs:

1. Network weights are defined using already-computed minor signs or xi zeros.
   That is a circular parametrization, not a proof.
2. The argument proves only `a_k > 0`, Hankel moment positivity, or finitely many
   Toeplitz minors. Those facts do not reach RH.
3. The construction works for every positive kernel. The non-PF control shows
   that a universal positive-kernel argument cannot be valid.
4. A finite network has no uniform limit controlling every rank and shift.
5. The proposed logarithmic form has poles away from mapped xi zeros, wrong
   residues, or uncancelled branch dependence from `sqrt(z)`.
6. The theta-Hydrotope cells cannot be embedded in one common carrier with a
   measure-preserving cancellation.
7. Any claimed sign depends on ordinary floating point when cancellation is
   comparable to precision.

## Next proof attempts

1. Interval-test the two-way Wronskians of the square-tail kernel through rank
   six; one wrong sign kills the corresponding sign-regularity conjecture.
2. Express rank-three permutation chambers over the common `6k`-dimensional
   simplex fiber and test whether an LGV-style crossing swap preserves the
   theta weight.
3. Recover finite planar-network face weights from the baseline slices, then
   test whether they admit formulas in local moment ratios rather than minors.
   Formula discovery is useful; using recovered positivity as proof is not.
4. Replace mpmath with Arb ball arithmetic before treating any unexplained
   finite sign pattern as a mathematical finding.
5. Compare the critical regime `k ~ r` with the known cubic tail wedge. A useful
   geometry must control the complementary regime, not merely reprove the tail.

## Primary sources

- O. Katkova, [Multiple positivity and the Riemann zeta-function](https://arxiv.org/abs/math/0505174).
- W. Michalowski, [An explicit uniform cubic wedge for consecutive Toeplitz minors of the Riemann xi coefficients](https://arxiv.org/abs/2607.16795).
- M. Griffin et al., [Jensen Polynomials for the Riemann Xi Function](https://arxiv.org/abs/1910.01227).
- A. Postnikov, [Total positivity, Grassmannians, and networks](https://arxiv.org/abs/math/0609764).
- N. Arkani-Hamed et al., [Surface Water Wave Scattering and the Hydrotope](https://arxiv.org/abs/2606.28280).
- D. H. J. Polymath, [Effective approximation of heat flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438).
- G. Csordas, T. Norfolk, and R. Varga, [The Riemann Hypothesis and the Turan inequalities](https://www.math.kent.edu/~varga/pub/paper_157.pdf).
- J. Nuttall, [Wronskians, Cumulants, and the Riemann Hypothesis](https://publish.uwo.ca/~jnuttall/det3CAEXP_2.pdf).
- J. Nuttall, [Cumulants, the Riemann Hypothesis, and Similar Problems](https://publish.uwo.ca/~jnuttall/cumanal_a.pdf).
- P. Planat and P. Sole, [Second-Level Concavity of the Riemann Xi Kernel](https://arxiv.org/abs/2608.19160).
