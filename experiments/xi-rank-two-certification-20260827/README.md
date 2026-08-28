# Rank-two xi certification attempt

## Follow-up

The effective conditional tail threshold is lowered from `10^9` to `10^5` in
[`../xi-rank-two-threshold-20260827/README.md`](../xi-rank-two-threshold-20260827/README.md)
by differentiating the curvature logarithm before applying Cauchy's estimate
and rebuilding the scalar saddle factorization with a smaller local window.

## Outcome

The full-theta tail transfer closes with enormous room **conditional on the
analytic bounds stated in Michałowski's arXiv paper**.  It does not yet produce
a global theorem.

The new exact reduction is that the real full-kernel saddle action obeys

\[
S'''(x)=-2\frac{K_\Phi'(u)+uK_\Phi''(u)}
                 {u^2K_\Phi'(u)^3},\qquad K_\Phi(u)=x.
\]

Writing `Phi=Phi_1(1+R)` and `g=log(1+R)`, one has

\[
K_\Phi=K-\frac u2g',\quad
\delta K'=-\frac12g'-\frac u2g'',\quad
\delta(K'+uK'')=-\frac12g'-\frac{3u}{2}g''-\frac{u^2}{2}g'''.
\]

Thus the transfer needs only three Cauchy bounds for `g`; no numerical
differentiation of a moving saddle is required.

Using the source-stated `sup |R| < 1.1e-14`, a Cauchy radius `rho=0.04`,
and the dominant-saddle inequalities `K'>4K` and `K'+uK''>4K`, the
relative correction to the dimensionless saddle coefficient is less than
`1e-9` (the actual certified enclosure is much smaller).  The previous
dominant coefficient bound therefore remains below `0.477`.

Combining this with the source-stated zero-free factorization gives,
conditionally,

\[
f'''(x)>0\qquad(x\ge 10^9).
\]

The third forward difference is a positive B-spline average of `f'''`, so

\[
Q_{2,k+1}<Q_{2,k}\qquad(k\ge 1{,}000{,}000{,}001)
\]

under those same hypotheses.

This is a checked scalar closure of a source-backed conditional argument, not
an internal theorem.  The interval `0 <= x < 10^9` remains open.

## Derivation of the transfer factor

Put

\[
H(u)=2\pi e^{4u}-\frac92-\frac6{2\pi e^{4u}-3},\qquad K=uH.
\]

On `u>=0.98`, `H'/H>4`, `K'>4K`, and `K''>0`.  With
`M=-log(1-1.1e-14)` and `rho=0.04`, Cauchy's estimate gives
`|g^(j)| <= j! M/rho^j`.  Consequently

\[
\begin{aligned}
\alpha&=\frac{|\delta K|}{K}
 \le \frac{M}{2\rho H},\\
\beta&=\frac{|\delta K'|}{K'}
 \le \frac{M}{8\rho uH}+\frac{M}{4\rho^2H},\\
\gamma&=\frac{|\delta(K'+uK'')|}{K'+uK''}
 \le \frac{M}{8\rho uH}+\frac{3M}{4\rho^2H}
      +\frac{3uM}{4\rho^3H}.
\end{aligned}
\]

All three right sides decrease for `u>=0.98`: each is a positive sum of
`1/H`, `1/(uH)`, and `u/H`, and
`H'/H>4>1/u`.  Therefore their endpoint values imply

\[
-x^2S'''(x)
\le U(0.98)\frac{(1+\alpha)^2(1+\gamma)}{(1-\beta)^3}<0.477,
\]

where

\[
U(u)=\frac{80000(3202u^2+2402u+197)}{(800u+197)^3}
\]

is the earlier dominant-saddle majorant.  Its derivative has numerator
`-160000(1280800u^2+1290806u-197)`, hence it is decreasing on this range.

For arbitrary real `x>=1e9`, center the source disk at the nearest integer
`n`.  A Cauchy circle of radius `n/25` about `x` remains inside
`|z-n|<=n/20`.  After multiplication by `x^2`, the curvature-logarithm and
relative-error contributions are bounded by

\[
\left(1+\frac1{2n}\right)^2
\frac{46875(\log(20n)+\pi/2)}n,
\qquad
\left(1+\frac1{2n}\right)^2\frac{1781.25}{n}.
\]

The positive Gamma contribution is larger than
`8x^2/(2x+1)^2`.  The Arb certificate leaves a coefficient margin greater
than `1.52` at the worst endpoint.

## Reproducibility boundary

The paper says its ancillary directory contains four certificate modules and
a 36-test suite.  On 2026-08-27, the arXiv v1 source archive downloaded from
`https://export.arxiv.org/e-print/2607.16795` contained only `main.tex` and
arXiv's `00README.json`; direct ancillary URLs returned HTTP 404.  Therefore
the source's tube, saddle, and factorization bounds could not be independently
replayed here.

This matters for claim level.  The algebra above and its scalar inequalities
are checked, but the argument remains conditional until those analytic
certificates are obtained or rebuilt.

Primary sources:

- Wojciech Michałowski, [An explicit uniform cubic wedge for consecutive
  Toeplitz minors of the Riemann xi coefficients](https://arxiv.org/abs/2607.16795).
- Jeremy J. F. Guo, [An inequality for coefficients of the real-rooted
  polynomials](https://arxiv.org/abs/2012.03530).  Its higher-order Turán
  inequalities are adjacent but do not imply the ratio monotonicity here.

## What would actually finish rank two

Two obligations remain:

1. obtain or independently rebuild the missing source certificates;
2. replace the enormous compact gap `0 <= x < 10^9` by a much lower effective
   saddle threshold plus interval quadrature on the remaining compact range.

The second item is now the substantive mathematical bottleneck.  Rechecking
the tiny theta perturbation is not.

## Reproduce

```bash
uv run --with python-flint \
  experiments/xi-rank-two-certification-20260827/certify_conditional_tail.py \
  --output experiments/xi-rank-two-certification-20260827/results/conditional-tail-v1.json
```
