# Lowering the conditional rank-two xi threshold

## Outcome

The effective tail threshold drops from `10^9` to `10^5`, conditional on the
same unreplayed analytic saddle lemmas used in the preceding attempt.

The directed 256-bit Arb certificate proves the scalar implication

\[
f'''(x)>0\quad(x\ge 100{,}000)
\]

and therefore

\[
Q_{2,k+1}<Q_{2,k}\quad(k\ge100{,}001),
\]

provided the source-stated kernel-tube, saddle-continuation, phase, and contour
lemmas hold with the bounds used below.  This is a factor `10,000` reduction
in the uncertified compact range.  It is not yet an unconditional global
rank-two theorem.

## The two changes that remove the billion-scale barrier

### 1. Differentiate curvature before using Cauchy

The earlier attempt bounded the third derivative of the curvature logarithm
by applying Cauchy's estimate directly to its absolute size.  That retained a
spurious `log k` and cost about `65/k^2` at `k=10^4`.

At the full saddle, let `K_Phi(u)=z`.  Stationarity gives the exact identity

\[
-\Psi_z''(u)=\frac{2K_\Phi'(u)}u.
\]

Hence, for `h(z)=log(-Psi_z''(u(z)))`,

\[
h'(z)=\frac{K_\Phi''}{K_\Phi'^2}-\frac1{uK_\Phi'}.
\]

The source phase bounds imply throughout the complex saddle disk

\[
|K_\Phi'|\ge2.955n,qquad |K_\Phi''|\le70.75n,
\]

so `|h'|<8.45/n`.  Applying Cauchy only to the remaining two derivatives
costs less than `0.053/k^2` at `10^5`.  A numerical diagnostic of the exact
dominant saddle gives only `8.7e-6/k^2` there, so the certified bound is still
very conservative.

### 2. Use a one-sigma local saddle window

The source factorization used `|y|<=2`, making its cubic Taylor error scale as
`8*sqrt(log(k)/k)`.  Taking `|y|<=1` reduces this by a factor of eight while
leaving manageable Gaussian tails.

At `k=10^5`, the directed components, normalized by the Gaussian main term,
give a relative factorization error below `0.223`:

- cubic replacement cost: approximately `0.09572` before normalization;
- omitted Gaussian tail: approximately `0.01660`;
- true local tail: approximately `0.03284`;
- global horizontal tail: approximately `0.00314`;
- radial connector: negligible.

The resulting Cauchy cost from `log(1+epsilon)` is below `0.237/k^2`.

## Final sign budget

After multiplying by `x^2`, the worst-endpoint bounds are

\[
\begin{array}{c|c}
\text{component} & \text{coefficient}\\ \hline
\text{positive Gamma term} & >1.99998\\
\text{full saddle action} & <0.477\\
\text{factorization error} & <0.237\\
\text{curvature logarithm} & <0.053
\end{array}
\]

leaving a certified coefficient margin greater than `1.23`.

All endpoint expressions worsen toward `10^5`: `log(k)/k`, the horizontal
tail, the connector, and the two Cauchy terms decrease beyond that point,
while the Gamma lower bound increases.  The full-saddle action bound is
uniform.

## Claim boundary

This result rebuilds and certifies the new scalar factorization and sign
budget.  It does not replay the source's underlying analytic continuation,
phase, and contour certificates, because the arXiv v1 source archive still
lacks the advertised ancillary modules.  Those lemmas are explicit
hypotheses here.

The remaining compact range is now

\[
0\le x<100{,}000.
\]

That is small enough for a serious saddle-centered interval quadrature or a
dyadic interval argument; it is no longer a billion-scale gap.

Primary source: Wojciech Michałowski,
[An explicit uniform cubic wedge for consecutive Toeplitz minors of the
Riemann xi coefficients](https://arxiv.org/abs/2607.16795).

## Reproduce

```bash
uv run --with python-flint \
  experiments/xi-rank-two-threshold-20260827/certify_threshold.py \
  --output experiments/xi-rank-two-threshold-20260827/results/threshold-v1.json
```
