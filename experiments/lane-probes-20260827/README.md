# Four RH Lane Probes — 2026-08-27

This experiment applies one bounded kill test to each of four candidate
research directions.  Run it from the repository root with:

```bash
uv run \
  --with-requirements experiments/lane-probes-20260827/requirements.txt \
  python experiments/lane-probes-20260827/probe.py \
  --output experiments/lane-probes-20260827/results/probes-v1.json
```

The checked output is
[`results/probes-v1.json`](results/probes-v1.json). It is an exploratory
artifact, not an RH claim.

The Toeplitz lead also has an independent 140-decimal, 18-theta-term rerun:

```bash
uv run \
  --with-requirements experiments/lane-probes-20260827/requirements.txt \
  python experiments/lane-probes-20260827/verify_toeplitz.py \
  --output experiments/lane-probes-20260827/results/toeplitz-independent-v1.json
```

## 1. Weil finite compression: window optimization is saturated

The first probe discretizes the published functional

\[
R(\psi)=\frac{\int\psi^2+\iint |u-v|\psi(u)\psi(v)\,du\,dv}
{(\int\psi)^2}.
\]

The numerical minimizer agrees with the normalized
`cos(sqrt(2) u)` Montgomery–Taylor window.  The resulting simple-zero
certificate is about `0.67250`; shortening the bandwidth also peaks at the
published endpoint.  Together with the published bandwidth-one ceiling below
`0.6819`, this kills window-only optimization as a source of a large advance.
New leverage requires wider-support pair correlation or rigorously available
higher trace moments.

Primary source: [Alpöge–Furman, arXiv:2608.13637](https://arxiv.org/abs/2608.13637).

## 2. Xi Toeplitz frontier: a candidate backward-propagation inequality

Let

\[
Q_{r,k}=\frac{D_{r,k}D_{r-2,k}}
{D_{r-1,k-1}D_{r-1,k+1}}
=\frac{D_{r-1,k}^2}{D_{r-1,k-1}D_{r-1,k+1}}-1.
\]

For xi coefficients through `a_54`, at 90 decimal digits, all sampled
`Q_(r,k)` for ranks 2–9 and shifts through 44 are positive and strictly
decrease with `k`.  If `Q_(r,k+1) < Q_(r,k)` held uniformly, the known positive
cubic tail would propagate backward and prove all consecutive minors positive.

This monotonicity is not generic total positivity: seeded degree-70
negative-real-root polynomial controls violate it at every tested rank 2–8.
That makes the xi pattern more specific, but the finite xi range is already
covered by known finite-order positivity.  The result is a candidate inequality,
not theorem evidence.  The independent 140-decimal rerun with 18 theta terms
reproduced every sign and the displayed terminal ratios.

Primary source: [Michałowski, arXiv:2607.16795](https://arxiv.org/abs/2607.16795).

## 3. De Bruijn–Newman heat flow: pointwise gap widening fails

Using the zero dynamics

\[
\gamma_j'(t)=2\sum_{k\ne j}'\frac{1}{\gamma_j-\gamma_k},
\]

the probe evaluates adjacent-gap derivatives at `t=0` from the first 100 zeta
zeros.  Among gaps 10–30, 12 shrink and 9 widen under forward flow.  Their signs
are unchanged when the symmetric sum is truncated at 40, 70, or 100 positive
zeros.  Pairwise repulsion therefore does not supply a pointwise no-collision
barrier; a viable argument needs a collective energy and a controlled infinite
tail.

Primary source: [Polymath, arXiv:1904.12438](https://arxiv.org/abs/1904.12438).

## 4. Nyman–Beurling: diagonal preconditioning is insufficient

The probe builds the exact finite Gram entries

\[
G_{k,l}=\frac{A(l/k)}{l}
\]

from the rational Vasyunin-sum formula, then compares raw and unit-diagonal
condition numbers through `N=320`.  Diagonal scaling improves the constant but
both condition numbers retain approximately quadratic growth.  Meanwhile
`log(N)` times the best squared approximation distance is nearly flat over the
larger samples, consistent with the familiar slow `1/log(N)` scale.

The next plausible preconditioner must exploit Mellin or block-Hankel structure,
not just normalize basis norms.

Primary sources:
[Alouges–Darses–Hillion, arXiv:2006.02953](https://arxiv.org/abs/2006.02953) and
[Báez-Duarte–Balazard–Landreau–Saias, arXiv:math/0306251](https://arxiv.org/abs/math/0306251).

## Claim boundary

- The xi and heat-flow calculations use arbitrary-precision numerics without
  interval enclosures.
- The Nyman–Beurling entries use exact finite formulas evaluated in double
  precision.
- The Weil probe reproduces published mathematics.
- None of the four probes proves RH, external novelty, or a new general theorem.
