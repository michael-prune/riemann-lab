# Positive geometry for xi: modular completion before a stable lift

Date: 2026-09-04. First Astra investigation, with a bounded independent Sol
mathematical review and a Luna audit of the previous work.

**Outcome:** a complete analytic obstruction to every finite, evenly folded
theta cutoff, a correction to the previous quartic counterexample, and a
precise constructive positive-geometry target. No RH proof or external
novelty claim. The propositions below have proofs and internal review; they
have not been admitted to the lab's protocol-defined Inner Core. The numerical
files are reproducible checks, not interval certificates of the propositions.

## What changed from the previous work

The earlier work represented xi-coefficient Toeplitz minors as Grassmannian
coordinates. Its square-tail carrier failed at rank five; its exact Toda
continuum lacked a uniform transfer to finite-rank xi. Those limitations
remain. See [the earlier record](../xi-uniform-rank-wedge-20260828/README.md).

One claimed limitation was false: the Fourier transform of `exp(-u^4)` has
only real zeros by Pólya's theorem, as reviewed in
[Newman–Wu, Section 2.2](https://arxiv.org/html/1901.06596v2#S2.SS2).
Its positive squared-variable coefficient sequence is PF-infinity. This is
different from translation total positivity of the density itself. The old
README, generated control, recorded JSON, and associated assertion are
corrected in this change. The supposed counterexample to the *combined*
first/second-level concavity condition is withdrawn; no conclusion about that
combined condition follows from it.

The present investigation asks for a geometry that **generates** xi from
positive local weights, giving zero-location control through stability.

## 1. Every finite folded theta cutoff has infinitely many nonreal zeros

Use the previous experiment's normalization. For real `u >= 0`, set

\[
a_n=\pi n^2,\qquad
f_n(u)=a_n(2a_ne^{9u}-3e^{5u})e^{-a_ne^{4u}},
\qquad \Phi(u)=\sum_{n\ge1}f_n(u).
\]

Each summand is positive on the half-line. The full theta identity makes
`Phi` even and analytic at zero, and

\[
2\int_0^\infty\Phi(u)\cos(tu)\,du
=\frac14\xi\!\left(\frac12+\frac{it}{2}\right).
\]

For a finite `N >= 1`, define `g_N(u)=sum_{n=1}^N f_n(u)` on the half-line
and reflect it evenly. Its Fourier transform is the even entire function

\[
F_N(z)=2\int_0^\infty g_N(u)\cos(zu)\,du.
\]

**Proposition.** For every finite `N >= 1`, `F_N` has finitely many real zeros
and infinitely many nonreal zeros.

**Proof.** Direct differentiation gives

\[
f_n'(0)=a_ne^{-a_n}(-8a_n^2+30a_n-15).
\]

The larger root of the quadratic is `(15+sqrt(105))/8 < 4`, so every term
with `n >= 2` has negative derivative. The exact modular identity gives
`Phi'(0)=0`. Termwise differentiation is valid, since the derivative series
is absolutely convergent with terms `O(n^6 exp(-pi n^2))`. Consequently

\[
\delta_N:=g_N'(0)=-\sum_{n>N}f_n'(0)>0.
\]

Thus the finite even reflection has a cusp which the full infinite sum
cancels exactly. All derivatives of `g_N` on the half-line are integrable
and decay superexponentially. Four integrations by parts yield, for real `t`,

\[
F_N(t)=-\frac{2\delta_N}{t^2}
+\frac{2g_N'''(0)}{t^4}
+\frac{2}{t^4}\int_0^\infty g_N^{(4)}(u)\cos(tu)\,du
=-\frac{2\delta_N}{t^2}+O(t^{-4}).
\]

It follows that `F_N(t)<0` for all sufficiently large real `|t|`. Because
`F_N(0)>0` and a nonzero entire function has isolated zeros, it has finitely
many real zeros.

For `|z| <= r`, positivity gives

\[
|F_N(z)|\le2\int_0^\infty g_N(u)e^{ru}\,du.
\]

The substitution `x=exp(4u)` bounds this integral by a constant times an
incomplete gamma integral with parameter `(r+9)/4`. Hence
`log M_N(r)=O(r log(r+2))`, so the entire function has order at most one.
If it had finitely many zeros in total, Hadamard factorization would give
`F_N(z)=exp(Az+B)P(z)` with `P` polynomial. Evenness forces `A=0`.
This contradicts real-axis Fourier decay and `F_N(0)>0`. Thus infinitely
many of its zeros are nonreal. QED.

This is compatible with `F_N -> xi(1/2+iz/2)/4` locally uniformly:
on a compact set, the positive tails are dominated by
`exp(Ru)Phi(u)`, which is integrable. Nonreal zeros of the approximants can
escape to infinity. No inference about nonreal zeros of the limiting xi is
being made.

**Implication for positive geometry.** A finite theta sum taken as an exact
carrier cannot have a stable positive-Grassmannian lift or a ferromagnetic
Lee–Yang realization with this Fourier transform. This holds however small
the omitted tail is. It does not invalidate finite-precision use of theta
truncation to approximate exact coefficients. It invalidates promoting that
finite truncation to an exact all-rank positive carrier.

The theorem concerns truncation in the theta summation index and integration
over the full half-line. A separate hard cutoff in `u` changes the endpoint
asymptotics and is not covered by this statement.

**Gaussian-tilt corollary.** For every real `lambda`, replacing `g_N(u)` by
`exp(lambda*u^2)g_N(u)` leaves its first derivative at zero equal to
`delta_N>0`. Its derivatives still decay superexponentially and its Fourier
transform still has order at most one: absorb the Gaussian factor into half
of the decaying `exp(-pi exp(4u))` factor. The same proof gives infinitely
many nonreal zeros for every Gaussian tilt. Thus the set of real-zero
Gaussian tilts for this finite folded cutoff is empty. Applying a finite
de Bruijn–Newman heat deformation cannot repair this particular approximation.

### Cancelling finitely many boundary derivatives does not solve this problem

Under the same decay and entire-order hypotheses, if the first nonzero odd
derivative at zero is `g^(2m+1)(0)`, repeated integration by parts gives

\[
2\int_0^\infty g(u)\cos(tu)\,du
=2(-1)^{m+1}\frac{g^{(2m+1)}(0)}{t^{2m+2}}
+O(t^{-2m-4}).
\]

The same eventual-sign argument applies. Cancelling finitely many odd
derivatives only moves the obstruction if a later one remains nonzero.
Analytically even approximants, or smooth approximants with every odd
derivative zero, escape this particular obstruction. Their stability still
needs a separate proof. Classical smooth Pólya approximants are therefore
not ruled out; see [the later Pólya-approximation account](https://arxiv.org/abs/1502.06844).

## 2. An individual theta summand also fails the Lee–Yang test

For a symmetric Lee–Yang law, the canonical product for its characteristic
function gives the necessary cumulant inequality

\[
\kappa_4=\mathbb E U^4-3(\mathbb E U^2)^2\le0.
\]

Normalize the positive even density `f_n(|u|)`, and put
`V=4a_n |U|`. Its unnormalized density relative to its value at zero is

\[
\frac{f_n(v/(4a_n))}{f_n(0)}
=e^{5v/(4a_n)}\frac{2a_ne^{v/a_n}-3}{2a_n-3}
  e^{-a_n(e^{v/a_n}-1)}\longrightarrow e^{-v}.
\]

For all sufficiently large `a_n` this is bounded by `2 exp(-v/2)`, so
dominated convergence applies to normalization and the first four moments.
The signed scaled variable approaches the centered Laplace distribution,
whose fourth cumulant is 12. Therefore

\[
\kappa_4(U)=\frac{12+o(1)}{(4\pi n^2)^4}>0
\quad\hbox{for all sufficiently large }n.
\]

This rules out realizing the individual large-index theta components as
ferromagnetic magnetization limits. Numerical quadrature already finds a
positive fourth cumulant at `n=2`; the all-large-index statement follows
from the analytic limit, not that computation.

Positive mixing is no repair. The Lee–Yang laws `delta_0` and the symmetric
two-point law on `{-1,1}` have a mixture with characteristic function
`p+(1-p)cos(z)`. For `p>1/2` its zeros include
`pi +/- i arcosh(p/(1-p))`. Disjoint Ising networks instead convolve their
magnetization laws and multiply their partition functions.

## 3. The constructive Grassmannian target

Define the entire squared-variable function by its power series,

\[
G(w)=\frac{\xi(1/2+\sqrt w/2)}{\xi(1/2)}
=\sum_{k\ge0}c_k w^k,\qquad c_0=1,\ c_k>0.
\]

The square root causes no branch ambiguity because xi is even around `1/2`.
Set

\[
J_N(w)=\sum_{k=0}^N\frac{(N)_k}{N^k}c_k w^k,
\qquad (N)_k=N(N-1)\cdots(N-k+1).
\]

These are the scaled Jensen polynomials for `gamma_k=k!c_k`, and converge
locally uniformly to `G`. The falling-factorial factors matter: ordinary
Taylor truncations are not generally real-rooted even when the limit is.

[Purbhoo's Theorem 1.1](https://arxiv.org/abs/1611.07548) says the homogeneous
multiaffine polynomial with the Plücker coordinates of a Grassmannian point
as coefficients is stable precisely when that point is totally nonnegative.
Consequently an explicit construction of points

\[
V_N\in\operatorname{Gr}_{\ge0}(N,2N)
\]

whose Plücker polynomials, specialized at
`(x_1,x_2,...,x_(2N-1),x_(2N))=(1,w,...,1,w)`, equal `J_N(w)` up to a
positive scalar, would prove RH. Stability survives the real specialization
and identification; a real stable univariate polynomial has real zeros,
positive coefficients exclude nonnegative zeros, and Hurwitz passes the
zero-free complement of the negative axis to `G`.

This is a construction target, not a new RH proof or novelty claim. Building
`V_N` from the unknown roots of `J_N`, or from minors whose positivity already
assumes RH, is circular. The missing ingredient is source-defined positive
edge weights plus an exact coefficient identity for every `N`.

### Why averaging positive geometries does not supply that ingredient

For `a>0`, the matrix with rows `(1,a,0,0)` and `(0,0,1,a)` represents a
totally nonnegative point of `Gr(2,4)`. Its nonzero minors are
`p13=1`, `p14=p23=a`, `p24=a^2`.
Average these minors over a nondegenerate positive random `a`. Their
Plücker defect becomes

\[
\overline p_{13}\overline p_{24}
-\overline p_{14}\overline p_{23}
-\overline p_{12}\overline p_{34}
=\operatorname{Var}(a)>0.
\]

The average leaves the Grassmannian. Its specialization
`1+2 E[a]w+E[a^2]w^2` has discriminant `-4 Var(a)<0`.
Thus positive theta weights cannot be inserted into a general
stability-preserving averaging step. A xi-specific identity must do more.

## 4. A more physical realization: positive orthogonal Grassmannian

[Galashin–Pylyavskyy](https://arxiv.org/abs/1807.03282) identify planar Ising
zero-field boundary correlation data with the totally nonnegative orthogonal
Grassmannian. This supplies genuine positive geometry, including an
interpretation of Kramers–Wannier duality. It does not itself identify the
full external-field partition function with a Plücker polynomial or with xi.

The concrete alternative target is a sequence of planar ferromagnetic
networks with `J_ij >= 0` and nonnegative magnetization weights `b_j,N`, such
that the laws of `M_N=sum_j b_j,N sigma_j` converge weakly to the normalized
full density `Phi(u)du` on the real line. Lee–Yang and the weak-closure theorem
of [Newman–Wu, Theorem 16](https://arxiv.org/html/1901.06596v2#S3.SS1)
would then prove RH. The density construction must be independent of RH.

The missing object is a **full field-dependent network construction** with
that limit. Matching boundary two-point correlations, enforcing a symmetry
of the Grassmannian chart, or mixing theta components does not establish it.

The strongest next target from this investigation is to perform modular
completion before constructing the positive geometry, and to establish
stability for the complete generating function. A proposed finite kernel
should first survive the odd-derivative obstruction; proposed components
should survive the fourth-cumulant obstruction; then the complete field
polynomial or its exact Plücker lift must be checked. These are discriminating
tests for a proposed construction, not an automatic research queue.

## Verification and limits

The finite-cutoff proof and its odd-derivative extension were independently
reviewed by Sol; Astra checked the integration-by-parts signs, gamma bound,
normalization, and the distinction between local uniform convergence and
zero control. The theta-component cumulant argument was derived by Sol and
checked by Astra. Neither review is external peer review.

`probe.py` checks the boundary derivative against automatic differentiation,
the incomplete-gamma integral formula against direct quadrature, finite-cutoff
real-axis asymptotics, and the component cumulants. The gamma formula used is

\[
\int_0^\infty f_n(u)e^{su}\,du
=\frac{a_n^{-(s+1)/4}}4
\left[2\Gamma\!\left(\frac{s+9}{4},a_n\right)
-3\Gamma\!\left(\frac{s+5}{4},a_n\right)\right].
\]

Run with Python and `mpmath==1.3.0`:

```sh
python probe.py --dps 80 --output results/checks-80dps.json
python probe.py --dps 120 --output results/checks-120dps.json
```

Quadrature uses explicit finite bounds, and all reported numerical values
remain nonrigorous approximations. The theorem proofs do not depend on
quadrature or on any finite list of zeros. The literature check was targeted;
external novelty of the cutoff proposition is unassessed. RH remains open.
