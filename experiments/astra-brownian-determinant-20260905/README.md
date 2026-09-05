# Xi as an average of positive determinants

Date: 2026-09-05. Astra continuation with two independently tasked Sol
mathematical investigators.

**Result:** an exact Brownian determinant representation identifies where a
positive construction for xi loses stability. Independent modes, a random
stopping interval, and the straightforward averaged operator do not repair
it. An exact rational sextic counterexample settles the earlier open question
about two simultaneous kernel-concavity conditions. An interval calculation
also excludes a specified class of Laguerre-potential factorizations of the
full theta kernel.

These are explicit identities, analytic arguments, and separately identified
computer checks. They do not establish RH, an all-order stable Grassmannian
lift, external novelty, or admission under the lab's Inner Core protocol.
The Brownian and Gamma representations use classical ingredients; their
availability alone is not an advance on RH.

## 1. A source-defined positive determinant identity

Retain the normalization of the
[previous investigation](../astra-positive-geometry-20260904/README.md):

\[
\Phi(u)=\sum_{n\ge1}\pi n^2(2\pi n^2e^{9u}-3e^{5u})e^{-\pi n^2e^{4u}},
\quad u\ge0,\qquad
d\mu(u)=\frac{\Phi(u)\,du}{\int_0^\infty\Phi(v)\,dv}.
\]

Write `U` for a nonnegative random variable with law `mu`. Then

\[
G(w)=\frac{\xi(1/2+\sqrt w/2)}{\xi(1/2)}
=\mathbb E\cosh(U\sqrt w)
=\sum_{k\ge0}\frac{\mathbb EU^{2k}}{(2k)!}w^k.
\]

Both expressions in `sqrt(w)` mean their even power series, so there is no
branch ambiguity. RH is equivalent to all zeros of `G` being negative real.

On `L²(0,u)`, let `B_u` be the Brownian covariance operator

\[
B_u(s,t)=\min(s,t).
\]

It is positive and trace class, with eigenvalues

\[
\lambda_j(B_u)=\frac{4u^2}{\pi^2(2j+1)^2},\qquad j=0,1,\ldots.
\]

Consequently

\[
\boxed{G(w)=\mathbb E_\mu\det(I+wB_U).}
\]

This identity also has a direct positive-minor proof. For ordered points
`0<t_1<...<t_k<u`, subtracting successive rows and columns gives

\[
\det[\min(t_i,t_j)]
=t_1(t_2-t_1)\cdots(t_k-t_{k-1}).
\]

Its simplex integral is `u^(2k)/(2k)!`, exactly the `k`th Fredholm
coefficient. Averaging therefore gives the coefficients of `G`. Exchanges
of expectation and expansion on compact `w` sets are justified by
`E cosh(U sqrt(|w|))<infinity`, using the superexponential theta tail.

Thus positive minors already exist **before** the common random length is
integrated out. What is absent is a theorem that the resulting average is
one stable determinant. A positive deterministic trace-class `L` satisfying
`det(I+wL)=G(w)`, constructed without assuming xi's zeros are real, would
settle RH. The displayed average does not provide that `L`.

### The obvious averaged operator gives the wrong second coefficient

Embed every `B_u` by zero extension into `L²(0,infinity)` and set
`S(t)=Pr(U>=t)`. The operator average is positive and trace class, with kernel

\[
K(s,t)=\min(s,t)S(\max(s,t)).
\]

Direct integration gives

\[
\operatorname{Tr}K=\frac{\mathbb EU^2}{2}=c_1,
\qquad
\operatorname{Tr}K^2=\frac23\int_0^\infty t^3S(t)^2\,dt.
\]

For the desired determinant, the second trace must instead be
`c_1²-2c_2`. Numerical quadrature gives

| Quantity | Value, rounded |
| --- | ---: |
| `c_1` | `0.00577624827885474` |
| `c_2` | `0.0000155208783618215` |
| required `Tr(L²)` | `0.00000232328745532936` |
| actual `Tr(K²)` | `0.0000103945079256276` |
| `[w²]det(I+wK)/c_2` | `0.739988282810335` |

The `S²` in this formula is essential: replacing it by `S` computes an
expectation of traces of squares, rather than the trace of the square of the
expectation. The mismatch is a numerical check of this particular candidate;
it is not a theorem excluding arbitrary positive operators.

## 2. Why three standard stability mechanisms do not close the gap

### Independent rank-one modes discard a required variance

On a fixed unit interval, write
`B_1=sum_j lambda_j e_j e_j*`. Replacing the common scale `U²` with
independent copies gives the positive rank-one sum

\[
A_{\rm ind}=\sum_j\lambda_jU_j^2e_je_j^*.
\]

Independence and determinant multi-affinity imply

\[
\mathbb E\det(I+wA_{\rm ind})
=\prod_j(1+w\lambda_j\mathbb EU^2)
=\cosh(\sqrt{\mathbb EU^2}\sqrt w).
\]

This is consistent with the real-rootedness mechanism in
[Marcus–Spielman–Srivastava](https://arxiv.org/abs/1306.3969). It matches
`c_1`, but its second coefficient is too small by exactly

\[
\boxed{c_2-[w^2]\mathbb E\det(I+wA_{\rm ind})
=\frac{\operatorname{Var}(U^2)}{24}>0.}
\]

For xi the deficit is approximately `0.00000996003766532611`.
General nonorthogonal covariance matrices remain possible, but they need a
source-defined rule reproducing the higher theta moments through mixed
discriminants. Independence by itself removes those moments.

### A stopping interval is not strongly Rayleigh

For two mesh points `r_1<r_2`, the occupancy indicators
`1_{U>=r_j}` have a prefix law with generating polynomial

\[
p_0+p_1z_1+p_2z_1z_2,\qquad p_0p_2>0.
\]

A nonnegative bivariate multiaffine polynomial
`a+bz_1+cz_2+dz_1z_2` is stable only if `bc>=ad`.
Here `c=0`, so the inequality fails. Since marginalization preserves
stability, a larger strongly Rayleigh ensemble cannot have this exact
occupancy marginal. See
[Borcea–Brändén–Liggett](https://arxiv.org/abs/0707.2340) for the general
framework. This argument concerns the nested occupancy variables, not every
possible representation of the scalar function `G`.

### Conditional Brownian determinants have no common interlacer

The zeros for length `u` are `-pi²(2j+1)²/(4u²)`. Their count in `[-R,0]`
is `u sqrt(R)/pi+O(1)`. For two distinct positive lengths, the count
difference is unbounded. Sequences with a common interlacer have bounded
count difference. Thus the exact continuum family does not support this
standard averaging argument either.

## 3. The Brownian-excursion source does not automatically supply exterior powers

[Biane–Pitman–Yor](https://arxiv.org/html/math/9912170v1) give

\[
\Sigma_2=\frac2{\pi^2}\sum_{n\ge1}\frac{\Gamma_{2,n}}{n^2},\quad
Y\overset d=\sqrt{\frac\pi2\Sigma_2}
\overset d=\sqrt{\frac2\pi}\max e_t,\quad
\mathbb EY^s=2\xi(s).
\]

Here the Gamma variables are independent with shape 2 and unit scale, and
`e_t` is a standard Brownian excursion. Under the `Y^(1/2)` size bias,
`V=(1/2)log Y` is symmetric by the functional equation and
`G(w)=E exp(V sqrt(w))`. The law of `|V|` is the `mu` above.

With `A=diag(2/(pi²n²))`, the Gamma sum has Laplace transform

\[
\mathbb Ee^{-\lambda\Sigma_2}
=\det(I+\lambda A)^{-2}
=\left(\frac{\sqrt{2\lambda}}{\sinh\sqrt{2\lambda}}\right)^2.
\]

This positive inverse determinant concerns the additive variable
`Sigma_2`. The target determinant concerns its logarithm after size bias.
The following explicit computation shows why the coordinate change matters.

### Two Gamma modes already lose the required zero location

Let `X=Gamma_2+q Gamma'_2`, `0<q<1`. Convolution and Gamma integration give

\[
M_q(z):=\mathbb EX^z
=\frac{\Gamma(z+1)}{(1-q)^3}
\left[(1-q)(z+1)-2q
+q^{z+2}\{2+(1-q)(z+1)\}\right],\qquad\Re z>-4,
\]

with removable singularities interpreted by continuity. For the actual
first two BPY weights, `q=1/4`, the bracket is one quarter of

\[
C(z)=3z+1+(3z+11)e^{-(z+2)\log4}.
\]

At `z=-2+iy`, `C(z)=0` iff

\[
(\log4)y-2\arctan(3y/5)=2\pi k.
\]

For `y>=0` the left side has strictly positive derivative, at least
`log4-6/5>0`, and tends to infinity. There is exactly one positive solution
for each integer `k>=1`. These are genuine Mellin zeros, since the Gamma
factor is finite and nonzero there. The solution `y=0` is different: it
cancels a Gamma pole, with

\[
M_{1/4}(-2)=\frac{16}{27}(5\log4-6)>0.
\]

Under `X^(1/4)` bias and `V=(1/4)log X`, the characteristic transform is
`M(1/4+it/4)/M(1/4)`. Its zeros therefore include

\[
t=\pm25.731184113960862\ldots+9i.
\]

They lie inside its analytic domain `Im(t)<17`. This finite model has only
a half-plane of analyticity, not an entire transform. Neither its positive
Gamma components nor their additive convolution supply Lee–Yang stability
after the logarithmic change of variable. Multiplying `X` by a constant
does not change these zeros.

## 4. A specified Laguerre-potential factorization is excluded

Set

\[
v(t)=-\log\frac{\Phi(\sqrt t)}{\Phi(0)}=\sum_{k\ge1}v_kt^k.
\]

One sufficient single-site Lee–Yang criterion asks that `b+v'` belong to
the Laguerre class for some `b>=0`. Such functions have nonnegative Taylor
coefficients. The criterion is due to Kozitsky; see
[Kozitsky–Pasurek, Proposition 8.2](https://bibos.math.uni-bielefeld.de/preprints/06-12-238.pdf)
and [Kozitsky (2003)](https://doi.org/10.1016/S0096-3003(02)00324-7).

An interval evaluation of the full theta derivatives, including a Cauchy
bound for the omitted tail, gives the following rounded values:

\[
(v_1,\ldots,v_6)\simeq
(37.4538098590,95.3432918410,-137.6662708212,
1549.9260979882,-15221.4359670945,150530.4870557805).
\]

The negative `v_3` excludes this criterion directly. A Gaussian tilt changes
only `v_1`, so it cannot fix this obstruction.

There is also a stronger, precisely scoped exclusion. Suppose

\[
\Phi(u)=e^{\lambda u^2}Q(u)\,C e^{-h(u^2)},\qquad
Q(u)=\prod_j(1+\alpha_j u^2),\quad
\alpha_j\ge0,\quad\sum_j\alpha_j<\infty,
\]

and try to place `b+h'` in the same Laguerre class. These are the usual
even imaginary-zero polynomial factors and their convergent products.
Writing `p_k=sum_j alpha_j^k`, the nonnegative coefficient conditions at
orders 3, 4, 5 of `h'` require

\[
p_4\le4v_4,\qquad p_5\ge-5v_5,\qquad p_6\le6v_6.
\]

But the nonnegative power sums satisfy `p_5²<=p_4 p_6`, whereas the interval
calculation establishes

\[
(-5v_5)^2-(4v_4)(6v_6)>192835692>0.
\]

This is a contradiction. Thus Gaussian tilts and these standard factors
cannot turn the exact theta potential into this particular sufficient
class. It does not exclude a different Lee–Yang construction.

For the tail certificate, on `|u|<=r=1/100` use
`Re exp(4u)>=beta=exp(-4r)cos(4r)`. The `n`th summand is bounded by

\[
a_n(2a_ne^{9r}+3e^{5r})e^{-\beta a_n},\quad a_n=\pi n^2.
\]

Beyond the retained terms, the ratio is bounded by
`(1+1/n)^4 exp(-beta pi(2n+1))`. Sum the resulting geometric majorant and
apply Cauchy's derivative bound. With 12 terms the twelfth-derivative tail
is less than `1.29e-183`. All coefficient operations and sign comparisons
are performed in interval arithmetic before decimal display rounding.

## 5. An exact sextic counterexample to the two-concavity shortcut

The previous investigation correctly withdrew the false negative control
`exp(-x^4)`: its Fourier transform has only real zeros. Here is a replacement
with an explicit parameter and an exact proof.

**Proposition.** The density proportional to `exp(-x^6-100x²)` has a Fourier
transform with a nonreal zero, although both of the following functions are
strictly log-concave for `t>0`:

\[
s(t)=e^{-t^3-100t},\qquad f(t)=s'(t)^2-s(t)s''(t).
\]

**Proof.** Directly, `f(t)=6t s(t)²`, so

\[
(\log s)''=-6t<0,\qquad(\log f)''=-12t-t^{-2}<0.
\]

To certify the existence of a zero off the real axis without numerical root finding, set
`Y=10X`, `epsilon=10^(-6)`, and use expectation under the reference density
`exp(-y²)/sqrt(pi)`. Its even moments are
`g_(2m)=(2m)!/(4^m m!)`. For `k=0,1,2,3` put

\[
p_k=\mathbb E_0[Y^{2k}e^{-\epsilon Y^6}],\quad
L_k=g_{2k}-\epsilon g_{2k+6},\quad
U_k=L_k+\frac{\epsilon^2}2g_{2k+12}.
\]

The elementary bounds `1-a<=exp(-a)<=1-a+a²/2` for `a>=0` give
`0<L_k<=p_k<=U_k`. The sixth cumulant has numerator

\[
\kappa_6(Y)=\frac{N}{p_0^3},\qquad
N=p_3p_0^2-15p_2p_1p_0+30p_1^3.
\]

Exact rational arithmetic yields

\[
N\le U_3U_0^2-15L_2L_1L_0+30U_1^3
=-\frac{60156783386913877734105236551761009}
{5368709120000000000000000000000000000000}
<-10^{-5}.
\]

Thus `kappa_6(X)<0`. Its even moment-generating function is entire of order
at most `6/5<2`, by a standard Young-inequality bound on `zx-x^6`.
If its zeros were all imaginary, pairing them in Hadamard factorization
would give `M(z)/M(0)=prod_j(1+z²/r_j²)`. This forces

\[
\kappa_6(X)=6![z^6]\log M(z)
=240\sum_jr_j^{-6}\ge0,
\]

a contradiction. The Fourier transform therefore has a nonreal zero. QED.

This shows that the two concavity assumptions alone do not give an
all-order stable lift. The proof does not depend on the earlier numerical
controls or on an approximate complex zero. Classical context for such
sextic failures is [Newman (1976)](https://sites.math.northwestern.edu/~auffing/papers/Newman.pdf).
No novelty is asserted for sextic Lee–Yang counterexamples.

## 6. A further test: exact functional symmetry is insufficient for one natural cutoff

Define

\[
A_n(s)=(\pi n^2)^{-s/2}\Gamma(s/2,\pi n^2),\qquad
\xi_N(s)=\frac12+\frac{s(s-1)}2
\sum_{n=1}^N[A_n(s)+A_n(1-s)].
\]

Each function is entire, has exact functional symmetry, and converges to
xi locally uniformly. Unlike the folded cutoffs in the previous note,
these have the symmetry built in. Nevertheless the `N=1` function
`xi_1(1/2+it/2)` has a **numerically located** off-axis zero at

\[
t\simeq45.8309315411216466283+9.33180405384913818388i.
\]

The two precision runs reproduce it with tiny residuals. A residual alone
is not an interval proof of a zero. This observation concerns this explicit
Mellin completion, not all finite models with modular structure.

## What remains open

The positive-geometry opportunity is to derive a coherent nonintersecting
path or mixed-discriminant construction **after** accounting for the theta
length correlations. It must reproduce `Var(U²)/24` at degree two and the
complete moment sequence at every degree. The constructions above do not
achieve that. Merely naming an operator with xi's presumed real zeros as
eigenvalues, or checking more finite minors, would not fill the gap.

The exact sextic theorem closes a proposed sufficient-condition shortcut.
The Brownian identity and failed constructions sharpen the remaining
question; they provide no evidence that the missing coupling exists.

## Reproduction and evidence boundaries

Use Python with the pinned dependency in `requirements.txt`:

```sh
python probe.py --dps 80 --output results/80dps.json
python probe.py --dps 120 --output results/120dps.json
```

The script checks the first three coefficients against direct completed-zeta
derivatives, the determinant identity at three positive arguments, and the
candidate operator coefficients. The quadrature uses 12 theta terms and
`0<=u<=2`; these floating-point integrals and the modular-cutoff root remain
numerical checks. The separate potential calculation includes an analytic
tail enclosure. The sextic certificate uses only exact rational arithmetic.
Rounded interval endpoints in JSON are display values, not replacement
certificates; the script checks the unrounded intervals.

Both recorded precision runs agree on the 35 displayed significant digits
of the compared values (apart from residuals). The interval and rational
certificates pass at both precisions. The repository's 23 contract tests
also pass; those tests are repository checks, not verification of RH.
