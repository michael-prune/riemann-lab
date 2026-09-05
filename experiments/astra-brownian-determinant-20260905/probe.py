#!/usr/bin/env python3
"""Reproduce the finite calculations in README.md; never infer RH from them.

The sextic control uses exact rational inequalities. Theta derivative bounds
use mpmath interval arithmetic and an analytic tail majorant. The quadrature
and approximate complex zeros are precision checks, not zero certificates.
"""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp


def sextic_certificate():
    # Gaussian reference density exp(-y*y)/sqrt(pi), epsilon=100**(-3).
    eps = F(1, 10**6)
    gm = lambda k: F(math.factorial(2*k), 4**k * math.factorial(k))
    p = [(gm(k)-eps*gm(k+3),
          gm(k)-eps*gm(k+3)+eps**2*gm(k+6)/2) for k in range(4)]
    assert all(a > 0 and a <= b for a, b in p)
    lo = p[3][0]*p[0][0]**2 - 15*p[2][1]*p[1][1]*p[0][1] + 30*p[1][0]**3
    hi = p[3][1]*p[0][1]**2 - 15*p[2][0]*p[1][0]*p[0][0] + 30*p[1][1]**3
    assert lo <= hi < 0
    # kappa_6(Y)=N/p_0^3; X=Y/10.
    k_lo, k_hi = lo/p[0][0]**3 / 10**6, hi/p[0][1]**3 / 10**6
    assert k_lo <= k_hi < 0
    return {'b': 100, 'epsilon': str(eps),
            'sixthCumulantNumeratorBounds': [str(lo), str(hi)],
            'sixthCumulantXBounds': [str(k_lo), str(k_hi)]}


def xi(s):
    return s*(s-1)*mp.power(mp.pi, -s/2)*mp.gamma(s/2)*mp.zeta(s)/2


def theta_term(n, u):
    a = mp.pi*n*n
    return a*(2*a*mp.exp(9*u)-3*mp.exp(5*u))*mp.exp(-a*mp.exp(4*u))


def theta_tail_integral(n, u):
    a = mp.pi*n*n
    x = a*mp.exp(4*u)
    return a**(-mp.mpf(1)/4)/4 * (
        2*mp.gammainc(mp.mpf(9)/4, x, mp.inf)
        - 3*mp.gammainc(mp.mpf(5)/4, x, mp.inf))


def bell_jet(h, count, ctx):
    out = [ctx.mpf(1)]
    for n in range(1, count+1):
        out.append(sum(math.comb(n-1, j-1)*h[j-1]*out[n-j]
                       for j in range(1, n+1)))
    return out


def interval_potential(nmax=12, order=12):
    iv = mp.iv
    pi = iv.pi
    derivatives = [iv.mpf(0) for _ in range(order+1)]
    for n in range(1, nmax+1):
        a = pi*n*n
        jets = [bell_jet([k*(j == 1)-a*4**j for j in range(1, order+1)],
                         order, iv) for k in (9, 5)]
        for j in range(order+1):
            derivatives[j] += a*iv.exp(-a)*(2*a*jets[0][j]-3*jets[1][j])
    # Cauchy on |u|<=r: Re exp(4u)>=exp(-4r)cos(4r)=beta.
    r = iv.mpf(1)/100
    beta = iv.exp(-4*r)*iv.cos(4*r)
    first = nmax+1
    a = pi*first*first
    q = (iv.mpf(first+1)/first)**4 * iv.exp(-beta*pi*(2*first+1))
    assert q.b < 1
    majorant = a*(2*a*iv.exp(9*r)+3*iv.exp(5*r))*iv.exp(-beta*a)/(1-q)
    tails = []
    for j in range(order+1):
        radius = (majorant*math.factorial(j)/r**j).b
        derivatives[j] += iv.mpf([-1, 1])*radius
        tails.append(mp.mpf(radius))
    # f(t)=Phi(sqrt(t))/Phi(0), and v(t)=-log f(t).
    f = [derivatives[2*k]/math.factorial(2*k)/derivatives[0]
         for k in range(order//2+1)]
    f[0] = iv.mpf(1)
    logarithm = [iv.mpf(0)]
    for n in range(1, len(f)):
        logarithm.append(f[n]-sum(k*logarithm[k]*f[n-k]
                                  for k in range(1, n))/n)
    v = [-x for x in logarithm]
    assert mp.mpf(v[3].b) < 0
    # Necessary power-sum inequalities for peeling off Q=prod(1+alpha_j*u^2).
    p4_upper, p5_lower, p6_upper = 4*v[4], -5*v[5], 6*v[6]
    gap = p5_lower**2-p4_upper*p6_upper
    assert gap.a > 192835692
    bounds = lambda x: [mp.nstr(mp.mpf(x.a), 35), mp.nstr(mp.mpf(x.b), 35)]
    # Decimal output is display-only; signs above are checked before rendering.
    return {'intervalArithmetic': 'mpmath.iv', 'thetaTerms': nmax,
            'cauchyRadius': '1/100',
            'derivative12TailBound': mp.nstr(tails[12], 12),
            'potentialCoefficientEnclosures_displayRounded': [bounds(x) for x in v[1:]],
            'powerSumContradictionGap_displayRounded': bounds(gap)}


def modular_cutoff(s, count):
    def a(n, z):
        t = mp.pi*n*n
        return t**(-z/2)*mp.gammainc(z/2, t, mp.inf)
    return mp.mpf(1)/2+s*(s-1)/2*mp.fsum(a(n, s)+a(n, 1-s)
                                         for n in range(1, count+1))


def run(dps):
    mp.mp.dps = dps
    mp.iv.dps = dps
    fmt = lambda x: mp.nstr(x, 35)
    nmax = 12
    phi = lambda u: mp.fsum(theta_term(n, u) for n in range(1, nmax+1))
    nodes = [0, mp.mpf('.125'), mp.mpf('.25'), mp.mpf('.5'), 1, 2]
    raw = [mp.quadgl(lambda u: u**k*phi(u), nodes) for k in (0, 2, 4, 6)]
    assert abs(raw[0]-xi(mp.mpf('.5'))/8) < mp.mpf('1e-50')
    m2, m4, m6 = (a/raw[0] for a in raw[1:])
    c1, c2, c3 = m2/2, m4/24, m6/720
    # Independent comparison against derivatives of completed zeta.
    direct = mp.taylor(lambda z: xi(mp.mpf('.5')+z/2)/xi(mp.mpf('.5')), 0, 6)
    coefficient_errors = [abs(c-direct[2*k]) for k, c in enumerate((c1, c2, c3), 1)]
    assert max(coefficient_errors) < mp.mpf('1e-50')
    survival = lambda t: mp.fsum(theta_tail_integral(n, t)
                                for n in range(1, nmax+1))/raw[0]
    trace_k2 = mp.mpf(2)/3*mp.quadgl(lambda t: t**3*survival(t)**2, nodes)
    naive_c2 = (c1**2-trace_k2)/2
    independent_c2 = m2**2/24
    assert 0 < independent_c2 < naive_c2 < c2
    # Direct quadrature of the determinant identity at three w values.
    determinant_checks = []
    for w in (1, 10, 100):
        integral = mp.quadgl(lambda u: mp.cosh(u*mp.sqrt(w))*phi(u), nodes)/raw[0]
        target = xi(mp.mpf('.5')+mp.sqrt(w)/2)/xi(mp.mpf('.5')) if w != 1 else mp.mpf('.5')/xi(mp.mpf('.5'))
        error = abs(integral-target)
        assert error < mp.mpf('1e-45')
        determinant_checks.append({'w': w, 'value': fmt(integral), 'absoluteError': fmt(error)})
    # Positive two-Gamma model: zeros are proved analytically in the note.
    bracket = lambda z: 3*z+1+(3*z+11)*mp.exp(-mp.log(4)*(z+2))
    gamma_roots = []
    for k in (1, 2):
        y = mp.findroot(lambda v: mp.log(4)*v-2*mp.atan(3*v/5)-2*mp.pi*k,
                        (4*k, 5*k+2))
        residual = abs(bracket(-2+1j*y))
        assert residual < mp.mpf(10)**(-dps+10)
        gamma_roots.append({'y': fmt(y), 'characteristicZeroRealPart': fmt(4*y),
                            'characteristicZeroImagPart': '9', 'bracketResidual': fmt(residual)})
    z = mp.findroot(lambda t: modular_cutoff(mp.mpf('.5')+1j*t/2, 1),
                    (mp.mpc(45, 9), mp.mpc(46, 10)))
    residual = abs(modular_cutoff(mp.mpf('.5')+1j*z/2, 1))
    assert abs(mp.im(z)) > 9 and residual < mp.mpf(10)**(-dps+10)
    return {'status': 'analytic_identities_and_separately_labelled_finite_checks_no_RH_proof',
            'decimalPrecision': dps, 'quadratureThetaTerms': nmax, 'quadratureUMax': 2,
            'coefficients': list(map(fmt, (c1, c2, c3))),
            'coefficientComparisonAbsoluteErrors': list(map(fmt, coefficient_errors)),
            'brownianDeterminantIdentityChecks': determinant_checks,
            'operatorCandidates': {'targetTraceSquare': fmt(c1**2-2*c2),
                'averagedOperatorTraceSquare': fmt(trace_k2),
                'averagedOperatorCoefficient2': fmt(naive_c2),
                'averagedOperatorCoefficient2Ratio': fmt(naive_c2/c2),
                'independentModesCoefficient2': fmt(independent_c2),
                'missingVarianceOver24': fmt(c2-independent_c2)},
            'twoGammaApproximateZeros': gamma_roots,
            'modularCompletionApproximateZero': {'realPart': fmt(mp.re(z)), 'imagPart': fmt(mp.im(z)),
                'residual': fmt(residual), 'status': 'numerical_not_interval_certified'},
            'thetaPotentialCertificate': interval_potential(),
            'sexticExactRationalCertificate': sextic_certificate()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dps', type=int, default=80)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 70:
        parser.error('use at least 70 decimal digits')
    result = run(args.dps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'decimalPrecision': args.dps}))
