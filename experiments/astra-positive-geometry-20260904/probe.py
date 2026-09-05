#!/usr/bin/env python3
"""Small numerical checks for the analytic arguments in README.md.

No zero census, interval certificate, or RH inference is made by this script.
The finite-cutoff obstruction is proved analytically in the accompanying note.
"""
import argparse
import json
from pathlib import Path

import mpmath as mp


def term(n, u):
    a = mp.pi * n * n
    return a * (2 * a * mp.exp(9*u) - 3 * mp.exp(5*u)) * mp.exp(-a*mp.exp(4*u))


def boundary_derivative(n):
    a = mp.pi * n * n
    return a * mp.exp(-a) * (-8*a*a + 30*a - 15)


def laplace_term(n, s):
    a = mp.pi * n * n
    return a**(-(s+1)/4) / 4 * (
        2*mp.gammainc((s+9)/4, a, mp.inf)
        - 3*mp.gammainc((s+5)/4, a, mp.inf)
    )


def fourier_cutoff(nmax, t):
    return mp.fsum(laplace_term(n, 1j*t) + laplace_term(n, -1j*t)
                   for n in range(1, nmax+1))


def component_cumulant(n):
    a = mp.pi * n * n
    # V=4*a*|U|; scale before integration to retain small-component accuracy.
    def weight(v):
        return (mp.exp(5*v/(4*a)) * (2*a*mp.exp(v/a)-3)/(2*a-3)
                * mp.exp(-a*mp.expm1(v/a)))
    raw = [mp.quad(lambda v: v**k * weight(v), [0, 1, 2, 4, 8, 16, 32, 64, 128])
           for k in (0, 2, 4)]
    kappa_v = raw[2]/raw[0] - 3*(raw[1]/raw[0])**2
    return kappa_v/(4*a)**4, kappa_v


def run(dps):
    mp.mp.dps = dps
    render = lambda x: mp.nstr(mp.re(x), 32)
    checks = []
    for n in (1, 2):
        explicit = boundary_derivative(n)
        automatic = mp.diff(lambda u: term(n, u), 0)
        error = abs((automatic-explicit)/explicit)
        assert error < mp.mpf('1e-40')
        checks.append({'check': 'boundary_derivative', 'n': n, 'relativeError': render(error)})
        for t in (0, 3):
            integral = 2*mp.quad(lambda u: term(n, u)*mp.cos(t*u),
                                [0, mp.mpf('.25'), mp.mpf('.5'), 1, 2])
            formula = fourier_cutoff(n, t) - (fourier_cutoff(n-1, t) if n>1 else 0)
            error = abs((integral-formula)/integral)
            assert error < mp.mpf('1e-40')
            checks.append({'check': 'incomplete_gamma_vs_quadrature', 'n': n,
                           't': t, 'relativeError': render(error)})
    cutoff_results = []
    for nmax in (1, 2, 4):
        direct = mp.fsum(boundary_derivative(n) for n in range(1, nmax+1))
        tail = -mp.fsum(boundary_derivative(n) for n in range(nmax+1, nmax+25))
        assert direct > 0 and tail > 0
        relative = abs((direct-tail)/tail)
        assert relative < mp.mpf('1e-30')
        samples = []
        for factor in (32, 64):
            t = factor*(nmax+1)**2
            value = mp.re(fourier_cutoff(nmax, t))
            asymptotic = -2*tail/t**2
            assert value < 0
            samples.append({'t': t, 'fourierValue': render(value),
                            'ratioToLeadingAsymptotic': render(value/asymptotic)})
        cutoff_results.append({'N': nmax, 'positiveBoundaryDerivative': render(tail),
                               'directVsTailRelativeError': render(relative),
                               'realAxisSamples': samples})
    components = []
    for n in (1, 2, 4, 16):
        cumulant, scaled = component_cumulant(n)
        if n >= 2:
            assert cumulant > 0
        components.append({'n': n, 'fourthCumulant': render(cumulant),
                           'scaledFourthCumulant_limit12': render(scaled)})
    return {'status': 'numerical_checks_of_analytic_propositions_not_interval_certified',
            'decimalPrecision': dps, 'checks': checks,
            'foldedThetaCutoffs': cutoff_results, 'individualThetaComponents': components}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dps', type=int, default=80)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.dps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'checks': len(result['checks']),
                      'cutoffs': len(result['foldedThetaCutoffs']),
                      'components': len(result['individualThetaComponents'])}))
