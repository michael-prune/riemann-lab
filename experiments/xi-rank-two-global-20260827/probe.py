#!/usr/bin/env python3
"""Probe the global rank-two xi exchange inequality.

The xi scan and the control are numerical.  The saddle budget isolates an
analytic proof obligation but does not certify the omitted full-kernel saddle
perturbation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.special import logsumexp, polygamma, roots_legendre


def text(value, digits=24):
    return mp.nstr(value, digits, strip_zeros=False)


def xi_dense_scan():
    node_count = 5000
    nodes, weights = roots_legendre(node_count)
    upper = 5.0
    u = (nodes + 1) * upper / 2
    log_u = np.log(u)
    log_weights = np.log(weights * upper / 2)
    e4 = np.exp(4 * u)
    theta_logs = []
    for n in range(1, 13):
        b = 2 * np.pi * n * n * e4
        theta_logs.append(
            np.log(np.pi * n * n) + 5 * u + np.log(b - 3) - b / 2
        )
    log_phi = logsumexp(np.asarray(theta_logs), axis=0)
    samples = np.unique(
        np.r_[np.linspace(0, 100, 2001), np.geomspace(100.1, 1.0e4, 500)]
    )
    minimum_margin = (np.inf, None)
    largest_ratio = (-np.inf, None)
    failures = []
    selected = {}
    requested = {0.0, 0.1, 0.5, 1.0, 5.0, 20.0, 100.0, 1000.0, 10000.0}
    for x in samples:
        log_density = log_weights + 2 * x * log_u + log_phi
        probabilities = np.exp(log_density - logsumexp(log_density))
        centered = log_u - np.dot(probabilities, log_u)
        kappa3 = np.dot(probabilities, centered**3)
        psi2 = float(polygamma(2, 2 * x + 1))
        margin = kappa3 - psi2
        ratio = kappa3 / psi2
        if margin < minimum_margin[0]:
            minimum_margin = (margin, x)
        if ratio > largest_ratio[0]:
            largest_ratio = (ratio, x)
        if margin <= 0:
            failures.append({"x": x, "margin": margin})
        if x in requested:
            selected[str(x)] = {
                "momentKappa3": f"{kappa3:.17g}",
                "gammaPsi2": f"{psi2:.17g}",
                "margin": f"{margin:.17g}",
                "ratio": f"{ratio:.17g}",
            }
    return {
        "status": "double_precision_gauss_legendre_not_interval_certified",
        "nodeCount": node_count,
        "thetaTerms": 12,
        "sampleCount": len(samples),
        "range": [0, 10000],
        "failureCount": len(failures),
        "failures": failures[:10],
        "smallestAbsoluteMargin": {
            "value": f"{minimum_margin[0]:.17g}",
            "x": minimum_margin[1],
        },
        "largestKappaToPsiRatio": {
            "value": f"{largest_ratio[0]:.17g}",
            "x": largest_ratio[1],
        },
        "selected": selected,
    }


def decreasing_log_concave_control():
    mp.mp.dps = 80

    def potential(u):
        return 2 * u if u <= 1 else 3 * u - 1

    raw = [
        mp.quad(
            lambda u, power=power: mp.exp(-potential(u)) * mp.log(u) ** power,
            [0, 1, mp.inf],
        )
        for power in range(4)
    ]
    moments = [value / raw[0] for value in raw]
    kappa3 = (
        moments[3]
        - 3 * moments[2] * moments[1]
        + 2 * moments[1] ** 3
    )
    psi2 = mp.polygamma(2, 1)
    return {
        "status": "high_precision_numeric_not_interval_certified",
        "decimalPrecision": mp.mp.dps,
        "density": "exp(-V(u)), V(u)=2u for 0<=u<=1 and V(u)=3u-1 for u>=1",
        "decreasing": True,
        "logConcave": True,
        "potentialConvex": True,
        "x": 0,
        "mass": text(raw[0], 60),
        "momentKappa3": text(kappa3, 60),
        "gammaPsi2": text(psi2, 60),
        "margin": text(kappa3 - psi2, 60),
        "normalizedMellinLogThirdDerivative": text(8 * (kappa3 - psi2), 60),
        "outcome": "violates_the_gamma_skew_comparison",
    }


def saddle_budget():
    mp.mp.dps = 60

    def K(u):
        b = 2 * mp.pi * mp.exp(4 * u)
        return u * (b - mp.mpf("4.5") - 6 / (b - 3))

    def coefficient(u):
        kp = mp.diff(K, u)
        kpp = mp.diff(K, u, 2)
        k = K(u)
        return 2 * k**2 * (kpp / (u * kp**3) + 1 / (u**2 * kp**2))

    table = []
    guesses = {
        16: mp.mpf("0.46"),
        100: mp.mpf("0.77"),
        1000: mp.mpf("1.22"),
        10**6: mp.mpf("2.74"),
        10**9: mp.mpf("4.35"),
    }
    for k, guess in guesses.items():
        u = mp.findroot(lambda value: K(value) - k, guess)
        table.append(
            {
                "k": k,
                "saddleU": text(u),
                "minusK2SaddleActionThird": text(coefficient(u)),
            }
        )

    u0 = mp.mpf("0.98")
    crude_main_bound = (
        80000 * (3202 * u0**2 + 2402 * u0 + 197) / (800 * u0 + 197) ** 3
    )
    k0 = mp.mpf(10) ** 9
    gamma_lower = 8 * k0**2 / (2 * k0 + 1) ** 2
    log_curvature_cauchy = (
        46875 * (mp.log(20 * k0) + mp.pi / 2) / k0
    )
    relative_error_cauchy = mp.mpf("1781.25") / k0
    unallocated = (
        gamma_lower
        - crude_main_bound
        - log_curvature_cauchy
        - relative_error_cauchy
    )
    return {
        "status": "analytic_budget_not_a_completed_proof",
        "dominantPhase": "Psi_z(u)=2z log(u)+5u+log(2pi exp(4u)-3)-pi exp(4u)",
        "saddleMap": "K(u)=u(2pi exp(4u)-9/2-6/(2pi exp(4u)-3))",
        "identity": "S0'''=-2(K''/(u K'^3)+1/(u^2 K'^2))",
        "table": table,
        "crudeMainBoundForUAtLeastPoint98": text(crude_main_bound),
        "crudeBoundDerivativeNumerator": "-160000(1280800u^2+1290806u-197)",
        "tailReferenceK": int(k0),
        "gammaLowerCoefficient": text(gamma_lower),
        "logCurvatureCauchyCoefficient": text(log_curvature_cauchy),
        "relativeErrorCauchyCoefficient": text(relative_error_cauchy),
        "unallocatedCoefficientMargin": text(unallocated),
        "missingProofObligation": "directed bound transferring the dominant-saddle estimate to the full theta saddle, including derivatives of log(1+R)",
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    result = {
        "schema": "riemann-lab.xi-rank-two-global-probe.v1",
        "target": "kappa_3,x(log u) > psi_2(2x+1), sufficient for global rank-two exchange monotonicity",
        "xiDenseScan": xi_dense_scan(),
        "decreasingLogConcaveControl": decreasing_log_concave_control(),
        "saddleBudget": saddle_budget(),
        "outcome": "xi_candidate_survives_but_generic_log_concavity_routes_are_closed",
        "claimBoundaries": [
            "The decreasing log-concave control is a numerical counterexample to a proposed general proof route, not to the xi inequality.",
            "The xi scan is finite, double precision, and non-interval.",
            "The saddle calculation leaves an explicit full-theta perturbation bound unproved and is therefore a proof skeleton, not an internal theorem.",
            "No global rank-two theorem or Riemann-Hypothesis result is claimed.",
        ],
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
