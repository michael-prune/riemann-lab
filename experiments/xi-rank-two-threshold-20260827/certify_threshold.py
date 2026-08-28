#!/usr/bin/env python3
"""Directed scalar certificate for lowering the conditional xi rank-two tail.

The certificate rebuilds the scalar part of the one-saddle factorization at
k >= 10^5, replaces the old curvature-log Cauchy estimate by a bound on its
first derivative, and closes the f''' sign budget.  The contour and analytic
continuation lemmas from arXiv:2607.16795 remain source-backed hypotheses
because the advertised ancillary certificate modules are unavailable.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import arb, ctx


ctx.prec = 256


def ball(value: arb, digits: int = 30) -> str:
    return value.str(digits, radius=True)


def certify() -> dict:
    pi = arb.pi()
    threshold = arb("1e5")
    theta_tail = arb("1.1e-14")
    saddle_floor = arb("0.98")
    u_drift = arb("0.015")
    z_radius = arb("0.04")
    contour_delta = arb(1) / 40
    window_y = arb(1)

    def dominant_k(u: arb) -> arb:
        b = 2 * pi * (4 * u).exp()
        return u * (b - arb(9) / 2 - arb(6) / (b - 3))

    # The source used k>=10^9 only to obtain u_k>1 in the scalar phase
    # estimates.  K is increasing and K(1)<10^5, so the same prerequisite
    # already holds at the new threshold.
    k_at_one = dominant_k(arb(1))
    assert k_at_one < threshold

    log_k = threshold.log()
    sigma_upper = (log_k / (4 * threshold)).sqrt()
    assert sigma_upper < contour_delta

    # Rebuilt one-saddle factorization with |y|<=1 instead of |y|<=2.
    # The normalized cubic remainder is bounded by
    # (100/6)*sqrt(u_k/k)*|y|^3 and u_k<log(k)/4.
    cubic_error = arb(100) / 6 * sigma_upper * window_y**3
    gaussian_lower = (2 * pi / arb("14.185")).sqrt()
    local_replacement = (cubic_error.exp() - 1) * (pi / 3).sqrt()
    gaussian_tail = (-3 * window_y**2).exp() / (3 * window_y)
    true_local_tail = (
        2 * (-arb("2.5") * window_y**2).exp() / (5 * window_y)
    )

    # The global horizontal envelope has b=(4/5)k/u_k^2.  Substituting
    # u_k<log(k)/4 gives the directed lower bound below.
    horizontal_b = arb("12.8") * threshold / log_k**2
    horizontal_tail = (
        threshold.sqrt()
        / (horizontal_b * contour_delta)
        * (-horizontal_b * contour_delta**2).exp()
    )
    connector = (
        arb("0.36")
        * 24
        / arb("0.665")
        * threshold.sqrt()
        * (-threshold).exp()
    )
    factorization_error = (
        local_replacement
        + gaussian_tail
        + true_local_tail
        + horizontal_tail
        + connector
    ) / gaussian_lower
    assert factorization_error < arb("0.223")
    factorization_log = -(1 - factorization_error).log()

    # Full-theta saddle-action transfer, inherited from the preceding
    # certificate and recomputed here.
    log_theta_remainder = -(1 - theta_tail).log()
    b_floor = 2 * pi * (4 * saddle_floor).exp()
    h_floor = b_floor - arb(9) / 2 - arb(6) / (b_floor - 3)
    alpha = log_theta_remainder / (2 * z_radius * h_floor)
    beta = (
        log_theta_remainder / (8 * z_radius * saddle_floor * h_floor)
        + log_theta_remainder / (4 * z_radius**2 * h_floor)
    )
    gamma = (
        log_theta_remainder / (8 * z_radius * saddle_floor * h_floor)
        + 3 * log_theta_remainder / (4 * z_radius**2 * h_floor)
        + 3
        * saddle_floor
        * log_theta_remainder
        / (4 * z_radius**3 * h_floor)
    )
    transfer_factor = (1 + alpha) ** 2 * (1 + gamma) / (1 - beta) ** 3
    dominant_saddle_bound = (
        arb(80000)
        * (arb(3202) * saddle_floor**2 + arb(2402) * saddle_floor + 197)
        / (arb(800) * saddle_floor + 197) ** 3
    )
    full_saddle_bound = dominant_saddle_bound * transfer_factor
    assert full_saddle_bound < arb("0.477")

    # Exact identity at the full saddle:
    #   -Psi'' = 2 K_phi'/u,
    #   d/dz log(-Psi'') = K_phi''/K_phi'^2 - 1/(u K_phi').
    # Source phase bounds give |K_phi'|>=2.955 n and
    # |K_phi''|<=70.75 n throughout the complex disk.
    k_prime_lower = 3 * (1 - u_drift)
    k_second_upper = 20 + 50 * (1 + u_drift)
    curvature_first_derivative_coefficient = (
        k_second_upper / k_prime_lower**2
        + 1 / ((1 - u_drift) * k_prime_lower)
    )
    assert curvature_first_derivative_coefficient < arb("8.45")

    # Center at the nearest integer n.  The radius .04n about x stays in the
    # source disk |z-n|<=.05n for n>=10^5.
    center_ratio = 1 + 1 / (2 * threshold)
    relative_error_cauchy = (
        center_ratio**2
        * 6
        * factorization_log
        / (z_radius**3 * threshold)
    )
    curvature_cauchy = (
        center_ratio**2
        * arb("8.45")
        / (z_radius**2 * threshold)
    )
    gamma_lower = 8 * threshold**2 / (2 * threshold + 1) ** 2
    final_margin = (
        gamma_lower
        - arb("0.477")
        - relative_error_cauchy
        - curvature_cauchy
    )

    assert relative_error_cauchy < arb("0.237")
    assert curvature_cauchy < arb("0.053")
    assert final_margin > arb("1.23")

    return {
        "schema": "riemann-lab.xi-rank-two-threshold-certificate.v1",
        "status": "checked_threshold_reduction_conditional_on_unreplayed_source_analytic_lemmas",
        "conditionalConclusion": "f'''(x)>0 for every real x>=100000",
        "conditionalDiscreteConsequence": "Q_(2,k+1)<Q_(2,k) for every integer k>=100001",
        "previousThreshold": 1_000_000_000,
        "newThreshold": 100_000,
        "improvementFactor": 10_000,
        "arbPrecisionBits": ctx.prec,
        "sourceHypothesesStillUnreplayed": [
            "kernel tube and full-saddle analytic continuation",
            "horizontal concavity and contour deformation",
            "uniform phase derivative bounds; their new-threshold scalar prerequisites are checked here",
        ],
        "factorization": {
            "windowY": ball(window_y),
            "KAtOne": ball(k_at_one),
            "sigmaUpperAtThreshold": ball(sigma_upper),
            "cubicError": ball(cubic_error),
            "gaussianLower": ball(gaussian_lower),
            "localReplacement": ball(local_replacement),
            "gaussianTail": ball(gaussian_tail),
            "trueLocalTail": ball(true_local_tail),
            "horizontalBLower": ball(horizontal_b),
            "horizontalTail": ball(horizontal_tail),
            "connector": ball(connector),
            "relativeError": ball(factorization_error),
            "logRelativeError": ball(factorization_log),
        },
        "signBudgetAtThreshold": {
            "gammaLower": ball(gamma_lower),
            "fullSaddleUpper": ball(full_saddle_bound),
            "relativeErrorCauchyUpper": ball(relative_error_cauchy),
            "curvatureFirstDerivativeCoefficient": ball(
                curvature_first_derivative_coefficient
            ),
            "curvatureCauchyUpper": ball(curvature_cauchy),
            "finalMarginLower": ball(final_margin),
        },
        "nonclaims": [
            "The source analytic lemmas remain unreplayed because the advertised ancillary modules are absent from the arXiv source archive.",
            "The compact range 0<=x<100000 is not interval-certified.",
            "No unconditional global rank-two theorem, rank-three theorem, RH implication, or novelty claim is made.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(certify(), indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
