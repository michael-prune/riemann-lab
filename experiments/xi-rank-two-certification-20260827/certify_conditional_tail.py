#!/usr/bin/env python3
"""Certify the scalar closure of a conditional rank-two xi tail bound.

This does not replay Michałowski's missing ancillary certificates.  It takes
the explicitly stated tube, saddle, and factorization bounds in arXiv:2607.16795
as hypotheses and checks, with Arb balls, that they imply f'''(x) > 0 for every
real x >= 10^9.  The analytic reductions and monotonicity arguments are
documented in README.md.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import arb, ctx


ctx.prec = 256


def interval_text(value: arb, digits: int = 30) -> str:
    return value.str(digits, radius=True)


def certify() -> dict:
    # Source-stated hypotheses.
    theta_tail = arb("1.1e-14")
    cauchy_radius = arb("0.04")
    saddle_floor = arb("0.98")
    tail_start = arb("1e9")

    # |log(1+R)| <= -log(1-|R|) on the zero-free tube.
    log_remainder = -(1 - theta_tail).log()

    pi = arb.pi()
    b = 2 * pi * (4 * saddle_floor).exp()
    h = b - arb(9) / 2 - arb(6) / (b - 3)

    # K_phi = K - u g'/2, g = log(1+R).  Cauchy and K'>4K,
    # K'+uK''>4K give uniform relative errors.  Every expression below is
    # decreasing for u>=0.98 because H'/H>4 and 1/u<4.
    alpha = log_remainder / (2 * cauchy_radius * h)
    beta = (
        log_remainder / (8 * cauchy_radius * saddle_floor * h)
        + log_remainder / (4 * cauchy_radius**2 * h)
    )
    gamma = (
        log_remainder / (8 * cauchy_radius * saddle_floor * h)
        + 3 * log_remainder / (4 * cauchy_radius**2 * h)
        + 3 * saddle_floor * log_remainder / (4 * cauchy_radius**3 * h)
    )
    transfer_factor = (1 + alpha) ** 2 * (1 + gamma) / (1 - beta) ** 3

    # Existing elementary dominant-saddle majorant at u=0.98.  Its exact
    # derivative numerator is negative on u>=0.98, so this is the endpoint
    # maximum.
    dominant_bound = (
        arb(80000)
        * (arb(3202) * saddle_floor**2 + arb(2402) * saddle_floor + 197)
        / (arb(800) * saddle_floor + 197) ** 3
    )
    full_saddle_bound = dominant_bound * transfer_factor

    # To cover every real x, center the source disk at the nearest integer n.
    # For n>=1e9, |x/n| <= 1+1/(2n), while the radius n/25 remains inside
    # |z-n|<=n/20.  Cauchy bounds the two analytic correction terms.
    center_ratio = 1 + 1 / (2 * tail_start)
    gamma_lower = 8 * tail_start**2 / (2 * tail_start + 1) ** 2
    curvature = (
        center_ratio**2
        * arb(46875)
        * ((20 * tail_start).log() + pi / 2)
        / tail_start
    )
    relative_error = center_ratio**2 * arb("1781.25") / tail_start
    margin = gamma_lower - full_saddle_bound - curvature - relative_error

    assert beta < 1
    assert transfer_factor < arb("1.000000001")
    assert full_saddle_bound < arb("0.477")
    assert margin > arb("1.52")

    return {
        "schema": "riemann-lab.xi-rank-two-conditional-tail-certificate.v1",
        "status": "checked_scalar_closure_conditional_on_unreplayed_source_lemmas",
        "conclusion": "f'''(x)>0 for every real x>=1e9, conditional on the source-stated analytic bounds",
        "discreteConsequence": "Q_(2,k+1)<Q_(2,k) for every integer k>=1000000001, under the same hypotheses",
        "sourceHypotheses": {
            "thetaRemainderSupremum": interval_text(theta_tail),
            "fullSaddleAnalyticity": "arXiv:2607.16795, Lemma 3.3",
            "uniformFactorization": "arXiv:2607.16795, Lemma 3.6",
            "factorizationRelativeError": "0.018",
            "curvatureSector": "6/sigma_n^2 <= |-Psi''| <= 20/sigma_n^2 in the right half-plane",
        },
        "arbPrecisionBits": ctx.prec,
        "constants": {
            "logRemainderBound": interval_text(log_remainder),
            "HAtPoint98": interval_text(h),
            "relativeKErrorAlpha": interval_text(alpha),
            "relativeKPrimeErrorBeta": interval_text(beta),
            "relativeNumeratorErrorGamma": interval_text(gamma),
            "fullThetaTransferFactor": interval_text(transfer_factor),
            "dominantSaddleCoefficientUpper": interval_text(dominant_bound),
            "fullSaddleCoefficientUpper": interval_text(full_saddle_bound),
            "gammaCoefficientLowerAt1e9": interval_text(gamma_lower),
            "curvatureCauchyCoefficientUpperAt1e9": interval_text(curvature),
            "relativeErrorCauchyCoefficientUpperAt1e9": interval_text(relative_error),
            "finalCoefficientMarginLowerAt1e9": interval_text(margin),
        },
        "nonclaims": [
            "The arXiv source archive did not contain the advertised certificate modules, so their bounds were not replayed.",
            "The compact range below the tail threshold is not certified.",
            "No global rank-two theorem, rank-three theorem, RH implication, or novelty claim is made.",
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
