#!/usr/bin/env python3
"""Exact and high-precision controls for the uniform-rank Toda lane."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path

import mpmath as mp


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_moment_mixture_control():
    points = [4, 5, 9, 11, 14]
    weights = [8, 9, 3, 1, 3]
    max_rank, max_shift = 5, 8
    size = max_rank + max_shift + 3
    sequence = [
        sum(Fraction(weight) * Fraction(point) ** k for point, weight in zip(points, weights))
        / factorial(2 * k)
        for k in range(size + 1)
    ]
    determinants = {(0, k): Fraction(1) for k in range(size + 1)}
    determinants.update({(1, k): sequence[k] for k in range(size + 1)})
    for rank in range(2, max_rank + 1):
        for shift in range(size - rank + 2):
            left = Fraction(0) if shift == 0 else determinants[rank - 1, shift - 1]
            determinants[rank, shift] = (
                determinants[rank - 1, shift] ** 2
                - left * determinants[rank - 1, shift + 1]
            ) / determinants[rank - 2, shift]

    def exchange(rank, shift):
        return (
            determinants[rank, shift]
            * determinants[rank - 2, shift]
            / (determinants[rank - 1, shift - 1] * determinants[rank - 1, shift + 1])
        )

    left, right = exchange(4, 3), exchange(4, 4)
    assert left > 0 and right > left
    return {
        "support": points,
        "weights": weights,
        "sequence": "a_k=sum_j w_j*t_j^k/(2k)!",
        "allInvolvedMinorsPositive": all(
            determinants[rank, shift] > 0
            for rank in range(2, 5)
            for shift in range(max(0, rank - 2), 6)
        ),
        "rank": 4,
        "leftShift": 3,
        "rightShift": 4,
        "leftExchange": fraction_text(left),
        "rightExchange": fraction_text(right),
        "relativeIncrease": fraction_text((right - left) / left),
        "conclusion": "positive moment structure plus even-factorial normalization does not imply exchange monotonicity",
    }


def generalized_gamma_precision_control(dps: int):
    mp.mp.dps = dps
    amplitude = mp.mpf(199) / 100
    rank, shift = 7, 255

    def log_a(index):
        return -mp.loggamma(amplitude * index + 1)

    slope = (log_a(shift + 1) - log_a(shift - 1)) / 2

    def normalized_coefficient(index):
        return mp.exp(log_a(index) - log_a(shift) - (index - shift) * slope)

    def determinant(size, center):
        return mp.det(
            mp.matrix(
                [
                    [
                        normalized_coefficient(center + column - row)
                        if center + column - row >= 0
                        else 0
                        for column in range(size)
                    ]
                    for row in range(size)
                ]
            )
        )

    previous_left = determinant(rank - 1, shift - 1)
    previous_center = determinant(rank - 1, shift)
    previous_right = determinant(rank - 1, shift + 1)
    condensation_ratio = previous_left * previous_right / previous_center**2
    value = determinant(rank, shift)
    assert condensation_ratio < 1 and value > 0
    return {
        "sequence": "a_k=1/Gamma(1.99*k+1), normalized by invariant c*lambda^k scaling",
        "rank": rank,
        "shift": shift,
        "decimalPrecision": dps,
        "condensationRatio": mp.nstr(condensation_ratio, 50),
        "determinantSign": "positive",
        "conclusion": "unscaled large-shift determinants can report a false negative; invariant local normalization is mandatory",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=600)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "schema": "riemann-lab.xi-uniform-rank-controls.v1",
        "status": "exact_rational_and_high_precision_controls",
        "positiveMomentMixture": exact_moment_mixture_control(),
        "nearbyGammaPrecisionControl": generalized_gamma_precision_control(args.dps),
        "smoothLogConcaveBoundary": {
            "kernel": "exp(-u^4)",
            "facts": [
                "strictly log-concave",
                "after s(t)=exp(-t^2), the second Laguerre expression is 2*exp(-2t^2) and has concave logarithm",
                "its Fourier transform has nonreal zeros, so the associated coefficient sequence is not PF-infinity",
            ],
            "conclusion": "log-concavity plus the known second-level kernel concavity cannot imply the all-rank exchange principle",
            "boundary": "The nonreal-zero fact is source-backed reasoning, not certified by this script.",
        },
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
