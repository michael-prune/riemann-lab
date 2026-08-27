#!/usr/bin/env python3
"""Focused numerical probe of xi Toeplitz exchange monotonicity.

This produces high-precision discovery evidence, not interval certificates.
The asymptotic statements documented beside it are analytic consequences of
published coefficient expansions and are not certified by this script.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import random
from functools import lru_cache
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[2]
SCAN_PATH = ROOT / "experiments" / "positive-grassmannian" / "scan.py"


def load_scan_module():
    spec = importlib.util.spec_from_file_location("positive_geometry_scan", SCAN_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {SCAN_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def decimal(value: mp.mpf, digits: int = 24) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def finite_stress(scan, max_rank: int, max_shift: int, max_coefficient: int, theta_terms: int):
    coefficients = scan.xi_coefficients(max_coefficient, theta_terms)
    sequence, _ = scan.normalized_coefficients(coefficients)

    @lru_cache(None)
    def determinant(rank: int, shift: int):
        return scan.determinant(sequence, rank, shift)

    records = []
    for rank in range(2, max_rank + 1):
        values = []
        for shift in range(max(1, rank - 2), max_shift + 1):
            denominator = determinant(rank - 1, shift - 1) * determinant(rank - 1, shift + 1)
            coordinate = determinant(rank, shift) * determinant(rank - 2, shift) / denominator
            values.append((shift, coordinate))
        relative_decreases = [
            (left - right) / max(abs(left), abs(right))
            for (_, left), (_, right) in zip(values, values[1:])
        ]
        records.append(
            {
                "rank": rank,
                "shiftRange": [values[0][0], values[-1][0]],
                "firstRatio": decimal(values[0][1]),
                "lastRatio": decimal(values[-1][1]),
                "allPositive": all(value > 0 for _, value in values),
                "strictlyDecreasing": all(value > 0 for value in relative_decreases),
                "minimumRelativeDecrease": decimal(min(relative_decreases), 16),
            }
        )
    return records


def tilted_log_moments(scan, x: mp.mpf, theta_terms: int):
    raw = [
        mp.quad(
            lambda u, power=power: (
                u ** (2 * x)
                * scan.classical_phi(u, theta_terms)
                * mp.log(u) ** power
            ),
            scan.INTEGRATION_BREAKS,
        )
        for power in range(4)
    ]
    moments = [value / raw[0] for value in raw]
    kappa_1 = moments[1]
    kappa_2 = moments[2] - moments[1] ** 2
    kappa_3 = moments[3] - 3 * moments[2] * moments[1] + 2 * moments[1] ** 3
    return {
        "logA": mp.log(raw[0]) - mp.loggamma(2 * x + 1),
        "f1": 2 * kappa_1 - 2 * mp.digamma(2 * x + 1),
        "f2": 4 * kappa_2 - 4 * mp.polygamma(1, 2 * x + 1),
        "f3": 8 * kappa_3 - 8 * mp.polygamma(2, 2 * x + 1),
        "momentKappa3": kappa_3,
        "gammaPsi2": mp.polygamma(2, 2 * x + 1),
    }


def continuous_curvature(scan, samples: list[int], theta_terms: int):
    cache = {}

    def stats(x: int):
        if x not in cache:
            cache[x] = tilted_log_moments(scan, mp.mpf(x), theta_terms)
        return cache[x]

    output = []
    for x in samples:
        left, center, right = stats(x - 1), stats(x), stats(x + 1)
        h = left["logA"] + right["logA"] - 2 * center["logA"]
        ratio = mp.exp(h)
        h1 = left["f1"] + right["f1"] - 2 * center["f1"]
        h2 = left["f2"] + right["f2"] - 2 * center["f2"]
        h3 = left["f3"] + right["f3"] - 2 * center["f3"]
        transform_1 = -ratio / (1 - ratio)
        transform_2 = -ratio / (1 - ratio) ** 2
        transform_3 = -ratio * (1 + ratio) / (1 - ratio) ** 3
        log_d2_third = (
            2 * center["f3"]
            + transform_3 * h1**3
            + 3 * transform_2 * h1 * h2
            + transform_1 * h3
        )
        output.append(
            {
                "x": x,
                "logAThirdDerivative": decimal(center["f3"]),
                "logD2ThirdDerivative": decimal(log_d2_third),
                "centeredSecondDifferenceLogA": decimal(h),
                "momentLogSkewness": decimal(center["momentKappa3"]),
                "gammaPolygamma2": decimal(center["gammaPsi2"]),
            }
        )
    return output


def polynomial_from_seeded_negative_roots(degree: int = 70):
    generator = random.Random(0)
    coefficients = [mp.mpf("1")]
    for _ in range(degree):
        root = mp.mpf(10) ** mp.mpf(generator.uniform(-2, 2))
        next_coefficients = [mp.mpf("0")] * (len(coefficients) + 1)
        for index, value in enumerate(coefficients):
            next_coefficients[index] += value
            next_coefficients[index + 1] += value / root
        coefficients = next_coefficients
    return coefficients


def control_violations(scan):
    sequence = polynomial_from_seeded_negative_roots()

    @lru_cache(None)
    def determinant(rank: int, shift: int):
        return scan.determinant(sequence, rank, shift)

    output = []
    for rank in range(2, 9):
        values = []
        for shift in range(max(1, rank - 2), 45):
            denominator = determinant(rank - 1, shift - 1) * determinant(rank - 1, shift + 1)
            values.append(
                (shift, determinant(rank, shift) * determinant(rank - 2, shift) / denominator)
            )
        for (left_shift, left), (right_shift, right) in zip(values, values[1:]):
            if right >= left:
                output.append(
                    {
                        "rank": rank,
                        "leftShift": left_shift,
                        "rightShift": right_shift,
                        "leftRatio": decimal(left),
                        "rightRatio": decimal(right),
                    }
                )
                break
    return output


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=130)
    parser.add_argument("--theta-terms", type=int, default=18)
    parser.add_argument("--max-rank", type=int, default=15)
    parser.add_argument("--max-shift", type=int, default=40)
    parser.add_argument("--max-coefficient", type=int, default=62)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    needed = args.max_shift + args.max_rank - 1
    if args.max_coefficient < needed:
        raise SystemExit(f"--max-coefficient must be at least {needed}")
    mp.mp.dps = args.dps
    scan = load_scan_module()
    result = {
        "schema": "riemann-lab.xi-exchange-monotonicity-probe.v1",
        "status": "high_precision_numeric_not_interval_certified",
        "candidate": "Q_(r,k+1) < Q_(r,k), where Q_(r,k)=D_(r,k)D_(r-2,k)/(D_(r-1,k-1)D_(r-1,k+1))",
        "parameters": {
            "decimalPrecision": args.dps,
            "thetaTerms": args.theta_terms,
            "maxCoefficient": args.max_coefficient,
            "maxRank": args.max_rank,
            "maxShift": args.max_shift,
        },
        "finiteStress": finite_stress(
            scan, args.max_rank, args.max_shift, args.max_coefficient, args.theta_terms
        ),
        "continuousCurvature": continuous_curvature(
            scan, [1, 2, 3, 5, 10, 20, 40, 80], args.theta_terms
        ),
        "pfInfinityControl": {
            "construction": "coefficients of product_(j=1)^70 (1+z/root_j), root_j=10^U_j, U_j uniform[-2,2] with seed 0",
            "allZerosNegativeReal": True,
            "firstMonotonicityViolations": control_violations(scan),
        },
        "claimBoundaries": [
            "The finite and continuous calculations are non-interval numerical evidence.",
            "Published all-orders coefficient asymptotics imply the rank-2 and rank-3 inequalities only for sufficiently large shift, with no effective threshold extracted here.",
            "No checked classical Turan, double-Turan, or PF-infinity theorem implies the candidate globally.",
            "Even a proof for every fixed rank in the tail would not resolve the RH-critical diagonal regime k comparable to r.",
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
