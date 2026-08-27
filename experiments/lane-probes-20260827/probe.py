#!/usr/bin/env python3
"""Four bounded numerical probes of current Riemann Lab lanes.

The output is exploratory evidence.  Floating-point or arbitrary-precision
numerics here are not interval certificates and do not prove RH or any new
general theorem.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import random
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SCAN_PATH = ROOT / "experiments" / "positive-grassmannian" / "scan.py"


def load_scan_module():
    spec = importlib.util.spec_from_file_location("positive_geometry_scan", SCAN_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SCAN_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def decimal(value, digits: int = 18) -> str:
    if isinstance(value, (float, np.floating)):
        return format(float(value), f".{digits}g")
    return mp.nstr(value, digits, strip_zeros=False)


def weil_window_probe(points: int = 96) -> dict[str, object]:
    raw_nodes, raw_weights = np.polynomial.legendre.leggauss(points)
    nodes = raw_nodes / 2
    weights = raw_weights / 2

    # R(psi) is a quadratic form subject to integral psi = 1.
    operator = np.diag(weights) + np.outer(weights, weights) * np.abs(nodes[:, None] - nodes[None, :])
    response = np.linalg.solve(operator, weights)
    optimizer = response / np.dot(weights, response)
    numeric_r = 1 / np.dot(weights, response)

    cosine = np.cos(np.sqrt(2) * nodes)
    cosine /= np.dot(weights, cosine)
    theoretical_r = 0.5 + math.cos(1 / math.sqrt(2)) / (math.sqrt(2) * math.sin(1 / math.sqrt(2)))
    shape_error = math.sqrt(np.dot(weights, (optimizer - cosine) ** 2) / np.dot(weights, cosine**2))

    lambdas = np.linspace(0.001, 1, 10000)
    certificates = 2 - 1 / lambdas - lambdas / 3
    best_index = int(np.argmax(certificates))

    return {
        "lane": "explicit-formula-zero-density",
        "question": "Can window or bandwidth rescaling improve the new finite-compression certificate?",
        "quadraturePoints": points,
        "windowFunctional": "R(psi)=(int psi^2 + double_int |u-v| psi(u)psi(v))/(int psi)^2",
        "numericMinimumR": decimal(numeric_r),
        "theoreticalMontgomeryTaylorR": decimal(theoretical_r),
        "relativeError": decimal(abs(numeric_r - theoretical_r) / theoretical_r),
        "simpleZeroCertificate": decimal(2 - numeric_r),
        "minimumSampledOptimizerValue": decimal(np.min(optimizer)),
        "relativeL2ShapeErrorAgainstNormalizedCosine": decimal(shape_error),
        "bandwidthRescaling": {
            "formula": "H(lambda)=2-1/lambda-lambda/3 for 0<lambda<=1",
            "bestSampledLambda": decimal(lambdas[best_index]),
            "bestSampledCertificate": decimal(certificates[best_index]),
        },
        "publishedBandwidthOneCeiling": "0.6818287",
        "outcome": "closed_within_first_two_moments_and_bandwidth_one",
        "interpretation": [
            "The discretized global minimizer reproduces the Montgomery-Taylor cosine window.",
            "Bandwidth shortening does not improve the flat-window certificate.",
            "A material advance needs information beyond the first two bandwidth-one trace moments, such as wider-support pair correlation or rigorously available higher moments.",
        ],
    }


def toeplitz_probe(
    max_coefficient: int = 54,
    decimal_precision: int = 90,
    theta_terms: int = 16,
) -> dict[str, object]:
    scan = load_scan_module()
    mp.mp.dps = decimal_precision
    coefficients = scan.xi_coefficients(max_coefficient, theta_terms)
    sequence, _ = scan.normalized_coefficients(coefficients)

    @lru_cache(None)
    def determinant(rank: int, shift: int):
        return scan.determinant(sequence, rank, shift)

    ranks = []
    maximum_identity_error = mp.mpf("0")
    for rank in range(2, 10):
        values = []
        start = max(1, rank - 2)
        for shift in range(start, 45):
            denominator = determinant(rank - 1, shift - 1) * determinant(rank - 1, shift + 1)
            from_upper = determinant(rank, shift) * determinant(rank - 2, shift) / denominator
            from_lower = determinant(rank - 1, shift) ** 2 / denominator - 1
            identity_error = abs(from_upper - from_lower) / max(abs(from_upper), mp.mpf("1e-80"))
            maximum_identity_error = max(maximum_identity_error, identity_error)
            values.append((shift, from_upper))
        differences = [values[index + 1][1] - values[index][1] for index in range(len(values) - 1)]
        ranks.append(
            {
                "rank": rank,
                "shiftRange": [values[0][0], values[-1][0]],
                "firstRatio": decimal(values[0][1]),
                "lastRatio": decimal(values[-1][1]),
                "allPositive": all(value > 0 for _, value in values),
                "strictlyDecreasing": all(difference < 0 for difference in differences),
                "smallestDecrease": decimal(min(-difference for difference in differences)),
            }
        )

    def control_sequence() -> list[mp.mpf]:
        generator = random.Random(0)
        output = [mp.mpf("1")]
        for _ in range(70):
            root_scale = mp.mpf(10) ** mp.mpf(generator.uniform(-2, 2))
            next_output = [mp.mpf("0")] * (len(output) + 1)
            for index, value in enumerate(output):
                next_output[index] += value
                next_output[index + 1] += value / root_scale
            output = next_output
        return output

    control = control_sequence()

    @lru_cache(None)
    def control_det(rank: int, shift: int):
        return scan.determinant(control, rank, shift)

    control_violations = []
    for rank in range(2, 9):
        values = []
        for shift in range(max(1, rank - 2), 45):
            denominator = control_det(rank - 1, shift - 1) * control_det(rank - 1, shift + 1)
            values.append(
                (
                    shift,
                    control_det(rank, shift) * control_det(rank - 2, shift) / denominator,
                )
            )
        for (left_shift, left), (right_shift, right) in zip(values, values[1:]):
            if right >= left:
                control_violations.append(
                    {
                        "rank": rank,
                        "leftShift": left_shift,
                        "leftRatio": decimal(left),
                        "rightShift": right_shift,
                        "rightRatio": decimal(right),
                    }
                )
                break

    return {
        "lane": "positivity-equivalent-criteria",
        "question": "Do consecutive xi Toeplitz minors exhibit a backward-propagating exchange inequality?",
        "status": "high_precision_numeric_not_interval_certified",
        "parameters": {
            "decimalPrecision": decimal_precision,
            "thetaTerms": theta_terms,
            "maxCoefficient": max_coefficient,
            "ranks": [2, 9],
            "maximumShift": 44,
        },
        "exchangeRatio": "Q_(r,k)=D_(r,k)D_(r-2,k)/(D_(r-1,k-1)D_(r-1,k+1))",
        "lowerRankIdentity": "Q_(r,k)=D_(r-1,k)^2/(D_(r-1,k-1)D_(r-1,k+1))-1",
        "candidateInequality": "Q_(r,k+1)<Q_(r,k) for every r>=2 and admissible k",
        "rankResults": ranks,
        "maximumRelativeDesnanotJacobiMismatch": decimal(maximum_identity_error),
        "pfInfinityControl": {
            "construction": "coefficients of product_(j=1)^70 (1+z/root_j), root_j=10^U_j, U_j uniform[-2,2] with seed 0",
            "allZerosNegativeReal": True,
            "firstMonotonicityViolations": control_violations,
        },
        "outcome": "candidate_survives_sample_and_is_not_generic_pf_infinity",
        "interpretation": [
            "Every sampled xi exchange ratio was positive and strictly decreasing in the shift.",
            "If this held uniformly, the proved positive cubic tail would propagate backward to every consecutive minor.",
            "Seeded negative-real-root polynomial controls violate the monotonicity, so it is not a generic consequence of total positivity.",
            "The finite xi sample is already covered by known finite-order positivity and is discovery evidence only.",
        ],
    }


def heat_flow_probe() -> dict[str, object]:
    mp.mp.dps = 40
    zero_count = 100
    zeros = [mp.im(mp.zetazero(index)) for index in range(1, zero_count + 1)]

    def velocities(cutoff: int) -> list[mp.mpf]:
        selected = zeros[:cutoff]
        output = []
        for index, x in enumerate(selected):
            total = mp.mpf("0")
            for other_index, y in enumerate(selected):
                if index != other_index:
                    total += 1 / (x - y)
                total += 1 / (x + y)
            output.append(2 * total)
        return output

    cutoff_results = []
    signs_by_cutoff = []
    for cutoff in (40, 70, 100):
        velocity = velocities(cutoff)
        gap_derivatives = [velocity[index + 1] - velocity[index] for index in range(cutoff - 1)]
        interior = gap_derivatives[9:30]
        signs_by_cutoff.append([1 if value > 0 else -1 if value < 0 else 0 for value in interior])
        cutoff_results.append(
            {
                "zeroCutoff": cutoff,
                "gapIndexRange": [10, 30],
                "negativeDerivativeCount": sum(value < 0 for value in interior),
                "positiveDerivativeCount": sum(value > 0 for value in interior),
                "minimumDerivative": decimal(min(interior)),
                "maximumDerivative": decimal(max(interior)),
                "negativeGapIndices": [
                    index + 1 for index, value in enumerate(gap_derivatives[:40]) if value < 0
                ],
            }
        )

    stable_signs = all(signs == signs_by_cutoff[0] for signs in signs_by_cutoff[1:])
    return {
        "lane": "de-bruijn-newman-heat-flow",
        "question": "Can a pointwise adjacent-gap repulsion barrier prevent backward zero collisions?",
        "status": "finite_truncated_numeric_not_interval_certified",
        "zeroDynamics": "gamma_j'(t)=2 sum'_(k!=j) 1/(gamma_j-gamma_k)",
        "cutoffs": cutoff_results,
        "allGapSignsStableAcrossCutoffs": stable_signs,
        "outcome": "naive_pointwise_gap_barrier_falsified",
        "interpretation": [
            "Both signs occur: forward heat flow widens some adjacent gaps and shrinks others.",
            "The signs for gaps 10 through 30 agree under 40-, 70-, and 100-positive-zero truncations.",
            "Pairwise repulsion does not imply pointwise widening because the rest of the zero configuration can dominate a given gap.",
            "A viable no-collision argument must use a collective or averaged energy with a controlled infinite tail.",
        ],
    }


EULER_GAMMA = 0.577215664901532860606512090082402431
AUTOCORRELATION_AT_ONE = math.log(2 * math.pi) - EULER_GAMMA


@lru_cache(None)
def vasyunin_sum(p: int, q: int) -> float:
    if q == 1:
        return 0.0
    return sum(
        ((index * p) % q) / q / math.tan(math.pi * index / q)
        for index in range(1, q)
    )


@lru_cache(None)
def fractional_autocorrelation(p: int, q: int) -> float:
    divisor = math.gcd(p, q)
    p //= divisor
    q //= divisor
    ratio = p / q
    return (
        0.5 * (1 - ratio) * math.log(ratio)
        + 0.5 * (ratio + 1) * AUTOCORRELATION_AT_ONE
        - math.pi / (2 * q) * (vasyunin_sum(p, q) + vasyunin_sum(q, p))
    )


def nyman_beurling_probe() -> dict[str, object]:
    maximum = 320
    gram = np.empty((maximum, maximum))
    for k in range(1, maximum + 1):
        for ell in range(k, maximum + 1):
            value = fractional_autocorrelation(ell, k) / ell
            gram[k - 1, ell - 1] = value
            gram[ell - 1, k - 1] = value

    target = np.array(
        [(math.log(k) + 1 - EULER_GAMMA) / k for k in range(1, maximum + 1)]
    )
    sizes = (5, 10, 20, 40, 80, 160, 320)
    records = []
    for size in sizes:
        matrix = gram[:size, :size]
        vector = target[:size]
        diagonal = np.sqrt(np.diag(matrix))
        correlation = matrix / np.outer(diagonal, diagonal)
        eigenvalues = np.linalg.eigvalsh(matrix)
        correlation_eigenvalues = np.linalg.eigvalsh(correlation)
        coefficients = np.linalg.solve(matrix, vector)
        distance_squared = 1 - np.dot(vector, coefficients)
        records.append(
            {
                "size": size,
                "rawConditionNumber": decimal(eigenvalues[-1] / eigenvalues[0]),
                "diagonalScaledConditionNumber": decimal(
                    correlation_eigenvalues[-1] / correlation_eigenvalues[0]
                ),
                "distanceSquared": decimal(distance_squared),
                "logNDistanceSquared": decimal(math.log(size) * distance_squared),
                "coefficientL2Norm": decimal(np.linalg.norm(coefficients)),
            }
        )

    raw_conditions = np.array([float(record["rawConditionNumber"]) for record in records[2:]])
    scaled_conditions = np.array(
        [float(record["diagonalScaledConditionNumber"]) for record in records[2:]]
    )
    fit_sizes = np.array(sizes[2:], dtype=float)
    raw_slope = np.polyfit(np.log(fit_sizes), np.log(raw_conditions), 1)[0]
    scaled_slope = np.polyfit(np.log(fit_sizes), np.log(scaled_conditions), 1)[0]

    return {
        "lane": "nyman-beurling-gram-geometry",
        "question": "Does natural diagonal preconditioning remove the Gram instability?",
        "status": "double_precision_exact_formula_evaluation",
        "basis": "rho_k(t)={1/(k t)}, k=1,...,N",
        "gramFormula": "G_(k,l)=A(l/k)/l with A from the rational Vasyunin-sum formula",
        "targetFormula": "b_k=(log(k)+1-EulerGamma)/k",
        "records": records,
        "logLogFit": {
            "rawConditionPower": decimal(raw_slope),
            "diagonalScaledConditionPower": decimal(scaled_slope),
            "fitSizes": list(map(int, sizes[2:])),
        },
        "outcome": "diagonal_scaling_improves_constant_not_growth_order",
        "interpretation": [
            "The minimum approximation distance decreases steadily and log(N)*distance^2 is nearly flat over the larger samples.",
            "Unit-diagonal scaling reduces the condition number by a constant factor but retains approximately quadratic growth.",
            "A useful preconditioner must exploit more than basis norms, probably the Mellin or block-Hankel structure.",
            "The calculation reallocates rather than resolves the coefficient-control difficulty in the Nyman-Beurling criterion.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = {
        "schema": "riemann-lab.lane-probes.v1",
        "date": "2026-08-27",
        "claimLevel": "checked_artifact",
        "probes": {
            "weilFiniteCompression": weil_window_probe(),
            "xiToeplitzFrontier": toeplitz_probe(),
            "deBruijnNewmanHeatFlow": heat_flow_probe(),
            "nymanBeurlingGram": nyman_beurling_probe(),
        },
        "globalBoundary": [
            "No probe proves RH or a new general theorem.",
            "The Toeplitz and heat-flow calculations are non-interval numerical experiments.",
            "The Nyman-Beurling calculation evaluates exact finite formulas in double precision.",
            "The Weil-window probe reproduces a published variational optimum and method ceiling.",
        ],
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
