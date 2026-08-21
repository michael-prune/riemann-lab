#!/usr/bin/env python3
"""Numerical probes for the Riemann-xi positive-Grassmannian program.

The calculation is deliberately evidence-producing, not proof-producing:
mpmath does not provide interval enclosures.  Every reported sign is therefore
labelled as a high-precision numerical observation.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import mpmath as mp


INTEGRATION_BREAKS = (0, 0.25, 0.5, 0.75, 1, 1.5, 2)


def classical_phi(u: mp.mpf, theta_terms: int) -> mp.mpf:
    """Classical positive theta kernel in the x/2 normalization."""
    e4 = mp.exp(4 * u)
    e5 = mp.exp(5 * u)
    e9 = mp.exp(9 * u)
    return mp.fsum(
        (
            2 * mp.pi**2 * n**4 * e9
            - 3 * mp.pi * n**2 * e5
        )
        * mp.exp(-mp.pi * n**2 * e4)
        for n in range(1, theta_terms + 1)
    )


def xi(s: mp.mpf | mp.mpc) -> mp.mpf | mp.mpc:
    if s == 0 or s == 1:
        return mp.mpf("0.5")
    return (
        mp.mpf("0.5")
        * s
        * (s - 1)
        * mp.power(mp.pi, -s / 2)
        * mp.gamma(s / 2)
        * mp.zeta(s)
    )


def xi_coefficients(max_index: int, theta_terms: int) -> list[mp.mpf]:
    """Return a_k where xi(1/2 + sqrt(z)/2)/8 = sum a_k z^k."""
    coefficients: list[mp.mpf] = []
    for k in range(max_index + 1):
        moment = mp.quad(
            lambda u: u ** (2 * k) * classical_phi(u, theta_terms),
            INTEGRATION_BREAKS,
        )
        coefficients.append(moment / mp.factorial(2 * k))
    return coefficients


def normalized_coefficients(coefficients: list[mp.mpf]) -> tuple[list[mp.mpf], mp.mpf]:
    """Apply a positive c*lambda^k scaling, which preserves every minor sign."""
    scale = coefficients[0] / coefficients[1]
    return [value / coefficients[0] * scale**k for k, value in enumerate(coefficients)], scale


def coefficient_at(sequence: list[mp.mpf], index: int) -> mp.mpf:
    if index < 0 or index >= len(sequence):
        return mp.mpf("0")
    return sequence[index]


def toeplitz_block(sequence: list[mp.mpf], rank: int, shift: int, columns: tuple[int, ...] | None = None) -> mp.matrix:
    selected = columns if columns is not None else tuple(range(rank))
    return mp.matrix(
        [
            [coefficient_at(sequence, shift + column - row) for column in selected]
            for row in range(rank)
        ]
    )


def determinant(sequence: list[mp.mpf], rank: int, shift: int) -> mp.mpf:
    if rank == 0:
        return mp.mpf("1")
    return mp.det(toeplitz_block(sequence, rank, shift))


def sign(value: mp.mpf) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def decimal(value: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def cancellation_margin(sequence: list[mp.mpf], rank: int, shift: int) -> mp.mpf:
    """det(A) / sum_sigma |term_sigma|, a diagnostic for cancellation."""
    if rank == 0:
        return mp.mpf("1")
    positive_and_negative_mass = mp.mpf("0")
    for permutation in itertools.permutations(range(rank)):
        term = mp.mpf("1")
        for row, column in enumerate(permutation):
            term *= coefficient_at(sequence, shift + column - row)
        positive_and_negative_mass += abs(term)
    value = determinant(sequence, rank, shift)
    return value / positive_and_negative_mass if positive_and_negative_mass else mp.mpf("0")


def consecutive_scan(sequence: list[mp.mpf], max_rank: int, max_shift: int) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for rank in range(1, max_rank + 1):
        for shift in range(0, max_shift + 1):
            value = determinant(sequence, rank, shift)
            records.append(
                {
                    "rank": rank,
                    "shift": shift,
                    "value": decimal(value),
                    "sign": sign(value),
                    "cancellationMargin": decimal(cancellation_margin(sequence, rank, shift), 18),
                }
            )
    return records


def plucker_scan(sequence: list[mp.mpf], rank: int, width: int, offsets: range) -> list[dict[str, object]]:
    scans: list[dict[str, object]] = []
    for offset in offsets:
        values: list[mp.mpf] = []
        negative_columns: list[list[int]] = []
        zero_columns: list[list[int]] = []
        for columns in itertools.combinations(range(width), rank):
            value = mp.det(toeplitz_block(sequence, rank, offset, columns))
            values.append(value)
            if value < 0:
                negative_columns.append(list(columns))
            elif value == 0:
                zero_columns.append(list(columns))
        scans.append(
            {
                "rank": rank,
                "width": width,
                "offset": offset,
                "pluckerCoordinateCount": len(values),
                "negativeCount": len(negative_columns),
                "zeroCount": len(zero_columns),
                "minimum": decimal(min(values)),
                "maximum": decimal(max(values)),
                "negativeColumnSets": negative_columns,
                "zeroColumnSets": zero_columns,
            }
        )
    return scans


def exchange_coordinates(sequence: list[mp.mpf], max_rank: int, max_shift: int) -> list[dict[str, object]]:
    """Positive coordinates from the Desnanot-Jacobi/Plucker exchange relation."""
    records: list[dict[str, object]] = []
    for rank in range(2, max_rank + 1):
        for shift in range(1, max_shift):
            numerator = determinant(sequence, rank, shift) * determinant(sequence, rank - 2, shift)
            denominator = determinant(sequence, rank - 1, shift - 1) * determinant(sequence, rank - 1, shift + 1)
            coordinate = numerator / denominator
            left = determinant(sequence, rank - 1, shift) ** 2
            residual = left - denominator - numerator
            records.append(
                {
                    "rank": rank,
                    "shift": shift,
                    "coordinate": decimal(coordinate, 24),
                    "sign": sign(coordinate),
                    "exchangeResidual": decimal(residual, 12),
                }
            )
    return records


def polynomial_from_negative_roots(roots: list[int]) -> list[mp.mpf]:
    coefficients = [mp.mpf("1")]
    for root in roots:
        next_coefficients = [mp.mpf("0")] * (len(coefficients) + 1)
        for index, value in enumerate(coefficients):
            next_coefficients[index] += value
            next_coefficients[index + 1] += value / root
        coefficients = next_coefficients
    return coefficients


def control_summary(max_rank: int, max_shift: int) -> dict[str, object]:
    pf_control = polynomial_from_negative_roots([1, 2, 4, 8, 16, 32, 64, 128])
    non_pf_control = [mp.mpf("1"), mp.mpf("1"), mp.mpf("10")]
    output: dict[str, object] = {}
    for name, sequence in (("negative-real-root-polynomial", pf_control), ("positive-coefficient-non-pf", non_pf_control)):
        records = consecutive_scan(sequence, min(max_rank, 5), min(max_shift, 7))
        output[name] = {
            "coefficients": [decimal(value, 18) for value in sequence],
            "negativeConsecutiveMinors": [
                {"rank": record["rank"], "shift": record["shift"], "value": record["value"]}
                for record in records
                if record["sign"] == "negative"
            ],
        }
    return output


def series_checks(coefficients: list[mp.mpf]) -> list[dict[str, str]]:
    checks: list[dict[str, str]] = []
    for z in (mp.mpf("0.25"), mp.mpf("1"), mp.mpf("4")):
        series_value = mp.fsum(value * z**k for k, value in enumerate(coefficients))
        direct_value = xi(mp.mpf("0.5") + mp.sqrt(z) / 2) / 8
        checks.append(
            {
                "z": decimal(z, 8),
                "series": decimal(series_value),
                "directCompletedXi": decimal(direct_value),
                "relativeError": decimal(abs(series_value - direct_value) / abs(direct_value), 12),
            }
        )
    return checks


def summarize_signs(records: list[dict[str, object]]) -> dict[str, int]:
    return {
        label: sum(record["sign"] == label for record in records)
        for label in ("positive", "negative", "zero")
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--theta-terms", type=int, default=14)
    parser.add_argument("--max-coefficient", type=int, default=20)
    parser.add_argument("--max-rank", type=int, default=6)
    parser.add_argument("--max-shift", type=int, default=12)
    parser.add_argument("--plucker-rank", type=int, default=4)
    parser.add_argument("--plucker-width", type=int, default=9)
    parser.add_argument("--plucker-max-offset", type=int, default=5)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    needed_index = max(
        args.max_shift + args.max_rank - 1,
        args.plucker_max_offset + args.plucker_width - 1,
    )
    if args.max_coefficient < needed_index:
        raise SystemExit(f"--max-coefficient must be at least {needed_index}")

    mp.mp.dps = args.dps
    coefficients = xi_coefficients(args.max_coefficient, args.theta_terms)
    normalized, scale = normalized_coefficients(coefficients)
    consecutive = consecutive_scan(normalized, args.max_rank, args.max_shift)
    exchange = exchange_coordinates(normalized, args.max_rank, args.max_shift)

    result = {
        "schema": "riemann-lab.positive-grassmannian-scan.v1",
        "status": "high_precision_numeric_not_interval_certified",
        "normalization": {
            "completedFunction": "xi(1/2 + sqrt(z)/2) / 8 = sum_k a_k z^k",
            "kernel": "Phi(u)=sum_n (2*pi^2*n^4*e^(9u)-3*pi*n^2*e^(5u))*e^(-pi*n^2*e^(4u))",
            "minorSignPreservingScale": "b_k=(a_k/a_0)*(a_0/a_1)^k",
            "scale": decimal(scale),
        },
        "parameters": {
            "decimalPrecision": args.dps,
            "thetaTerms": args.theta_terms,
            "integrationBreaks": list(INTEGRATION_BREAKS),
            "maxCoefficient": args.max_coefficient,
            "maxRank": args.max_rank,
            "maxShift": args.max_shift,
        },
        "coefficients": [
            {"index": k, "a": decimal(value), "normalized": decimal(normalized[k])}
            for k, value in enumerate(coefficients)
        ],
        "seriesChecks": series_checks(coefficients),
        "consecutiveMinorSummary": summarize_signs(consecutive),
        "consecutiveMinors": consecutive,
        "pluckerSlices": plucker_scan(
            normalized,
            args.plucker_rank,
            args.plucker_width,
            range(args.plucker_max_offset + 1),
        ),
        "exchangeCoordinateSummary": summarize_signs(exchange),
        "exchangeCoordinates": exchange,
        "controls": control_summary(args.max_rank, args.max_shift),
        "interpretationBoundary": [
            "Positive finite scans do not prove RH.",
            "The signs are not certified because mpmath does not use interval arithmetic.",
            "The positive-coefficient control must fail, showing that coefficient positivity alone is insufficient.",
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
