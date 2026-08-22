#!/usr/bin/env python3
"""Test the square-tail kernel at coalescing points through rank six.

Two numerical methods are recorded:

1. derivative Wronskians of K at a logarithmic/linear sweep of t values;
2. actual nearby-point minors, computed without numerical differentiation.

The calculation uses arbitrary-precision quadrature, not interval arithmetic.
It is therefore a reproducible numerical falsification, not a formal proof.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


SWEEP_POINTS = (
    "0.000001",
    "0.000003",
    "0.00001",
    "0.00003",
    "0.0001",
    "0.0003",
    "0.001",
    "0.003",
    "0.01",
    "0.02",
    "0.03",
    "0.05",
    "0.07",
    "0.1",
    "0.15",
    "0.2",
    "0.3",
    "0.5",
)
COUNTEREXAMPLE_T = "0.001"
COUNTEREXAMPLE_STEPS = ("0.0001", "0.00005", "0.00002", "0.00001")


def classical_phi(u: mp.mpf, theta_terms: int) -> mp.mpf:
    e4 = mp.exp(4 * u)
    e5 = mp.exp(5 * u)
    e9 = mp.exp(9 * u)
    return mp.fsum(
        (2 * mp.pi**2 * n**4 * e9 - 3 * mp.pi * n**2 * e5)
        * mp.exp(-mp.pi * n**2 * e4)
        for n in range(1, theta_terms + 1)
    )


def square_kernel(t: mp.mpf, theta_terms: int) -> mp.mpf:
    return classical_phi(mp.sqrt(t), theta_terms)


def square_tail(t: mp.mpf, theta_terms: int) -> mp.mpf:
    lower = mp.sqrt(t)
    breaks = [lower]
    breaks.extend(
        point
        for point in map(mp.mpf, ("0.05", "0.1", "0.2", "0.3", "0.5", "0.75", "1", "1.5", "2"))
        if point > lower
    )
    return mp.quadgl(lambda u: u * classical_phi(u, theta_terms), breaks)


def expected_multiplier(rank: int) -> int:
    return -1 if (rank * (rank - 1) // 2) % 2 else 1


def decimal(value: mp.mpf, digits: int = 40) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def derivative_wronskians(t: mp.mpf, theta_terms: int) -> list[dict[str, object]]:
    s = lambda value: square_kernel(value, theta_terms)
    derivatives = [square_tail(t, theta_terms)]
    derivatives.extend(-mp.diff(s, t, order - 1) / 2 for order in range(1, 11))

    records = []
    for rank in (4, 5, 6):
        wronskian = mp.det(
            mp.matrix(
                [[derivatives[row + column] for column in range(rank)] for row in range(rank)]
            )
        )
        oriented = expected_multiplier(rank) * wronskian
        records.append(
            {
                "rank": rank,
                "expectedWronskianSign": "positive" if expected_multiplier(rank) > 0 else "negative",
                "wronskian": decimal(wronskian),
                "orientedWronskian": decimal(oriented),
                "passes": oriented > 0,
            }
        )
    return records


def nearby_tail_values(t: mp.mpf, h: mp.mpf, maximum_index: int, theta_terms: int) -> list[mp.mpf]:
    """Use K(t+d)=K(t)-(1/2) integral_t^(t+d) s(v)dv.

    This gives closely spaced K values from one shared base integral and avoids
    subtracting independently computed quadratures.
    """
    base = square_tail(t, theta_terms)
    values = [base]
    for index in range(1, maximum_index + 1):
        endpoint = t + index * h
        local_mass = mp.quadgl(lambda v: square_kernel(v, theta_terms), (t, endpoint)) / 2
        values.append(base - local_mass)
    return values


def direct_minor_records(t: mp.mpf, h: mp.mpf, theta_terms: int) -> list[dict[str, object]]:
    values = nearby_tail_values(t, h, 10, theta_terms)
    records = []
    for rank in (4, 5, 6):
        determinant = mp.det(
            mp.matrix([[values[row + column] for column in range(rank)] for row in range(rank)])
        )
        oriented = expected_multiplier(rank) * determinant
        scale = h ** (rank * (rank - 1))
        records.append(
            {
                "rank": rank,
                "determinant": decimal(determinant),
                "orientedDeterminant": decimal(oriented),
                "orientedDividedByCoalescenceScale": decimal(oriented / scale),
                "passes": oriented > 0,
            }
        )
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--derivative-dps", type=int, default=120)
    parser.add_argument("--direct-dps", type=int, default=220)
    parser.add_argument("--theta-terms", type=int, default=16)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    derivative_sweep = []
    with mp.workdps(args.derivative_dps):
        for rendered_t in SWEEP_POINTS:
            t = mp.mpf(rendered_t)
            derivative_sweep.append(
                {"t": rendered_t, "ranks": derivative_wronskians(t, args.theta_terms)}
            )

    direct_checks = []
    with mp.workdps(args.direct_dps):
        t = mp.mpf(COUNTEREXAMPLE_T)
        for rendered_h in COUNTEREXAMPLE_STEPS:
            h = mp.mpf(rendered_h)
            direct_checks.append(
                {
                    "t": COUNTEREXAMPLE_T,
                    "h": rendered_h,
                    "xPoints": [decimal(index * h, 12) for index in range(6)],
                    "yPoints": [decimal(t + index * h, 12) for index in range(6)],
                    "ranks": direct_minor_records(t, h, args.theta_terms),
                }
            )

    rank_summary = {}
    for rank in (4, 5, 6):
        derivative_failures = [
            record["t"]
            for record in derivative_sweep
            if not next(item for item in record["ranks"] if item["rank"] == rank)["passes"]
        ]
        direct_failures = [
            record["h"]
            for record in direct_checks
            if not next(item for item in record["ranks"] if item["rank"] == rank)["passes"]
        ]
        rank_summary[str(rank)] = {
            "derivativeFailurePoints": derivative_failures,
            "directCounterexampleSteps": direct_failures,
            "passesFiniteSweep": not derivative_failures,
            "passesDirectChecksAtCounterexample": not direct_failures,
        }

    result = {
        "schema": "riemann-lab.square-tail-wronskian-test.v1",
        "status": "high_precision_numeric_not_interval_certified",
        "kernel": "K(t)=integral_(sqrt(t))^infinity u*Phi(u)du",
        "requiredSigns": {"rank4": "positive", "rank5": "positive", "rank6": "negative"},
        "parameters": {
            "derivativeDecimalPrecision": args.derivative_dps,
            "directDecimalPrecision": args.direct_dps,
            "thetaTerms": args.theta_terms,
            "sweepRange": [SWEEP_POINTS[0], SWEEP_POINTS[-1]],
        },
        "rankSummary": rank_summary,
        "derivativeSweep": derivative_sweep,
        "directNearbyPointChecks": direct_checks,
        "conclusion": [
            "Rank four had the required sign at every sampled Wronskian and direct-minor point; this is not a global proof.",
            "Rank five has the wrong sign near t=0.001 in both the derivative Wronskian and actual nearby-point minors.",
            "The square-tail kernel therefore fails the numerical RR5 test and cannot be RR6, because RR6 includes all lower ranks.",
            "The isolated rank-six determinant had its expected sign at the tested points, but that does not repair the rank-five failure.",
        ],
        "interpretationBoundary": [
            "mpmath does not provide interval enclosures.",
            "A formal disproof should replace the counterexample calculation with certified ball arithmetic and theta-tail bounds.",
            "No finite sweep can establish the global rank-four sign on t>=0.",
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
