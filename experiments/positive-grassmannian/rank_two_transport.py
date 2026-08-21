#!/usr/bin/env python3
"""Rank-two transport and square-tail sign-regularity probes.

This is a high-precision numerical experiment, not interval arithmetic and not
a proof.  It tests two sharply separated questions:

1. Does the signed rank-two integrand have a positive fiber at fixed uv?
2. Does the CNV square-tail kernel K(x+y) exhibit reverse sign-regularity on a
   small ordered grid?
"""

from __future__ import annotations

import argparse
import itertools
import json
from functools import lru_cache
from pathlib import Path

import mpmath as mp


GRID = ("0.001", "0.003", "0.01", "0.03", "0.1", "0.3")
FIBER_CASES = (
    (1, "0.001"),
    (1, "0.01"),
    (2, "0.01"),
    (2, "0.03"),
    (3, "0.01"),
    (3, "0.03"),
    (5, "0.03"),
    (5, "0.1"),
)


def classical_phi(u: mp.mpf, theta_terms: int) -> mp.mpf:
    """Positive Riemann theta kernel in the lab's normalization."""
    e4 = mp.exp(4 * u)
    e5 = mp.exp(5 * u)
    e9 = mp.exp(9 * u)
    return mp.fsum(
        (2 * mp.pi**2 * n**4 * e9 - 3 * mp.pi * n**2 * e5)
        * mp.exp(-mp.pi * n**2 * e4)
        for n in range(1, theta_terms + 1)
    )


def decimal(value: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def sign(value: mp.mpf) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def reverse_sign(rank: int) -> int:
    return -1 if (rank * (rank - 1) // 2) % 2 else 1


def fixed_product_fiber(k: int, x: mp.mpf, theta_terms: int) -> mp.mpf:
    """Inner y-integral after x=uv and y=log(u/v)."""
    c_k = mp.mpf(2 * k * (2 * k - 1)) / ((2 * k + 1) * (2 * k + 2))
    root_x = mp.sqrt(x)

    def integrand(y: mp.mpf) -> mp.mpf:
        u = root_x * mp.exp(y / 2)
        v = root_x * mp.exp(-y / 2)
        weight = classical_phi(u, theta_terms) * classical_phi(v, theta_terms)
        return weight * (1 - c_k * mp.cosh(2 * y))

    # The integrand is even.  At y=12 one theta factor is already far below
    # the working precision for every committed x value.
    return 2 * mp.quadgl(integrand, (0, 1, 2, 4, 6, 8, 10, 12))


def square_tail_kernel(t: mp.mpf, theta_terms: int) -> mp.mpf:
    """K(t)=integral_{sqrt(t)}^infinity u Phi(u) du."""
    lower = mp.sqrt(t)
    breaks = [lower]
    breaks.extend(
        point
        for point in map(mp.mpf, ("0.25", "0.5", "0.75", "1", "1.5", "2"))
        if point > lower
    )
    # Phi(2) is smaller than exp(-9000), so truncating at 2 is invisible at
    # the committed precision and avoids a costly infinite-interval transform.
    return mp.quadgl(lambda u: u * classical_phi(u, theta_terms), breaks)


def sign_regular_scan(theta_terms: int) -> dict[str, object]:
    grid = tuple(mp.mpf(value) for value in GRID)

    @lru_cache(maxsize=None)
    def kernel_at(value: str) -> mp.mpf:
        return square_tail_kernel(mp.mpf(value), theta_terms)

    matrix = mp.matrix(
        [[kernel_at(decimal(x + y, 50)) for y in grid] for x in grid]
    )
    rank_records: list[dict[str, object]] = []
    for rank in range(1, len(grid) + 1):
        expected = reverse_sign(rank)
        oriented_values: list[mp.mpf] = []
        bad: list[dict[str, object]] = []
        for rows in itertools.combinations(range(len(grid)), rank):
            for columns in itertools.combinations(range(len(grid)), rank):
                minor = mp.det(mp.matrix([[matrix[i, j] for j in columns] for i in rows]))
                oriented = expected * minor
                oriented_values.append(oriented)
                if oriented <= 0:
                    bad.append(
                        {
                            "rows": list(rows),
                            "columns": list(columns),
                            "determinant": decimal(minor),
                        }
                    )
        rank_records.append(
            {
                "rank": rank,
                "expectedSign": "positive" if expected > 0 else "negative",
                "minorCount": len(oriented_values),
                "badSignCount": len(bad),
                "minimumOrientedMinor": decimal(min(oriented_values)),
                "badMinors": bad,
            }
        )

    return {
        "grid": list(GRID),
        "kernelMatrix": [
            [decimal(matrix[i, j]) for j in range(len(grid))]
            for i in range(len(grid))
        ],
        "ranks": rank_records,
        "totalMinorCount": sum(record["minorCount"] for record in rank_records),
        "totalBadSignCount": sum(record["badSignCount"] for record in rank_records),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--theta-terms", type=int, default=14)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mp.mp.dps = args.dps

    fibers = []
    for k, rendered_x in FIBER_CASES:
        value = fixed_product_fiber(k, mp.mpf(rendered_x), args.theta_terms)
        fibers.append({"k": k, "x": rendered_x, "value": decimal(value), "sign": sign(value)})

    result = {
        "schema": "riemann-lab.rank-two-transport.v1",
        "status": "high_precision_numeric_not_interval_certified",
        "parameters": {
            "decimalPrecision": args.dps,
            "thetaTerms": args.theta_terms,
            "fiberCutoff": 12,
            "squareTailUpperCutoff": 2,
        },
        "identities": {
            "squareTail": "K(t)=integral_(sqrt(t))^infinity u*Phi(u)du=(1/2)integral_t^infinity Phi(sqrt(v))dv",
            "fixedProductFiber": "integral_R Phi(sqrt(x)e^(y/2))*Phi(sqrt(x)e^(-y/2))*(1-c_k*cosh(2y))dy",
            "cK": "c_k=(2k)(2k-1)/((2k+1)(2k+2))",
        },
        "fixedProductFibers": fibers,
        "fixedProductConclusion": (
            "Both signs occur. A positivity proof cannot pair chambers while preserving x=uv fiber-by-fiber."
        ),
        "squareTailSignRegularity": sign_regular_scan(args.theta_terms),
        "interpretationBoundary": [
            "The finite grid is not evidence that all minors at all points have the predicted sign.",
            "mpmath does not provide interval enclosures.",
            "Even global sign-regularity of this transformed kernel would still require a non-circular transfer to every xi Toeplitz minor.",
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
