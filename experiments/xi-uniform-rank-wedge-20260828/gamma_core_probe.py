#!/usr/bin/env python3
"""Exact-Gamma core convergence to the rank/shift continuum model."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
from scipy.integrate import solve_ivp


def text(value, digits: int = 24) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def lambda_of_t(t: mp.mpf, amplitude: mp.mpf) -> mp.mpf:
    y = 1 - mp.exp(-t)
    a = 1 - 1 / amplitude
    b = 1 + 1 / amplitude
    value = mp.hyp2f1(a, b, 2, y)
    derivative = a * b / 2 * mp.hyp2f1(a + 1, b + 1, 3, y)
    u = y * value
    u_t = (1 - y) * (value + y * derivative)
    return amplitude * u_t / u


def continuum_values(lambdas, amplitude):
    start = 500.0
    amplitude_float = float(amplitude)
    initial = amplitude_float / start - 0.5 / start**2

    def rhs(lam, value):
        t = value[0]
        reciprocal = 0.0 if t > 700 else 1.0 / math.expm1(t)
        return [-amplitude_float / (lam * lam + reciprocal)]

    solution = solve_ivp(
        rhs,
        (start, float(min(lambdas))),
        [initial],
        rtol=1e-12,
        atol=1e-14,
        dense_output=True,
        max_step=0.02,
    )
    assert solution.success
    return {str(lam): mp.mpf(float(solution.sol(float(lam))[0])) for lam in lambdas}


def log_determinant_table(max_rank: int, max_index: int):
    rows = {
        0: [mp.mpf("0") for _ in range(max_index + 1)],
        1: [-mp.loggamma(2 * shift + 1) for shift in range(max_index + 1)],
    }
    for rank in range(2, max_rank + 1):
        row = []
        for shift in range(max_index - rank + 2):
            previous = rows[rank - 1]
            ratio = (
                mp.mpf("0")
                if shift == 0
                else mp.exp(previous[shift - 1] + previous[shift + 1] - 2 * previous[shift])
            )
            if not 0 <= ratio < 1:
                raise ArithmeticError(f"precision loss at rank={rank}, shift={shift}")
            row.append(
                2 * previous[shift]
                - rows[rank - 2][shift]
                + mp.log1p(-ratio)
            )
        rows[rank] = row
    return rows


def tau(rows, rank: int, shift: int) -> mp.mpf:
    row = rows[rank]
    return (
        2 * row[shift] - row[shift + 1]
        if shift == 0
        else 2 * row[shift] - row[shift - 1] - row[shift + 1]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=350)
    parser.add_argument("--max-rank", type=int, default=120)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    lambdas = [mp.mpf(value) for value in ("0.25", "0.5", "0.75", "1", "1.5", "2", "3")]
    continuum = continuum_values(lambdas, mp.mpf(2))
    max_index = 4 * args.max_rank + 3
    rows = log_determinant_table(args.max_rank, max_index)
    records = []
    monotonicity = []
    for rank in (10, 20, 50, 80, args.max_rank):
        if rank > args.max_rank:
            continue
        first_shift = max(1, rank - 1)
        values = [tau(rows, rank, shift) for shift in range(first_shift, 3 * rank + 2)]
        differences = [left - right for left, right in zip(values, values[1:])]
        monotonicity.append(
            {
                "rank": rank,
                "shiftRange": [first_shift, 3 * rank + 1],
                "strictlyDecreasing": all(value > 0 for value in differences),
                "minimumDecrease": text(min(differences)),
            }
        )
        for lam in lambdas:
            shift = int(mp.nint(lam * rank))
            value = tau(rows, rank, shift)
            limit = continuum[str(lam)]
            records.append(
                {
                    "rank": rank,
                    "lambda": text(lam, 6),
                    "tau": text(value),
                    "continuumT": text(limit),
                    "relativeError": text(abs(value - limit) / limit),
                }
            )

    family = []
    for amplitude in map(mp.mpf, ("1", "1.25", "1.5", "2")):
        values = continuum_values(lambdas, amplitude)
        samples = []
        for lam in lambdas:
            value = values[str(lam)]
            denominator = lam**2 + 1 / mp.expm1(value)
            samples.append(
                {
                    "lambda": text(lam, 6),
                    "T": text(value),
                    "derivative": text(-amplitude / denominator),
                    "normalizedU": text(lam * value / amplitude),
                }
            )
        family.append({"amplitude": text(amplitude, 6), "samples": samples})

    result = {
        "schema": "riemann-lab.xi-uniform-rank-gamma-core.v1",
        "status": "high_precision_gamma_model_numeric_not_interval_certified",
        "gammaSequence": "a_k=1/Gamma(2k+1)",
        "continuumParametrization": "y=1-exp(-T), u=y*2F1(1-1/A,1+1/A;2;y), lambda=A*u_T/u",
        "hypergeometricIdentitySpotCheck": {
            "A": "2",
            "T": "1",
            "lambda": text(lambda_of_t(mp.mpf(1), mp.mpf(2))),
        },
        "parameters": {
            "decimalPrecision": args.dps,
            "maxRank": args.max_rank,
            "maxCoefficient": max_index,
        },
        "continuumFamily": family,
        "gammaMonotonicity": monotonicity,
        "gammaConvergence": records,
        "allGammaMonotone": all(item["strictlyDecreasing"] for item in monotonicity),
        "boundary": "The hypergeometric identity is analytic; the finite Gamma condensation and convergence rates reported here are high-precision diagnostics, not interval certificates.",
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
