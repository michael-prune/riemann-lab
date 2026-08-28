#!/usr/bin/env python3
"""High-precision diagonal xi/Toda probe against the continuum limit."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import mpmath as mp
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parents[2]
SCAN_PATH = ROOT / "experiments" / "positive-grassmannian" / "scan.py"


def load_scan():
    spec = importlib.util.spec_from_file_location("positive_geometry_scan", SCAN_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {SCAN_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def text(value, digits: int = 22) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def log_determinant_table(log_a: list[mp.mpf], max_rank: int):
    max_index = len(log_a) - 1
    negative_infinity = mp.ninf
    table = {
        0: [mp.mpf("0") for _ in range(max_index + 1)],
        1: log_a,
    }
    for rank in range(2, max_rank + 1):
        row = []
        for shift in range(max_index - rank + 2):
            center = table[rank - 1][shift]
            left = negative_infinity if shift == 0 else table[rank - 1][shift - 1]
            right = table[rank - 1][shift + 1]
            ratio = mp.mpf("0") if left == negative_infinity else mp.exp(left + right - 2 * center)
            if not 0 <= ratio < 1:
                raise ArithmeticError(f"nonpositive condensation pivot at rank={rank}, shift={shift}")
            row.append(2 * center - table[rank - 2][shift] + mp.log1p(-ratio))
        table[rank] = row
    return table


def tau(table, rank: int, shift: int) -> mp.mpf:
    center = table[rank][shift]
    left = mp.ninf if shift == 0 else table[rank][shift - 1]
    right = table[rank][shift + 1]
    return 2 * center - right if left == mp.ninf else 2 * center - left - right


def continuum_solution(lambdas: list[float]):
    start = 500.0
    t0 = 2 / start - 1 / (2 * start**2) + 5 / (12 * start**3)

    def rhs(lam, value):
        t = value[0]
        reciprocal = 0.0 if t > 700 else 1.0 / mp.expm1(t)
        return [-2.0 / (lam * lam + float(reciprocal))]

    end = min(lambdas)
    solution = solve_ivp(
        rhs,
        (start, end),
        [t0],
        rtol=1e-12,
        atol=1e-14,
        dense_output=True,
        max_step=0.02,
    )
    assert solution.success
    return {lam: float(solution.sol(lam)[0]) for lam in lambdas}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=120)
    parser.add_argument("--theta-terms", type=int, default=16)
    parser.add_argument("--max-rank", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    scan = load_scan()
    lambdas = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    max_shift = int(max(lambdas) * args.max_rank) + 2
    max_coefficient = max_shift + args.max_rank
    coefficients = scan.xi_coefficients(max_coefficient, args.theta_terms)
    normalized, _ = scan.normalized_coefficients(coefficients)
    table = log_determinant_table([mp.log(value) for value in normalized], args.max_rank)
    continuum = continuum_solution(lambdas)

    records = []
    monotonicity = []
    gamma_one_barrier = []
    for rank in (5, 10, 15, args.max_rank):
        if rank > args.max_rank:
            continue
        first_shift = max(1, rank - 1)
        values = [tau(table, rank, shift) for shift in range(first_shift, 3 * rank + 2)]
        decreases = [left - right for left, right in zip(values, values[1:])]
        monotonicity.append(
            {
                "rank": rank,
                "shiftRange": [first_shift, 3 * rank + 1],
                "strictlyDecreasing": all(value > 0 for value in decreases),
                "minimumDecrease": text(min(decreases)),
            }
        )
        barrier_records = []
        for shift in range(first_shift, 3 * rank + 1):
            value = tau(table, rank, shift)
            next_value = tau(table, rank, shift + 1)
            lower = mp.log1p(mp.mpf(rank) / shift)
            next_lower = mp.log1p(mp.mpf(rank) / (shift + 1))
            barrier_records.append(
                {
                    "shift": shift,
                    "curvatureMargin": value - lower,
                    "decreaseMargin": (value - next_value) - (lower - next_lower),
                }
            )
        gamma_one_barrier.append(
            {
                "rank": rank,
                "shiftRange": [first_shift, 3 * rank],
                "curvatureDominates": all(item["curvatureMargin"] > 0 for item in barrier_records),
                "decreaseDominates": all(item["decreaseMargin"] > 0 for item in barrier_records),
                "minimumCurvatureMargin": text(min(item["curvatureMargin"] for item in barrier_records)),
                "minimumDecreaseMargin": text(min(item["decreaseMargin"] for item in barrier_records)),
            }
        )
        for lam in lambdas:
            shift = max(1, round(lam * rank))
            effective_lambda = shift / rank
            value = tau(table, rank, shift)
            limit = continuum[lam] if effective_lambda == lam else continuum_solution([effective_lambda])[effective_lambda]
            records.append(
                {
                    "rank": rank,
                    "shift": shift,
                    "lambda": effective_lambda,
                    "tau": text(value),
                    "normalizedU": text(mp.mpf(shift) * value / (2 * rank)),
                    "continuumT": text(limit),
                    "relativeContinuumError": text(abs(value - limit) / limit),
                }
            )

    result = {
        "schema": "riemann-lab.xi-uniform-rank-diagonal-probe.v1",
        "status": "high_precision_numeric_not_interval_certified",
        "parameters": {
            "decimalPrecision": args.dps,
            "thetaTerms": args.theta_terms,
            "maxRank": args.max_rank,
            "maxCoefficient": max_coefficient,
        },
        "continuumEquation": "T'=-2/(lambda^2+1/(exp(T)-1)), T~2/lambda",
        "monotonicity": monotonicity,
        "gammaOneBarrier": {
            "model": "a_k=1/k!, tau_(r,k)=log(1+r/k)",
            "records": gamma_one_barrier,
            "allCurvaturesDominate": all(item["curvatureDominates"] for item in gamma_one_barrier),
            "allDecreasesDominate": all(item["decreaseDominates"] for item in gamma_one_barrier),
        },
        "diagonalComparison": records,
        "allMonotone": all(item["strictlyDecreasing"] for item in monotonicity),
        "boundary": "Coefficient quadrature, condensation, and continuum integration are not interval-certified; convergence is diagnostic only.",
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
