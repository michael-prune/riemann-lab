#!/usr/bin/env python3
"""Exact formal Toda expansion for the rank-normalized xi comparison model.

The base model has f'''(x)=2/x^2 exactly.  We propagate
tau_(r+1)=2 tau_r-tau_(r-1)+L log(1-exp(-tau_r))
as a rational power series in z=1/x.  This isolates universal rank growth;
it is discovery algebra, not a uniform remainder theorem for xi.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import comb
from pathlib import Path

import sympy as sp


def add(left, right):
    return [a + b for a, b in zip(left, right)]


def scale(value, factor):
    return [factor * item for item in value]


def multiply(left, right):
    size = len(left)
    output = [Fraction(0) for _ in range(size)]
    for i, a in enumerate(left):
        for j, b in enumerate(right[: size - i]):
            output[i + j] += a * b
    return output


def log_one_plus(value):
    size = len(value)
    output = [Fraction(0) for _ in range(size)]
    power = [Fraction(0) for _ in range(size)]
    power[0] = Fraction(1)
    for degree in range(1, size):
        power = multiply(power, value)
        output = add(output, scale(power, Fraction((-1) ** (degree + 1), degree)))
    return output


def shift(value, sign):
    """Substitute z/(1+sign*z), corresponding to x -> x+sign."""
    size = len(value)
    output = [Fraction(0) for _ in range(size)]
    output[0] = value[0]
    for degree in range(1, size):
        for extra in range(size - degree):
            output[degree + extra] += (
                value[degree]
                * Fraction(comb(degree + extra - 1, extra))
                * ((-sign) ** extra)
            )
    return output


def spatial_laplacian(value):
    return add(scale(value, 2), scale(add(shift(value, -1), shift(value, 1)), -1))


def b_series(tau):
    """B(t)=log((1-exp(-t))/t) through the available order."""
    size = len(tau)
    output = scale(tau, Fraction(-1, 2))
    bernoulli = {
        2: Fraction(1, 6),
        4: Fraction(-1, 30),
        6: Fraction(1, 42),
        8: Fraction(-1, 30),
        10: Fraction(5, 66),
        12: Fraction(-691, 2730),
        14: Fraction(7, 6),
        16: Fraction(-3617, 510),
    }
    power = [Fraction(0) for _ in range(size)]
    power[0] = Fraction(1)
    powers = {0: power}
    for degree in range(1, size):
        power = multiply(power, tau)
        powers[degree] = power
    for degree, number in bernoulli.items():
        if degree >= size:
            continue
        factorial = 1
        for item in range(2, degree + 1):
            factorial *= item
        coefficient = number / (degree * factorial)
        output = add(output, scale(powers[degree], coefficient))
    return output


def laplacian_j(tau):
    size = len(tau)
    leading = tau[1]
    assert leading > 0 and tau[0] == 0
    normalized_tail = [Fraction(0) for _ in range(size)]
    for degree in range(1, size):
        normalized_tail[degree] = tau[degree + 1] / leading if degree + 1 < size else 0
    regular = add(log_one_plus(normalized_tail), b_series(tau))
    output = spatial_laplacian(regular)
    # L log(z)=log(1-z^2)=-sum_(m>=1) z^(2m)/m.
    for degree in range(2, size, 2):
        output[degree] -= Fraction(1, degree // 2)
    return output


def tau_table(max_rank: int, order: int):
    size = order + 1
    zero = [Fraction(0) for _ in range(size)]
    base = zero.copy()
    for m in range(1, (order + 2) // 2):
        degree = 2 * m - 1
        if degree <= order:
            base[degree] = Fraction(2, m * (2 * m - 1))
    table = [zero, base]
    for rank in range(1, max_rank):
        table.append(
            add(
                add(scale(table[rank], 2), scale(table[rank - 1], -1)),
                laplacian_j(table[rank]),
            )
        )
    return table


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rank", type=int, default=24)
    parser.add_argument("--order", type=int, default=10)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    table = tau_table(args.max_rank, args.order)
    r = sp.symbols("r")
    coefficient_laws = []
    for power in range(2, args.order + 2):
        tau_degree = power - 1
        points = []
        for rank in range(1, args.max_rank + 1):
            value = tau_degree * table[rank][tau_degree]
            points.append((rank, sp.Rational(value.numerator, value.denominator)))
        training_count = min(args.max_rank, args.order + 2)
        polynomial = sp.factor(sp.interpolate(points[:training_count], r))
        # Ranks beyond the interpolation set are exact held-out checks.
        assert all(sp.simplify(polynomial.subs(r, rank) - value) == 0 for rank, value in points)
        coefficient_laws.append(
            {
                "inversePower": power,
                "polynomialInRank": str(polynomial),
                "degreeInRank": int(sp.Poly(polynomial, r).degree()),
                "rankValues": [fraction_text(Fraction(value)) for _, value in points[:8]],
            }
        )

    result = {
        "schema": "riemann-lab.xi-uniform-rank-formal-toda.v1",
        "status": "exact_rational_formal_comparison_model_not_uniform_xi_theorem",
        "model": "f'''(x)=2/x^2 exactly; tau_1=L f; exact Toda curvature recursion",
        "maxRankChecked": args.max_rank,
        "interpolationRanks": min(args.max_rank, args.order + 2),
        "heldOutRanks": max(0, args.max_rank - min(args.max_rank, args.order + 2)),
        "seriesOrder": args.order,
        "monotonicityExpansion": "-tau_r'(x)=sum c_(m,r)/x^m",
        "coefficientLaws": coefficient_laws,
        "boundary": "Interpolation identities are exact on the generated ranks but do not bound the truncated remainder or xi-model error uniformly in rank.",
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
