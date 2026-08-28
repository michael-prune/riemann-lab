#!/usr/bin/env python3
"""Certify a conditional effective rank-three xi exchange tail.

The scalar argument composes the conditional rank-two saddle certificate with
an analytic derivative hierarchy for tau=2f-f(x-1)-f(x+1).  The underlying
source saddle/contour lemmas remain unreplayed hypotheses.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx


ctx.prec = 256


def ball(value: arb, digits: int = 30) -> str:
    return value.str(digits, radius=True)


def exact_control() -> dict:
    q = [
        Fraction(11, 100),
        Fraction(7, 50),
        Fraction(6, 25),
        Fraction(23, 50),
        Fraction(11, 20),
        Fraction(73, 100),
        Fraction(22, 25),
    ]
    a = [Fraction(1), Fraction(1)]
    for value in q:
        a.append(value * a[-1] * a[-1] / a[-2])

    def b(index: int) -> Fraction:
        return a[index] * a[index] - a[index - 1] * a[index + 1]

    def q3(index: int) -> Fraction:
        return b(index) * b(index) / (b(index - 1) * b(index + 1)) - 1

    q2 = [1 / value - 1 for value in q]
    q3_values = [q3(index) for index in range(2, 6)]
    assert all(left > right for left, right in zip(q2, q2[1:]))
    assert q3_values[3] > q3_values[2]

    def fraction(value: Fraction) -> str:
        return f"{value.numerator}/{value.denominator}"

    return {
        "qSequence": [fraction(value) for value in q],
        "qDefinition": "q_k=a_(k-1)a_(k+1)/a_k^2",
        "rankTwoExchangeRatios": [fraction(value) for value in q2],
        "rankTwoStrictlyDecreasing": True,
        "rankThreeShifts": [2, 3, 4, 5],
        "rankThreeExchangeRatios": [fraction(value) for value in q3_values],
        "rankThreeViolation": "Q_(3,5)>Q_(3,4)",
        "conclusion": "decreasing rank-two exchange ratios do not formally imply rank-three monotonicity",
    }


def certify() -> dict:
    # The rank-two threshold certificate supplies x^2 f'''(x)>1.23 for
    # real x>=10^5, conditional on the source analytic lemmas.
    rank_two_floor = arb("1.23")
    analytic_base = arb("1e5")
    rank_three_threshold = arb("1e6")

    # Rebuild a complex |f'''| bound on |z-n|<=.03n from the same source
    # factorization.  Nested Cauchy radius is .02n.
    nested_radius = arb(1) / 50
    saddle_second = arb("0.645")
    curvature_first = arb("8.45")
    factorization_log = arb("0.252052")
    gamma_real_part = arb("1.9")

    saddle_third_coefficient = saddle_second / nested_radius
    curvature_third_coefficient = (
        curvature_first / nested_radius**2 / analytic_base
    )
    factorization_third_coefficient = (
        6 * factorization_log / nested_radius**3 / analytic_base
    )
    gamma_third_coefficient = (
        8 / gamma_real_part**2
        + 16 / (gamma_real_part**3 * analytic_base)
    )
    complex_f3_coefficient = (
        saddle_third_coefficient
        + curvature_third_coefficient
        + factorization_third_coefficient
        + gamma_third_coefficient
    )
    assert complex_f3_coefficient < arb(37)

    # Cauchy applied to f''' on the nested disk.
    f3_coefficient = arb(37)
    computed_f4 = f3_coefficient / nested_radius
    computed_f5 = 2 * f3_coefficient / nested_radius**2
    assert computed_f4 < arb(1851)
    assert computed_f5 < arb(185001)
    # These are exact rational simplifications of 37/(1/50) and
    # 2*37/(1/50)^2; the preceding outward checks guard the implementation.
    f4_coefficient = arb(1850)
    f5_coefficient = arb(185000)

    # At integer n, the source-backed Turan curvature window gives
    # tau_n>1/(2n).  Moving at most 1/2 using |tau'|<=37/n^2 leaves
    # tau(x)>.499/n for n>=10^5.
    tau_floor_exact = arb("0.5") - f3_coefficient / (2 * analytic_base)
    tau_floor = arb("0.499")
    assert tau_floor_exact > tau_floor

    # For J(t)=log(1-exp(-t)) and t>0:
    # J'<=1/t, |J''|<=1/t^2, J'''<=2/t^3.
    # tau derivatives are bounded by the corresponding central averages of
    # f''', f'''', and f'''''.
    correction_constant = (
        2 * f3_coefficient**3 / tau_floor**3
        + 3 * f3_coefficient * f4_coefficient / tau_floor**2
        + f5_coefficient / tau_floor
    )

    center_ratio = 1 + 1 / (2 * rank_three_threshold)
    correction_coefficient = (
        correction_constant * center_ratio**2 / rank_three_threshold
    )
    inherited_coefficient = 2 * rank_two_floor
    final_margin = inherited_coefficient - correction_coefficient

    assert correction_coefficient < arb("2.011")
    assert final_margin > arb("0.449")

    return {
        "schema": "riemann-lab.xi-rank-three-kill-test-certificate.v1",
        "status": "checked_effective_rank_three_tail_conditional_on_unreplayed_source_analytic_lemmas",
        "conditionalContinuousConclusion": "(log D_2(x))'''>0 for every real x>=1000000",
        "conditionalDiscreteConsequence": "Q_(3,k+1)<Q_(3,k) for every integer k>=1000001",
        "arbPrecisionBits": ctx.prec,
        "exactGenericControl": exact_control(),
        "derivativeHierarchy": {
            "complexF3CoefficientUpper": ball(complex_f3_coefficient),
            "roundedF3CoefficientUpper": ball(f3_coefficient),
            "F4CoefficientUpper": ball(f4_coefficient),
            "F5CoefficientUpper": ball(f5_coefficient),
            "tauFloorExactAtAnalyticBase": ball(tau_floor_exact),
            "tauFloorUsed": ball(tau_floor),
        },
        "rankThreeBudgetAtThreshold": {
            "inheritedTwoF3CoefficientLower": ball(inherited_coefficient),
            "nonlinearCorrectionConstant": ball(correction_constant),
            "nonlinearCorrectionCoefficientUpper": ball(correction_coefficient),
            "finalCoefficientMarginLower": ball(final_margin),
        },
        "mechanism": [
            "g=2f+J(tau), tau=2f-f(x-1)-f(x+1), J(t)=log(1-exp(-t))",
            "the inherited term 2f''' is order x^-2",
            "the nonlinear J(tau) correction is bounded by order x^-3",
            "the lift uses the same saddle derivative hierarchy as rank two, but is not a generic consequence of rank-two monotonicity",
        ],
        "nonclaims": [
            "The source saddle and contour certificates remain unavailable and unreplayed.",
            "The interval below 1000000 is supported numerically, not interval-certified.",
            "This fixed-rank tail mechanism does not control ranks growing with k or the RH-critical k comparable to r regime.",
            "No global rank-three theorem, all-rank theorem, RH implication, or novelty claim is made.",
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
