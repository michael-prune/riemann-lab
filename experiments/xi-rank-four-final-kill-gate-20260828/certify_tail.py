#!/usr/bin/env python3
"""Directed rank-four xi tail budget and rank-growth diagnostic.

The effective statement is conditional on the same unreplayed source analytic
lemmas as the rank-two and rank-three certificates.  Arb checks only the
scalar implications assembled here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import arb, ctx


ctx.prec = 256


def ball(value: arb, digits: int = 30) -> str:
    return value.str(digits, radius=True)


def certify() -> dict:
    rank_two_floor = arb("1.23")
    f3, f4, f5 = arb(37), arb(1850), arb(185000)
    cauchy_radius = arb(1) / 50
    f6 = 6 * f3 / cauchy_radius**3
    f7 = 24 * f3 / cauchy_radius**4
    tau_f_floor = arb("0.499")

    # For J(t)=log(1-exp(-t)), t>0, use
    # |J^(m)(t)| <= (m-1)!/t^m for m=1,...,5.
    h2 = f3**2 / tau_f_floor**2 + f4 / tau_f_floor
    h3 = (
        2 * f3**3 / tau_f_floor**3
        + 3 * f3 * f4 / tau_f_floor**2
        + f5 / tau_f_floor
    )
    h4 = (
        6 * f3**4 / tau_f_floor**4
        + 12 * f3**2 * f4 / tau_f_floor**3
        + 3 * f4**2 / tau_f_floor**2
        + 4 * f3 * f5 / tau_f_floor**2
        + f6 / tau_f_floor
    )
    h5 = (
        24 * f3**5 / tau_f_floor**5
        + 60 * f3**3 * f4 / tau_f_floor**4
        + 30 * f3 * f4**2 / tau_f_floor**3
        + 20 * f3**2 * f5 / tau_f_floor**3
        + 10 * f4 * f5 / tau_f_floor**2
        + 5 * f3 * f6 / tau_f_floor**2
        + f7 / tau_f_floor
    )

    # g=log D_2=2f+J(tau_f).
    # H^(m) has one extra inverse power of n relative to f^(m), so convert
    # every term to the coefficient of g^(m)=O(n^(1-m)) at n>=10^6.
    g3 = 2 * f3 + h3 / arb("1e6")
    g4 = 2 * f4 + h4 / arb("1e6")
    g5 = 2 * f5 + h5 / arb("1e6")

    # tau_g=2g-g(x-1)-g(x+1).  At integer n,
    # tau_g >= 1/n-h2/n^2.  Move by at most 1/2 with |tau_g'|<=|g'''|.
    curvature_base = arb("1e6")
    tau_g_floor_exact = 1 - h2 / curvature_base - g3 / (2 * curvature_base)
    tau_g_floor = arb("0.98")
    assert tau_g_floor_exact > tau_g_floor

    rank_four_correction = (
        2 * g3**3 / tau_g_floor**3
        + 3 * g3 * g4 / tau_g_floor**2
        + g5 / tau_g_floor
    )
    inherited_correction = 2 * h3
    total_correction = inherited_correction + rank_four_correction

    # F_3''' = 3f''' + 2J(tau_f)''' + J(tau_g)'''.
    # Convert n^-3 errors to an x^-2 coefficient at the nearest integer.
    threshold = arb("2e6")
    center_ratio = 1 + 1 / (2 * threshold)
    correction_coefficient = total_correction * center_ratio**2 / threshold
    inherited_coefficient = 3 * rank_two_floor
    margin = inherited_coefficient - correction_coefficient
    assert margin > arb("0.36")

    # Exact formal rank recurrence through the first universal correction:
    # c_(r+1)=2c_r-c_(r-1), c_0=0,c_1=2;
    # e_(r+1)=2e_r-e_(r-1)-2, e_0=e_1=0.
    ranks = []
    c_prev, c_now = 0, 2
    e_prev, e_now = 0, 0
    for rank in range(1, 13):
        if rank > 1:
            c_prev, c_now = c_now, 2 * c_now - c_prev
            e_prev, e_now = e_now, 2 * e_now - e_prev - 2
        assert c_now == 2 * rank
        assert e_now == -rank * (rank - 1)
        ranks.append(
            {
                "rank": rank,
                "xMinus2Coefficient": c_now,
                "universalXMinus3RankCorrection": e_now,
            }
        )

    return {
        "schema": "riemann-lab.xi-rank-four-final-kill-gate.v1",
        "status": "checked_rank_four_scalar_tail_conditional_on_unreplayed_source_analytic_lemmas",
        "conditionalContinuousConclusion": "(log D_3(x))'''>0 for every real x>=2000000",
        "conditionalDiscreteConsequence": "Q_(4,k+1)<Q_(4,k) for every integer k>=2000001",
        "arbPrecisionBits": ctx.prec,
        "exactRecursion": "F_r=2F_(r-1)-F_(r-2)+J(2F_(r-1)-F_(r-1)(x-1)-F_(r-1)(x+1))",
        "derivativeBudget": {
            "f3": ball(f3),
            "f4": ball(f4),
            "f5": ball(f5),
            "f6": ball(f6),
            "f7": ball(f7),
            "JTauFSecondCorrection": ball(h2),
            "JTauFThirdCorrection": ball(h3),
            "JTauFFourthCorrection": ball(h4),
            "JTauFFifthCorrection": ball(h5),
            "g3AtCurvatureBase": ball(g3),
            "g4": ball(g4),
            "g5": ball(g5),
            "tauGFloorExactAtCurvatureBase": ball(tau_g_floor_exact),
            "tauGFloorUsed": ball(tau_g_floor),
        },
        "rankFourBudgetAtThreshold": {
            "inheritedThreeF3CoefficientLower": ball(inherited_coefficient),
            "rankThreeCorrectionConstantTwice": ball(inherited_correction),
            "rankFourNonlinearCorrectionConstant": ball(rank_four_correction),
            "totalCorrectionCoefficientUpper": ball(correction_coefficient),
            "finalCoefficientMarginLower": ball(margin),
        },
        "formalRankScaling": {
            "statement": "F_r'''(x)=2r/x^2+(r*d_1-r(r-1))/x^3+lower terms for each fixed r, if f'''=2/x^2+d_1/x^3+lower terms",
            "universalRecurrence": ranks,
            "interpretation": "the first rank-dependent correction is quadratic, not factorial; at x=lambda*r its ratio to the leading term tends 1/(2*lambda)",
            "boundary": "fixed-r asymptotic algebra is not a remainder estimate uniform in r",
        },
        "decision": {
            "fixedRankFour": "passes with a conservative threshold only twice the rank-three threshold",
            "rawCauchyConstants": "grow quickly before normalization and are not themselves a uniform proof",
            "normalizedMechanism": "survives because the exact first two asymptotic orders have polynomial rank growth",
            "finalGate": "not killed mathematically; no theorem-level uniform-r bridge has been established",
        },
        "nonclaims": [
            "The source saddle and contour certificates remain unavailable and unreplayed.",
            "The derivative budget is deliberately conservative and does not estimate the true rank-four onset.",
            "The formal rank scaling is not uniform in rank and cannot be substituted for an all-rank inequality.",
            "No global rank-four theorem, all-rank theorem, RH implication, or novelty claim is made.",
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
