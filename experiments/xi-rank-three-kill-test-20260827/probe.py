#!/usr/bin/env python3
"""Saddle-centered high-precision rank-three xi diagnostic."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


def text(value: mp.mpf, digits: int = 24) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def dominant_k(u: mp.mpf) -> mp.mpf:
    b = 2 * mp.pi * mp.exp(4 * u)
    return u * (b - mp.mpf("4.5") - 6 / (b - 3))


def phi(u: mp.mpf, terms: int) -> mp.mpf:
    e4 = mp.exp(4 * u)
    e5 = mp.exp(5 * u)
    e9 = mp.exp(9 * u)
    return mp.fsum(
        (2 * mp.pi**2 * n**4 * e9 - 3 * mp.pi * n**2 * e5)
        * mp.exp(-mp.pi * n * n * e4)
        for n in range(1, terms + 1)
    )


def tilted_stats(x: mp.mpf, terms: int) -> dict[str, mp.mpf]:
    saddle = mp.findroot(lambda u: dominant_k(u) - x, mp.log(max(x, 16)) / 4)
    sigma = mp.sqrt(saddle / x)
    lower = max(mp.mpf("1e-8"), saddle - 12 * sigma)
    upper = saddle + 12 * sigma
    shift = 2 * x * mp.log(saddle) + mp.log(phi(saddle, terms))
    log_saddle = mp.log(saddle)

    def integrand(u: mp.mpf, power: int) -> mp.mpf:
        centered = mp.log(u) - log_saddle
        return (
            mp.exp(2 * x * mp.log(u) + mp.log(phi(u, terms)) - shift)
            * centered**power
        )

    raw = [
        mp.quad(lambda u, power=power: integrand(u, power), [lower, saddle, upper])
        for power in range(4)
    ]
    moments = [value / raw[0] for value in raw]
    variance = moments[2] - moments[1] ** 2
    third = (
        moments[3]
        - 3 * moments[2] * moments[1]
        + 2 * moments[1] ** 3
    )
    return {
        "logA": shift + mp.log(raw[0]) - mp.loggamma(2 * x + 1),
        "f1": 2 * (log_saddle + moments[1]) - 2 * mp.digamma(2 * x + 1),
        "f2": 4 * variance - 4 * mp.polygamma(1, 2 * x + 1),
        "f3": 8 * third - 8 * mp.polygamma(2, 2 * x + 1),
        "saddle": saddle,
    }


def rank_three_sample(x: int, terms: int) -> dict[str, object]:
    left, center, right = [
        tilted_stats(mp.mpf(x + offset), terms) for offset in (-1, 0, 1)
    ]
    h = left["logA"] + right["logA"] - 2 * center["logA"]
    ratio = mp.exp(h)
    h1 = left["f1"] + right["f1"] - 2 * center["f1"]
    h2 = left["f2"] + right["f2"] - 2 * center["f2"]
    h3 = left["f3"] + right["f3"] - 2 * center["f3"]
    correction_terms = [
        -ratio * (1 + ratio) / (1 - ratio) ** 3 * h1**3,
        -3 * ratio / (1 - ratio) ** 2 * h1 * h2,
        -ratio / (1 - ratio) * h3,
    ]
    inherited = 2 * center["f3"]
    correction = mp.fsum(correction_terms)
    value = inherited + correction
    return {
        "x": x,
        "saddle": text(center["saddle"]),
        "tau": text(-h),
        "fThird": text(center["f3"]),
        "inheritedTwoFThird": text(inherited),
        "nonlinearCorrection": text(correction),
        "logD2Third": text(value),
        "xSquaredLogD2Third": text(mp.mpf(x) ** 2 * value),
        "positive": value > 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--theta-terms", type=int, default=12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    samples = [100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 2000000]
    records = [rank_three_sample(x, args.theta_terms) for x in samples]
    result = {
        "schema": "riemann-lab.xi-rank-three-wide-probe.v1",
        "status": "high_precision_saddle_centered_numeric_not_interval_certified",
        "decimalPrecision": args.dps,
        "thetaTerms": args.theta_terms,
        "integrationWindow": "dominant saddle plus or minus 12 sqrt(u_s/x)",
        "samples": records,
        "allPositive": all(record["positive"] for record in records),
        "boundary": "The integration window, theta cutoff, and quadrature are not interval-certified.",
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
