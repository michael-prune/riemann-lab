#!/usr/bin/env python3
"""Saddle-centered high-precision rank-four xi diagnostic."""

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
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    return mp.fsum(
        (2 * mp.pi**2 * n**4 * e9 - 3 * mp.pi * n**2 * e5)
        * mp.exp(-mp.pi * n * n * e4)
        for n in range(1, terms + 1)
    )


def tilted_stats(x: mp.mpf, terms: int) -> list[mp.mpf]:
    saddle = mp.findroot(lambda u: dominant_k(u) - x, mp.log(max(x, 16)) / 4)
    sigma = mp.sqrt(saddle / x)
    lower, upper = max(mp.mpf("1e-8"), saddle - 12 * sigma), saddle + 12 * sigma
    shift = 2 * x * mp.log(saddle) + mp.log(phi(saddle, terms))
    log_saddle = mp.log(saddle)

    def integrand(u: mp.mpf, power: int) -> mp.mpf:
        centered = mp.log(u) - log_saddle
        return mp.exp(2 * x * mp.log(u) + mp.log(phi(u, terms)) - shift) * centered**power

    raw = [
        mp.quad(lambda u, power=power: integrand(u, power), [lower, saddle, upper])
        for power in range(4)
    ]
    moments = [value / raw[0] for value in raw]
    variance = moments[2] - moments[1] ** 2
    third = moments[3] - 3 * moments[2] * moments[1] + 2 * moments[1] ** 3
    return [
        shift + mp.log(raw[0]) - mp.loggamma(2 * x + 1),
        2 * (log_saddle + moments[1]) - 2 * mp.digamma(2 * x + 1),
        4 * variance - 4 * mp.polygamma(1, 2 * x + 1),
        8 * third - 8 * mp.polygamma(2, 2 * x + 1),
    ]


def lift(left: list[mp.mpf], center: list[mp.mpf], right: list[mp.mpf]) -> list[mp.mpf]:
    h = [left[index] + right[index] - 2 * center[index] for index in range(4)]
    ratio = mp.exp(h[0])
    transforms = [
        mp.log(1 - ratio),
        -ratio / (1 - ratio),
        -ratio / (1 - ratio) ** 2,
        -ratio * (1 + ratio) / (1 - ratio) ** 3,
    ]
    return [
        2 * center[0] + transforms[0],
        2 * center[1] + transforms[1] * h[1],
        2 * center[2] + transforms[2] * h[1] ** 2 + transforms[1] * h[2],
        2 * center[3]
        + transforms[3] * h[1] ** 3
        + 3 * transforms[2] * h[1] * h[2]
        + transforms[1] * h[3],
    ]


def sample(x: int, terms: int) -> dict[str, object]:
    f = {offset: tilted_stats(mp.mpf(x + offset), terms) for offset in range(-2, 3)}
    g = {offset: lift(f[offset - 1], f[offset], f[offset + 1]) for offset in (-1, 0, 1)}
    p = lift(g[-1], g[0], g[1])
    # lift(g) computes 2g+J(tau_g); Desnanot-Jacobi requires 2g-f+J(tau_g).
    p = [p[index] - f[0][index] for index in range(4)]
    return {
        "x": x,
        "fThird": text(f[0][3]),
        "logD2Third": text(g[0][3]),
        "logD3Third": text(p[3]),
        "xSquaredLogD3Third": text(mp.mpf(x) ** 2 * p[3]),
        "positive": p[3] > 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--theta-terms", type=int, default=12)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    records = [sample(x, args.theta_terms) for x in (100, 300, 1000, 3000, 10000, 100000, 1000000)]
    result = {
        "schema": "riemann-lab.xi-rank-four-wide-probe.v1",
        "status": "high_precision_saddle_centered_numeric_not_interval_certified",
        "decimalPrecision": args.dps,
        "thetaTerms": args.theta_terms,
        "samples": records,
        "allPositive": all(record["positive"] for record in records),
        "boundary": "The integration window, theta cutoff, quadrature, and signs are not interval-certified.",
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
