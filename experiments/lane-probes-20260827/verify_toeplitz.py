#!/usr/bin/env python3
"""Independent higher-precision rerun of the xi exchange-ratio probe."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path


PROBE_PATH = Path(__file__).with_name("probe.py")


def load_probe_module():
    spec = importlib.util.spec_from_file_location("lane_probe", PROBE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {PROBE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    probe = load_probe_module()
    full = probe.toeplitz_probe(
        max_coefficient=54,
        decimal_precision=140,
        theta_terms=18,
    )
    result = {
        "schema": "riemann-lab.xi-exchange-ratio-independent-check.v1",
        "status": full["status"],
        "parameters": full["parameters"],
        "candidateInequality": full["candidateInequality"],
        "rankResults": [
            {
                "rank": record["rank"],
                "allPositive": record["allPositive"],
                "strictlyDecreasing": record["strictlyDecreasing"],
                "lastRatio": record["lastRatio"],
            }
            for record in full["rankResults"]
        ],
        "boundary": "Independent higher-precision numerical agreement is not an interval certificate or a theorem.",
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
