import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const load = (relative) =>
  JSON.parse(fs.readFileSync(new URL(relative, import.meta.url), "utf8"));
const midpoint = (value) => Number(value.match(/^\[?([^ ]+)/)?.[1]);

const tail = load(
  "../experiments/xi-rank-four-final-kill-gate-20260828/results/tail-v1.json",
);
const probe = load(
  "../experiments/xi-rank-four-final-kill-gate-20260828/results/wide-probe-v1.json",
);

test("rank-four gate records a positive conditional tail and exact rank scaling", () => {
  assert.equal(tail.schema, "riemann-lab.xi-rank-four-final-kill-gate.v1");
  assert.match(tail.conditionalDiscreteConsequence, /2000001/);
  assert.ok(midpoint(tail.rankFourBudgetAtThreshold.finalCoefficientMarginLower) > 0.36);
  assert.equal(tail.formalRankScaling.universalRecurrence.at(-1).xMinus2Coefficient, 24);
  assert.equal(tail.formalRankScaling.universalRecurrence.at(-1).universalXMinus3RankCorrection, -132);
  assert.match(tail.formalRankScaling.boundary, /not a remainder estimate uniform in r/i);
});

test("wide rank-four diagnostic stays positive and preserves its numerical boundary", () => {
  assert.equal(probe.schema, "riemann-lab.xi-rank-four-wide-probe.v1");
  assert.equal(probe.allPositive, true);
  assert.equal(probe.samples.at(-1).x, 1_000_000);
  assert.match(probe.boundary, /not interval-certified/i);
});
