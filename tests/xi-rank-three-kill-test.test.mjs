import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const load = (relative) =>
  JSON.parse(fs.readFileSync(new URL(relative, import.meta.url), "utf8"));
const midpoint = (value) => Number(value.match(/^\[([^ ]+)/)?.[1]);

const tail = load(
  "../experiments/xi-rank-three-kill-test-20260827/results/tail-v1.json",
);
const probe = load(
  "../experiments/xi-rank-three-kill-test-20260827/results/wide-probe-v1.json",
);

test("rank-three kill test records the effective tail and exact obstruction", () => {
  assert.equal(
    tail.schema,
    "riemann-lab.xi-rank-three-kill-test-certificate.v1",
  );
  assert.match(tail.conditionalDiscreteConsequence, /1000001/);
  assert.ok(
    midpoint(tail.rankThreeBudgetAtThreshold.finalCoefficientMarginLower) > 0.449,
  );
  assert.equal(tail.exactGenericControl.rankTwoStrictlyDecreasing, true);
  assert.equal(tail.exactGenericControl.rankThreeViolation, "Q_(3,5)>Q_(3,4)");
  assert.ok(tail.nonclaims.some((line) => /fixed-rank/i.test(line)));
});

test("wide rank-three xi diagnostic stays positive with a numerical boundary", () => {
  assert.equal(probe.schema, "riemann-lab.xi-rank-three-wide-probe.v1");
  assert.equal(probe.allPositive, true);
  assert.equal(probe.samples.at(-1).x, 2_000_000);
  assert.match(probe.boundary, /not interval-certified/i);
});
