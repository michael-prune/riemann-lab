import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const result = JSON.parse(
  fs.readFileSync(
    new URL("../experiments/xi-exchange-monotonicity-20260827/results/probe-v1.json", import.meta.url),
    "utf8",
  ),
);

test("exchange-monotonicity evidence retains its claim boundaries", () => {
  assert.equal(result.schema, "riemann-lab.xi-exchange-monotonicity-probe.v1");
  assert.equal(result.status, "high_precision_numeric_not_interval_certified");
  assert.equal(result.parameters.maxRank, 15);
  assert.ok(result.finiteStress.every((record) => record.allPositive));
  assert.ok(result.finiteStress.every((record) => record.strictlyDecreasing));
  assert.ok(result.continuousCurvature.every((record) => Number(record.logAThirdDerivative) > 0));
  assert.ok(result.continuousCurvature.every((record) => Number(record.logD2ThirdDerivative) > 0));
  assert.equal(result.pfInfinityControl.allZerosNegativeReal, true);
  assert.equal(result.pfInfinityControl.firstMonotonicityViolations.length, 7);
  assert.ok(result.claimBoundaries.some((line) => line.includes("non-interval")));
  assert.ok(result.claimBoundaries.some((line) => line.includes("sufficiently large shift")));
});
