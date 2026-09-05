import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const load = (relative) =>
  JSON.parse(fs.readFileSync(new URL(relative, import.meta.url), "utf8"));

const formal = load(
  "../experiments/xi-uniform-rank-wedge-20260828/results/formal-toda-v1.json",
);
const gamma = load(
  "../experiments/xi-uniform-rank-wedge-20260828/results/gamma-core-v1.json",
);
const diagonal = load(
  "../experiments/xi-uniform-rank-wedge-20260828/results/diagonal-v1.json",
);
const controls = load(
  "../experiments/xi-uniform-rank-wedge-20260828/results/controls-v1.json",
);

test("formal Toda coefficients have the rank degree predicted by diagonal scaling", () => {
  assert.equal(formal.schema, "riemann-lab.xi-uniform-rank-formal-toda.v1");
  assert.equal(formal.heldOutRanks, 12);
  for (const coefficient of formal.coefficientLaws) {
    assert.equal(coefficient.degreeInRank, coefficient.inversePower - 1);
  }
  assert.match(formal.boundary, /do not bound the truncated remainder/i);
});

test("Gamma core converges to the continuum and remains monotone on the solid range", () => {
  assert.equal(gamma.schema, "riemann-lab.xi-uniform-rank-gamma-core.v1");
  assert.equal(gamma.allGammaMonotone, true);
  const rank120 = gamma.gammaConvergence.filter((item) => item.rank === 120 && Number(item.lambda) >= 1);
  assert.ok(rank120.every((item) => Number(item.relativeError) < 2e-6));
});

test("xi diagonal probe dominates the explicit Gamma-one barrier", () => {
  assert.equal(diagonal.schema, "riemann-lab.xi-uniform-rank-diagonal-probe.v1");
  assert.equal(diagonal.allMonotone, true);
  assert.equal(diagonal.gammaOneBarrier.allCurvaturesDominate, true);
  assert.equal(diagonal.gammaOneBarrier.allDecreasesDominate, true);
  assert.match(diagonal.boundary, /not interval-certified/i);
});

test("moment control remains valid and the false quartic counterexample is withdrawn", () => {
  assert.equal(controls.schema, "riemann-lab.xi-uniform-rank-controls.v1");
  assert.equal(controls.positiveMomentMixture.allInvolvedMinorsPositive, true);
  assert.match(controls.positiveMomentMixture.conclusion, /does not imply exchange monotonicity/i);
  assert.match(controls.smoothLogConcaveBoundary.conclusion, /withdrawn counterexample/i);
  assert.match(controls.smoothLogConcaveBoundary.facts.join(" "), /only real zeros/i);
  assert.match(controls.smoothLogConcaveBoundary.source, /1901\.06596/);
  assert.equal(controls.nearbyGammaPrecisionControl.determinantSign, "positive");
});
