import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const root = new URL("../", import.meta.url);
const json = async (file) => JSON.parse(await readFile(new URL(file, root), "utf8"));

test("four-lane probe records bounded positive and negative results", async () => {
  const result = await json("experiments/lane-probes-20260827/results/probes-v1.json");
  assert.equal(result.schema, "riemann-lab.lane-probes.v1");

  const weil = result.probes.weilFiniteCompression;
  assert.equal(weil.outcome, "closed_within_first_two_moments_and_bandwidth_one");
  assert.ok(Number(weil.relativeError) < 5e-5);
  assert.ok(Number(weil.minimumSampledOptimizerValue) > 0);

  const toeplitz = result.probes.xiToeplitzFrontier;
  assert.equal(toeplitz.rankResults.length, 8);
  assert.ok(toeplitz.rankResults.every((record) => record.allPositive));
  assert.ok(toeplitz.rankResults.every((record) => record.strictlyDecreasing));
  assert.equal(toeplitz.pfInfinityControl.firstMonotonicityViolations.length, 7);

  const heat = result.probes.deBruijnNewmanHeatFlow;
  assert.equal(heat.outcome, "naive_pointwise_gap_barrier_falsified");
  assert.equal(heat.allGapSignsStableAcrossCutoffs, true);
  assert.ok(heat.cutoffs.every((record) => record.negativeDerivativeCount > 0));
  assert.ok(heat.cutoffs.every((record) => record.positiveDerivativeCount > 0));

  const nyman = result.probes.nymanBeurlingGram;
  assert.equal(nyman.outcome, "diagonal_scaling_improves_constant_not_growth_order");
  assert.ok(Number(nyman.logLogFit.rawConditionPower) > 1.5);
  assert.ok(Number(nyman.logLogFit.diagonalScaledConditionPower) > 1.5);
  assert.ok(
    Number(nyman.records.at(-1).distanceSquared) < Number(nyman.records[0].distanceSquared),
  );

  assert.ok(result.globalBoundary.some((line) => /No probe proves RH/i.test(line)));
});

test("xi exchange-ratio lead survives the independent precision rerun", async () => {
  const result = await json(
    "experiments/lane-probes-20260827/results/toeplitz-independent-v1.json",
  );
  assert.equal(result.schema, "riemann-lab.xi-exchange-ratio-independent-check.v1");
  assert.equal(result.parameters.decimalPrecision, 140);
  assert.equal(result.parameters.thetaTerms, 18);
  assert.ok(result.rankResults.every((record) => record.allPositive));
  assert.ok(result.rankResults.every((record) => record.strictlyDecreasing));
  assert.match(result.boundary, /not.*theorem/i);
});
