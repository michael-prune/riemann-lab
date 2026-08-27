import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const root = new URL("../", import.meta.url);
const json = async (file) => JSON.parse(await readFile(new URL(file, root), "utf8"));

test("positive-geometry baseline records finite scans and a falsifying control", async () => {
  const result = await json("experiments/positive-grassmannian/results/baseline-v1.json");
  assert.equal(result.schema, "riemann-lab.positive-grassmannian-scan.v1");
  assert.equal(result.status, "high_precision_numeric_not_interval_certified");
  assert.equal(result.parameters.decimalPrecision, 100);
  assert.equal(result.parameters.thetaTerms, 14);
  assert.equal(result.consecutiveMinorSummary.negative, 0);
  assert.equal(result.exchangeCoordinateSummary.negative, 0);
  assert.equal(result.pluckerSlices.reduce((sum, slice) => sum + slice.pluckerCoordinateCount, 0), 756);
  assert.ok(
    result.controls["positive-coefficient-non-pf"].negativeConsecutiveMinors.some(
      (minor) => minor.rank === 2 && minor.shift === 1 && minor.value.startsWith("-9."),
    ),
  );
  assert.ok(result.interpretationBoundary.some((line) => /do not prove RH/i.test(line)));
});

test("positive-geometry lane anchors the relevant primary literature", async () => {
  const known = await json("references/known-work.json");
  for (const id of [
    "katkova-multiple-positivity",
    "xi-toeplitz-cubic-wedge",
    "postnikov-positive-grassmannian",
    "hydrotope",
    "cnv-turan-1986",
    "nuttall-rh-rr3",
    "nuttall-cumulants",
    "xi-kernel-second-level-concavity",
  ]) {
    assert.ok(known.records.some((record) => record.id === id), `missing ${id}`);
  }
});

test("rank-two transport records the fiber obstruction and finite sign-regular scan", async () => {
  const result = await json("experiments/positive-grassmannian/results/rank-two-transport-v1.json");
  assert.equal(result.schema, "riemann-lab.rank-two-transport.v1");
  assert.equal(result.status, "high_precision_numeric_not_interval_certified");
  assert.ok(result.fixedProductFibers.some((fiber) => fiber.sign === "positive"));
  assert.ok(result.fixedProductFibers.some((fiber) => fiber.sign === "negative"));
  assert.equal(result.squareTailSignRegularity.totalMinorCount, 923);
  assert.equal(result.squareTailSignRegularity.totalBadSignCount, 0);
  assert.deepEqual(
    result.squareTailSignRegularity.ranks.map((record) => record.minorCount),
    [36, 225, 400, 225, 36, 1],
  );
  assert.ok(result.interpretationBoundary.some((line) => /not.*all minors/i.test(line)));
});

test("coalescing-point test falsifies square-tail rank five", async () => {
  const result = await json("experiments/positive-grassmannian/results/square-tail-wronskian-v1.json");
  assert.equal(result.schema, "riemann-lab.square-tail-wronskian-test.v1");
  assert.equal(result.status, "high_precision_numeric_not_interval_certified");
  assert.equal(result.rankSummary["4"].derivativeFailurePoints.length, 0);
  assert.ok(result.rankSummary["5"].derivativeFailurePoints.includes("0.001"));
  assert.deepEqual(
    result.rankSummary["5"].directCounterexampleSteps,
    ["0.0001", "0.00005", "0.00002", "0.00001"],
  );
  assert.equal(result.rankSummary["6"].directCounterexampleSteps.length, 0);
  assert.ok(result.conclusion.some((line) => /cannot be RR6/i.test(line)));
  assert.ok(result.interpretationBoundary.some((line) => /interval/i.test(line)));
});
