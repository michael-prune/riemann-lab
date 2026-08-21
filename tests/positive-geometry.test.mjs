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
  ]) {
    assert.ok(known.records.some((record) => record.id === id), `missing ${id}`);
  }
});
