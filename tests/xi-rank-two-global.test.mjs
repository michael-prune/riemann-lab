import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const result = JSON.parse(
  fs.readFileSync(
    new URL("../experiments/xi-rank-two-global-20260827/results/probe-v1.json", import.meta.url),
    "utf8",
  ),
);

test("global rank-two probe separates xi evidence from the generic no-go", () => {
  assert.equal(result.schema, "riemann-lab.xi-rank-two-global-probe.v1");
  assert.equal(result.xiDenseScan.failureCount, 0);
  assert.equal(
    result.decreasingLogConcaveControl.outcome,
    "violates_the_gamma_skew_comparison",
  );
  assert.ok(Number(result.decreasingLogConcaveControl.margin) < 0);
  assert.ok(Number(result.saddleBudget.crudeMainBoundForUAtLeastPoint98) < 0.477);
  assert.ok(Number(result.saddleBudget.unallocatedCoefficientMargin) > 1.5);
  assert.match(result.saddleBudget.missingProofObligation, /full theta saddle/);
  assert.ok(result.claimBoundaries.some((line) => line.includes("not an internal theorem")));
});
