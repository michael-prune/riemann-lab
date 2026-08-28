import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const result = JSON.parse(
  fs.readFileSync(
    new URL(
      "../experiments/xi-rank-two-threshold-20260827/results/threshold-v1.json",
      import.meta.url,
    ),
    "utf8",
  ),
);

const midpoint = (value) => Number(value.match(/^\[([^ ]+)/)?.[1]);

test("rank-two threshold reduction remains conditional and bounded", () => {
  assert.equal(
    result.schema,
    "riemann-lab.xi-rank-two-threshold-certificate.v1",
  );
  assert.equal(result.previousThreshold / result.newThreshold, 10_000);
  assert.ok(midpoint(result.factorization.relativeError) < 0.223);
  assert.ok(
    midpoint(result.signBudgetAtThreshold.relativeErrorCauchyUpper) < 0.237,
  );
  assert.ok(midpoint(result.signBudgetAtThreshold.curvatureCauchyUpper) < 0.053);
  assert.ok(midpoint(result.signBudgetAtThreshold.finalMarginLower) > 1.23);
  assert.match(result.conditionalDiscreteConsequence, /100001/);
  assert.ok(result.nonclaims.some((line) => /unreplayed/i.test(line)));
  assert.ok(result.nonclaims.some((line) => /0<=x<100000/i.test(line)));
  assert.ok(result.nonclaims.some((line) => /No unconditional global/i.test(line)));
});
