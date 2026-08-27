import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const result = JSON.parse(
  fs.readFileSync(
    new URL(
      "../experiments/xi-rank-two-certification-20260827/results/conditional-tail-v1.json",
      import.meta.url,
    ),
    "utf8",
  ),
);

const midpoint = (ball) => Number(ball.match(/^\[([^ ]+)/)?.[1]);

test("rank-two tail certificate preserves its conditional boundary", () => {
  assert.equal(
    result.schema,
    "riemann-lab.xi-rank-two-conditional-tail-certificate.v1",
  );
  assert.equal(
    result.status,
    "checked_scalar_closure_conditional_on_unreplayed_source_lemmas",
  );
  assert.ok(midpoint(result.constants.fullThetaTransferFactor) < 1.000000001);
  assert.ok(midpoint(result.constants.fullSaddleCoefficientUpper) < 0.477);
  assert.ok(midpoint(result.constants.finalCoefficientMarginLowerAt1e9) > 1.52);
  assert.match(result.discreteConsequence, /1000000001/);
  assert.ok(result.nonclaims.some((line) => /not contain.*certificate modules/i.test(line)));
  assert.ok(result.nonclaims.some((line) => /compact range.*not certified/i.test(line)));
  assert.ok(result.nonclaims.some((line) => /No global rank-two theorem/i.test(line)));
});
