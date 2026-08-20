import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { readdir, readFile } from "node:fs/promises";
import test from "node:test";

const root = new URL("../", import.meta.url);
const json = async (file) => JSON.parse(await readFile(new URL(file, root), "utf8"));
const sha = (value) => createHash("sha256").update(value).digest("hex");

test("schemas are strict and public indexes resolve", async () => {
  for (const file of ["submission","core-decision","review-task","review-result","verification-manifest"]) {
    const schema = await json(`schemas/${file}.v1.schema.json`);
    assert.equal(schema.additionalProperties, false);
  }
  const lanes = await json("lanes/index.json");
  assert.equal(lanes.lanes.length, 8);
  for (const lane of lanes.lanes) assert.match(await readFile(new URL(`lanes/${lane.readme}`, root), "utf8"), /#/);
  for (const index of ["outer-core/index.json", "inner-core/index.json"]) assert.ok((await json(index)).layer.endsWith("core"));
});

test("valid submission has no contributor-controlled decision fields", async () => {
  const packet = await json("examples/valid-literature-map.json");
  assert.equal(packet.schema, "riemann-lab.submission.v1");
  assert.equal(packet.labId, "riemann-lab");
  assert.equal(Object.hasOwn(packet, "decision"), false);
  assert.equal(Object.hasOwn(packet, "verifier"), false);
  assert.match(packet.evidence.artifacts[0].sha256, /^[a-f0-9]{64}$/);
  assert.notEqual(sha(JSON.stringify(packet)), "");
});

test("invalid example demonstrates extra-field rejection", async () => {
  const packet = await json("examples/invalid-extra-field.json");
  const schema = await json("schemas/submission.v1.schema.json");
  assert.equal(schema.additionalProperties, false);
  assert.equal(Object.hasOwn(schema.properties, "decision"), false);
  assert.ok(Object.hasOwn(packet, "decision"));
});

test("finding IDs and known-work anchors are stable", async () => {
  const files = await readdir(new URL("findings/", root));
  for (const file of files) assert.match((await json(`findings/${file}`)).id, /^RL-(INNER|OUTER)-\d{4}$/);
  const known = await json("references/known-work.json");
  assert.ok(known.records.some((r) => r.id === "zeta-23" && /not RH/i.test(r.note)));
});

test("Lean manifest pins the controlled statement and build command", async () => {
  const manifest = await json("verification-manifest.v1.json");
  assert.equal(manifest.proofSystem, "lean4");
  assert.equal(manifest.buildCommand, "lake build");
  assert.match(manifest.mathlibRevision, /c44e0c8ee63ca166450922a373c7409c5d26b00b/);
});
