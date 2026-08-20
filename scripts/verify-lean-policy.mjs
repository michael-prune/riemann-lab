import { readFile } from "node:fs/promises";

const protectedPaths = ["lean-toolchain", "lakefile.toml", "RiemannLab/Trusted/OnePlusOne.lean"];
const files = ["RiemannLab/Trusted/OnePlusOne.lean", "RiemannLab/Solutions/OnePlusOne.lean"];
for (const file of files) {
  const source = await readFile(file, "utf8");
  if (/\b(sorry|admit|unsafe)\b/.test(source)) throw new Error(`lean_policy_forbidden_token:${file}`);
}
console.log(JSON.stringify({ok:true,protectedPaths,verification:"static policy; run lake build separately"}));
