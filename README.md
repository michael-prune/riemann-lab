# Riemann Lab

Can humans and AI turn a famous open problem into a stream of claims that can
actually be checked?

Riemann Lab is a public constitution and research-state repository for a
Labarchy research lab. It organizes bounded, evidence-bearing work around the
Riemann Hypothesis (RH); it does **not** claim that RH is solved.

- RH remains open.
- Checking finitely many zeros is not a proof of RH.
- A Lean-checked lemma can still formalize the wrong mathematical statement.
- Inner Core means internally verified under a pinned protocol, not peer review,
  novelty, Clay-prize eligibility, or community acceptance.
- The recent result proving that more than two thirds of zeta zeros are simple
  and on the critical line is an important neighboring milestone, not RH.

Start with [LAB_PROTOCOL.md](LAB_PROTOCOL.md),
[CLAIM_LEVELS.md](CLAIM_LEVELS.md), and [lab.json](lab.json). The eight initial
lanes are indexed in [lanes/index.json](lanes/index.json). Submissions must
validate against [schemas/submission.v1.schema.json](schemas/submission.v1.schema.json).

## Lean starter

The root Lean project is a small pinned Lean 4 project for controlled future contributions.
It deliberately separates a trusted statement from a solution. Run `lake build`
after installing the toolchain in `lean-toolchain`; `npm test` verifies the
repository contracts and deterministic verifier fixtures.

## Public links

- Labarchy landing page: `https://lablab-en5.pages.dev/riemann-lab.html`
- Outer Core and Inner Core are public Labarchy views; private workspaces are
  account scoped and never enumerated here.
