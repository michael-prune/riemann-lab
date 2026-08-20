# Riemann Lab Protocol

## Purpose

Riemann Lab turns research contributions into small, inspectable artifacts.
It studies RH and adjacent zeta/L-function questions without presenting a
community consensus, novelty claim, or prize claim as an internal decision.

## Layers

1. **Workspace** is private and account scoped.
2. **Submission inbox** holds immutable structured packets.
3. **Outer Core** is public provisional memory. Its records carry claim level,
   assumptions, evidence, nonclaims, failure mode, verification state, and Core
   decision.
4. **Inner Core** is public internal knowledge admitted only after the named
   deterministic gates and a bounded Core semantic decision.

An Inner Core record is neither external peer review nor a claim of novelty,
acceptance, prize eligibility, or a proof of RH.

## Core policy

The bouncer first applies deterministic schema, identity, lane, hash,
duplicate, and artifact checks. Only then may a bounded low-cost model classify
the packet. Raw model text is retained separately from normalized packets and
is never executable policy. Decisions are immutable and cached by submission,
policy, and model hashes. A circuit breaker and paid-call kill switch yield
`review_pending`, not invented progress.

Any claim to prove/refute full RH, establish external novelty, or qualify for
the Clay prize must receive `human_attention_required`; it cannot be promoted
automatically.

## Lean Inner-Core gate

Lean contributions are accepted only through the controlled tree. The gate
requires pinned Lean/Mathlib, an allowlisted import set, owned build command,
no changed trusted files or workflow/toolchain/configuration, no `sorry` or
`admit`, no unsafe declarations or disallowed axioms, built declaration,
restricted `#print axioms`, trusted-statement comparison, source and artifact
hashes, a second clean replay where practical, semantic audit, and explicit
nonclaims. `lake build` alone is insufficient.

## Lane proposals

“Suggest a lane” creates a structured Outer-Core proposal, never a lane. It
must name the missing object, why existing lanes fail to cover it, known-work
anchors, verification mode, expected overlap, and a closure/kill condition.
