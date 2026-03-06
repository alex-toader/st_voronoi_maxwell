# TODOs & Open Directions — Paper 3

**Date:** Mar 2026

---

## Infrastructure gaps

| ID | What | Needed for | Status |
|----|------|-----------|--------|
| I1 | Non-cubic box builder (L_vec with Lx ≠ Ly ≠ Lz) | Anisotropic cell test (§2) — verify G = Vol·I on non-cubic periodic domain | open |

---

## Test consolidations (reviewer suggestions)

Tests to add to existing files. Each reviewer suggestion evaluated and assigned.

### File 1 (§2 metric identity)

| ID | Test | What it verifies | Status |
|----|------|-----------------|--------|
| T1 | test_divergence_theorem_term3 | Term3=0 (face centroid offset cancellation) — needs dual face geometry | open (complex) |
| T2 | test_T2_per_vertex_flux_closure | Per-vertex flux Σ_{e∋v} sign·⋆₁·Δx = 0 (Term 1 = 0) | DONE |
| T3 | — | Rotation invariance — redundant with R1 random (15 meshes, no symmetry) | skip |
| T4 | test_T4_tight_frame | {√⋆₁ · Δx_e} tight frame, F†F = Vol·I, σ all equal | DONE |
| T5 | test_T5_admissible_positivity | Cone radius: Kelvin 1.35×, C15 5.95×, WP 2.59× | DONE |
| T6 | — | Subdomain additivity — trivial (linearity of sum) | skip |
| T7 | test_T7_domain_scaling | G → s³·I under V→sV, L→sL (dimensional correctness) | DONE |
| T8 | test_anisotropic_cell | G = Vol·I on non-cubic box (needs I1 infrastructure) | blocked by I1 |

---

## Skeleton writing — proof gaps

Items identified from critical review of skeleton_v1.tex physics logic.

### MAJOR (proof incomplete without these)

| ID | What | Section | Notes |
|----|------|---------|-------|
| L1 | Prove discrete curl identity ∂(h₂†S₂ d₁(k) h)/∂k = iε | §3.2 | DONE — Lemma + proof added: Stokes integral, parallelogram areas, H=Vol·I contraction. |
| L3 | Prove Schur complement = u_perp†S₂u_perp | §3.2 | DONE — Lemma + proof added: u = u_im + u_perp decomposition, B·D⁻¹·B† cancels u_im part. |
| L4 | Expand H = Vol·I proof | §2.2 | DONE — Dual complex argument with dual perpendicularity (Voronoi edge ⊥ Delaunay face). |

### MODERATE (presentation/logic)

| ID | What | Section | Notes |
|----|------|---------|-------|
| P1 | Fix attribution of exactness role | §3.2 | DONE — Added NOTE clarifying real use is at k≠0 for gauge/physical separation. |
| P2 | Sketch z≠4 argument | §3.2 | DONE — Expanded NOTE: d₂d₁=0 ⇒ u ∈ ker(d₂), im(d₂†) enters at O(k²) only, c²=1 for general Voronoi. |
| P3 | Expand converse proof sketch | §3.3 | DONE — Full proof: 6-direction argument with explicit (k,α) pairs extracting all G,H components. |

### MINOR

| ID | What | Section | Notes |
|----|------|---------|-------|
| M1 | Move direction dependence line | §6 | DONE — Removed from §6 (belongs in §4). |
| M2 | Fix abstract random mesh count | Abstract | DONE — Changed to "5 random". |

---

## Skeleton critical review — Mar 2026

Issues found from deep scientific review of skeleton_v1.tex. Two rounds.

### Round 1 — ERRORS

| ID | What | Section | Notes |
|----|------|---------|-------|
| E1 | V-E+F = C ≠ 0; complement dim = C+2 not 3 | §3.2 | DONE — Rewrote with correct general Euler relation. |
| E2 | Converse equation for general G,H | §3.3 | RETRACTED — equation IS correct (H⁻¹ cancels). Added derivation note. |

### Round 1 — HIDDEN ASSUMPTIONS

| ID | What | Section | Notes |
|----|------|---------|-------|
| A1 | H¹(k≠0) = 0 assumed not proved | §3.1 | DONE — Restated as numerical fact. |
| A2 | Gauge decoupling silent in L3 | §3.2 | DONE — Added explicit block. |
| A3 | Optical gap stays open at small k | §3.2 | DONE — Added gap values + continuity. |

### Round 1 — FLOW

| ID | What | Section | Notes |
|----|------|---------|-------|
| F1 | z≠4 NOTE should be main argument | §3.2 | DONE — Merged into general argument. |
| F2 | h_α M-orthonormality unstated | §3.1 | DONE — Added explicit derivation. |
| F3 | Girth → trace mapping imprecise | §6 | DONE — Clarified walk length 2n. |

### Round 2

| ID | What | Section | Notes |
|----|------|---------|-------|
| C1 | "Both necessary" — B necessity numerical only | §1 | DONE — Changed to "sufficient" + qualified. |
| C2 | Term 3 cancellation is pair-cell not per-face | §2.1 | DONE — Rewrote proof sketch with correct mechanism. |
| C3 | Dual perpendicularity wording | §2.2 | DONE — Explicit: ẽ ⊥ e* because ẽ lies in ⊥-bisector. |
| C4 | "face centroid = edge midpoint" wrong for dual | §2.2 | DONE — Removed, replaced with pair-cell cancellation. |
| C7 | Δd₁ undefined | §3.2 | DONE — Wrote [d₁(k)−d₁(0)]·harm. |
| C8 | d₂d₁=0 assumed not attributed | §3.2 | DONE — Added as ingredient (3), attributed to ∂²=0 on chain complex. |
| C9 | Index β 3× in L1 contraction | §3.2 L1 | DONE — Fixed with β' dummy index, shows H_{ββ'} before H=Vol·I. |
| C11 | Converse "general G,H" requires exactness | §3.3 | DONE — Added "Assuming exactness (ingredient 3)". |
| C12 | Off-diag extraction incomplete | §3.3 | DONE — Expanded with explicit (k×α)ᵀH(k×α) algebra. |
| C13 | Random c_std ref R39→R38 | §4 | DONE — Fixed. |
| C14 | Killer table Kelvin-only not labeled | §5 | DONE — Added label + note about full paper. |
| C15 | "Holonomy" undefined in girth argument | §6 | DONE — Added: product of Bloch phases around cycle. |
| C16 | "(C) combinatorial (z=4)" wrong | §7 | DONE — Changed to topological β₁(T³)=3. |
| C17 | Perpendicularity conflates G and H | §7 | DONE — Separated: face∥edge for G, edge⊥Delaunay for H. |

---

## Open directions (for paper discussion / future work)

| ID | Direction | Priority |
|----|-----------|----------|
| D1 | Tight frame interpretation — edge vectors form tight frame, connects to harmonic analysis | high — §2 framing |
| D2 | Homogenization interpretation — G = Vol·I = trivial homogenized metric | high — §7 discussion |
| D3 | d-dimensional generalization — same identity in ℝ^d | medium — §7 |
| D4 | Regge calculus connection — Voronoi DEC metric vs Regge metric | low — future |
| D5 | Topological perturbations (Pachner moves) — how fast does c degrade | medium — future |
| D6 | Generalization beyond Voronoi — c²=1 proof uses (A) exactness, (B) dim H¹=3, (C) G=H=Vol·I, (D) curl identity. Voronoi satisfies all four but may not be the only class. Circumcentric DEC complexes (Delaunay dual) are candidates, but 3D circumcentric is not guaranteed (obtuse tetrahedra). Need to check if any non-Voronoi complex satisfies (C). Unclear if this is a true extension or just Voronoi in disguise. | low — investigate later |

---
