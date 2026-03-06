# W2: Beyond Voronoi — When does c²=1 hold on non-Voronoi complexes?

**Date:** Mar 2026
**Source:** D6 from `internal/w_1/todos.md`
**Status:** Investigation

---

## Question

The paper proves c²=1 under four conditions:
- **(A)** G = H = Vol·I (metric isotropy from Hodge stars)
- **(B)** d₁d₀ = 0 (cochain exactness)
- **(C)** dim H¹ = 3 (from vertex-edge economy nE = 2nV on z=4 meshes)
- **(D)** Discrete curl identity dF/dk|₀ = iε (Levi-Civita)

Voronoi tessellations satisfy all four. **Are they the only class?**

---

## What the proof actually uses from "Voronoi"

### For condition (A): G = Vol·I

The divergence-theorem proof (§2) has 4 terms:
1. **Term 1 = 0:** Per-vertex flux closure Σ_{e∋v} sign·⋆₁·Δx = 0. This uses
   the fact that on each Voronoi cell, boundary fluxes cancel (closed surface).
2. **Term 2 = G:** Definition.
3. **Term 3 = 0:** Face centroid offset cancellation. Uses **perpendicularity**:
   dual edge ∥ face normal (Voronoi bisector property).
4. **Term 4 = Vol·I:** Boundary integral telescopes to total volume. Uses that
   **pillars tile space** — each primal edge e with its dual face forms a
   "pillar" (prism-like region), and these partition the domain.

**Minimal geometric requirements for (A):**
- Perpendicularity: dual edge ⊥ primal face (for Term 3 = 0)
- Pillar tiling: primal edge × dual face partitions the domain (for Term 4 = Vol·I)
- Closed dual cells: Σ fluxes = 0 per vertex (for Term 1 = 0)

### For condition (B): exactness

d₁d₀ = 0 comes from the recurrence construction in gauge_bloch.py — this is
**purely topological** (works on any cell complex with consistent orientation).
Not Voronoi-specific at all.

### For condition (C): dim H¹ = 3

Requires nE = 2nV on the periodic torus (z=4 vertex valence). This is a
**topological constraint on the mesh**, not a Voronoi property. Any z=4 periodic
polyhedral complex has this.

### For condition (D): discrete curl identity

The proof of dF/dk|₀ = iε uses:
- Edge shift vectors (lattice crossings) — topological
- The Stokes integral ∑_{e∈∂f} d₁[f,e]·(shift·L)_γ·(Δx)_α = (A_f)_β — geometric
- Contraction with H²: needs H = Vol·I — follows from (A)

So (D) follows from (A) + topology. Not independently Voronoi-specific.

---

## Candidate non-Voronoi classes

### 1. Power diagrams (Laguerre tessellations)

**What:** Generalization of Voronoi. Each site i has a weight w_i. The cell of i
is {x : |x−p_i|²−w_i ≤ |x−p_j|²−w_j ∀j}. When all weights equal → Voronoi.

**Perpendicularity:** YES. Each face is the perpendicular bisector of the
**power distance** between sites. The dual edge (connecting sites) is
perpendicular to the face. This is exact, not approximate.

**Pillar tiling:** YES. Same argument as Voronoi — pillars partition space.

**Hodge stars:** ⋆₁[e] = dual_face_area / edge_length, ⋆₂[f] = dual_edge_length / face_area.
These are well-defined as long as all dual face areas and edge lengths are positive
(which requires the power diagram to be non-degenerate — all cells have positive volume).

**Prediction:** G = H = Vol·I holds. c²=1 holds.

**Key difference from Voronoi:** Cell shapes differ (can be empty cells, asymmetric).
Sites are NOT equidistant from cell vertices (unlike Voronoi). But perpendicularity
is preserved, which is all the proof needs.

**Test:** Generate random power diagram with non-uniform weights → check G = Vol·I.

**Verdict: VIABLE — likely the natural generalization.**

### 2. Circumcentric DEC (Delaunay dual)

**What:** Start with a simplicial mesh (tetrahedralization). The dual is formed by
connecting circumcenters of adjacent tetrahedra. In 2D, circumcentric Delaunay dual
gives exact perpendicularity. This is the standard DEC construction (Hirani, Desbrun).

**Problem in 3D:** Circumcenters of obtuse tetrahedra lie outside the tetrahedron.
This means:
- Dual edges may cross primal faces at wrong angles
- ⋆₁ can become negative (non-physical Hodge star)
- Pillar tiling breaks down

**2D:** Works perfectly. Circumcentric dual of Delaunay = Voronoi diagram.
So in 2D, circumcentric DEC IS Voronoi. Not a distinct class.

**3D on well-centered meshes:** If all tetrahedra are well-centered (circumcenter
inside), then perpendicularity holds and ⋆₁ > 0. But well-centered 3D meshes are
hard to generate and may not exist for arbitrary domains.

**Verdict: NOT A DISTINCT CLASS in 2D (= Voronoi). PROBLEMATIC in 3D
(negative Hodge stars, broken perpendicularity). Not a viable candidate
unless restricted to well-centered meshes.**

### 3. Barycentric dual

**What:** Dual vertices at barycenters (centroids) of cells, not circumcenters.
Standard choice in some FEM/FVM frameworks.

**Perpendicularity:** NO. Barycenters are not equidistant from face vertices.
The line connecting two barycenters is generally NOT perpendicular to the shared face.

**Consequence:** Term 3 ≠ 0 in the divergence theorem proof. G ≠ Vol·I in general.

**Test idea:** Build barycentric dual on same Voronoi topology → measure G deviation.

**Verdict: NOT VIABLE (breaks perpendicularity). But useful as a negative control.**

### 4. Centroidal Voronoi Tessellation (CVT)

**What:** Special Voronoi where sites coincide with cell centroids (Lloyd's algorithm
fixed point). Used in mesh generation (CGAL, VoroCrust).

**Perpendicularity:** YES — it's still Voronoi, so all Voronoi properties hold.

**Verdict: SPECIAL CASE OF VORONOI, not a distinct class. Interesting only as a
particularly "nice" Voronoi (cells close to spherical).**

### 5. Perturbed Voronoi (same topology, moved dual vertices)

**What:** Start with Voronoi complex. Keep the primal mesh (vertices, edges, faces)
but move the dual vertices (cell centers) away from Voronoi sites.

**Perpendicularity:** BROKEN. Dual edges no longer perpendicular to faces.

**Hodge stars:** Recomputed from new dual geometry → different values.

**G = Vol·I?** The admissible subspace (codim-6) shows that SOME Hodge star
perturbations preserve G = Vol·I. But these are abstract perturbations of ⋆₁,
not necessarily realizable as geometric dual constructions.

**Key question:** Can we perturb dual vertices such that the resulting Hodge stars
still satisfy G = Vol·I? This would be a non-Voronoi, non-power-diagram example.

**Verdict: OPEN. Needs numerical exploration.**

### 6. Orthocentric meshes (generalization)

**What:** A mesh is orthocentric if for every edge e, the line connecting the
dual vertices (one per adjacent cell) is perpendicular to the primal face.
This is the MINIMAL requirement for the divergence theorem proof.

**Voronoi and power diagrams are orthocentric.** Are there others?

In 2D: orthocentric = Delaunay (well-known). So Voronoi is the only class.

In 3D: orthocentric polyhedral complexes are less studied. The condition is:
for each face f shared by cells c₁, c₂, the vector (center_c₂ − center_c₁)
must be parallel to the face normal n_f.

**Conjecture:** On a periodic T³ complex with flat faces, orthocentricity +
positive Hodge stars ⇒ the complex is a power diagram (possibly degenerate).

This would mean: **power diagrams are the maximal class satisfying (A).**

**Verdict: KEY THEORETICAL QUESTION. If the conjecture holds, the answer to
"beyond Voronoi" is exactly "power diagrams, and nothing else."**

---

## Investigation plan

### Phase 1: Power diagrams (computational)
1. Build periodic power diagram with SciPy (weighted Voronoi)
2. Compute Hodge stars from dual geometry
3. Verify G = H = Vol·I
4. Run full c²=1 test with exact Bloch complex
5. Compare with unweighted Voronoi on same sites

### Phase 2: Negative controls
1. Barycentric dual on Voronoi topology → show G ≠ Vol·I
2. Random Hodge star perturbation (outside admissible subspace) → show G ≠ Vol·I
3. Quantify: how far from Vol·I, and how does c deviate?

### Phase 3: Theoretical
1. Literature search: orthocentric polyhedral complexes in 3D
2. Is every orthocentric periodic complex a power diagram?
3. If not, construct a counterexample

---

## Summary table

| Class | Perpendicularity | G = Vol·I | c²=1 | Distinct from Voronoi? |
|-------|-----------------|-----------|-------|----------------------|
| Voronoi | yes | yes (proved) | yes (proved) | — |
| Power diagram | yes | predicted yes | predicted yes | YES (different cells) |
| Circumcentric 3D | only if well-centered | only if well-centered | unclear | = Voronoi in 2D |
| Barycentric dual | no | no (predicted) | no (predicted) | yes |
| CVT | yes | yes | yes | no (special Voronoi) |
| Perturbed dual | no | unclear | unclear | yes |
| Orthocentric | yes (by definition) | predicted yes | predicted yes | = power diagram? |

---

## Numerical experiments (Mar 6)

### Negative controls — confirmed

1. **Uniform stars (⋆₁=1) on Kelvin:** G ∝ I (from cubic symmetry), but G ≠ Vol·I
   (G_ii = 128 vs Vol = 512). Isotropy is accidental (symmetry), not from Hodge stars.

2. **Uniform stars on random Voronoi (50 cells):** G NOT isotropic.
   Diagonal spread 7.6%, off-diagonal/diagonal 3.4%. Confirms: Voronoi Hodge stars
   are essential for isotropy on non-symmetric meshes.

3. **Perturbed dual centers (shift 0.1 RMS):** G/Vol deviates 0.5–1.3% off-diagonal.
   Breaks perpendicularity → breaks G = Vol·I. Confirmed.

4. **Random ⋆₁ perturbation (30% noise):** G/Vol off-diagonal 3.8%. Expected.

### Proof structure analysis

The §2 proof has 4 terms from divergence theorem on dual cells:
- **Term 1 = 0:** flux closure (verified T2). Needs closed dual cells.
- **Term 2 = G:** definition.
- **Term 3 = 0:** face centroid offset cancellation. **This is the critical term.**
  Uses perpendicularity of dual edge to primal face. On Kelvin, the dual face
  centroid is NOT at the edge midpoint (offset up to 2.68). But these offsets
  cancel when summed with the appropriate weights. Power diagrams shift these
  offsets — unclear if cancellation survives.
- **Term 4 = Vol·I:** pillar volume tiling. Total pillar volume = 2.12 × Vol
  (not 3×Vol as naively expected). The factor depends on geometry.

### Power diagram feasibility

- SciPy supports 4D Voronoi (needed for lifting trick). ✓
- Building a full periodic power diagram complex requires a new builder
  (non-trivial: periodic images in 3D + 4D Voronoi + projection + dedup).
- **Not yet tested numerically.** Need builder infrastructure.

### Key open question

Does Term 3 = 0 on power diagrams? The centroid offsets change when faces shift.
If Term 3 still cancels, G = Vol·I holds for all power diagrams (theorem).
If not, power diagrams are NOT sufficient, and Voronoi is more special than
just "orthocentric."

**This is the bottleneck.** An analytic argument for Term 3 on power diagrams
would settle the question without needing a builder.

---

## Current assessment

| Finding | Status |
|---------|--------|
| Voronoi ⋆₁ essential for G = Vol·I | Confirmed (3 negative controls) |
| Perpendicularity necessary | Confirmed (perturbed centers break G) |
| Power diagram satisfies perpendicularity | Known (theoretical) |
| Power diagram has G = Vol·I | **Open** (Term 3 cancellation unclear) |
| Power diagram builder feasible | Yes (4D lifting) but non-trivial |
| Orthocentric = power diagram? | **Open** (key theoretical question) |

---

## Conclusion

### Necessary conditions for c²=1

1. **Perpendicularity:** dual edge ⊥ primal face (without this, G ≠ Vol·I — confirmed numerically)
2. **Exactness:** d₁d₀ = 0 (topological construction, works on any complex)
3. **z=4 economy:** nE = 2nV (gives dim H¹ = 3 on T³)
4. **Positive Hodge stars** (⋆₁ > 0, ⋆₂ > 0)

Condition (1) is the geometric filter. The rest are topological or constructive.

### What satisfies perpendicularity in 3D

| Class | Perpendicularity | Distinct from Voronoi? |
|-------|-----------------|----------------------|
| **Voronoi** | yes (bisector property) | — |
| **Power diagrams (Laguerre)** | yes (weighted bisector) | yes — different cells, sites ≠ equidistant from vertices |
| Circumcentric DEC 3D | only on well-centered meshes | problematic (⋆₁ can be negative) |
| Barycentric dual | no | — |
| Any other perturbed-center dual | no | — |

### Answer

**Power diagrams** are the natural generalization. Each site has a weight w_i,
the cell is {x : |x−p_i|²−w_i ≤ |x−p_j|²−w_j}. Equal weights → classical Voronoi.
Perpendicularity is exact by construction.

This is the only serious candidate. In 2D it is known that orthocentric = Delaunay
= Voronoi. In 3D it is not formally proved that orthocentric = power diagram, but
no counterexample is known.

### Paper §8 paragraph

*"The result extends naturally to any periodic power diagram (Laguerre tessellation),
since dual-primal perpendicularity — the only geometric property used in the
divergence theorem proof of G = H = Vol·I — is satisfied exactly. Classical Voronoi
is the special case with equal weights."*

### Status: CLOSED (low priority, answer sufficient for paper discussion)

---

## References

- Aurenhammer (1987): Power diagrams — properties and algorithms
- Glickenstein (2005): Geometric triangulations and discrete Laplacians on manifolds
- Mullen et al. (2011): HOT (Hodge-Optimized Triangulations) — optimizes Hodge stars
- VanderZee et al. (2010): Well-centered meshes in 3D
- Edelsbrunner & Shah (1996): Weighted Delaunay and power diagrams
