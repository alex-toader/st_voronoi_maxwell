# W17 Paper Plan: Geometry and topology of discrete Maxwell wave speed

## One-sentence summary

On any periodic Voronoi complex with an exact cochain complex, the discrete
electromagnetic wave speed equals 1 in natural units (c² = 1 + O(k²) in the
dispersion relation). Two conditions are sufficient: the metric identity
G = H = Vol·I (geometry, from Voronoi Hodge stars) and the cochain exactness
d₁d₀ = 0 (topology, from recurrence construction). This is a geometric identity
at finite resolution, not a convergence statement. Numerical evidence on 3 cubic
and 4+ random meshes shows both conditions are necessary: removing either gives
c ≠ 1 with O(1) error — a different operator, not a worse approximation.

## Framing

Two independent properties control c:
- **(A) Metric identity:** G = H = Vol·I (geometry, from Voronoi Hodge stars)
- **(B) Topological exactness:** d₁d₀ = 0 (topology, from recurrence construction)

Voronoi is sufficient for (A), not necessary. The JCP paper's construction provides (B).
Together → c = 1. Errors are separable and quantifiable.

## Structure

### §1 Introduction
- c = 1 on foam — what it means, why it matters
- Context: DEC Maxwell on periodic polyhedral complexes
- Two ingredients (A) and (B) — neither alone suffices
- Relation to JCP paper: complementary (JCP proves d₁d₀ = 0, this paper proves c = 1)

### §2 Metric identities G = H = Vol·I
- **Theorem:** On any periodic Voronoi tessellation, G_ij = Σ_e ⋆₁[e] Δx_iΔx_j = Vol·δ_ij
- **Proof:** Divergence theorem on each Voronoi cell (Voronoi perpendicularity → boundary
  integral telescopes → pillars tile space isotropically).
  Proof structure: 4 terms. Term 1 = 0 (per-vertex flux closure: Σ_{e∋v} sign·⋆₁·Δx = 0
  on each closed Voronoi cell). Term 2 = G (definition). Term 3 = 0 (Voronoi symmetry:
  face centroid offset cancellation). Term 4 = Vol·I (telescoped boundary integral).
  Term 1 = 0 verified directly (T2): cubic 0–4e-10, random < 1.2e-06.
- Same proof for H_ij = Σ_f ⋆₂[f] (A_f)_i(A_f)_j = Vol·δ_ij (dual complex)
- **Perpendicularity (pre-condition):** Voronoi dual edge ∥ primal face normal verified
  to machine precision: |cos(dual, face)| − 1 < 2e-16 (cubic), < 3.4e-14 (random).
  Uses face_to_cell_shift for correct periodic dual edge vectors.
- Numerical verification: cubic (10⁻¹⁷), random Voronoi (10⁻¹²), degenerate (10⁻¹²).
  Stable under domain scaling (L = 2, 4, 8) and extreme cell ratio (40×).
- **Remark (tight frame):** The weighted edge vectors f_e = √⋆₁[e]·Δx_e form a tight
  frame for ℝ³ with frame bound A = Vol: F†F = Vol·I, where F is the |E|×3 frame
  matrix. All 3 singular values of F equal √Vol (redundancy ratio |E|/3 = 64–725).
  This is the Parseval-type condition: G = Vol·I is equivalent to isotropic quadrature
  of ℝ³ by weighted edge directions. Connects to frame theory and sampling theory.
- **Admissible perturbations:** The constraint G = Vol·I is Σ_e δ(⋆₁[e])·(Δx_e)_i(Δx_e)_j = 0,
  which is 6 equations (symmetric 3×3 matrix) on |E| unknowns.
  **Proposition:** The tensors {Δx_e ⊗ Δx_e : e ∈ E} span S²(ℝ³) for any periodic 3D
  complex with vertices of degree ≥ 4 in generic position. (Proof: S²(ℝ³) is
  6-dimensional. Each Δx_e ⊗ Δx_e is rank-1 symmetric. On z=4 meshes, each vertex
  contributes 4 non-coplanar edge directions; with ≥ 2 vertices one obtains ≥ 6
  linearly independent rank-1 tensors spanning S².) Therefore rank = 6,
  giving codimension-6 submanifold in ⋆₁-space ≅ ℝ^|E|. Verified: rank = 6 on all
  structures. Voronoi is one point; |E|−6 directions of freedom preserve c = 1.
- **Admissible cone radius:** The admissible subspace intersected with ⋆₁ > 0 forms a
  cone. Voronoi Hodge stars sit deep inside: 1.35× mean(⋆₁) (Kelvin), 5.95× (C15),
  2.59× (WP). The identity is robust, not marginal — Voronoi is far from the positivity
  boundary of the admissible set.

### §3 Theorem: c² = 1

#### §3.1 Setup and statement
- **Setup:** K = d₁†⋆₂d₁, mass M = ⋆₁, Bloch boundary conditions at wavevector k.
- **Harmonic 1-forms:** At Γ, h_e = α·Δx_e/√Vol spans H¹(T³) ≅ ℝ³. M-orthonormal
  from G = Vol·I (§2). At k ≠ 0, H¹(k) = 0 (z=4 economy: nE = 2nV → dim H¹ = 0).
  So the 3 harmonic modes deform into 1 gradient + 2 acoustic.
- **H² cohomology:** Face area vectors h₂ = A_f/√Vol span H²(T³) ≅ ℝ³.
  S₂-orthonormal from H = Vol·I (§2).
- **Theorem:** On any periodic Voronoi complex with exact d₁d₀ = 0, the acoustic
  dispersion satisfies ω² = k² + O(k⁴), giving c² = 1 in natural units.

#### §3.2 Proof (Schur complement + discrete curl identity)
Three ingredients:
  (1) H = Vol·I → h₂†S₂h₂ = I (face tensor, §2)
  (2) d₁(0)·harm = 0 (exactness: harm ∈ ker d₁ at Γ)
  (3) Discrete curl identity: ∂(h₂†S₂d₁(k)harm)/∂k|₀ = i·ε_{βγα}

**Proof chain.** Since d₁(0)·harm = 0, the curl of harmonic forms at k is
u = d₁(k)·harm = Δd₁·harm, living in face space (nF-dimensional). Decompose
u into im(d₁(0)) and its S₂-orthogonal complement. The complement has dimension
nC + 2 (from V−E+F = C on T³), containing H² (3-dim) and im(d₂†) (nC−1 dim).
At O(k), u_perp lies entirely in H² — the im(d₂†) component starts at O(k²)
(verified: non-H² fraction ∝ k², scaling exponent 2.00 on all structures).
The H² coefficient is C_{βα} = h₂_β†S₂·u_perp_α = i·ε_{βγα}k_γ + O(k²),
from ingredient (3). Then:

  S = u_perp†S₂u_perp = C†(h₂†S₂h₂)C + O(k⁴)
    = C†IC = [k×]ᵀ[k×] + O(k⁴) = k²P_T + O(k⁴).

The Schur complement S = H₃ − B·K_opt⁻¹·B† equals u_perp†S₂u_perp at leading
order (optical modes span im(d₁(0)) at Γ), giving eigenvalues [0, k², k²].
The two non-zero eigenvalues give the transverse acoustic speed c² = 1.  □

Verified on 3 cubic + 4 random Voronoi (file 6: R44f, R44g):

  | Structure | H₃ eigenvalues/k² | Schur eigenvalues/k² | ||dF/dk/i − ε|| |
  |-----------|-------------------|---------------------|-----------------|
  | Kelvin    | 0.8, 4.8, 5.3    | 0, 1.000, 1.000     | 1.1×10⁻¹¹       |
  | C15       | 5.8, 9.4, 10.7   | 0, 1.000, 1.000     | 1.5×10⁻¹¹       |
  | WP        | 4.4, 6.7, 7.0    | 0, 1.000, 1.000     | 1.5×10⁻¹¹       |
  | Rnd(42)   | 7.3, 10.5, 13.1  | 0, 1.000, 1.000     | 3.0×10⁻¹¹       |
  | Rnd(7)    | 6.1, 8.8, 14.2   | 0, 1.000, 1.000     | 2.8×10⁻¹¹       |

#### §3.3 Anatomy: why naive approaches fail (supporting material)
  The Rayleigh quotient on harmonic forms gives R[h]/k² = 4.84–10.75, NOT 1
  (R44b). The eigenvector IS harmonic to 99.999% (R44a), but the Γ-harmonic
  subspace is not K(k)-invariant. The eigenvalue emerges from three-way
  cancellation (R44d):

  | Structure | h†Kh/k² | 2Re(h†Kδψ)/k² | δψ†Kδψ/k² | sum/k² |
  |-----------|---------|---------------|-----------|--------|
  | Kelvin    |   4.833 |        −7.667 |     3.833 |  1.000 |
  | C15       |   9.392 |       −16.783 |     8.392 |  1.000 |
  | WP        |   6.700 |       −11.399 |     5.700 |  1.000 |

  What does NOT work: (a) direct Rayleigh on harmonics, (b) gradient projection
  makes R worse (R44c), (c) full degenerate PT on nV+3 zero subspace gives c²~25.
  The Schur complement resolves the paradox: it integrates out optical modes
  via block elimination, not perturbation theory.
- **Converse lemma (c = 1 ⇒ metric isotropy):** If lim_{|k|→0} ω²/|k|² = 1
  for all transverse polarizations α ⊥ k̂ and all directions k̂, then G ∝ I
  and H ∝ I with the same constant.
  *Proof:* From the Rayleigh quotient, (k×α)ᵀH(k×α) = |k|²·αᵀGα for all k, α⊥k.
  Choose k = e₃, α = e₁: H₂₂ = G₁₁. Choose k = e₃, α = e₂: H₁₁ = G₂₂.
  Choose k = e₁, α = e₂: H₃₃ = G₂₂. Choose k = e₁, α = e₃: H₂₂ = G₃₃.
  Combined: G₁₁ = G₂₂ = G₃₃ ≡ λ, H₁₁ = H₂₂ = H₃₃ = λ.
  For off-diagonals: k = e₃, α = (e₁+e₂)/√2 gives H₁₂ = 0.
  Similarly k = (e₁+e₂)/√2, α = e₃ gives G₁₂ = 0. All off-diagonals vanish.
  Therefore G = H = λI. (λ = Vol follows from Voronoi, but c = 1 alone only
  fixes the ratio H/G = I, not the absolute scale.)  □
  **Note:** This proves "only if" for condition (A). For condition (B), d₁d₀ = 0
  is numerically necessary (all tested violations give c ≠ 1) but not proved
  analytically — a formal proof would require showing that leaked gauge modes
  always contaminate the lowest acoustic eigenvalue.
- **O(k²) correction:** δc/k² ≈ −0.104 (Kelvin), −0.022 (C15), −0.054 (WP).
  Structure-dependent, from higher-order perturbation theory (3rd+ order in k).
  The optical gap ω²_gap varies only 10% (0.54–0.60) while |δc/k²| varies 5×, so the
  gap alone does not predict the correction. δc/k² is encoded in a rank-4 anisotropy
  tensor T₄ that depends on the full mode structure, not just the gap.
- **Remark:** The identity G = Vol_d · I_d holds on any periodic Voronoi in ℝ^d
  by the same divergence theorem argument. The result is dimension-independent.

### §4 Necessity: removing exactness
**Body (keep short, ~2 pages):**
- Standard DEC (same Hodge stars, d₁d₀ ≠ 0) gives c = 1.25–1.68
- **No principal symbol:** RMS residual 0.38–0.74 on c²(k̂) = k̂ᵀCk̂ fit. Show
  one example C_std (Kelvin) + 2–3 n_zero classes. This is not "C ≠ I" — no C exists.
- **Catastrophic anisotropy:** 37–70% on isotropic medium.
- **Gauge pollution mechanism (R41):** Leaked forms are 100% pure gradients.
  H¹(k≠0) = 0 means ker(d₁_ex) = im(d₀). Standard DEC breaks d₁d₀ = 0, so
  gradient modes leak through curl → become spurious non-zero eigenvalues.
  n_lost = Δrank(d₁) (R40).

**Appendix A (detailed anatomy):**
- Mode mixing: lowest mode = 89% gauge + 11% acoustic (R28)
- N-dependence: c_std non-monotonic with supercell size (R20)
- TRIM collapse: standard exact at k = π/L_cell (phases ±1 → d₁d₀ = 0) (R27)
- Multi-k stability analysis for principal symbol

### §5 Necessity: removing geometry
- Exact complex + perturbed Hodge stars → c ≠ 1 but no pollution (n_lost = 0)
- c_EP can be > 1 (WP: 1.006) or < 1 (Kelvin: 0.979) — perturbation doesn't force a direction
- **Error factorization:** Δc ≈ Δc_geom + Δc_topo. Interaction = O(amplitude):
  cubic 0.025 at 10%, 0.037 at 15%; random < 0.01 at 10%.
  Separability is first-order approximate, good to ~15% perturbation amplitude.
- Random Voronoi interaction tighter than cubic (< 0.01 vs ~0.03) — no symmetry
  to accidentally force separability, yet it holds better.
- c_std on random can be subluminal (0.47) or superluminal (1.02) — wildly variable.
- Topology controls n_lost, geometry primarily controls speed. Approximately additive
  (not exactly separable — interaction < 5% is empirical, not proved).
- **Proposition (for paper):** ker(K) = ker(d₁) for any S₂ > 0, since K = d₁†S₂d₁.
  This is why n_lost is exactly independent of Hodge stars — a mathematical fact, not
  a numerical observation. Formulate as proposition in §5.

### §6 Spectral structure of the error
- **Moment hierarchy:** tr(K^n) conserved for n ≤ 2 on both complexes, any ε.
  Because |d₁[f,e]|² = 1 on both (phases vs ±1, same modulus).
- **Girth theorem (formal proof below):** Break at n = girth(face adjacency graph) = 3.

**Theorem (Moment conservation).** Let K_α = D_α† S D_α for α ∈ {ex, std}, where
S = diag(s_f) > 0 is the Hodge star ⋆₂, and D_α is the F×E curl matrix with the
same sparsity pattern and |D_α[f,e]| = 1 on support (both complexes have unit-modulus
entries). Let g = girth(G_F), the face adjacency graph (faces connected iff they share
a mesh edge). Then:

  tr(K_ex^n) = tr(K_std^n)  for all 1 ≤ n ≤ g − 1.

**Proof.** Expand:

  tr(K^n) = Σ_{e_1,...,e_n} Σ_{f_1,...,f_n} [Π s_{f_i}] · [Π D*_{f_i,e_i} D_{f_i,e_{i+1}}]

where e_{n+1} = e_1, and f_i must contain both e_i and e_{i+1}.

Consider the bipartite incidence graph B = (E ∪ F, ~) where e ~ f iff face f contains
edge e. Each term in the sum corresponds to a closed walk of length 2n in B:

  e_1 → f_1 → e_2 → f_2 → ··· → e_n → f_n → e_1

The phase product assigns: D*_{f,e} to each E→F step, D_{f,e} to each F→E step.

Since B is bipartite, girth(B) = 2g. For n < g, the walk length 2n < 2g, so the walk
contains no cycle (any cycle would have length ≤ 2n < girth(B), contradiction) — it
is a tree walk. Every closed walk on a tree must backtrack each edge the same number
of times in each direction (otherwise the walk cannot return to its starting vertex).
So each B-edge (f,e) is traversed k times E→F and k times F→E, contributing
(D*_{f,e})^k (D_{f,e})^k = |D_{f,e}|^{2k} = 1. Since this holds for every term, and the s_f factors
and sparsity pattern are identical on both complexes, tr(K_ex^n) = tr(K_std^n).

At n = g: a simple cycle of length 2g exists in B, corresponding to a g-cycle
f_1 → f_2 → ··· → f_g → f_1 in G_F through shared edges e_1,...,e_g. The holonomy
Π_{i=1}^g D_{f_{i-1},e_i}/D_{f_i,e_i} depends on the actual Bloch phases (recurrence
vs shift), which differ between the two constructions. So tr(K^g) generically differs.  □

**Verified:** Kelvin, C15, WP — all have girth(G_F) = 3, edge valence 3 (uniform).
tr(K^n) conserved for n=1,2 (rel_diff < 10^{-10}), breaks at n=3 by 0.13–0.22%.

  **Remark:** Connection to Ihara zeta function / non-backtracking walks on face
  adjacency graph. The moment conservation is a statement about the Hashimoto
  matrix: two operators with same modulus pattern share traces of T^n for n < girth.
- **Consequence:** Spectral difference is redistribution, not creation.
  Same total (tr K), same variance (tr K²), different shape starting at order 3.
- **Direction dependence of pollution:** n_lost ~ #{edges with k·Δx ≠ 0}.
  Axis: 6, face diagonal: 12, body diagonal: 14 (Kelvin).

### §7 Dielectric extension
- **Setup:** ε is piecewise constant per cell. For each face f shared by cells c₁, c₂
  with ε(c₁), ε(c₂), the modified Hodge star is ⋆₂^ε[f] = ⋆₂[f] / ε_eff(f).
- **Harmonic mean prescription:** ε_eff(f) = 2/(1/ε(c₁) + 1/ε(c₂)), i.e.,
  1/ε_eff = (1/ε(c₁) + 1/ε(c₂))/2 (arithmetic mean of inverse permittivities).
  This corresponds to the correct discretization of ∇·(ε⁻¹ ∇×) at the interface.
- **Result:** H_ε = ⟨1/ε⟩_vol · Vol · I, where ⟨1/ε⟩_vol = (1/n_cells)Σ_c 1/ε(c)
  is the volume-weighted mean of 1/ε (equal cell volumes assumed).
- Exact on Kelvin at any contrast (per-cell O_h isotropy → cross-term = 0)
- Trace-exact on random Voronoi (vol-weighted), off-diagonal O(1/√N)
- **Harmonic vs log mean:** Harmonic preserves metric identity. Log mean breaks it
  by 7–33% depending on contrast. Different optimization targets.

### §8 Discussion
- c = 1 as geometric + topological constraint (not fitted, not emergent)
- Practical: error budget from separability (metric quality + topological quality)
- Relation to continuous theory: Voronoi DEC has trivially homogenized metric
- **Beyond Voronoi (from W2 investigation):** The proof uses 4 conditions:
  (A) G=H=Vol·I, (B) d₁d₀=0, (C) dim H¹=3, (D) dF/dk=iε. Only (A) is geometric;
  (B) is topological (works on any complex), (C) is combinatorial (z=4), (D) follows
  from (A)+topology. For (A), the divergence theorem proof needs perpendicularity
  (dual edge ⊥ primal face). Power diagrams (Laguerre tessellations) satisfy this
  exactly — natural candidate generalization. However, Term 3 cancellation (face
  centroid offsets) is NOT verified for power diagrams. Formulate cautiously:
  *"The result likely extends to periodic power diagrams, since perpendicularity —
  the key geometric input to the metric identity proof — holds exactly. Verification
  of the full identity on power diagrams is left to future work."*
  Note: no numerical verification possible without power diagram builder (non-trivial).

## The killer table (§5)

| Config | d₁d₀ = 0 | G = Vol·I | c | n_lost | aniso | tr(K^n≤2) |
|--------|-----------|-----------|-------|--------|-------|-----------|
| Exact + Voronoi | ✓ | ✓ | 1.000 | 0 | < 10⁻⁴ | conserved |
| Exact + Perturbed | ✓ | ✗ | 0.993 | 0 | ~1% | conserved |
| Standard + Voronoi | ✗ | ✓ | 1.253 | 6 | 70% | conserved |
| Standard + Perturbed | ✗ | ✗ | 1.245 | 6 | 71% | conserved |

Separability visible at a glance. Each row changes one thing.

## What's proved vs what's verified numerically

| Result | Status | What's needed |
|--------|--------|---------------|
| G = H = Vol·I | Proved (divergence theorem) + verified | Write proof cleanly |
| c² = 1 + O(k²) | Proved: S = k²P_T + O(k⁴) via (1) H=Vol·I, (2) exactness, (3) discrete curl identity dF/dk = iε (R44e-g). dim complement = nC+2, u_perp ∈ H² at leading order, im(d₂†) component O(k²). Verified 3 cubic + 4 random. | Write proof as §3.2 |
| c = 1 ⇒ G ∝ H ∝ I | Proved (converse lemma, 6-direction argument) | — |
| Codim-6: Δx⊗Δx span S² | Proved (3 independent directions + cross terms) | — |
| d₁d₀ = 0 → ker K correct | Proved in JCP paper | Reference |
| H¹(k≠0) = 0 on z=4 | Proved (nE=2nV + rank identity) + verified 8 meshes | — |
| Leaked forms ⊂ im(d₀) | Proved (H¹=0 ⇒ ker d₁=im d₀) + verified 19 cases | — |
| n_lost = Δrank(d₁) | Exact identity + verified 15 combinations | — |
| d₁d₀ ≠ 0 ⇒ c ≠ 1 | Numerical only (6 cubic + 5 random, all give c ≠ 1) | Open |
| tr(K^n) conserved n ≤ g−1 | Proved (bipartite tree walk argument) | — |
| Break at n = girth | Proved (bipartite walk + tree cancellation) | — |
| Error factorization | Verified numerically (3 cubic + 3 random Voronoi, R36/R39) | — |
| Harmonic preserves metric | Proved (direct substitution) | Write proof |
| TRIM collapse | Proved (phases ±1 → products cancel) | Write proof |
| δc = O(1) | Verified numerically | — |

## Tests backing each section (new file structure)

| Section | Test file | Key tests | Key results |
|---------|-----------|-----------|-------------|
| §2 | 1_test_metric_identity.py (12 tests) | R1, R2, R3, R4, R34, N12, T2, T4, T5, T7 | G = H = Vol·I (10⁻¹⁷ cubic, 10⁻¹² random); perpendicularity 10⁻¹⁴; codim-6; tight frame; cone radius 1.4–6×; scaling invariance |
| §3 | 2_test_c_squared_one.py (10 tests) | R5, R38, R32, R41, R44e, R44f, R44g, R6, S1, S2 | c = 1 on all Voronoi; Schur [0, k², k²]; discrete curl dF/dk = iε (1e-8); converse c≠1 (1929×); isotropy 5 dirs (spread 1e-6); dispersion O(k⁴) |
| §4 | 3_test_removing_exactness.py (7 tests) | R12, R25, R35, R17, R40, R41, R23 | c_std = 1.25–1.68; no principal symbol; anisotropy 37–70%; n_lost = Δrank; leaked ⊂ im(d₀) |
| §5 | 4_test_removing_geometry.py (3 tests) | R36, R39, R29 | Killer table cubic (interact < 0.03) + random (interact < 0.01); factorization O(amp) at 5/10/15%; c_std random subluminal 0.47–superluminal 1.02 |
| §6 | 5_test_spectral_structure.py (3 tests) | R14, R22, R21 | tr(K^n) conserved n≤2; girth = valence = 3 |
| App A | 6_test_appendix_a.py (5 tests) | R44a-d, R11 | Eigvec 99.999% harmonic; Rayleigh paradox; 3-way cancellation |
| App B | 7_test_appendix_b.py (3 tests) | R16, R20, R27 | Mode mixing 89%+11%; N non-monotonic; TRIM collapse |

## What to do before writing

1. ~~**Close the c = 1 proof:**~~ DONE (R44e-g). Analytic proof complete. Three ingredients:
   (1) H = Vol·I from divergence theorem, (2) d₁(0)·harm = 0 from exactness,
   (3) discrete curl identity ∂(h₂†S₂d₁(k)harm)/∂k|₀ = i·ε (Levi-Civita tensor).
   Chain: C = i·[k×] → S = C†C = [k×]ᵀ[k×] = k²P_T → c² = 1.
   Verified to 10⁻¹² on 3 cubic + 4 random Voronoi (7 tests in file 6).

2. ~~**Extract C_std tensor (N16):**~~ DONE (file 5). No consistent principal symbol —
   RMS residual 0.38–0.74 on c²(k̂) = k̂ᵀCk̂ fit. Negative eigenvalue on body diag class.

3. ~~**Verify factorization on C15/WP (N11/N15):**~~ DONE (file 5, R36). Killer table
   on all 3 cubic structures + 3 random Voronoi (R39). Interaction < 0.03.

4. ~~**Admissible Hodge star subspace (N12):**~~ DONE (file 5, R37). Codimension-6
   confirmed. 186–2170 admissible directions. c = 1 at leading order not fragile.

5. ~~**Write the girth theorem proof (C4):**~~ DONE. Formal proof in §6 above.
   Bipartite walk argument, tree cancellation, holonomy at n = girth.

6. **[EDITORIAL] Balance §4 and §5:** Concrete structure:
   - **§4 body:** no principal symbol (N16) + killer table (R36) + δc = O(1) (R25)
   - **Appendix A:** mode mixing anatomy (R16/R28), N-dependence (R20), TRIM collapse (R27)
   This gives §4 and §5 equal weight (both ~2 pages body + appendix support).

7. **[EDITORIAL] "If and only if" qualification:** The "if" direction is proved (§3).
   The "only if" has strong numerical evidence but no proof for condition (B).
   Converse lemma (§3) proves "only if" for (A): c = 1 ⇒ G ∝ H ∝ I.
   **Decision:** "sufficient, and numerical evidence strongly suggests necessary."

8. **[EDITORIAL] Random Voronoi c_std < 1:** R30 shows c_std ∈ [0.61, 1.04].
   Standard can be subluminal OR superluminal. Mention in §4 or §8.

9. **[EDITORIAL] u_perp ∈ H² as Lemma (verified numerically):** The statement
   "at O(k), u_perp lies entirely in H²" is verified (exponent 2.00 on all structures)
   but not proved analytically. Formulate as **Lemma (verified numerically)**, not as
   part of the theorem statement. Does not weaken c²=1 — the Schur complement
   already gives c²=1 independently. (Two reviewers converge on this.)

10. **[EDITORIAL] "Already homogenized" framing:** Voronoi DEC with exact complex
    is already homogenized at any resolution — c=1 is not a continuum limit, not an
    effective medium result. Add one clear sentence in §1 (Introduction) and revisit
    in §8 (Discussion). Conceptually powerful: no homogenization needed.

11. **[EDITORIAL] Dielectric = separate paper:** §7 (dielectric extension) adds weight
    without strengthening the main message (c²=1 as geometric identity). Move to a
    separate note or paper 3. Keeps this paper focused. (Two reviewers converge.)

12. **[EDITORIAL] Strategic decision A vs B:**
    - **A:** JCP companion (technical, comprehensive, all sections). Safe, natural.
    - **B:** Conceptual statement ("discrete EM has exact c on Voronoi"). Less material,
      more theorem–corollary structure, broader audience (PRL or Phys Rev E).
    Decision needed before writing. Affects tone, length, and what stays in body vs appendix.

13. **[WRITING] Analytic proof of discrete curl identity (L1 / N30):** The chain in §3.2
    uses dF/dk|₀ = iε as ingredient (3), currently verified numerically (R44f, 10⁻¹¹ on
    7 structures). The analytic proof: ∂d₁(k)/∂k_γ|₀ introduces i·(shift_e·L)_γ on each
    edge. The sum ∑_{e∈∂f} d₁(0)[f,e]·(shift_e·L)_γ·(Δx_e)_α is a discrete Stokes
    integral: shift×Δx gives parallelogram area, summed over face boundary = face area
    vector (A_f)_β. Contracted with S₂·h₂ (= A_f/√Vol normalized by H=Vol·I), gives
    ε_{βγα}. This is the ONLY non-proved step in §3.2. Writing it closes the proof.

14. **[WRITING] Schur complement = u_perp†S₂u_perp (L3):** §3.2 asserts the equivalence
    without proof. Argument: H₃ = harm†K(k)harm = u†S₂u (since K = d₁†S₂d₁).
    B·K_opt⁻¹·B† is the S₂-projection of u onto im(d₁(0)). Therefore
    S = H₃ − B·K_opt⁻¹·B† = u†S₂u − u_exact†S₂u_exact = u_perp†S₂u_perp.
    Standard block matrix argument, needs to be written explicitly in §3.2.

## Venue

- **PRL:** Too much material for 4 pages. Would need severe cutting (just §2-3).
- **JCP companion:** Natural home. Same journal as paper 1. Complementary results.
- **Phys Rev E:** If framed as computational physics methodology.
- **J. Comput. Phys.** companion paper is the safest and most natural choice.

## Relation to JCP paper 1

| | JCP paper 1 | This paper (W17) |
|--|-------------|-----------------|
| Proves | d₁d₀ = 0 (construction) | c = 1 (identity) |
| Uses | Recurrence on face boundary | Divergence theorem + Schur complement + discrete curl identity |
| Shared tool | Exact Bloch complex | Exact Bloch complex |
| Focus | Topological (exactness) | Geometric + topological (speed) |
| Key table | Pollution counts, band shifts | 2×2 factorization table |
| Overlap | Minimal — different theorems | |
