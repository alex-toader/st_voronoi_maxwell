# Tests Map — Paper 3 (st_voronoi_maxwell)

**Date:** Mar 2026

Each test file maps to a paper section. Tests live in ONE file only.

---

## 1_test_metric_identity.py (§2) — DONE, 12 tests

| Test | ID | What it verifies |
|------|----|-----------------|
| test_R1_edge_tensor_cubic | R1 | G = Vol·I on Kelvin, C15, WP |
| test_R1_edge_tensor_random | R1 | G = Vol·I on random Voronoi |
| test_R2_face_tensor_cubic | R2 | H = Vol·I on Kelvin, C15, WP |
| test_R2_face_tensor_random | R2 | H = Vol·I on random Voronoi |
| test_R3_perpendicularity | R3 | ⋆₁ > 0, ℓ > 0 + dual edge ∥ face normal |
| test_R4_trace_identity | R4 | tr(G) = 3·Vol |
| test_R34_degenerate_voronoi | R34 | G = H = Vol·I under cell ratio 40× |
| test_N12_admissible_subspace | R37 | Codim-6: rank = 6, null_dim = |E|−6 |
| test_T2_per_vertex_flux_closure | T2 | Per-vertex Σ sign·⋆₁·Δx = 0 (Term 1 = 0) |
| test_T4_tight_frame | T4 | {√⋆₁·Δx} tight frame, F†F = Vol·I, σ equal |
| test_T5_admissible_positivity | T5 | Cone radius: 1.4–6× mean(⋆₁) |
| test_T7_domain_scaling | T7 | G → s³·I under V→sV, L→sL (dimensional) |

## 2_test_c_squared_one.py (§3) — DONE, 10 tests

| Test | ID | What it verifies |
|------|----|-----------------|
| test_R5_gauge_speed | R5 | c = 1 on Kelvin, C15, WP |
| test_R38_random_voronoi | R38 | c = 1 on 5 random Voronoi seeds |
| test_R32_scaling_invariance | R32 | c = 1 independent of box scale |
| test_R41_H1_vanishes | R41 | H¹(k≠0) = 0 on z=4 (cubic + random), d₁d₀ = 0 |
| test_R44e_schur_complement | R44e | Schur eigenvalues [0, k², k²] |
| test_R44f_discrete_curl | R44f | dF/dk = iε on 6 structures (tol 1e-8) |
| test_R44g_full_proof_chain | R44g | G=H=Vol·I + rank(d₁(0))=nV−2 + u_perp_T∈H²(99.99%) → c²=1 |
| test_R6_converse | R6 | Non-Voronoi stars → c ≠ 1 (1929× worse) |
| test_S1_multi_direction | S1 | c = 1 along 5 k-directions (isotropy) |
| test_S2_dispersion_convergence | S2 | δc/k² → const (dispersion O(k⁴)) |

## 3_test_removing_exactness.py (§4) — DONE, 7 tests

| Test | ID | What it verifies |
|------|----|-----------------|
| test_I6_standard_speed | R12 | c_std = 1.25–1.68 on cubic |
| test_R25_delta_c_O1 | R25 | δc = O(1), not O(k²) |
| test_N16_no_principal_symbol | R35 | No consistent c²(k̂) = k̂ᵀCk̂ (RMS 0.39–0.84) |
| test_R17_anisotropy | R17 | Standard anisotropy 37–70% |
| test_N19_nlost_identity | R40 | n_lost = Δrank(d₁) on 3×3 combinations |
| test_N13_leaked_gradients | R41 | Leaked forms ⊂ im(d₀) (grad_frac ≥ 0.885) |
| test_R23_pollution_direction | R23 | n_lost: axis < face diag < body diag |

## 4_test_removing_geometry.py (§5) — DONE, 3 tests

| Test | ID | What it verifies |
|------|----|-----------------|
| test_N15_killer_table_cubic | R36 | 2×2 table on Kelvin, C15, WP (interaction < 0.05) |
| test_R39_killer_table_random | R39 | 2×2 table on 3 random Voronoi (interaction < 0.10) |
| test_M1_error_factorization | R29 | Δc_geom + Δc_topo separable at 5/10/15% (interaction < 0.05) |

## 5_test_spectral_structure.py (§6) — DONE, 2 tests

| Test | ID | What it verifies |
|------|----|-----------------|
| test_R22_moment_hierarchy | R22 (includes R14) | tr(K^n) conserved n=1,2 (< 1e-10), breaks n=3 (0.4–0.6%) |
| test_R21_girth_valence | R21 | Face adjacency girth = 3, edge valence = 3 uniform |

**Note:** Theorem uses girth(B) = 2×girth(face_adj) = 6, where B = bipartite incidence graph (E∪F).
Implicit assumption: no two faces share two edges (trivial on Voronoi — convex faces).

## 6_test_appendix_a.py (App A — naive PT) — DONE, 5 tests

| Test | ID | What it verifies | Source |
|------|----|-----------------|--------|
| test_R44a_eigvec_harmonic | R44a | Eigenvector 99.999% harmonic | old/6 |
| test_R44b_rayleigh_paradox | R44b | R[h]/k² = 4.8–10.7 yet ω²/k² = 1 | old/6 |
| test_R44c_projection_worse | R44c | Gradient projection increases R | old/6 |
| test_R44d_cancellation | R44d | h†Kh + cross + δψ†Kδψ = k² | old/6 |
| test_R11_plane_wave_overlap | R11 | Overlap → 1 as k → 0 | old/1 |

## 7_test_appendix_b.py (App B — std anatomy) — DONE, 3 tests

| Test | ID | What it verifies | Source |
|------|----|-----------------|--------|
| test_R16_mode_mixing | R16 | Lowest mode = 89% gauge + 11% acoustic | old/3 |
| test_R20_N_dependence | R20 | c_std non-monotonic with N | old/3 |
| test_R27_TRIM_collapse | R27 | Standard = exact at k = π/L_cell | old/4 |

---

**Total:** 7 files, 42 tests.
