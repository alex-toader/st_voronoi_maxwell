# Advancements Map — st_voronoi_maxwell

**Date:** Mar 2026

---

## A. Metric identity: G = H = Vol·I (§2)

| ID | Result | Status | Test |
|----|--------|--------|------|
| R1 | G = Vol·I on cubic (Kelvin, C15, WP) | Done | file 1 |
| R1 | G = Vol·I on random Voronoi | Done | file 1 |
| R2 | H = Vol·I on cubic | Done | file 1 |
| R2 | H = Vol·I on random Voronoi | Done | file 1 |
| R3 | Voronoi edge ⊥ dual face | Done | file 1 |
| R4 | tr(G) = tr(H) = 3·Vol | Done | file 1 |
| R7 | Dual edge ⊥ primal face | Done | file 1 |
| R8 | Per-vertex flux closure (∫ n dS = 0) | Done | file 1 |
| R34 | Degenerate Voronoi: G = H = Vol·I survives | Done | file 1 |
| R37 | Admissible subspace codim-6 (N12) | Done | file 5 |

## B. Theorem c² = 1 (§3)

| ID | Result | Status | Test |
|----|--------|--------|------|
| R5 | c = 1 spectral on Kelvin, C15, WP | Done | file 1 |
| R5 | c = 1 isotropic to < 10⁻⁴ | Done | file 1 |
| R38 | c = 1 on 5 random Voronoi seeds | Done | file 5 |
| R32 | c = 1 independent of box scale | Done | file 1 |
| R11 | Acoustic eigenvector ≈ plane wave (overlap → 1) | Done | file 1 |
| R44a | Eigenvector 99.999% harmonic | Done | file 6 |
| R44b | Rayleigh paradox: R[h]/k² = 4.8–10.7, yet ω²/k² = 1 | Done | file 6 |
| R44c | Gradient projection makes Rayleigh worse | Done | file 6 |
| R44d | Three-way cancellation: h†Kh + cross + δψ†Kδψ = k² | Done | file 6 |
| R44e | Schur complement eigenvalues [0, k², k²] | Done | file 6 |
| R44f | Discrete curl identity dF/dk = iε (10⁻¹¹, 7 structures) | Done | file 6 |
| R44g | Full proof chain: G=H=Vol·I + exactness + curl → c²=1 | Done | file 6 |
| R43 | 1st + 2nd order PT cancellation anatomy | Done | file 5 |

## C. Necessity — removing exactness (§4)

| ID | Result | Status | Test |
|----|--------|--------|------|
| R12 | c_std ≠ 1 (1.25–1.68 on cubic) | Done | file 2 |
| R13 | Plane wave overlap → 0.06 on standard | Done | file 2 |
| R16 | Lowest std mode: 89% gauge + 11% acoustic | Done | file 3 |
| R17 | Standard anisotropy 37–70% | Done | file 3 |
| R23 | n_lost depends on k-direction | Done | file 4 |
| R25 | δc = O(1), not O(k²) | Done | file 4 |
| R27 | TRIM collapse: standard = exact at k = π/L_cell | Done | file 4 |
| R28 | Band contamination profile | Done | file 4 |
| R35 | No consistent principal symbol (N16) | Done | file 5 |
| R40 | n_lost = Δrank(d₁) exact identity (N19) | Done | file 5 |
| R41 | Leaked forms ⊂ im(d₀), H¹(k≠0) = 0 (N13) | Done | file 5 |
| R30 | c_std on random Voronoi: [0.61, 1.04] | Done | file 4 |
| R20 | c_std non-monotonic with N | Done | file 3 |

## D. Necessity — removing geometry (§5)

| ID | Result | Status | Test |
|----|--------|--------|------|
| R6 | Perturbed Hodge stars → c ≠ 1 | Done | file 1 |
| R29 | Error factorization separable (M1) | Done | file 4 |
| R36 | Killer table on 3 cubic (N15) | Done | file 5 |
| R39 | Killer table on 3 random Voronoi | Done | file 5 |

## E. Spectral structure (§6)

| ID | Result | Status | Test |
|----|--------|--------|------|
| R14 | tr(K) conserved on exact/standard (I8) | Done | file 2 |
| R22 | Moment hierarchy: tr(K^n) conserved for n < girth | Done | file 3 |
| R21 | Face adjacency girth = 3 on all Voronoi | Done | file 3 |

## F. Dielectric (§7 — future paper)

| ID | Result | Status | Test |
|----|--------|--------|------|
| R15 | H_ε = ⟨1/ε⟩·Vol·I (I11) | Done | file 2 |
| R19 | Harmonic mean preserves metric, log breaks it (M2) | Done | file 3 |
| R24 | Dielectric breaks G isotropy in layering direction | Done | file 4 |

## G. Converse and auxiliary

| ID | Result | Status | Test |
|----|--------|--------|------|
| R33 | V − E + F = C on T³ (Euler) | Done | file 1 |
| R9 | Diamond tiling | Done | file 1 |
| R10 | Kahan summation parity | Done | file 1 |

---

## Open directions (not yet investigated)

| ID | Direction | Priority |
|----|-----------|----------|
| N29 | c² = 1 in 2D and 4D (dimension-independent?) | High — paper 4 |
| N30 | Analytic proof of dF/dk = iε via Stokes | High — needed for paper 3 |
| N31 | im(d₂†) component starts at O(k²) — analytic proof | Medium |
| N25 | "Only if" for condition (B): d₁d₀ ≠ 0 ⇒ c ≠ 1 | Open problem |
| N20 | Ihara zeta / non-backtracking walk connection | Low |
| N24 | Non-Voronoi mesh with exact d₁d₀ = 0 | Medium |
| N26 | p-forms in d dimensions | Future |
| C8 | Girth ≠ 3 test (mesh with girth 4) | Optional |
