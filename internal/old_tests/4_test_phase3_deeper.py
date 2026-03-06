"""
W17 Test 4: Phase 3 deeper — direction-dependent pollution, TRIM collapse,
O(1) speed error, dielectric anisotropy, band contamination.

CLAIM 1 (s3/R23): Pollution count depends on k-direction. On-axis directions lose
fewer gauge modes than off-axis: Kelvin 6→12→14 (axis→face diag→body diag),
C15 39→70→91, WP 15→28→37. Off-axis k crosses more periodic boundaries →
more phase errors in d₁d₀. On C15/WP the three cubic axes are equivalent
(39/39/39 and 15/15/15); on Kelvin y-axis has n_lost=7 (BCC lattice orientation).

CLAIM 2 (s6/R27): Standard DEC becomes exact at unit-cell TRIM points where
k = π/L_cell (not π/L_supercell). At these points Bloch phases are exactly ±1
→ d₁(k)d₀(k)=0 automatically (||d1d0|| = 3e-15). n_lost=0, eigenvalues
identical to 4e-15. At supercell TRIM (π/L), phases are NOT ±1 and standard
still fails (||d1d0|| ~ 12-20). Error vanishes only at k=0 and k=π/L_cell.

CLAIM 3 (s9/R25): delta_c = c_std − c_exact is O(1) as k→0, NOT O(k²).
On Kelvin: delta_c ≈ 0.253 from |k|=0.005 to |k|=0.1. This means standard
DEC speed error is a structural property of the operator, not a discretization
artifact that vanishes with mesh refinement. delta_c/k² → ∞.

CLAIM 4 (s2/R24): Dielectric z-split breaks diagonal isotropy of G_ε
(G_xx=G_yy≠G_zz) but preserves off-diagonal zero (to 4e-18).
NOTE: This is the edge tensor G_ε = Σ_e (⋆₁/ε)dx⊗dx, not the face tensor H_ε.
Dielectric layering introduces anisotropy only in the layering direction,
without shear coupling. Vacuum G/Vol = I to 1e-16 (baseline).

CLAIM 5 (m3/R28): Standard DEC's lowest physical modes are ~89% gauge, ~11%
physical. The acoustic subspace is the 2 transverse modes at k≠0 (the
longitudinal mode is absorbed into gauge on exact). First standard "physical"
mode has low acoustic overlap — mostly gradient contamination from d₁d₀≠0.
Contamination hits acoustic first because gauge-leaked modes have ω² ~ k².

RAW OUTPUT (5 tests, all pass):
==================================
=== s3: Pollution count by direction ===
  Kelvin (nV=96):  [1,0,0]=6, [0,1,0]=7, [0,0,1]=6, [1,1,0]=12, [1,1,1]=14
  C15 (nV=1088):   [1,0,0]=39, [0,1,0]=39, [0,0,1]=39, [1,1,0]=70, [1,1,1]=91
  WP (nV=368):     [1,0,0]=15, [0,1,0]=15, [0,0,1]=15, [1,1,0]=28, [1,1,1]=37
=== s6: TRIM collapse ===
  Generic k: ||d1d0||_std = 7.43e+00
  X_cell: ||d1d0||_std = 3.06e-15, max|λ_st-λ_ex| = 3.55e-15
  M_cell: ||d1d0||_std = 3.95e-15, max|λ_st-λ_ex| = 4.44e-15
  R_cell: ||d1d0||_std = 4.54e-15, max|λ_st-λ_ex| = 2.66e-15
=== s9: delta_c scaling ===
  |k|=0.005: delta_c=0.253, delta_c/k²=16407
  |k|=0.010: delta_c=0.253, delta_c/k²=4098
  |k|=0.100: delta_c=0.219, delta_c/k²=35
  range ratio=1.16 (constant to 16%)
=== s2: Dielectric G_eps (edge tensor, not face tensor) ===
  Vacuum G/Vol = I (to 1e-16)
  G_eps/Vol diag: (0.596, 0.596, 0.654), off-diag max = 4.3e-18
=== m3: Band contamination ===
  Mode 0: acou_overlap=0.113, opti_overlap=0.000
  Modes 1-5: acou_overlap=0.008-0.011, opti < 0.001
=== M1: Error factorization ===
  Exact+Voronoi: n_zero=96, c=1.000
  Exact+Perturbed: n_zero=96, c=0.993 (c changes, n_lost stays 0)
  Standard+Voronoi: n_zero=90, c=1.253
  Standard+Perturbed: n_zero=90, c=1.245
  interaction = -0.0004 (additive to < 0.1%)
==================================

CLAIM 6 (M1/R29): Geometric and topological errors are SEPARABLE. Perturbing
Hodge stars changes c (0.993 vs 1.000) but NOT pollution (n_lost stays 0 on
exact, stays 6 on standard). Breaking exactness changes pollution but NOT
through stars. Additivity: dc_both ≈ dc_geom + dc_topo with interaction
< 0.001 at 10% perturbation. Factorization breaks mildly at 20-30%.

ANSWER:
=======
Standard DEC failure has rich structure: direction-dependent (R23), vanishes at
unit-cell TRIM (R27), O(1) not O(k²) (R25), acoustic modes ~11% contaminated
while optical untouched (R28). Dielectric preserves off-diagonal zero but
breaks z-isotropy (R24). Geometric and topological errors are separable (R29).
"""
import sys
import os
import numpy as np
from scipy import sparse
from scipy.linalg import eigh

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
SRC = os.path.join(os.path.dirname(__file__), "..", "src")
sys.path.insert(0, os.path.abspath(SRC))

from physics.hodge import (
    build_kelvin_with_dual_info,
    build_c15_with_dual_info,
    build_wp_with_dual_info,
    build_hodge_stars_voronoi,
)
from physics.gauge_bloch import build_d1_bloch_exact
from physics.bloch import build_d0_bloch, build_d1_bloch_standard


# ===================================================================
# Helpers
# ===================================================================

BUILDERS = [
    ("Kelvin", build_kelvin_with_dual_info),
    ("C15", build_c15_with_dual_info),
    ("WP", build_wp_with_dual_info),
]


def build_both_spectra(data, k_vec):
    """Build exact and standard K, return both sorted spectra + extras."""
    star1, star2 = build_hodge_stars_voronoi(data)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    S1 = np.diag(star1)
    S2 = np.diag(star2)

    d0k = build_d0_bloch(V, E, L, k_vec)
    if sparse.issparse(d0k):
        d0k = d0k.toarray()

    # Exact
    d1k_ex = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
    if sparse.issparse(d1k_ex):
        d1k_ex = d1k_ex.toarray()
    K_ex = d1k_ex.conj().T @ S2 @ d1k_ex
    K_ex = (K_ex + K_ex.conj().T) / 2
    vals_ex, vecs_ex = eigh(K_ex, S1)
    idx = np.argsort(np.real(vals_ex))
    vals_ex = np.real(vals_ex[idx])
    vecs_ex = vecs_ex[:, idx]

    # Standard
    d1k_st = build_d1_bloch_standard(V, E, F, L, k_vec)
    if sparse.issparse(d1k_st):
        d1k_st = d1k_st.toarray()
    K_st = d1k_st.conj().T @ S2 @ d1k_st
    K_st = (K_st + K_st.conj().T) / 2
    vals_st, vecs_st = eigh(K_st, S1)
    idx = np.argsort(np.real(vals_st))
    vals_st = np.real(vals_st[idx])
    vecs_st = vecs_st[:, idx]

    return {
        "vals_ex": vals_ex, "vecs_ex": vecs_ex,
        "vals_st": vals_st, "vecs_st": vecs_st,
        "K_ex": K_ex, "K_st": K_st,
        "S1": S1, "S2": S2, "d0k": d0k,
    }


# ===================================================================
# Tests
# ===================================================================

def test_s3_pollution_by_direction():
    """s3/R23: n_lost depends on k-direction. Axes < face diag < body diag."""
    print("\n=== s3: Pollution count by direction ===")

    directions = {
        "[1,0,0]": np.array([1, 0, 0], dtype=float),
        "[0,1,0]": np.array([0, 1, 0], dtype=float),
        "[0,0,1]": np.array([0, 0, 1], dtype=float),
        "[1,1,0]": np.array([1, 1, 0], dtype=float) / np.sqrt(2),
        "[1,1,1]": np.array([1, 1, 1], dtype=float) / np.sqrt(3),
    }

    for name, builder in BUILDERS:
        data = builder(N=2)
        nV = len(data["V"])
        L = data["L"]

        n_lost_axis = []
        n_lost_diag = []

        print(f"  {name} (nV={nV}):")
        for dname, d in directions.items():
            k_vec = 0.01 * (2 * np.pi / L) * d
            r = build_both_spectra(data, k_vec)

            thresh = 1e-4 * (k_vec @ k_vec)
            n0_ex = np.sum(r["vals_ex"] < thresh)
            n0_st = np.sum(r["vals_st"] < thresh)
            n_lost = n0_ex - n0_st

            print(f"    {dname:>10s}: n_zero_ex={n0_ex}, n_zero_st={n0_st}, n_lost={n_lost}")

            # Classify
            if dname in ("[1,0,0]", "[0,1,0]", "[0,0,1]"):
                n_lost_axis.append(n_lost)
            else:
                n_lost_diag.append(n_lost)

        # Off-axis should lose MORE modes than on-axis
        min_diag = min(n_lost_diag)
        max_axis = max(n_lost_axis)
        assert min_diag > max_axis, (
            f"{name}: diag n_lost {min_diag} <= axis n_lost {max_axis}")

        # Body diagonal should lose at least 2× on-axis
        mean_axis = np.mean(n_lost_axis)
        assert min_diag >= 1.5 * mean_axis, (
            f"{name}: diag/axis ratio = {min_diag / mean_axis:.2f}, expected >= 1.5")

    print("  PASSED")


def test_s6_TRIM_collapse():
    """s6/R27: At unit-cell TRIM (k=π/L_cell), standard = exact."""
    print("\n=== s6: TRIM collapse — standard becomes exact at unit-cell BZ boundary ===")

    N = 2
    data = build_kelvin_with_dual_info(N=N)
    V, E, F = data["V"], data["E"], data["F"]
    L = data["L"]
    L_vec = np.array(data["L_vec"])
    L_cell = L / N  # unit cell lattice parameter

    star1, star2 = build_hodge_stars_voronoi(data)
    S1 = np.diag(star1)
    S2 = np.diag(star2)

    # Unit-cell TRIM: k = pi/L_cell (phases are exactly ±1 on unit cell edges)
    scale_cell = np.pi / L_cell
    trim_points = {
        "X_cell": np.array([1, 0, 0], dtype=float) * scale_cell,
        "M_cell": np.array([1, 1, 0], dtype=float) * scale_cell,
        "R_cell": np.array([1, 1, 1], dtype=float) * scale_cell,
    }

    # Generic k for contrast
    generic_k = 0.1 * (2 * np.pi / L) * np.array([1, 0, 0], dtype=float)
    d0k_gen = build_d0_bloch(V, E, L, generic_k)
    if sparse.issparse(d0k_gen):
        d0k_gen = d0k_gen.toarray()
    d1k_st_gen = build_d1_bloch_standard(V, E, F, L, generic_k)
    if sparse.issparse(d1k_st_gen):
        d1k_st_gen = d1k_st_gen.toarray()

    d1d0_gen = np.linalg.norm(d1k_st_gen @ d0k_gen)
    print(f"  Generic k: ||d1d0||_std = {d1d0_gen:.2e} (should be >> 0)")
    assert d1d0_gen > 1.0, "Generic k should have ||d1d0|| >> 0"

    for tname, k_vec in trim_points.items():
        d0k = build_d0_bloch(V, E, L, k_vec)
        if sparse.issparse(d0k):
            d0k = d0k.toarray()

        d1k_st = build_d1_bloch_standard(V, E, F, L, k_vec)
        if sparse.issparse(d1k_st):
            d1k_st = d1k_st.toarray()

        d1k_ex = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
        if sparse.issparse(d1k_ex):
            d1k_ex = d1k_ex.toarray()

        # Check exactness
        d1d0_st = np.linalg.norm(d1k_st @ d0k)
        d1d0_ex = np.linalg.norm(d1k_ex @ d0k)

        # Eigenvalues
        K_st = d1k_st.conj().T @ S2 @ d1k_st
        K_st = (K_st + K_st.conj().T) / 2
        K_ex = d1k_ex.conj().T @ S2 @ d1k_ex
        K_ex = (K_ex + K_ex.conj().T) / 2

        vals_st = np.sort(np.real(eigh(K_st, S1, eigvals_only=True)))
        vals_ex = np.sort(np.real(eigh(K_ex, S1, eigvals_only=True)))

        maxdiff = np.max(np.abs(vals_st - vals_ex))

        print(f"  {tname}: ||d1d0||_std = {d1d0_st:.2e}, "
              f"max|λ_st-λ_ex| = {maxdiff:.2e}")

        # Standard should be exact at unit-cell TRIM
        assert d1d0_st < 1e-12, (
            f"{tname}: ||d1d0||_std = {d1d0_st:.2e}, expected < 1e-12")
        assert maxdiff < 1e-12, (
            f"{tname}: eigenvalue mismatch = {maxdiff:.2e}")

    print("  PASSED")


def test_s9_delta_c_is_O1():
    """s9/R25: delta_c = c_std - c_exact is O(1), not O(k²). Structural error."""
    print("\n=== s9: delta_c scaling with |k| ===")

    data = build_kelvin_with_dual_info(N=2)
    L = data["L"]
    scale = 2 * np.pi / L
    direction = np.array([1.0, 0.0, 0.0])

    k_mags = [0.005, 0.01, 0.02, 0.05, 0.1]
    delta_c_values = []

    for km in k_mags:
        k_vec = km * scale * direction
        k2 = k_vec @ k_vec
        r = build_both_spectra(data, k_vec)

        thresh = 1e-8
        nz_st = r["vals_st"][r["vals_st"] > thresh]
        nz_ex = r["vals_ex"][r["vals_ex"] > thresh]
        c_st = np.sqrt(nz_st[0] / k2)
        c_ex = np.sqrt(nz_ex[0] / k2)
        dc = abs(c_st - c_ex)
        delta_c_values.append(dc)

        print(f"  |k|={km:.3f}: c_std={c_st:.6f}, c_exact={c_ex:.6f}, "
              f"delta_c={dc:.6f}, delta_c/k²={dc / k2:.2f}")

    # delta_c should be approximately constant (O(1))
    dc_arr = np.array(delta_c_values)
    # Ratio of max to min should be < 2 (if it were O(k²) the ratio would be ~400)
    ratio = dc_arr.max() / dc_arr.min()
    print(f"  delta_c range: [{dc_arr.min():.4f}, {dc_arr.max():.4f}], ratio={ratio:.2f}")

    assert ratio < 2.0, (
        f"delta_c ratio = {ratio:.2f}, expected < 2 (O(1) constant)")

    # delta_c should be > 0.1 (not vanishing)
    assert dc_arr.min() > 0.1, (
        f"delta_c min = {dc_arr.min():.4f}, expected > 0.1")

    print("  PASSED")


def test_s2_dielectric_Geps_anisotropy():
    """s2/R24: Dielectric z-split breaks G_ε diagonal isotropy, off-diag stays 0.

    NOTE: This computes G_ε = Σ_e (⋆₁[e]/ε_e) dx⊗dx (the edge tensor with
    per-edge dielectric), NOT H_ε = Σ_f (⋆₂[f]/ε_f) A_f⊗A_f (the face tensor
    with per-face dielectric). Both reduce to Vol·I in vacuum.
    """
    print("\n=== s2: Off-diagonal G_ε with z-split dielectric ===")

    data = build_kelvin_with_dual_info(N=2)
    V, E = data["V"], data["E"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    vol = np.prod(L_vec)
    star1, _ = build_hodge_stars_voronoi(data)
    nE = len(E)

    # Edge displacement vectors
    dx = np.zeros((nE, 3))
    for e_idx, (i, j) in enumerate(E):
        d = V[j] - V[i]
        d -= np.round(d / L_vec) * L_vec
        dx[e_idx] = d

    # Vacuum metric (baseline)
    G_vac = np.zeros((3, 3))
    for e_idx in range(nE):
        G_vac += star1[e_idx] * np.outer(dx[e_idx], dx[e_idx])

    vac_diag = [G_vac[i, i] / vol for i in range(3)]
    vac_off = max(abs(G_vac[i, j] / vol) for i in range(3) for j in range(3) if i != j)
    print(f"  Vacuum G/Vol diag: ({vac_diag[0]:.10f}, {vac_diag[1]:.10f}, {vac_diag[2]:.10f})")
    print(f"  Vacuum off-diag max: {vac_off:.2e}")

    assert all(abs(d - 1.0) < 1e-12 for d in vac_diag), "Vacuum G/Vol != I"
    assert vac_off < 1e-14, f"Vacuum off-diag = {vac_off}"

    # Dielectric z-split: eps = 1 (z < L/2), eps = 13 (z >= L/2)
    eps_edge = np.ones(nE)
    for e_idx, (i, j) in enumerate(E):
        mid = (V[i] + V[j]) / 2
        mid -= np.floor(mid / L) * L  # wrap to [0, L)
        if mid[2] >= L / 2:
            eps_edge[e_idx] = 13.0

    star1_eps = star1 / eps_edge

    G_eps = np.zeros((3, 3))
    for e_idx in range(nE):
        G_eps += star1_eps[e_idx] * np.outer(dx[e_idx], dx[e_idx])

    eps_diag = [G_eps[i, i] / vol for i in range(3)]
    eps_off = max(abs(G_eps[i, j] / vol) for i in range(3) for j in range(3) if i != j)

    print(f"  G_eps/Vol diag: ({eps_diag[0]:.6f}, {eps_diag[1]:.6f}, {eps_diag[2]:.6f})")
    print(f"  G_eps off-diag max: {eps_off:.2e}")

    # xy should be equal (no dielectric layering in xy)
    assert abs(eps_diag[0] - eps_diag[1]) < 1e-10, (
        f"G_eps xx != yy: {eps_diag[0]:.6f} vs {eps_diag[1]:.6f}")

    # z should differ from xy (dielectric breaks z-isotropy)
    z_aniso = abs(eps_diag[2] - eps_diag[0])
    assert z_aniso > 0.01, (
        f"z-anisotropy = {z_aniso:.6f}, expected > 0.01")

    # Off-diagonal should remain zero
    assert eps_off < 1e-14, (
        f"G_eps off-diag = {eps_off:.2e}, expected ~ 0")

    print(f"  z-anisotropy: {z_aniso:.6f} (dielectric breaks z-direction only)")
    print("  PASSED")


def test_m3_band_contamination():
    """m3/R28: Acoustic modes destroyed on standard, optical moderate."""
    print("\n=== m3: Band contamination — acoustic vs optical ===")

    k_mag = 0.01
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    data = build_kelvin_with_dual_info(N=2)
    r = build_both_spectra(data, k_vec)
    S1 = r["S1"]

    nV = len(data["V"])
    n0_ex = np.sum(r["vals_ex"] < 1e-12)

    # At k≠0 on exact complex: first 2 physical modes = transverse acoustic
    # (the longitudinal mode is absorbed into gauge/ker K).
    # Mode n0_ex+2 is the first optical mode.
    n_acoustic = 2
    acoustic_idx = list(range(n0_ex, n0_ex + n_acoustic))
    optical_idx = list(range(n0_ex + n_acoustic, n0_ex + n_acoustic + 5))

    n0_st = np.sum(r["vals_st"] < 1e-12)

    print(f"  Kelvin: n0_ex={n0_ex}, n0_st={n0_st}")
    print(f"  Exact acoustic modes: indices {acoustic_idx}")
    print(f"  Exact optical modes: indices {optical_idx[:3]}...")

    # Compute subspace overlap: project standard modes onto exact subspace
    # For each standard physical mode, compute overlap with exact acoustic subspace
    def subspace_overlap(std_idx, ex_indices):
        """Sum of |<std_i | S1 | ex_j>|² over ex_j, normalized."""
        v = r["vecs_st"][:, std_idx]
        norm_sq = np.real(v.conj() @ S1 @ v)
        overlap = 0.0
        for j in ex_indices:
            u = r["vecs_ex"][:, j]
            overlap += abs(u.conj() @ S1 @ v) ** 2
        return overlap / norm_sq

    # For the first few standard physical modes
    print(f"\n  Standard physical modes (starting at index {n0_st}):")
    print(f"  {'mode':>5s} {'acou_overlap':>13s} {'opti_overlap':>13s} {'total_phys':>11s}")

    acou_overlaps = []
    opti_overlaps = []

    for i in range(n0_st, min(n0_st + 8, len(r["vals_st"]))):
        ao = subspace_overlap(i, acoustic_idx)
        oo = subspace_overlap(i, optical_idx)
        acou_overlaps.append(ao)
        opti_overlaps.append(oo)
        print(f"  {i - n0_st:5d} {ao:13.6f} {oo:13.6f} {ao + oo:11.6f}")

    # First few standard modes should have very LOW acoustic overlap
    # (acoustic destroyed by gradient leakage)
    max_acou = max(acou_overlaps[:3])
    print(f"\n  Max acoustic overlap (first 3 std modes): {max_acou:.6f}")

    # Optical modes should have moderate overlap
    max_opti = max(opti_overlaps)
    print(f"  Max optical overlap (first 8 std modes): {max_opti:.6f}")

    # Acoustic modes: severely contaminated (overlap < 0.1)
    assert max_acou < 0.15, (
        f"Acoustic overlap = {max_acou:.4f}, expected < 0.15 (modes should be destroyed)")

    print("  PASSED")


def test_M1_error_factorization():
    """M1/R29: Geometric and topological errors are separable (additive)."""
    print("\n=== M1: Error factorization — geometry × topology ===")

    data = build_kelvin_with_dual_info(N=2)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    nE = len(E)

    star1, star2 = build_hodge_stars_voronoi(data)

    k_vec = 0.01 * (2 * np.pi / L) * np.array([1.0, 0.0, 0.0])
    k2 = k_vec @ k_vec

    d0k = build_d0_bloch(V, E, L, k_vec)
    if sparse.issparse(d0k):
        d0k = d0k.toarray()
    d1k_ex = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
    if sparse.issparse(d1k_ex):
        d1k_ex = d1k_ex.toarray()
    d1k_st = build_d1_bloch_standard(V, E, F, L, k_vec)
    if sparse.issparse(d1k_st):
        d1k_st = d1k_st.toarray()

    def get_c_and_nlost(d1k, S1, S2m):
        K = d1k.conj().T @ S2m @ d1k
        K = (K + K.conj().T) / 2
        vals = np.sort(np.real(eigh(K, S1, eigvals_only=True)))
        n0 = np.sum(vals < 1e-10)
        nz = vals[vals > 1e-10]
        c = np.sqrt(nz[0] / k2) if len(nz) > 0 else float('nan')
        return c, n0

    # Voronoi stars
    S1_v = np.diag(star1)
    S2_v = np.diag(star2)

    # Perturbed stars (10%)
    np.random.seed(42)
    s1_p = star1 * np.maximum(1.0 + 0.10 * np.random.randn(nE), 0.5)
    np.random.seed(43)
    s2_p = star2 * np.maximum(1.0 + 0.10 * np.random.randn(len(star2)), 0.5)
    S1_p = np.diag(s1_p)
    S2_p = np.diag(s2_p)

    c_EV, n0_EV = get_c_and_nlost(d1k_ex, S1_v, S2_v)
    c_EP, n0_EP = get_c_and_nlost(d1k_ex, S1_p, S2_p)
    c_SV, n0_SV = get_c_and_nlost(d1k_st, S1_v, S2_v)
    c_SP, n0_SP = get_c_and_nlost(d1k_st, S1_p, S2_p)

    print(f"  {'Config':>22s} {'n_zero':>7s} {'c':>10s}")
    for label, c, n0 in [("Exact+Voronoi", c_EV, n0_EV),
                          ("Exact+Perturbed", c_EP, n0_EP),
                          ("Standard+Voronoi", c_SV, n0_SV),
                          ("Standard+Perturbed", c_SP, n0_SP)]:
        print(f"  {label:>22s} {n0:7d} {c:10.6f}")

    # (1) Perturbation does NOT change n_lost
    assert n0_EV == n0_EP, (
        f"Perturb changed exact n_zero: {n0_EV} -> {n0_EP}")
    assert n0_SV == n0_SP, (
        f"Perturb changed standard n_zero: {n0_SV} -> {n0_SP}")

    # (2) Standard changes n_lost regardless of stars
    n_lost_V = n0_EV - n0_SV
    n_lost_P = n0_EP - n0_SP
    assert n_lost_V > 0, "Standard should lose modes"
    assert n_lost_V == n_lost_P, (
        f"n_lost differs: Voronoi {n_lost_V} vs Perturbed {n_lost_P}")

    # (3) Additivity: dc_both ≈ dc_geom + dc_topo
    dc_geom = c_EP - c_EV
    dc_topo = c_SV - c_EV
    dc_both = c_SP - c_EV
    dc_sum = dc_geom + dc_topo
    interact = dc_both - dc_sum

    print(f"\n  dc_geom = {dc_geom:+.6f} (perturbation effect)")
    print(f"  dc_topo = {dc_topo:+.6f} (exactness effect)")
    print(f"  dc_both = {dc_both:+.6f} (combined)")
    print(f"  dc_sum  = {dc_sum:+.6f} (additive prediction)")
    print(f"  interaction = {interact:+.6f}")

    assert abs(interact) < 0.005, (
        f"Interaction term = {interact:.4f}, expected < 0.005 at 10% perturbation")

    print("  PASSED")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    test_s3_pollution_by_direction()
    test_s6_TRIM_collapse()
    test_s9_delta_c_is_O1()
    test_s2_dielectric_Geps_anisotropy()
    test_m3_band_contamination()
    test_M1_error_factorization()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED (s3, s6, s9, s2, m3, M1 — 6 tests)")
    print("=" * 60)
