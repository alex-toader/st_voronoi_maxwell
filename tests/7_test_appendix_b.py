"""
Appendix B: Standard DEC anatomy — mode mixing and special points.

Supports: Appendix B (supporting material for §4)

CLAIM: The standard DEC failure (§4) has additional structure beyond c ≠ 1.
(1) The lowest nonzero standard mode is ~89% gauge + ~11% acoustic — a mixture,
not a pure gradient. (2) c_std depends on supercell size N non-monotonically
(1.25, 1.58, 1.39 for N=2,3,4). (3) At unit-cell TRIM points (k = π/L_cell),
standard DEC becomes exact: d₁d₀ = 0 automatically because Bloch phases are ±1.

RAW OUTPUT (3 tests, all pass):
=================================
=== R16: mode anatomy (Kelvin) ===
  Mode 0: omega/k=1.253, gauge_frac=0.886, phys_frac=0.114 (mixing)
  Modes 1-5: gauge_frac ≥ 0.988 (nearly pure gauge)
  Modes 6-7: gauge_frac ≥ 0.998 (pure gauge)
=== R20: c_standard vs N (Kelvin) ===
  N=2: 1.253, N=3: 1.576, N=4: 1.392 (non-monotonic, spread 0.32)
=== R27: TRIM collapse (Kelvin N=2) ===
  Generic k: ||d1d0||_std = 7.43e+00
  X_cell: ||d1d0||_std = 3.1e-15, max|λ_st-λ_ex| = 3.0e-15
  M_cell: ||d1d0||_std = 4.0e-15, max|λ_st-λ_ex| = 3.1e-15
  R_cell: ||d1d0||_std = 4.5e-15, max|λ_st-λ_ex| = 2.4e-15
ALL TESTS PASSED (App B — 3 tests)

ANSWER:
=======
Standard DEC pathologies go beyond c ≠ 1: mode mixing destroys the gauge/physical
separation, c_std is non-monotonic in N (no convergence), and exactness is restored
only at the isolated TRIM points where phases become ±1.
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
    build_hodge_stars_voronoi,
)
from physics.gauge_bloch import build_d1_bloch_exact
from physics.bloch import build_d0_bloch, build_d1_bloch_standard


# ===================================================================
# Helpers
# ===================================================================

def build_both_spectra(data, k_vec):
    """Build exact and standard K, return sorted spectra + eigenvectors."""
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
        "S1": S1, "S2": S2, "d0k": d0k,
    }


# ===================================================================
# Tests
# ===================================================================

def test_R16_mode_mixing():
    """R16: Lowest standard nonzero mode is ~89% gauge + ~11% acoustic.

    Project each standard physical mode onto the exact gauge subspace
    (Γ zero modes lifted to k). The lowest mode is a mixture — not a
    pure gradient. Modes 1-5 are >95% gauge. This shows standard DEC
    doesn't just shift speeds: it destroys the gauge/physical boundary.
    """
    print("\n=== R16: mode anatomy — gauge/physical decomposition ===")

    k_vec = 0.005 * np.array([1.0, 0.0, 0.0])
    k_mag = np.linalg.norm(k_vec)

    data = build_kelvin_with_dual_info(N=2)
    r = build_both_spectra(data, k_vec)
    S1 = r["S1"]

    n0_ex = np.sum(r["vals_ex"] < 1e-12)
    n0_st = np.sum(r["vals_st"] < 1e-12)
    n_lost = n0_ex - n0_st

    print(f"  Kelvin: n0_ex={n0_ex}, n0_st={n0_st}, n_lost={n_lost}")
    print(f"  {'mode':>5s} {'omega/k':>10s} {'gauge_frac':>11s} {'phys_frac':>10s}")

    gauge_frac_0 = None
    for i in range(n0_st, min(n0_st + 8, len(r["vals_st"]))):
        v = r["vecs_st"][:, i]
        omega_k = np.sqrt(max(r["vals_st"][i], 0)) / k_mag

        # Project onto exact gauge subspace
        gauge_overlap = sum(
            abs(r["vecs_ex"][:, j].conj() @ S1 @ v) ** 2
            for j in range(n0_ex)
        )
        norm_sq = np.real(v.conj() @ S1 @ v)
        gf = gauge_overlap / norm_sq
        pf = 1.0 - gf

        print(f"  {i - n0_st:5d} {omega_k:10.4f} {gf:11.6f} {pf:10.6f}")

        if i == n0_st:
            gauge_frac_0 = gf
        elif i < n0_st + n_lost:
            assert gf > 0.95, (
                f"Mode {i - n0_st} gauge_frac = {gf:.4f}, expected > 0.95")

    # Lowest mode: should be 70-98% gauge (mixture, not pure)
    assert 0.7 < gauge_frac_0 < 0.98, (
        f"Mode 0 gauge_frac = {gauge_frac_0:.4f}, expected 0.7-0.98 (mixture)")
    # Physical fraction should be 2-30%
    phys_frac_0 = 1.0 - gauge_frac_0
    assert 0.02 < phys_frac_0 < 0.30, (
        f"Mode 0 phys_frac = {phys_frac_0:.4f}, expected 0.02-0.30")

    print("  PASSED")


def test_R20_N_dependence():
    """R20: c_std depends on supercell size N non-monotonically.

    On exact DEC, c = 1 for all N (physical, mesh-independent). On standard,
    c_std varies: 1.25 (N=2), 1.58 (N=3), 1.39 (N=4). Non-monotonic means
    there's no convergence — the error is structural, not a finite-size artifact.
    """
    print("\n=== R20: c_standard vs supercell size N ===")

    k_mag = 0.005
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    c_values = []
    for N in [2, 3, 4]:
        data = build_kelvin_with_dual_info(N=N)
        star1, star2 = build_hodge_stars_voronoi(data)
        V, E, F = data["V"], data["E"], data["F"]
        L = data["L"]
        S1 = np.diag(star1)
        S2 = np.diag(star2)

        d1k_st = build_d1_bloch_standard(V, E, F, L, k_vec)
        if sparse.issparse(d1k_st):
            d1k_st = d1k_st.toarray()
        K_st = d1k_st.conj().T @ S2 @ d1k_st
        K_st = (K_st + K_st.conj().T) / 2
        vals_st = np.sort(np.real(eigh(K_st, S1, eigvals_only=True)))

        nz = vals_st[vals_st > 1e-10]
        c_st = np.sqrt(nz[0]) / k_mag
        c_values.append(c_st)
        print(f"  N={N}: nV={len(V):5d}, c_std = {c_st:.6f}")

    # c_std should vary with N (not constant)
    spread = max(c_values) - min(c_values)
    assert spread > 0.1, (
        f"c_std spread = {spread:.4f}, expected > 0.1 (N-dependent)")

    # Should NOT be monotonic
    monotone_up = all(c_values[i] <= c_values[i + 1] for i in range(len(c_values) - 1))
    monotone_down = all(c_values[i] >= c_values[i + 1] for i in range(len(c_values) - 1))
    assert not (monotone_up or monotone_down), (
        f"c_std is monotonic in N (unexpected): {c_values}")

    print(f"  Spread = {spread:.4f} (N-dependent, non-monotonic)")
    print("  PASSED")


def test_R27_TRIM_collapse():
    """R27: At unit-cell TRIM (k = π/L_cell), standard becomes exact.

    At TRIM points, Bloch phases e^{ik·Δx} are exactly ±1 on unit-cell edges.
    This makes d₁(k)d₀(k) = 0 automatically — the standard DEC inherits
    exactness at these isolated k-points. ||d₁d₀|| drops from O(1) to ~1e-15,
    and eigenvalues match to machine precision.
    """
    print("\n=== R27: TRIM collapse — standard becomes exact ===")

    N = 2
    data = build_kelvin_with_dual_info(N=N)
    V, E, F = data["V"], data["E"], data["F"]
    L = data["L"]
    L_vec = np.array(data["L_vec"])
    L_cell = L / N

    star1, star2 = build_hodge_stars_voronoi(data)
    S1 = np.diag(star1)
    S2 = np.diag(star2)

    # Unit-cell TRIM: k = π/L_cell
    scale = np.pi / L_cell
    trim_points = {
        "X_cell": np.array([1, 0, 0], dtype=float) * scale,
        "M_cell": np.array([1, 1, 0], dtype=float) * scale,
        "R_cell": np.array([1, 1, 1], dtype=float) * scale,
    }

    # Generic k for contrast
    generic_k = 0.1 * (2 * np.pi / L) * np.array([1, 0, 0], dtype=float)
    d0g = build_d0_bloch(V, E, L, generic_k)
    if sparse.issparse(d0g): d0g = d0g.toarray()
    d1g_st = build_d1_bloch_standard(V, E, F, L, generic_k)
    if sparse.issparse(d1g_st): d1g_st = d1g_st.toarray()
    d1d0_gen = np.linalg.norm(d1g_st @ d0g)
    print(f"  Generic k: ||d1d0||_std = {d1d0_gen:.2e} (should be >> 0)")
    assert d1d0_gen > 1.0, "Generic k should have ||d1d0|| >> 0"

    for tname, k_vec in trim_points.items():
        d0k = build_d0_bloch(V, E, L, k_vec)
        if sparse.issparse(d0k): d0k = d0k.toarray()

        d1k_st = build_d1_bloch_standard(V, E, F, L, k_vec)
        if sparse.issparse(d1k_st): d1k_st = d1k_st.toarray()

        d1k_ex = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
        if sparse.issparse(d1k_ex): d1k_ex = d1k_ex.toarray()

        # Check exactness
        d1d0_st = np.linalg.norm(d1k_st @ d0k)

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

        assert d1d0_st < 1e-12, (
            f"{tname}: ||d1d0||_std = {d1d0_st:.2e}, expected < 1e-12")
        assert maxdiff < 1e-12, (
            f"{tname}: eigenvalue mismatch = {maxdiff:.2e}")

    print("  PASSED")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    test_R16_mode_mixing()
    test_R20_N_dependence()
    test_R27_TRIM_collapse()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED (App B — 3 tests)")
    print("=" * 60)
