"""
Paper §4: Removing exactness — what breaks when d₁d₀ ≠ 0.

Supports: §4 (Standard DEC failure modes)

CLAIM: On the SAME Voronoi complex with the SAME Hodge stars, replacing the exact
d₁ (which satisfies d₁d₀ = 0) with the standard DEC d₁ (which does NOT) breaks
c = 1. The failure has rich structure: c_std = 1.25–1.68, direction-dependent,
O(1) not O(k²), no consistent principal symbol, and n_lost = Δrank(d₁).

RAW OUTPUT (7 tests, all pass):
=================================
=== R12: c_standard on cubic ===
  Kelvin: c_exact=1.000, c_std=1.253, n_lost=6, ||d₁d₀||_ex=4e-16, ||d₁d₀||_std=0.50
  C15: c_exact=1.000, c_std=1.685, n_lost=39, ||d₁d₀||_ex=8e-16, ||d₁d₀||_std=1.1
  WP: c_exact=1.000, c_std=1.482, n_lost=15, ||d₁d₀||_ex=4e-16, ||d₁d₀||_std=0.69
=== R25: δc = O(1) ===
  δc ≈ 0.253 constant (ratio 1.16 across 20× k range), δc/k² → ∞
=== R35: No principal symbol ===
  Kelvin: C eigs [0.14,0.60,1.32], RMS=0.390
  C15: C eigs [0.68,1.93,2.09], RMS=0.841
  WP: C eigs [0.49,1.53,1.64], RMS=0.632
=== R17: Standard anisotropy ===
  Kelvin 70.2%, C15 42.8%, WP 36.6% (exact < 0.001%)
=== R40: n_lost = Δrank(d₁) ===
  Identity exact on 3 structures × 3 directions. rank(d₁_ex) = nV always.
=== R41: Leaked forms = pure gradients ===
  grad_frac ≥ 0.885 (Kelvin), ≥ 0.980 (C15), ≥ 0.992 (WP)
=== R23: Pollution by direction ===
  Kelvin: axis 6–7, [110]=12, [111]=14. C15: 39→70→91. WP: 15→28→37.
ALL TESTS PASSED (§4 removing exactness — 7 tests)

ANSWER:
=======
Standard DEC (d₁d₀ ≠ 0) breaks c = 1 on the same mesh with the same Hodge stars.
c_std = 1.25–1.68 (structure-dependent). The error is O(1) not O(k²), with no
consistent principal symbol (RMS 0.39–0.84), catastrophic anisotropy (37–70%),
and direction-dependent pollution (n_lost = 6–91). The mechanism: n_lost = Δrank(d₁),
leaked forms are pure gradients (≥ 88.5%).
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
    build_foam_with_dual_info,
    build_hodge_stars_voronoi,
)
from physics.gauge_bloch import build_d1_bloch_exact
from physics.bloch import build_d0_bloch, build_d1_bloch_standard

BUILDERS = [
    ("Kelvin", build_kelvin_with_dual_info),
    ("C15", build_c15_with_dual_info),
    ("WP", build_wp_with_dual_info),
]


# ===================================================================
# Helpers
# ===================================================================

def build_both(data, k_vec):
    """Build exact and standard operators at k, return spectra + extras."""
    star1, star2 = build_hodge_stars_voronoi(data)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    S1 = np.diag(star1)
    S2 = np.diag(star2)

    d0k = build_d0_bloch(V, E, L, k_vec)
    if sparse.issparse(d0k):
        d0k = d0k.toarray()

    d1k_ex = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
    if sparse.issparse(d1k_ex):
        d1k_ex = d1k_ex.toarray()

    d1k_st = build_d1_bloch_standard(V, E, F, L, k_vec)
    if sparse.issparse(d1k_st):
        d1k_st = d1k_st.toarray()

    def spectrum(d1k):
        K = d1k.conj().T @ S2 @ d1k
        K = (K + K.conj().T) / 2
        vals, vecs = eigh(K, S1)
        idx = np.argsort(np.real(vals))
        return np.real(vals[idx]), vecs[:, idx], K

    vals_ex, vecs_ex, K_ex = spectrum(d1k_ex)
    vals_st, vecs_st, K_st = spectrum(d1k_st)

    return {
        "vals_ex": vals_ex, "vecs_ex": vecs_ex, "K_ex": K_ex,
        "vals_st": vals_st, "vecs_st": vecs_st, "K_st": K_st,
        "S1": S1, "S2": S2, "d0k": d0k,
        "d1k_ex": d1k_ex, "d1k_st": d1k_st,
    }


# ===================================================================
# Tests
# ===================================================================

def test_I6_standard_speed():
    """R12: c_std = 1.25–1.68 on cubic structures (standard DEC gives c ≠ 1).

    Same Hodge stars, same mesh. Only d₁ differs: exact (d₁d₀=0) vs standard (d₁d₀≠0).
    Standard DEC loses gauge modes (n_lost > 0), lowest physical mode has c > 1.
    """
    print("\n=== R12: c_standard on cubic ===")

    k_mag = 0.005
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    for name, builder in BUILDERS:
        data = builder(N=2)
        r = build_both(data, k_vec)

        # Verify exactness violation directly
        d1d0_ex = np.linalg.norm(r["d1k_ex"] @ r["d0k"])
        d1d0_st = np.linalg.norm(r["d1k_st"] @ r["d0k"])
        assert d1d0_ex < 1e-12, f"{name}: exact d₁d₀ = {d1d0_ex:.2e}"
        assert d1d0_st > 0.01, f"{name}: standard d₁d₀ = {d1d0_st:.2e}, expected >> 0"

        n0_ex = int(np.sum(r["vals_ex"] < 1e-12))
        n0_st = int(np.sum(r["vals_st"] < 1e-12))
        n_lost = n0_ex - n0_st

        nz_ex = r["vals_ex"][r["vals_ex"] > 1e-10]
        nz_st = r["vals_st"][r["vals_st"] > 1e-10]
        c_ex = np.sqrt(nz_ex[0]) / k_mag
        c_st = np.sqrt(nz_st[0]) / k_mag

        print(f"  {name:8s}: c_exact={c_ex:.6f}, c_std={c_st:.4f}, n_lost={n_lost}, "
              f"||d₁d₀||_ex={d1d0_ex:.1e}, ||d₁d₀||_std={d1d0_st:.1e}")

        assert abs(c_ex - 1.0) < 1e-4, f"{name}: c_exact = {c_ex}"
        assert c_st > 1.1, f"{name}: c_std = {c_st}, expected > 1.1"
        assert n_lost > 0, f"{name}: no lost gauge modes"

    print("  PASSED")


def test_R25_delta_c_O1():
    """R25: δc = c_std − c_exact is O(1), not O(k²). Structural error.

    If δc were O(k²) it would vanish in the continuum limit. Instead δc ≈ 0.25
    constant across 20× range in k. This means standard DEC has a structural
    speed error that does NOT improve with mesh refinement.
    """
    print("\n=== R25: δc = O(1) ===")

    data = build_kelvin_with_dual_info(N=2)
    L = data["L"]
    scale = 2 * np.pi / L

    k_mags = [0.005, 0.01, 0.02, 0.05, 0.1]
    delta_c = []

    for km in k_mags:
        k_vec = km * scale * np.array([1.0, 0.0, 0.0])
        k2 = k_vec @ k_vec
        r = build_both(data, k_vec)

        nz_ex = r["vals_ex"][r["vals_ex"] > 1e-8]
        nz_st = r["vals_st"][r["vals_st"] > 1e-8]
        c_ex = np.sqrt(nz_ex[0] / k2)
        c_st = np.sqrt(nz_st[0] / k2)
        dc = abs(c_st - c_ex)
        delta_c.append(dc)
        print(f"  |k|={km:.3f}: δc={dc:.4f}, δc/k²={dc/k2:.0f}")

    dc_arr = np.array(delta_c)
    ratio = dc_arr.max() / dc_arr.min()
    print(f"  Range ratio = {ratio:.2f} (O(1) → ratio < 2)")

    assert ratio < 2.0, f"δc ratio = {ratio:.2f}, expected < 2"
    assert dc_arr.min() > 0.1, f"δc min = {dc_arr.min():.4f}, expected > 0.1"
    print("  PASSED")


def test_N16_no_principal_symbol():
    """R35: Standard DEC has no quadratic principal symbol c²(k̂) = k̂ᵀCk̂.

    Probe 17 directions, fit c²_std to quadratic form. RMS residual is huge
    (0.3–0.8) because n_lost jumps by direction class → "lowest physical mode"
    changes identity. c²(k̂) depends on direction non-quadratically — no single
    constant tensor C describes the standard operator.
    """
    print("\n=== R35: No principal symbol ===")

    # 17 probe directions: 3 axes + 6 face diags + 4 body diags + 4 generic
    dirs_raw = []
    for i in range(3):
        d = np.zeros(3); d[i] = 1.0; dirs_raw.append(d)
    for i in range(3):
        for j in range(i + 1, 3):
            for s in [+1, -1]:
                d = np.zeros(3); d[i] = 1.0; d[j] = s
                dirs_raw.append(d / np.linalg.norm(d))
    for sx in [+1, -1]:
        for sy in [+1, -1]:
            d = np.array([1.0, sx, sy])
            dirs_raw.append(d / np.linalg.norm(d))
    np.random.seed(123)
    for _ in range(4):
        d = np.random.randn(3)
        dirs_raw.append(d / np.linalg.norm(d))

    k_mag = 0.01

    for name, builder in BUILDERS:
        data = builder(N=2)
        c2_vals = []
        k_hat_list = []

        for d in dirs_raw:
            k_vec = k_mag * d
            r = build_both(data, k_vec)
            nz = r["vals_st"][r["vals_st"] > 1e-10]
            c2 = nz[0] / (k_mag**2)
            c2_vals.append(c2)
            k_hat_list.append(d)

        # Fit c² = k̂ᵀ C k̂ via least squares
        # Design matrix: [kx², ky², kz², 2kxky, 2kxkz, 2kykz]
        A_fit = np.zeros((len(dirs_raw), 6))
        for i, kh in enumerate(k_hat_list):
            A_fit[i] = [kh[0]**2, kh[1]**2, kh[2]**2,
                        2*kh[0]*kh[1], 2*kh[0]*kh[2], 2*kh[1]*kh[2]]
        c2_arr = np.array(c2_vals)
        coeffs, res, _, _ = np.linalg.lstsq(A_fit, c2_arr, rcond=None)
        c2_pred = A_fit @ coeffs
        rms = np.sqrt(np.mean((c2_arr - c2_pred)**2))

        C_mat = np.array([
            [coeffs[0], coeffs[3], coeffs[4]],
            [coeffs[3], coeffs[1], coeffs[5]],
            [coeffs[4], coeffs[5], coeffs[2]],
        ])
        C_eigs = np.sort(np.linalg.eigvalsh(C_mat))

        print(f"  {name:8s}: C eigs [{C_eigs[0]:.2f},{C_eigs[1]:.2f},{C_eigs[2]:.2f}], "
              f"RMS={rms:.3f}")

        assert rms > 0.1, f"{name}: RMS = {rms:.3f}, fit too good for 'no symbol'"

    print("  PASSED")


def test_R17_anisotropy():
    """R17: Standard DEC anisotropy 37–70% on cubic structures.

    c_standard varies wildly with k-direction (4 directions tested).
    Exact DEC: isotropic to < 0.01%. The anisotropy comes from how Bloch phases
    interact with the face-edge topology when d₁d₀ ≠ 0.
    """
    print("\n=== R17: Standard anisotropy ===")

    k_mag = 0.005
    dirs = [
        np.array([1, 0, 0], dtype=float),
        np.array([0, 1, 0], dtype=float),
        np.array([0, 0, 1], dtype=float),
        np.array([1, 1, 1], dtype=float) / np.sqrt(3),
    ]

    for name, builder in BUILDERS:
        data = builder(N=2)
        speeds_st = []
        speeds_ex = []

        for d in dirs:
            k_vec = k_mag * d
            r = build_both(data, k_vec)
            nz_st = r["vals_st"][r["vals_st"] > 1e-10]
            nz_ex = r["vals_ex"][r["vals_ex"] > 1e-10]
            speeds_st.append(np.sqrt(nz_st[0]) / k_mag)
            speeds_ex.append(np.sqrt(nz_ex[0]) / k_mag)

        aniso_st = (max(speeds_st) - min(speeds_st)) / np.mean(speeds_st) * 100
        aniso_ex = (max(speeds_ex) - min(speeds_ex)) / np.mean(speeds_ex) * 100
        print(f"  {name:8s}: std {aniso_st:.1f}% (c={min(speeds_st):.3f}–{max(speeds_st):.3f}), "
              f"exact {aniso_ex:.4f}%")

        assert aniso_st > 30, f"{name}: std anisotropy = {aniso_st:.1f}%"
        assert aniso_ex < 0.1, f"{name}: exact anisotropy = {aniso_ex:.4f}%"

    print("  PASSED")


def test_N19_nlost_identity():
    """R40: n_lost = rank(d₁_std) − rank(d₁_ex) — exact algebraic identity.

    On exact complex: rank(d₁_ex) = nV (exactness constrains kernel to im(d₀)).
    Standard releases this constraint → higher rank. The difference counts
    exactly how many gauge modes leak into the physical spectrum.
    """
    print("\n=== R40: n_lost = Δrank(d₁) ===")

    k_mag = 0.01
    dirs = [
        ("[1,0,0]", np.array([1, 0, 0], dtype=float)),
        ("[1,1,0]", np.array([1, 1, 0], dtype=float) / np.sqrt(2)),
        ("[1,1,1]", np.array([1, 1, 1], dtype=float) / np.sqrt(3)),
    ]

    for name, builder in BUILDERS:
        data = builder(N=2)
        nV = len(data["V"])
        L = data["L"]

        for dname, d in dirs:
            k_vec = k_mag * (2 * np.pi / L) * d
            r = build_both(data, k_vec)

            rk_ex = np.linalg.matrix_rank(r["d1k_ex"], tol=1e-8)
            rk_st = np.linalg.matrix_rank(r["d1k_st"], tol=1e-8)
            delta_rank = rk_st - rk_ex

            n0_ex = int(np.sum(r["vals_ex"] < 1e-10))
            n0_st = int(np.sum(r["vals_st"] < 1e-10))
            n_lost = n0_ex - n0_st

            print(f"  {name:6s} {dname}: rank_ex={rk_ex}, rank_st={rk_st}, "
                  f"Δrank={delta_rank}, n_lost={n_lost}")

            assert rk_ex == nV, f"{name} {dname}: rank(d₁_ex) = {rk_ex} ≠ nV = {nV}"
            assert delta_rank == n_lost, (
                f"{name} {dname}: Δrank = {delta_rank} ≠ n_lost = {n_lost}")
            assert n_lost > 0, f"{name} {dname}: no lost modes"

    print("  PASSED")


def test_N13_leaked_gradients():
    """R41: Leaked forms ⊂ im(d₀) — gauge pollution is pure gradient leakage.

    The n_lost modes that are in ker(d₁_ex) but not in ker(d₁_std) are predominantly
    gradients. At finite k, mode mixing with physical modes reduces the gradient
    fraction slightly (Kelvin: 88.5%, C15: 98.0%, WP: 99.2%). In the k → 0 limit,
    the leaked modes become pure gradients (frac → 1).
    """
    print("\n=== R41: Leaked forms are pure gradients ===")

    k_mag = 0.01
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    for name, builder in BUILDERS:
        data = builder(N=2)
        r = build_both(data, k_vec)
        S1 = r["S1"]
        d0k = r["d0k"]

        # Gradient projector: P_grad = d₀(d₀†S₁d₀)⁻¹d₀†S₁
        G0 = d0k.conj().T @ S1 @ d0k
        G0inv = np.linalg.inv(G0)
        P_grad = d0k @ G0inv @ d0k.conj().T @ S1

        n0_ex = int(np.sum(r["vals_ex"] < 1e-12))
        n0_st = int(np.sum(r["vals_st"] < 1e-12))
        n_lost = n0_ex - n0_st

        # Leaked modes: first n_lost nonzero modes on standard
        min_frac = 1.0
        for i in range(n0_st, n0_st + n_lost):
            v = r["vecs_st"][:, i]
            pv = P_grad @ v
            gf = np.real(pv.conj() @ S1 @ pv) / np.real(v.conj() @ S1 @ v)
            min_frac = min(min_frac, gf)

        # Physical modes on exact: should be pure curl (0% gradient)
        max_phys_grad = 0.0
        for i in range(n0_ex, n0_ex + 2):
            v = r["vecs_ex"][:, i]
            pv = P_grad @ v
            gf = np.real(pv.conj() @ S1 @ pv) / np.real(v.conj() @ S1 @ v)
            max_phys_grad = max(max_phys_grad, gf)

        print(f"  {name:8s}: n_lost={n_lost}, leaked grad_frac ≥ {min_frac:.4f}, "
              f"exact phys grad_frac ≤ {max_phys_grad:.1e}")

        assert min_frac > 0.85, f"{name}: leaked mode not gradient: frac = {min_frac:.4f}"
        assert max_phys_grad < 0.01, f"{name}: exact mode has gradient: frac = {max_phys_grad:.4f}"

    print("  PASSED")


def test_R23_pollution_direction():
    """R23: n_lost depends on k-direction. On-axis < face diagonal < body diagonal.

    Off-axis k crosses more periodic boundaries → more phase errors in d₁d₀.
    """
    print("\n=== R23: Pollution by direction ===")

    directions = {
        "[1,0,0]": np.array([1, 0, 0], dtype=float),
        "[0,1,0]": np.array([0, 1, 0], dtype=float),
        "[0,0,1]": np.array([0, 0, 1], dtype=float),
        "[1,1,0]": np.array([1, 1, 0], dtype=float) / np.sqrt(2),
        "[1,1,1]": np.array([1, 1, 1], dtype=float) / np.sqrt(3),
    }

    for name, builder in BUILDERS:
        data = builder(N=2)
        L = data["L"]

        n_lost_axis = []
        n_lost_diag = []

        results = []
        for dname, d in directions.items():
            k_vec = 0.01 * (2 * np.pi / L) * d
            r = build_both(data, k_vec)
            thresh = 1e-3 * (k_vec @ k_vec)
            n0_ex = int(np.sum(r["vals_ex"] < thresh))
            n0_st = int(np.sum(r["vals_st"] < thresh))
            n_lost = n0_ex - n0_st
            results.append((dname, n_lost))

            if dname in ("[1,0,0]", "[0,1,0]", "[0,0,1]"):
                n_lost_axis.append(n_lost)
            else:
                n_lost_diag.append(n_lost)

        nlost_str = ", ".join(f"{d}={n}" for d, n in results)
        print(f"  {name:8s}: {nlost_str}")

        assert min(n_lost_diag) > max(n_lost_axis), (
            f"{name}: diag n_lost {min(n_lost_diag)} ≤ axis {max(n_lost_axis)}")

    print("  PASSED")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    test_I6_standard_speed()
    test_R25_delta_c_O1()
    test_N16_no_principal_symbol()
    test_R17_anisotropy()
    test_N19_nlost_identity()
    test_N13_leaked_gradients()
    test_R23_pollution_direction()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED (§4 removing exactness — 7 tests)")
    print("=" * 60)
