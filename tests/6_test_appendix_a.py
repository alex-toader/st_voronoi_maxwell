"""
Appendix A: Why naive perturbation theory fails.

Supports: Appendix A (supporting material for §3)

CLAIM: The acoustic eigenvector at k≠0 is >99.99% harmonic, yet the harmonic
Rayleigh quotient R[h]/k² = 4.8–10.7 (NOT 1). The eigenvalue ω²/k² = 1 emerges
from three-way cancellation: h†Kh + 2Re(h†Kδψ) + δψ†Kδψ = k². Projecting out
the gradient component makes R WORSE (not better). The Schur complement resolves
the paradox via block elimination, not perturbation theory.

RAW OUTPUT (5 tests, all pass):
=================================
=== R44a: eigenvector is harmonic ===
  Kelvin: harm=0.999995  grad=0.000001  opt=0.000004
  C15:    harm=0.999995  grad=0.000000  opt=0.000005
  WP:     harm=0.999995  grad=0.000000  opt=0.000005
=== R44b: Rayleigh paradox ===
  Kelvin: H3_trans/k² = [4.840, 5.333]  ω²/k² = 1.000
  C15:    H3_trans/k² = [9.392, 10.747]  ω²/k² = 1.000
  WP:     H3_trans/k² = [6.705, 6.968]  ω²/k² = 1.000
=== R44c: projection makes Rayleigh worse ===
  Kelvin: R[h]/k²=4.840  R[h_perp]/k²=4.849  (increased)
  C15:    R[h]/k²=9.392  R[h_perp]/k²=9.412  (increased)
  WP:     R[h]/k²=6.705  R[h_perp]/k²=6.774  (increased)
=== R44d: cross-term anatomy ===
  Kelvin: [4.8, -7.7, 3.8] → 1.000
  C15:    [9.4, -16.8, 8.4] → 1.000
  WP:     [6.7, -11.4, 5.7] → 1.000
=== R11: plane wave overlap (Kelvin) ===
  k=0.1: 0.947, k=0.05: 0.987, k=0.01: 0.999, k=0.005: 1.000, k=0.001: 1.000
ALL TESTS PASSED (App A — 5 tests)

ANSWER:
=======
Naive perturbation theory fails because the harmonic subspace is not K(k)-invariant.
The correction δψ is tiny in norm (~0.001%) but contains optical components amplified
by K. The Schur complement integrates out optical modes exactly, giving c² = 1 without
perturbative expansion.
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
from physics.bloch import build_d0_bloch

BUILDERS = [
    ("Kelvin", build_kelvin_with_dual_info),
    ("C15", build_c15_with_dual_info),
    ("WP", build_wp_with_dual_info),
]


# ===================================================================
# Helpers
# ===================================================================

def build_all(data, k_vec):
    """Build Hodge stars, operators at Gamma and at k. Return dict."""
    star1, star2 = build_hodge_stars_voronoi(data)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    S1 = np.diag(star1)
    S2 = np.diag(star2)
    nE = len(E)
    nV = len(V)
    Vol = L**3

    # Edge vectors (periodic wrapping)
    dx = np.zeros((nE, 3))
    for ie, (v1, v2) in enumerate(E):
        dv = np.array(V[v2]) - np.array(V[v1])
        dv -= np.round(dv / L) * L
        dx[ie] = dv
    harm = dx / np.sqrt(Vol)  # M-orthonormal harmonic forms

    # Gamma operators
    k0 = np.zeros(3)
    d0g = build_d0_bloch(V, E, L, k0)
    if sparse.issparse(d0g): d0g = d0g.toarray()
    d1g = build_d1_bloch_exact(V, E, F, k0, L_vec, d0g)
    if sparse.issparse(d1g): d1g = d1g.toarray()
    Kg = d1g.conj().T @ S2 @ d1g
    Kg = (Kg + Kg.conj().T) / 2
    vals_g, vecs_g = eigh(Kg, S1)
    n_zero = int(np.sum(np.real(vals_g) < 1e-8))

    # Operators at k
    d0k = build_d0_bloch(V, E, L, k_vec)
    if sparse.issparse(d0k): d0k = d0k.toarray()
    d1k = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
    if sparse.issparse(d1k): d1k = d1k.toarray()
    Kk = d1k.conj().T @ S2 @ d1k
    Kk = (Kk + Kk.conj().T) / 2
    vals_k, vecs_k = eigh(Kk, S1)

    return {
        "S1": S1, "S2": S2, "nE": nE, "nV": nV, "Vol": Vol,
        "harm": harm, "d0k": d0k,
        "Kg": Kg, "vals_g": np.real(vals_g), "vecs_g": vecs_g, "n_zero": n_zero,
        "Kk": Kk, "vals_k": np.real(vals_k), "vecs_k": vecs_k,
    }


# ===================================================================
# Tests
# ===================================================================

def test_R44a_eigvec_harmonic():
    """R44a: Acoustic eigenvector at k≠0 is >99.99% harmonic.

    Decompose ψ(k) into Γ-eigenbasis: harm_frac = |projection onto 3 harmonic
    forms|² / ||ψ||². The eigenvector IS harmonic — but harmonicity doesn't
    determine the eigenvalue (that's the paradox R44b addresses).
    """
    print("\n=== R44a: eigenvector is harmonic ===")

    k_vec = 0.001 * np.array([1., 0., 0.])
    results = []

    for name, builder in BUILDERS:
        b = build_all(builder(N=2), k_vec)
        S1, harm = b["S1"], b["harm"]
        vecs_g, n_zero = b["vecs_g"], b["n_zero"]
        vals_k, vecs_k = b["vals_k"], b["vecs_k"]

        nz_idx = np.where(vals_k > 1e-10)[0]
        fracs = []
        for i in range(2):  # 2 transverse acoustic
            psi = vecs_k[:, nz_idx[i]]
            pn2 = np.real(psi.conj() @ S1 @ psi)

            hc = harm.T @ S1 @ psi
            hn = np.sum(np.abs(hc)**2)

            Z = vecs_g[:, :n_zero]
            zc = Z.conj().T @ S1 @ psi
            zn = np.sum(np.abs(zc)**2)

            fracs.append((hn / pn2, (zn - hn) / pn2, (pn2 - zn) / pn2))

        hf = np.mean([f[0] for f in fracs])
        gf = np.mean([f[1] for f in fracs])
        of = np.mean([f[2] for f in fracs])
        print(f"  {name:8s}: harm={hf:.6f}  grad={gf:.6f}  opt={of:.6f}")
        results.append((name, hf))

    for name, hf in results:
        assert hf > 0.9999, f"{name}: harm_frac={hf:.6f}"

    print("  PASSED")


def test_R44b_rayleigh_paradox():
    """R44b: R[h_transverse]/k² = 4.8–10.7, yet ω²/k² = 1.

    The 3×3 harmonic block H₃ = harm†K(k)harm has eigenvalues all ≫ k².
    The transverse Rayleigh quotient R/k² > 2 on every structure.
    Yet the exact eigenvalue ω²/k² = 1.000. This is the paradox that
    motivates the Schur complement approach.
    """
    print("\n=== R44b: Rayleigh paradox ===")

    k_small = 0.001
    k_vec = k_small * np.array([1., 0., 0.])
    k2 = k_small**2

    for name, builder in BUILDERS:
        b = build_all(builder(N=2), k_vec)
        H3 = np.real(b["harm"].conj().T @ b["Kk"] @ b["harm"])
        h3_eigs = np.sort(np.linalg.eigvalsh(H3))
        h3_trans = h3_eigs[1:]  # 2 transverse

        actual = b["vals_k"][b["vals_k"] > 1e-10][0]

        print(f"  {name:8s}: H3_trans/k² = [{h3_trans[0]/k2:.3f}, {h3_trans[1]/k2:.3f}]"
              f"  ω²/k² = {actual/k2:.3f}")

        for j, ht in enumerate(h3_trans):
            assert ht / k2 > 2.0, f"{name}: H3_trans[{j}]/k² = {ht/k2:.3f}"
        assert abs(actual / k2 - 1.0) < 1e-4, f"{name}: ω²/k² = {actual/k2:.6f}"

    print("  PASSED")


def test_R44c_projection_worse():
    """R44c: Projecting h onto M⊥(im d₀(k)) increases R.

    K kills gradients (K·grad = 0), so h†Kh = h_perp†Kh_perp (numerator
    unchanged). But ||h_perp||² < ||h||² (denominator shrinks). Therefore
    R[h_perp] ≥ R[h]. Gradient removal makes the Rayleigh quotient WORSE.
    """
    print("\n=== R44c: projection makes Rayleigh worse ===")

    k_small = 0.001
    k_vec = k_small * np.array([1., 0., 0.])
    k2 = k_small**2

    for name, builder in BUILDERS:
        b = build_all(builder(N=2), k_vec)
        S1, harm, Kk, d0k = b["S1"], b["harm"], b["Kk"], b["d0k"]

        MG = d0k.conj().T @ S1 @ d0k
        MG_inv = np.linalg.inv(MG)

        H3 = np.real(harm.conj().T @ Kk @ harm)
        h3_evecs = np.linalg.eigh(H3)[1]
        harm_eig = harm @ h3_evecs

        for j in [1, 2]:  # transverse
            h = harm_eig[:, j]
            grad_coeff = MG_inv @ (d0k.conj().T @ S1 @ h)
            h_perp = h - d0k @ grad_coeff

            R_h = np.real(h.conj() @ Kk @ h) / np.real(h.conj() @ S1 @ h)
            R_perp = np.real(h_perp.conj() @ Kk @ h_perp) / np.real(h_perp.conj() @ S1 @ h_perp)

            if j == 1:
                print(f"  {name:8s}: R[h]/k²={R_h/k2:.3f}  R[h_perp]/k²={R_perp/k2:.3f}"
                      f"  ({'increased' if R_perp > R_h else 'decreased'})")

            assert R_perp >= R_h - 1e-10, (
                f"{name} mode {j}: R_perp < R_h")

    print("  PASSED")


def test_R44d_cancellation():
    """R44d: Three-way cancellation h†Kh + 2Re(h†Kδψ) + δψ†Kδψ = k².

    δψ is tiny (~0.001% of ||ψ||) but contains optical components with large
    eigenvalues (ω²_opt ~ 0.5), so K amplifies it. All three terms are O(k²)
    and none is negligible. The cancellation [~5, ~-8, ~4] → 1 is the mechanism.
    """
    print("\n=== R44d: cross-term anatomy ===")

    k_small = 0.001
    k_vec = k_small * np.array([1., 0., 0.])
    k2 = k_small**2

    for name, builder in BUILDERS:
        b = build_all(builder(N=2), k_vec)
        S1, harm, Kk = b["S1"], b["harm"], b["Kk"]
        vals_k, vecs_k = b["vals_k"], b["vecs_k"]

        nz_idx = np.where(vals_k > 1e-10)[0]
        psi = vecs_k[:, nz_idx[0]]

        hc = harm.T @ S1 @ psi
        h_proj = harm @ hc
        delta = psi - h_proj

        hKh = np.real(h_proj.conj() @ Kk @ h_proj)
        cross = 2 * np.real(h_proj.conj() @ Kk @ delta)
        dKd = np.real(delta.conj() @ Kk @ delta)
        pKp = np.real(psi.conj() @ Kk @ psi)
        pMp = np.real(psi.conj() @ S1 @ psi)

        print(f"  {name:8s}: h†Kh/k²={hKh/k2:.3f}  2Re(h†Kδψ)/k²={cross/k2:.3f}"
              f"  δψ†Kδψ/k²={dKd/k2:.3f}  sum/k²={(hKh+cross+dKd)/pMp/k2:.3f}")

        assert hKh / k2 > 2.0, f"{name}: h†Kh too small"
        assert cross / k2 < -2.0, f"{name}: cross term not large negative"
        assert dKd / k2 > 0.5, f"{name}: δψ†Kδψ too small"
        total = (hKh + cross + dKd) / pMp
        assert abs(total / k2 - 1.0) < 1e-3, f"{name}: total/k² = {total/k2:.6f}"

    print("  PASSED")


def test_R11_plane_wave_overlap():
    """R11: Acoustic eigenspace coincides with plane-wave subspace as k → 0.

    Subspace overlap Tr(P_eig·P_pw)/2 between the 2 transverse eigenvectors
    and the 2 transverse plane waves (pol·Δx, M₁-normalized). Overlap → 1
    monotonically as k → 0, confirming eigenfunctions become plane waves.
    """
    print("\n=== R11: plane wave overlap ===")

    data = build_kelvin_with_dual_info(N=2)
    star1, star2 = build_hodge_stars_voronoi(data)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    S1 = np.diag(star1)
    S2 = np.diag(star2)
    nE = len(E)

    def plane_wave_trial(pol):
        a = np.zeros(nE)
        for e_idx, (i, j) in enumerate(E):
            dx = V[j] - V[i]
            dx -= np.round(dx / L_vec) * L_vec
            a[e_idx] = np.dot(pol, dx)
        norm = np.sqrt(a @ S1 @ a)
        return a / norm

    trial_y = plane_wave_trial(np.array([0, 1, 0]))
    trial_z = plane_wave_trial(np.array([0, 0, 1]))

    print(f"  {'k':>8s} {'overlap':>10s}")
    overlaps = []
    for k_mag in [0.1, 0.05, 0.01, 0.005, 0.001]:
        k_vec = k_mag * np.array([1.0, 0.0, 0.0])
        d0k = build_d0_bloch(V, E, L, k_vec)
        d1k = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
        if sparse.issparse(d1k): d1k = d1k.toarray()
        K = d1k.conj().T @ S2 @ d1k
        K = (K + K.conj().T) / 2

        vals, vecs = eigh(K, S1)
        idx = np.argsort(np.real(vals))
        vals, vecs = np.real(vals[idx]), vecs[:, idx]

        nz_idx = np.where(vals > 1e-10)[0][:2]
        ov = 0.0
        for vi_idx in nz_idx:
            vi = vecs[:, vi_idx]
            for tj in [trial_y, trial_z]:
                ov += abs(tj @ S1 @ vi)**2
        overlaps.append((k_mag, ov / 2))
        print(f"  {k_mag:8.4f} {ov/2:10.6f}")

    assert overlaps[-1][1] > 0.999, f"Overlap at k={overlaps[-1][0]}: {overlaps[-1][1]:.6f}"
    for i in range(1, len(overlaps)):
        assert overlaps[i][1] >= overlaps[i-1][1] - 0.01, "Overlap not monotonic"

    print("  PASSED")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    test_R44a_eigvec_harmonic()
    test_R44b_rayleigh_paradox()
    test_R44c_projection_worse()
    test_R44d_cancellation()
    test_R11_plane_wave_overlap()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED (App A — 5 tests)")
    print("=" * 60)
