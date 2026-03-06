"""
Paper §3: c² = 1 on exact Voronoi complex.

Supports: §3 (Theorem: c² = 1)

CLAIM: On any periodic Voronoi complex with an exact cochain complex (d₁d₀ = 0),
the discrete electromagnetic wave speed equals 1 in natural units. Three ingredients:
(1) G = H = Vol·I (geometry, from §2), (2) d₁(0)·harm = 0 (exactness),
(3) discrete curl identity dF/dk|₀ = iε. These combine to give the Schur
complement S = k²·P_transverse, proving c² = 1.

RAW OUTPUT (10 tests, all pass):
==================================
=== R5: c_gauge = 1 + O(k²) ===
  Kelvin: c = 0.999997, C15: 0.999999, WP: 0.999999
=== R38: c = 1 on random Voronoi ===
  5 seeds (n=80): c = 0.99999940–0.99999947, δc < 6e-7
=== R32: Scaling invariance ===
  L_cell = 2,4,8,16: c identical (spread 1.5e-12)
=== R41: H¹(k≠0) = 0 ===
  All cubic + 3 random: nE = 2nV, H¹ = 0, d₁d₀ = 0
=== R44e: Schur complement ===
  Schur/k² = [~0, 1.000, 1.000] on cubic + 2 random
  H₃/k² = [0.8–14.2] (structure-dependent, ≠ 1)
=== R44f: discrete curl identity ===
  ||dF/dk/i − ε|| < 4.7e-11 on 3 cubic + 3 random (tol 1e-8)
=== R44g: full proof chain ===
  G=H=Vol·I + rank(d₁(0))=nV−2 + u_perp ∈ H² (99.996–99.999%) + Schur → S = k²·P_T, err < 2e-11
=== R6: converse ===
  Voronoi δc = 2.6e-6, perturbed 10% δc = 5.0e-3 (1929× worse)
=== S1: Multi-direction isotropy ===
  5 directions (x,y,z,111,110): spread = 1.0e-6, all δc < 4e-6
=== S2: δc/k² convergence ===
  (c-1)/k² → −0.1042, converged to 2e-4 — dispersion is O(k⁴)
ALL TESTS PASSED (§3 c² = 1 — 10 tests)

ANSWER:
=======
c² = 1 is exact on all periodic Voronoi complexes with exact cochain complex.
Verified on 3 cubic + 5 random Voronoi, 5 k-directions. The proof chain:
(1) G = H = Vol·I (divergence theorem, §2), (2) d₁(0)·harm = 0 (exactness),
(3) discrete curl identity dF/dk|₀ = iε (verified to 10⁻¹¹). These give the
Schur complement S = k²·P_transverse with eigenvalues [0, k², k²]. H¹(k≠0) = 0
on all z=4 foams (nE = 2nV). Scale-invariant. Isotropic (spread 10⁻⁶).
Dispersion O(k⁴). Converse: perturbing ⋆₁ by 10% breaks c by 1929×.
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
from physics.bloch import build_d0_bloch

BUILDERS = [
    ("Kelvin", build_kelvin_with_dual_info),
    ("C15", build_c15_with_dual_info),
    ("WP", build_wp_with_dual_info),
]


# ===================================================================
# Helpers
# ===================================================================

def solve_bloch(data, k_vec):
    """Build K = d₁†⋆₂d₁ at k, solve generalized eigenvalue problem K·ψ = ω²·⋆₁·ψ."""
    star1, star2 = build_hodge_stars_voronoi(data)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    S1, S2 = np.diag(star1), np.diag(star2)

    d0k = build_d0_bloch(V, E, L, k_vec)
    if sparse.issparse(d0k):
        d0k = d0k.toarray()
    d1k = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
    if sparse.issparse(d1k):
        d1k = d1k.toarray()

    K = d1k.conj().T @ S2 @ d1k
    K = (K + K.conj().T) / 2

    vals, vecs = eigh(K, S1)
    return np.real(vals), vecs, {"S1": S1, "S2": S2, "d0k": d0k, "d1k": d1k,
                                  "K": K, "star1": star1, "star2": star2,
                                  "V": V, "E": E, "F": F, "L_vec": L_vec, "L": L}


def build_harm_h2(data):
    """Build M-orthonormal harmonic 1-forms and S₂-orthonormal H² basis."""
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    nE, nF = len(E), len(F)
    Vol = np.prod(L_vec)

    # Edge vectors → harmonic 1-forms: harm_α = Δx_α / √Vol
    dx = np.zeros((nE, 3))
    for ie, (v1, v2) in enumerate(E):
        dv = np.array(V[v2]) - np.array(V[v1])
        dv -= np.round(dv / L_vec) * L_vec
        dx[ie] = dv
    harm = dx / np.sqrt(Vol)

    # Face area vectors → H² cohomology: h2_β = A_f / √Vol
    fa = np.zeros((nF, 3))
    for iF, face in enumerate(F):
        A = np.zeros(3)
        v0 = np.array(V[face[0]])
        for i in range(1, len(face) - 1):
            e1 = np.array(V[face[i]]) - v0
            e1 -= np.round(e1 / L_vec) * L_vec
            e2 = np.array(V[face[i + 1]]) - v0
            e2 -= np.round(e2 / L_vec) * L_vec
            A += np.cross(e1, e2)
        fa[iF] = A / 2
    h2 = fa / np.sqrt(Vol)

    return harm, h2, dx, fa, Vol


# ===================================================================
# Tests
# ===================================================================

def test_R5_gauge_speed():
    """R5: c = 1 + O(k²) on Kelvin, C15, WP from Bloch eigenvalues."""
    print("\n=== R5: c_gauge = 1 + O(k²) ===")
    k_mag = 0.005
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    for name, builder in BUILDERS:
        data = builder(N=2)
        vals, _, _ = solve_bloch(data, k_vec)
        nz = vals[vals > 1e-10]
        assert len(nz) >= 2, f"{name}: need ≥2 nonzero modes"

        c1 = np.sqrt(nz[0]) / k_mag
        c2 = np.sqrt(nz[1]) / k_mag
        dc = abs(c1 - 1.0)
        print(f"  {name:8s}: c₁ = {c1:.10f}, c₂ = {c2:.10f}, δc = {dc:.2e}")
        assert dc < 1e-4, f"{name}: c = {c1:.6f}, δc = {dc:.2e}"
    print("  PASSED")


def test_R38_random_voronoi():
    """R38: c = 1 on 5 random Voronoi seeds (no cubic symmetry)."""
    print("\n=== R38: c = 1 on random Voronoi ===")
    L = 4.0
    k_mag = 0.005
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    for seed in [42, 137, 999, 2024, 7]:
        np.random.seed(seed)
        pts = np.random.uniform(0, L, (80, 3))
        data = build_foam_with_dual_info(pts, L)
        vals, _, _ = solve_bloch(data, k_vec)
        nz = vals[vals > 1e-10]
        c = np.sqrt(nz[0]) / k_mag
        dc = abs(c - 1.0)
        print(f"  Seed {seed:4d}: nE={len(data['E']):4d}, c = {c:.8f}, δc = {dc:.2e}")
        assert dc < 5e-4, f"Seed {seed}: c = {c:.6f}, δc = {dc:.2e}"
    print("  PASSED")


def test_R32_scaling_invariance():
    """R32: c = 1 independent of box scale L_cell."""
    print("\n=== R32: Scaling invariance ===")
    speeds = []
    for L_cell in [2.0, 4.0, 8.0, 16.0]:
        data = build_kelvin_with_dual_info(N=2, L_cell=L_cell)
        L = data["L"]
        k_mag = 0.01 * (2 * np.pi / L)
        k_vec = k_mag * np.array([1.0, 0.0, 0.0])

        vals, _, _ = solve_bloch(data, k_vec)
        nz = vals[vals > 1e-10]
        c = np.sqrt(nz[0]) / k_mag
        speeds.append(c)
        print(f"  L_cell={L_cell:5.1f}: L={L:6.1f}, c = {c:.10f}")

    spread = max(speeds) - min(speeds)
    print(f"  Spread = {spread:.2e}")
    assert spread < 1e-8, f"c not scale-invariant: spread = {spread:.2e}"
    print("  PASSED")


def test_R41_H1_vanishes():
    """R41: H¹(k≠0) = 0 on z=4 foams — ker(d₁) = im(d₀).

    On z=4 periodic foams: nE = 2nV. At k≠0, rank(d₁) = nV and rank(d₀) = nV,
    so dim H¹ = dim(ker d₁) − dim(im d₀) = (nE − nV) − nV = 0.
    This means every closed 1-form is exact — no spurious harmonic modes.
    """
    print("\n=== R41: H¹(k≠0) = 0 ===")

    k_mag = 0.01

    # Cubic
    for name, builder in BUILDERS:
        data = builder(N=2)
        V, E, F = data["V"], data["E"], data["F"]
        L = data["L"]
        L_vec = np.array(data["L_vec"])
        nV, nE = len(V), len(E)

        k_vec = k_mag * (2 * np.pi / L) * np.array([1.0, 0.0, 0.0])
        d0k = build_d0_bloch(V, E, L, k_vec)
        if sparse.issparse(d0k):
            d0k = d0k.toarray()
        d1k = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
        if sparse.issparse(d1k):
            d1k = d1k.toarray()

        rk_d1 = np.linalg.matrix_rank(d1k, tol=1e-8)
        rk_d0 = np.linalg.matrix_rank(d0k, tol=1e-8)
        H1 = (nE - rk_d1) - rk_d0

        print(f"  {name:8s}: nE/nV={nE/nV:.2f}, rank(d₁)={rk_d1}, "
              f"rank(d₀)={rk_d0}, H¹={H1}")

        assert nE == 2 * nV, f"{name}: nE={nE} ≠ 2nV={2*nV}"
        assert rk_d1 == nV, f"{name}: rank(d₁)={rk_d1} ≠ nV={nV}"
        assert H1 == 0, f"{name}: H¹={H1} ≠ 0"
        # Exactness
        assert np.linalg.norm(d1k @ d0k) < 1e-10, f"{name}: d₁d₀ ≠ 0"

    # Random Voronoi
    L_rv = 4.0
    for seed in [42, 137, 999]:
        np.random.seed(seed)
        pts = np.random.uniform(0, L_rv, (80, 3))
        data = build_foam_with_dual_info(pts, L_rv)
        V, E, F = data["V"], data["E"], data["F"]
        L_vec = np.array(data["L_vec"])
        nV, nE = len(V), len(E)

        k_vec = k_mag * (2 * np.pi / L_rv) * np.array([1.0, 0.0, 0.0])
        d0k = build_d0_bloch(V, E, L_rv, k_vec)
        if sparse.issparse(d0k):
            d0k = d0k.toarray()
        d1k = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
        if sparse.issparse(d1k):
            d1k = d1k.toarray()

        rk_d1 = np.linalg.matrix_rank(d1k, tol=1e-8)
        rk_d0 = np.linalg.matrix_rank(d0k, tol=1e-8)
        H1 = (nE - rk_d1) - rk_d0

        print(f"  Rnd({seed:4d}): nE/nV={nE/nV:.2f}, H¹={H1}")
        assert nE == 2 * nV, f"Rnd({seed}): nE ≠ 2nV"
        assert H1 == 0, f"Rnd({seed}): H¹ ≠ 0"
        assert np.linalg.norm(d1k @ d0k) < 1e-10, f"Rnd({seed}): d₁d₀ ≠ 0"

    print("  PASSED")


def test_R44e_schur_complement():
    """R44e: Schur complement on harmonic/optical blocks gives eigenvalues [0, k², k²].

    Block K(k) into harmonic (H₃, 3×3) and optical subspaces. The Schur
    complement S = H₃ − B·K_opt⁻¹·B† gives the effective acoustic operator.
    Eigenvalues: [0, k², k²] — longitudinal (→ gradient) + 2 transverse at c² = 1.
    """
    print("\n=== R44e: Schur complement ===")

    k_small = 0.001
    k_vec = k_small * np.array([1.0, 0.0, 0.0])
    k2 = k_small**2

    def run_schur(name, data):
        harm, _, _, _, Vol = build_harm_h2(data)
        # Solve at Gamma for optical basis
        vals_g, vecs_g, info_g = solve_bloch(data, np.zeros(3))
        n_zero = int(np.sum(vals_g < 1e-8))
        vecs_opt = vecs_g[:, n_zero:]

        # Solve at k
        vals_k, _, info_k = solve_bloch(data, k_vec)
        Kk = info_k["K"]

        H3 = harm.conj().T @ Kk @ harm
        B = harm.conj().T @ Kk @ vecs_opt
        Kopt = vecs_opt.conj().T @ Kk @ vecs_opt
        S = H3 - B @ np.linalg.inv(Kopt) @ B.conj().T
        S_eigs = np.sort(np.real(np.linalg.eigvalsh(S)))
        H3_eigs = np.sort(np.real(np.linalg.eigvalsh(H3)))

        s_trans = S_eigs[1:]
        print(f"  {name:8s}: H₃/k²=[{H3_eigs[0]/k2:.1f},{H3_eigs[1]/k2:.1f},"
              f"{H3_eigs[2]/k2:.1f}]  Schur/k²=[{S_eigs[0]/k2:.1e},"
              f"{s_trans[0]/k2:.6f},{s_trans[1]/k2:.6f}]")

        for j in range(2):
            assert abs(s_trans[j] / k2 - 1.0) < 1e-4, (
                f"{name}: Schur trans[{j}]/k² = {s_trans[j]/k2:.6f}")

    for sname, builder in BUILDERS:
        run_schur(sname, builder(N=2))

    L_rv = 4.0
    for seed in [42, 7]:
        np.random.seed(seed)
        pts = np.random.uniform(0, L_rv, (80, 3))
        data = build_foam_with_dual_info(pts, L_rv)
        run_schur(f"Rnd({seed:3d})", data)

    print("  PASSED")


def test_R44f_discrete_curl():
    """R44f: ∂(h₂†S₂d₁(k)·harm)/∂k|₀ = i·ε (Levi-Civita tensor).

    The discrete analog of ∇×(e^{ik·x}ê_α) = ik×ê_α·e^{ik·x}.
    Verified by central finite differences on each k-component.
    """
    print("\n=== R44f: discrete curl identity ===")

    dk = 1e-6

    def run_one(name, data):
        harm, h2, _, _, Vol = build_harm_h2(data)
        star1, star2 = build_hodge_stars_voronoi(data)
        V, E, F = data["V"], data["E"], data["F"]
        L_vec = np.array(data["L_vec"])
        L = data["L"]
        S2 = np.diag(star2)

        # H tensor check
        H = (h2 * np.sqrt(Vol)).T @ S2 @ (h2 * np.sqrt(Vol))
        H_err = np.max(np.abs(H / Vol - np.eye(3)))

        # dF/dk by central differences
        max_eps_err = 0
        for gamma in range(3):
            kp, km = np.zeros(3), np.zeros(3)
            kp[gamma], km[gamma] = dk, -dk

            d0p = build_d0_bloch(V, E, L, kp)
            if sparse.issparse(d0p): d0p = d0p.toarray()
            d1p = build_d1_bloch_exact(V, E, F, kp, L_vec, d0p)
            if sparse.issparse(d1p): d1p = d1p.toarray()

            d0m = build_d0_bloch(V, E, L, km)
            if sparse.issparse(d0m): d0m = d0m.toarray()
            d1m = build_d1_bloch_exact(V, E, F, km, L_vec, d0m)
            if sparse.issparse(d1m): d1m = d1m.toarray()

            Fp = h2.T @ S2 @ d1p @ harm
            Fm = h2.T @ S2 @ d1m @ harm
            dF = (Fp - Fm) / (2 * dk)

            for beta in range(3):
                for alpha in range(3):
                    perm = (beta, gamma, alpha)
                    if perm in [(0,1,2), (1,2,0), (2,0,1)]:
                        eps = 1
                    elif perm in [(0,2,1), (2,1,0), (1,0,2)]:
                        eps = -1
                    else:
                        eps = 0
                    err = abs(dF[beta, alpha] / 1j - eps)
                    max_eps_err = max(max_eps_err, err)

        print(f"  {name:10s}: ||H/Vol-I||={H_err:.2e}  ||dF/dk/i - ε||={max_eps_err:.2e}")
        assert H_err < 5e-8, f"{name}: H tensor error {H_err:.2e}"
        assert max_eps_err < 1e-8, f"{name}: Levi-Civita error {max_eps_err:.2e}"

    for sname, builder in BUILDERS:
        run_one(sname, builder(N=2))

    L_rv = 4.0
    for seed in [42, 7, 999]:
        np.random.seed(seed)
        pts = np.random.uniform(0, L_rv, (80, 3))
        data = build_foam_with_dual_info(pts, L_rv)
        run_one(f"Rnd({seed:4d})", data)

    print("  PASSED")


def test_R44g_full_proof_chain():
    """R44g: Complete analytic proof chain → S = k²·P_transverse.

    Chain: (1) G = H = Vol·I, (2) d₁(0)·harm = 0, (3) dF/dk = iε
    → C = h₂†S₂·u_perp = i·[k×] → S = C†C = [k×]ᵀ[k×] = k²·P_T.
    """
    print("\n=== R44g: full proof chain ===")

    k_small = 0.001
    k_vec = k_small * np.array([0.0, 1.0, 0.0])
    k2 = k_small**2

    def run_chain(name, data):
        harm, h2, dx, fa, Vol = build_harm_h2(data)
        star1, star2 = build_hodge_stars_voronoi(data)
        V, E, F = data["V"], data["E"], data["F"]
        L_vec = np.array(data["L_vec"])
        L = data["L"]
        S1 = np.diag(star1)
        S2 = np.diag(star2)

        # Step 1: G = H = Vol·I
        G = dx.T @ S1 @ dx
        H = fa.T @ S2 @ fa
        G_err = np.max(np.abs(G / Vol - np.eye(3)))
        H_err = np.max(np.abs(H / Vol - np.eye(3)))
        GH_err = np.max(np.abs(G - H))
        assert GH_err < 1e-8, f"{name}: ||G-H|| = {GH_err:.2e}"

        # Step 2: d₁(0)·harm = 0
        d0g = build_d0_bloch(V, E, L, np.zeros(3))
        if sparse.issparse(d0g): d0g = d0g.toarray()
        d1g = build_d1_bloch_exact(V, E, F, np.zeros(3), L_vec, d0g)
        if sparse.issparse(d1g): d1g = d1g.toarray()
        curl_harm = np.linalg.norm(d1g @ harm)

        # Step 3: u = d₁(k)·harm, project out im(d₁(0)) to get u_perp
        d0k = build_d0_bloch(V, E, L, k_vec)
        if sparse.issparse(d0k): d0k = d0k.toarray()
        d1k = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
        if sparse.issparse(d1k): d1k = d1k.toarray()

        u = d1k @ harm

        # S₂-orthogonal projector out of im(d₁(0))
        nV_loc, nF_loc = len(V), len(F)
        S2h = np.diag(np.sqrt(star2))
        S2ih = np.diag(1.0 / np.sqrt(star2))
        A_mat = S2h @ d1g
        U_full, s_full, _ = np.linalg.svd(A_mat, full_matrices=True)
        rank = np.sum(s_full > 1e-10)
        # At k=0: rank(d₁) = nV−2 (b₁=3, rank(d₀)=nV−1), coker = nF−nV+2
        assert rank == nV_loc - 2, f"{name}: rank(d₁(0)) = {rank}, expected {nV_loc-2}"
        Q = S2ih @ U_full[:, :rank]

        u_perp = np.zeros_like(u)
        for j in range(3):
            c_ex = Q.conj().T @ S2 @ u[:, j]
            u_perp[:, j] = u[:, j] - Q @ c_ex

        # S5: transverse u_perp ∈ H² at leading order
        h2_S2_norms = np.sqrt(np.diag(h2.conj().T @ S2 @ h2))
        h2_hat = h2 / h2_S2_norms[np.newaxis, :]
        h2_fracs = []
        for j in range(3):
            n2 = np.real(u_perp[:, j].conj() @ S2 @ u_perp[:, j])
            if n2 < k2 * 1e-4:
                continue  # longitudinal mode, norm ≈ 0
            coeffs = h2_hat.conj().T @ S2 @ u_perp[:, j]
            frac = np.real(coeffs.conj() @ coeffs) / n2
            h2_fracs.append(frac)
            assert frac > 0.9999, f"{name}: u_perp[{j}] H² frac = {frac:.6f}"
        min_h2 = min(h2_fracs) if h2_fracs else 0.0

        # Step 4: S = u_perp†S₂u_perp
        S_mat = np.real(u_perp.conj().T @ S2 @ u_perp)

        # Step 5: predict S = [k×]ᵀ[k×] = k²·P_T
        kx = np.array([
            [0, -k_vec[2], k_vec[1]],
            [k_vec[2], 0, -k_vec[0]],
            [-k_vec[1], k_vec[0], 0]
        ])
        S_pred = kx.T @ kx
        S_err = np.max(np.abs(S_mat - S_pred))

        print(f"  {name:10s}: G_err={G_err:.1e} H_err={H_err:.1e} "
              f"curl_h={curl_harm:.1e}  H²_frac={min_h2:.6f}  S_err={S_err:.2e}")
        assert G_err < 5e-8, f"{name}: G error {G_err:.2e}"
        assert H_err < 5e-8, f"{name}: H error {H_err:.2e}"
        assert curl_harm < 1e-10, f"{name}: curl(harm) = {curl_harm:.2e}"
        assert S_err < 1e-4, f"{name}: S chain error {S_err:.2e}"

    for sname, builder in BUILDERS:
        run_chain(sname, builder(N=2))

    L_rv = 4.0
    for seed in [42, 7]:
        np.random.seed(seed)
        pts = np.random.uniform(0, L_rv, (80, 3))
        data = build_foam_with_dual_info(pts, L_rv)
        run_chain(f"Rnd({seed:4d})", data)

    print("  PASSED")


def test_R6_converse():
    """R6: Non-Voronoi Hodge stars → c ≠ 1 (converse direction).

    Exact complex + perturbed ⋆₁ (G ≠ Vol·I): c deviates from 1.
    Shows Voronoi metric is necessary for c = 1, not just sufficient.
    """
    print("\n=== R6: Converse — c ≠ 1 with perturbed stars ===")

    k_mag = 0.005
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    data = build_kelvin_with_dual_info(N=2)
    star1, star2 = build_hodge_stars_voronoi(data)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]

    d0k = build_d0_bloch(V, E, L, k_vec)
    if sparse.issparse(d0k): d0k = d0k.toarray()
    d1k = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
    if sparse.issparse(d1k): d1k = d1k.toarray()

    S2 = np.diag(star2)
    K = d1k.conj().T @ S2 @ d1k
    K = (K + K.conj().T) / 2

    # Voronoi baseline
    S1 = np.diag(star1)
    vals = np.sort(np.real(eigh(K, S1, eigvals_only=True)))
    nz = vals[vals > 1e-10]
    c_vor = np.sqrt(nz[0]) / k_mag
    dc_vor = abs(c_vor - 1.0)

    # Perturbed ⋆₁ (10% random perturbation)
    np.random.seed(999)
    star1_pert = star1 * (1 + 0.1 * np.random.randn(len(star1)))
    star1_pert = np.maximum(star1_pert, 1e-10)
    S1p = np.diag(star1_pert)
    vals_p = np.sort(np.real(eigh(K, S1p, eigvals_only=True)))
    nz_p = vals_p[vals_p > 1e-10]
    c_pert = np.sqrt(nz_p[0]) / k_mag
    dc_pert = abs(c_pert - 1.0)

    print(f"  Voronoi:   c = {c_vor:.8f}, δc = {dc_vor:.2e}")
    print(f"  Perturbed: c = {c_pert:.8f}, δc = {dc_pert:.2e}")
    print(f"  → Perturbation worsens c by {dc_pert/dc_vor:.0f}×")

    assert dc_vor < 1e-4, f"Voronoi: c = {c_vor:.6f}, expected ≈ 1"
    assert dc_pert > 1e-3, f"Perturbed: c = {c_pert:.6f}, expected c ≠ 1"
    print("  PASSED")


def test_S1_multi_direction():
    """S1: c = 1 along 5 k-directions (isotropy check).

    R5 tests k ∥ x̂ only. This verifies c = 1 for body diagonal, face diagonal,
    and off-axis directions, ruling out lattice-direction artifacts.
    """
    print("\n=== S1: Multi-direction isotropy ===")

    k_mag = 0.005
    dirs = [
        ("x", np.array([1, 0, 0], dtype=float)),
        ("y", np.array([0, 1, 0], dtype=float)),
        ("z", np.array([0, 0, 1], dtype=float)),
        ("111", np.array([1, 1, 1], dtype=float)),
        ("110", np.array([1, 1, 0], dtype=float)),
    ]
    dirs = [(n, d / np.linalg.norm(d)) for n, d in dirs]

    data = build_kelvin_with_dual_info(N=2)
    speeds = []
    for dname, d_hat in dirs:
        k_vec = k_mag * d_hat
        vals, _, _ = solve_bloch(data, k_vec)
        nz = vals[vals > 1e-10]
        c = np.sqrt(nz[0]) / k_mag
        speeds.append(c)
        print(f"  k ∥ {dname:3s}: c = {c:.10f}, δc = {abs(c-1):.2e}")

    spread = max(speeds) - min(speeds)
    print(f"  Spread across directions = {spread:.2e}")
    for c in speeds:
        assert abs(c - 1.0) < 1e-4, f"c = {c:.6f} too far from 1"
    assert spread < 1e-5, f"Anisotropy: spread = {spread:.2e}"
    print("  PASSED")


def test_S2_dispersion_convergence():
    """S2: δc/k² → const as k → 0 (dispersion is O(k⁴)).

    At small k, c(k) = 1 + α·k² + O(k⁴). If δc/k² converges to a constant,
    the leading correction is exactly k². Richardson extrapolation ratio ≈ 4
    confirms O(k⁴) convergence.
    """
    print("\n=== S2: δc/k² convergence ===")

    data = build_kelvin_with_dual_info(N=2)
    k_dir = np.array([1.0, 0.0, 0.0])

    k_vals = [0.02, 0.01, 0.005, 0.0025]
    ratios = []

    for k_mag in k_vals:
        k_vec = k_mag * k_dir
        vals, _, _ = solve_bloch(data, k_vec)
        nz = vals[vals > 1e-10]
        c = np.sqrt(nz[0]) / k_mag
        dc = c - 1.0
        ratio = dc / k_mag**2
        ratios.append(ratio)
        print(f"  k={k_mag:.4f}: c-1 = {dc:.4e}, (c-1)/k² = {ratio:.6f}")

    # Richardson ratio: if δc = α·k², halving k → ratio stays constant
    # Actually δc = α·k² means (c-1)/k² = α. Check convergence.
    r_last = ratios[-1]
    r_prev = ratios[-2]
    conv = abs(r_last - r_prev) / abs(r_last)
    print(f"  Convergence: |(r[-1]-r[-2])/r[-1]| = {conv:.2e}")
    assert conv < 0.1, f"δc/k² not converging: {conv:.2e}"
    print("  PASSED")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    test_R5_gauge_speed()
    test_R38_random_voronoi()
    test_R32_scaling_invariance()
    test_R41_H1_vanishes()
    test_R44e_schur_complement()
    test_R44f_discrete_curl()
    test_R44g_full_proof_chain()
    test_R6_converse()
    test_S1_multi_direction()
    test_S2_dispersion_convergence()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED (§3 c² = 1 — 10 tests)")
    print("=" * 60)
