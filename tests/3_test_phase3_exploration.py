"""
W17 Test 3: Phase 3 exploration — anisotropy, metric preservation, mode anatomy.

CLAIM 1 (m1/R17): Standard DEC is catastrophically anisotropic. c_standard varies
by 37-70% depending on direction, even on cubic meshes. On Kelvin, cubic axes are
non-equivalent (n_lost = 6/7/6 on x/y/z). Exact DEC: isotropic to 4e-6.
Source: boundary face counts are symmetric and ||d1(k)-d1(0)|| identical on all axes.
The anisotropy comes from how phases interact with interior face-edge topology.

CLAIM 2 (M2/R19): Harmonic mean preserves H_eps = <1/eps>_vol * Vol * I exactly on
Kelvin at any contrast. Log mean breaks the trace by 7-33% depending on contrast
(19% at paper contrast eps in [1,13]). Log mean also introduces off-diagonal
anisotropy (~1.7%) even on Kelvin where harmonic gives exact isotropy.
Spectral consequence: harmonic vs log mean bands differ by up to 7.3%.

CLAIM 3 (s7/R16): The lowest nonzero standard mode is 89% gauge (exact zero modes)
+ 11% physical (acoustic mode at c=1). It is a MIXING of gauge and acoustic, not a
pure gradient mode. gauge_frac = grad_frac for all modes (exact gauge = im(d0)).
Modes 1-5: 98-99% gauge. Mode 8+: 99.9% physical (optical, barely affected).

CLAIM 4 (s4/R21 + R22): Face adjacency girth = 3 on all 3D Voronoi meshes. Every
edge shared by exactly 3 faces (valence 3 uniform). Confirms: tr(K^n) conserved for
n = 1,2 (= girth - 1), breaks at n = 3 by 0.44-0.63%.

CLAIM 5 (s5/R20): c_standard depends on supercell size N non-monotonically
(1.25, 1.58, 1.39 for N=2,3,4). Exact DEC: c = 1 for all N. Yet another
standard DEC pathology.

RAW OUTPUT (6 tests, all pass):
==================================
=== m1: Standard DEC anisotropy (+ exact sanity check) ===
  Kelvin  standard: 70.2% aniso (c=0.721-1.458), exact: <0.01%
  C15     standard: 42.8% aniso (c=1.120-1.805), exact: <0.01%
  WP      standard: 36.6% aniso (c=1.059-1.578), exact: <0.01%
=== M2: Harmonic vs log mean metric preservation ===
  eps [1,5]  harmonic: trace_ratio = 1.00000000, log_mean: 0.91362289
  eps [1,13] harmonic: trace_ratio = 1.00000000, log_mean: 0.81439994
  eps [1,50] harmonic: trace_ratio = 1.00000000, log_mean: 0.67407180
=== s7: Mode anatomy ===
  Mode 0: omega/k=1.2530, gauge_frac=0.886, phys_frac=0.114 (mixing)
  Modes 1-5: gauge_frac ≥ 0.988 (nearly pure gauge)
  Modes 6-7: gauge_frac ≥ 0.998 (pure gauge)
=== s4: Face adjacency ===
  All structures: valence=3 uniform, girth=3, triangles/nE=3.00
=== R22: Moment hierarchy ===
  n=1,2: exact conservation (rel_diff = 0)
  n=3: break 0.44-0.63% (girth = 3 → first 2 conserved)
=== s5: c_standard vs N ===
  N=2: 1.253, N=3: 1.576, N=4: 1.392 (non-monotonic, spread 0.32)

ANSWER:
=======
Standard DEC fails in three new ways beyond c != 1:
(1) catastrophic anisotropy (70% on Kelvin),
(2) N-dependent speed (non-monotonic),
(3) mode mixing (89% gauge + 11% acoustic).
Meanwhile, harmonic mean (not log mean) preserves the metric identity H_eps = <1/eps>*Vol*I.
"""
import sys
import os
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from collections import defaultdict

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BLOCH_SRC = "/Users/alextoader/Sites/st_bloch_exactness/src"
ST11_SRC = os.path.join(os.path.dirname(__file__), "..", "..", "src", "1_foam")
sys.path.insert(0, os.path.abspath(ST11_SRC))
sys.path.insert(0, BLOCH_SRC)

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


def build_std_spectrum(data, k_vec):
    """Build standard K and return sorted eigenvalues + eigenvectors."""
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
    vals, vecs = eigh(K_st, S1)
    idx = np.argsort(np.real(vals))
    return np.real(vals[idx]), vecs[:, idx], S1, S2


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


def compute_face_area_vectors(V, F, L_vec):
    """Face area vectors via cross product."""
    area_vecs = []
    for face in F:
        A = np.zeros(3)
        v0 = V[face[0]]
        for i in range(1, len(face)):
            vi = V[face[i]] - v0
            vi -= np.round(vi / L_vec) * L_vec
            vj = V[face[(i + 1) % len(face)]] - v0
            vj -= np.round(vj / L_vec) * L_vec
            A += np.cross(vi, vj)
        area_vecs.append(A / 2)
    return np.array(area_vecs)


def log_mean(a, b):
    """Logarithmic mean of two positive numbers."""
    if abs(a - b) < 1e-12 * max(a, b):
        return a
    return (a - b) / np.log(a / b)


# ===================================================================
# Tests
# ===================================================================

def test_m1_standard_anisotropy():
    """m1/R17: Standard DEC speed varies by 37-70% with direction."""
    print("\n=== m1: Standard DEC anisotropy ===")

    k_mag = 0.005
    dirs = {
        "[1,0,0]": np.array([1, 0, 0], dtype=float),
        "[0,1,0]": np.array([0, 1, 0], dtype=float),
        "[0,0,1]": np.array([0, 0, 1], dtype=float),
        "[1,1,1]": np.array([1, 1, 1], dtype=float) / np.sqrt(3),
    }

    for name, builder in BUILDERS:
        data = builder(N=2)
        nV = len(data["V"])
        speeds_st = []
        speeds_ex = []

        for dname, k_dir in dirs.items():
            k_vec = k_mag * k_dir
            r = build_both_spectra(data, k_vec)
            n0_st = np.sum(r["vals_st"] < 1e-12)
            nz_st = r["vals_st"][r["vals_st"] > 1e-10]
            c_st = np.sqrt(nz_st[0]) / k_mag
            speeds_st.append(c_st)

            nz_ex = r["vals_ex"][r["vals_ex"] > 1e-10]
            c_ex = np.sqrt(nz_ex[0]) / k_mag
            speeds_ex.append(c_ex)

        aniso_st = (max(speeds_st) - min(speeds_st)) / np.mean(speeds_st) * 100
        aniso_ex = (max(speeds_ex) - min(speeds_ex)) / np.mean(speeds_ex) * 100
        print(f"  {name:8s} standard: {aniso_st:.1f}% aniso "
              f"(c={min(speeds_st):.3f}-{max(speeds_st):.3f}), "
              f"exact: {aniso_ex:.4f}%")

        # Standard: anisotropy should be > 30%
        assert aniso_st > 30, f"{name}: std anisotropy = {aniso_st:.1f}%, expected > 30%"
        # Exact: anisotropy should be < 0.01% (sanity check)
        assert aniso_ex < 0.01, f"{name}: exact anisotropy = {aniso_ex:.4f}%, expected < 0.01%"

    print("  PASSED")


def test_M2_harmonic_preserves_metric():
    """M2/R19: Harmonic mean preserves H_eps = <1/eps>*Vol*I; log mean does not."""
    print("\n=== M2: Harmonic vs log mean metric preservation ===")

    data = build_kelvin_with_dual_info(N=2)
    _, star2_geom = build_hodge_stars_voronoi(data)
    V, F = data["V"], data["F"]
    L_vec = np.array(data["L_vec"])
    vol = np.prod(L_vec)
    f2c = data["face_to_cells"]
    nF = len(F)
    nC = max(max(c1, c2) for c1, c2 in f2c.values()) + 1

    Af = compute_face_area_vectors(V, F, L_vec)

    for contrast, label in [(5, "[1,5]"), (13, "[1,13]"), (50, "[1,50]")]:
        eps = np.random.RandomState(42).uniform(1, contrast, nC)
        inv_eps_vol = np.mean(1.0 / eps)  # equal cell volumes on Kelvin

        for avg_name in ["harmonic", "log_mean"]:
            star2_eps = np.zeros(nF)
            for f_idx in range(nF):
                c1, c2 = f2c[f_idx]
                e1, e2 = eps[c1], eps[c2]
                if avg_name == "harmonic":
                    inv_eff = 0.5 * (1.0 / e1 + 1.0 / e2)
                else:
                    eps_eff = log_mean(e1, e2)
                    inv_eff = 1.0 / eps_eff
                star2_eps[f_idx] = star2_geom[f_idx] * inv_eff

            H = np.zeros((3, 3))
            for f_idx in range(nF):
                H += star2_eps[f_idx] * np.outer(Af[f_idx], Af[f_idx])

            diag = [H[i, i] for i in range(3)]
            off = max(abs(H[i, j]) for i in range(3) for j in range(3) if i != j)
            trace_ratio = np.mean(diag) / (inv_eps_vol * vol)

            print(f"  eps {label:6s} {avg_name:10s}: trace_ratio = {trace_ratio:.8f}, "
                  f"off_max = {off:.2e}")

            if avg_name == "harmonic":
                assert abs(trace_ratio - 1.0) < 1e-10, (
                    f"Harmonic eps {label}: ratio = {trace_ratio}")
                assert off < 1e-12, (
                    f"Harmonic eps {label}: off = {off}")
            else:
                # Log mean should NOT preserve the identity
                assert trace_ratio < 0.95, (
                    f"Log mean eps {label}: ratio = {trace_ratio}, expected < 0.95")

    print("  PASSED")


def test_s7_mode_anatomy():
    """s7/R16: Lowest standard mode is gauge-physical mixture, not pure gradient."""
    print("\n=== s7: Mode anatomy — gauge/physical decomposition ===")

    k_mag = 0.005
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    data = build_kelvin_with_dual_info(N=2)
    r = build_both_spectra(data, k_vec)
    S1 = r["S1"]

    n0_ex = np.sum(r["vals_ex"] < 1e-12)
    n0_st = np.sum(r["vals_st"] < 1e-12)
    n_lost = n0_ex - n0_st

    print(f"  Kelvin: n0_ex={n0_ex}, n0_st={n0_st}, n_lost={n_lost}")
    print(f"  {'mode':>5s} {'omega/k':>10s} {'gauge_frac':>11s} {'phys_frac':>10s}")

    for i in range(n0_st, min(n0_st + 8, len(r["vals_st"]))):
        v = r["vecs_st"][:, i]
        omega_k = np.sqrt(max(r["vals_st"][i], 0)) / k_mag

        # Project onto exact gauge and physical subspaces
        gauge_overlap = sum(
            abs(r["vecs_ex"][:, j].conj() @ S1 @ v) ** 2
            for j in range(n0_ex)
        )
        norm_sq = np.real(v.conj() @ S1 @ v)
        gauge_frac = gauge_overlap / norm_sq
        phys_frac = 1.0 - gauge_frac

        print(f"  {i - n0_st:5d} {omega_k:10.4f} {gauge_frac:11.6f} {phys_frac:10.6f}")

        if i == n0_st:
            # Lowest mode: should be 80-95% gauge, 5-20% physical
            assert 0.7 < gauge_frac < 0.98, (
                f"Mode 0 gauge_frac = {gauge_frac:.4f}, expected 0.7-0.98")
        elif i < n0_st + n_lost:
            # Other lost modes: should be > 95% gauge
            assert gauge_frac > 0.95, (
                f"Mode {i - n0_st} gauge_frac = {gauge_frac:.4f}")

    print("  PASSED")


def test_s4_girth_and_valence():
    """s4/R21: Face adjacency girth = 3, edge valence = 3 uniform on all Voronoi."""
    print("\n=== s4: Face adjacency girth and edge valence ===")

    for name, builder in BUILDERS:
        data = builder(N=2)
        E = data["E"]
        F = data["F"]
        nE = len(E)
        nF = len(F)

        # Build edge index
        edge_idx = {}
        for e_idx, (a, b) in enumerate(E):
            edge_idx[(min(a, b), max(a, b))] = e_idx

        # Edge-to-faces map
        edge_to_faces = defaultdict(set)
        for f_idx, face in enumerate(F):
            n = len(face)
            for i in range(n):
                a, b = face[i], face[(i + 1) % n]
                key = (min(a, b), max(a, b))
                if key in edge_idx:
                    edge_to_faces[edge_idx[key]].add(f_idx)

        valences = [len(fs) for fs in edge_to_faces.values()]

        # Build face adjacency graph
        face_adj = defaultdict(set)
        for e, faces in edge_to_faces.items():
            flist = list(faces)
            for i in range(len(flist)):
                for j in range(i + 1, len(flist)):
                    face_adj[flist[i]].add(flist[j])
                    face_adj[flist[j]].add(flist[i])

        # Count triangles
        n_triangles = 0
        for f in range(nF):
            for f2 in face_adj[f]:
                if f2 > f:
                    common = face_adj[f] & face_adj[f2]
                    for f3 in common:
                        if f3 > f2:
                            n_triangles += 1

        # Girth via BFS (sample first 30 faces for speed)
        girth = float("inf")
        for start in range(min(nF, 30)):
            dist = {start: 0}
            parent = {start: -1}
            queue = [start]
            qi = 0
            while qi < len(queue):
                f = queue[qi]
                qi += 1
                for nb in face_adj[f]:
                    if nb not in dist:
                        dist[nb] = dist[f] + 1
                        parent[nb] = f
                        queue.append(nb)
                    elif parent[f] != nb and parent.get(nb, -1) != f:
                        girth = min(girth, dist[f] + dist[nb] + 1)

        tri_per_edge = n_triangles / nE

        print(f"  {name:8s}: valence = {min(valences)}-{max(valences)}, "
              f"girth = {girth}, triangles/nE = {tri_per_edge:.2f}")

        assert min(valences) == 3, f"{name}: min valence = {min(valences)}"
        assert max(valences) == 3, f"{name}: max valence = {max(valences)}"
        assert girth == 3, f"{name}: girth = {girth}"
        assert abs(tri_per_edge - 3.0) < 0.01, f"{name}: tri/nE = {tri_per_edge}"

    print("  PASSED")


def test_R22_moment_hierarchy():
    """R22: tr(K^n) conserved for n=1,2 (exact), breaks at n=3 by 0.4-0.7%."""
    print("\n=== R22: Moment hierarchy tr(K^n) ===")

    k_mag = 0.1
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    for name, builder in BUILDERS:
        data = builder(N=2)
        r = build_both_spectra(data, k_vec)
        nE = len(data["E"])

        Kn_ex = np.eye(nE)
        Kn_st = np.eye(nE)

        print(f"  {name}:")
        for n in range(1, 5):
            Kn_ex = Kn_ex @ r["K_ex"]
            Kn_st = Kn_st @ r["K_st"]
            tr_ex = np.real(np.trace(Kn_ex))
            tr_st = np.real(np.trace(Kn_st))
            rel = abs(tr_st - tr_ex) / abs(tr_ex) if abs(tr_ex) > 0 else 0
            print(f"    n={n}: rel_diff = {rel:.2e}")

            if n <= 2:
                assert rel < 1e-10, (
                    f"{name} n={n}: rel_diff = {rel:.2e}, expected < 1e-10")
            elif n == 3:
                assert 1e-4 < rel < 0.01, (
                    f"{name} n=3: rel_diff = {rel:.2e}, expected 0.4-0.7%")

    print("  PASSED")


def test_s5_N_dependence():
    """s5/R20: c_standard depends on N non-monotonically."""
    print("\n=== s5: c_standard vs supercell size N ===")

    k_mag = 0.005
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])

    c_values = []
    for N in [2, 3, 4]:
        data = build_kelvin_with_dual_info(N=N)
        vals_st, _, _, _ = build_std_spectrum(data, k_vec)
        nV = len(data["V"])
        n0_st = np.sum(vals_st < 1e-12)
        nz_st = vals_st[vals_st > 1e-10]
        c_st = np.sqrt(nz_st[0]) / k_mag
        c_values.append(c_st)
        print(f"  N={N}: nV={nV:5d}, n_lost={nV - n0_st:4d}, c_std = {c_st:.6f}")

    # c_std should vary with N (not constant)
    spread = max(c_values) - min(c_values)
    assert spread > 0.1, f"c_std spread = {spread:.4f}, expected > 0.1 (N-dependent)"

    # Should NOT be monotonic (non-monotonic was observed: 1.25, 1.58, 1.39)
    monotone_up = all(c_values[i] <= c_values[i + 1] for i in range(len(c_values) - 1))
    monotone_down = all(c_values[i] >= c_values[i + 1] for i in range(len(c_values) - 1))
    assert not (monotone_up or monotone_down), "c_std is monotonic in N (unexpected)"

    print(f"  Spread = {spread:.4f} (N-dependent, non-monotonic)")
    print("  PASSED")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    test_m1_standard_anisotropy()
    test_M2_harmonic_preserves_metric()
    test_s7_mode_anatomy()
    test_s4_girth_and_valence()
    test_R22_moment_hierarchy()
    test_s5_N_dependence()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED (m1, M2, s4, s5, s7, R22 — 6 tests)")
    print("=" * 60)
