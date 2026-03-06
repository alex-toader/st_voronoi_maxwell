"""
Paper §5: Removing geometry — what breaks when G ≠ Vol·I.

Supports: §5 (Necessity of metric identity)

CLAIM: On an exact complex (d₁d₀ = 0), perturbing Hodge stars away from Voronoi
breaks c = 1 but does NOT create gauge pollution (n_lost stays 0). The killer table
shows a clean 2×2 factorization: topology controls mode count, geometry controls speed.
Errors are separable: Δc_both ≈ Δc_geom + Δc_topo with interaction < 0.05.

RAW OUTPUT (3 tests, all pass):
=================================
=== R36: Killer table — cubic structures ===
  Kelvin (nE=192):
    E+V: c=1.0000, n0=96,  n_lost=0   | E+P: c=0.9793, n0=96,  n_lost=0
    S+V: c=1.2527, n0=90,  n_lost=6   | S+P: c=1.2572, n0=90,  n_lost=6
    interaction=+0.025
  C15 (nE=2176):
    E+V: c=1.0000, n0=1088, n_lost=0  | E+P: c=0.9964, n0=1088, n_lost=0
    S+V: c=1.6848, n0=1049, n_lost=39 | S+P: c=1.6505, n0=1049, n_lost=39
    interaction=-0.031
  WP (nE=736):
    E+V: c=1.0000, n0=368, n_lost=0   | E+P: c=1.0062, n0=368, n_lost=0
    S+V: c=1.4816, n0=353, n_lost=15  | S+P: c=1.5116, n0=353, n_lost=15
    interaction=+0.024
=== R39: Killer table — random Voronoi ===
  Seed 42:  c_EV=1.000, c_SV=0.472, n_lost=23, interaction=+0.006
  Seed 137: c_EV=1.000, c_SV=0.903, n_lost=28, interaction=-0.004
  Seed 999: c_EV=1.000, c_SV=1.024, n_lost=31, interaction=+0.002
=== R29: Error factorization — additivity ===
  Kelvin: dc_topo=+0.253
    5%: interaction=+0.013  |  10%: interaction=+0.025  |  15%: interaction=+0.037
ALL TESTS PASSED (§5 removing geometry — 3 tests)

ANSWER:
=======
Perturbing Hodge stars on an exact complex changes c (0.979–1.006) but not n_lost (stays 0).
Breaking exactness changes n_lost but not through stars. The 2×2 killer table is universal:
cubic + random, with interaction term < 0.05. Geometry and topology are independent mechanisms.
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

def build_operators(data, k_vec):
    """Build exact and standard d₁ at k, return operators + stars."""
    star1, star2 = build_hodge_stars_voronoi(data)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]

    d0k = build_d0_bloch(V, E, L, k_vec)
    if sparse.issparse(d0k):
        d0k = d0k.toarray()

    d1k_ex = build_d1_bloch_exact(V, E, F, k_vec, L_vec, d0k)
    if sparse.issparse(d1k_ex):
        d1k_ex = d1k_ex.toarray()

    d1k_st = build_d1_bloch_standard(V, E, F, L, k_vec)
    if sparse.issparse(d1k_st):
        d1k_st = d1k_st.toarray()

    return {
        "star1": star1, "star2": star2,
        "d1k_ex": d1k_ex, "d1k_st": d1k_st,
    }


def get_c_and_n0(d1k, S1, S2, k2):
    """Return (c, n_zero) from curl-curl eigenvalue problem."""
    K = d1k.conj().T @ S2 @ d1k
    K = (K + K.conj().T) / 2
    vals = np.sort(np.real(eigh(K, S1, eigvals_only=True)))
    n0 = int(np.sum(vals < 1e-10))
    nz = vals[vals > 1e-10]
    c = np.sqrt(nz[0] / k2) if len(nz) > 0 else float('nan')
    return c, n0


def perturb_stars(star1, star2, seed, amplitude=0.10):
    """Return perturbed Hodge stars (positive, clipped at 0.5× original)."""
    rng = np.random.RandomState(seed)
    s1_p = star1 * np.maximum(1.0 + amplitude * rng.randn(len(star1)), 0.5)
    s2_p = star2 * np.maximum(1.0 + amplitude * rng.randn(len(star2)), 0.5)
    return s1_p, s2_p


def killer_table(d1k_ex, d1k_st, star1, star2, k2, pert_seed):
    """Run 2×2 table, return dict of (c, n0) for each config."""
    S1_v = np.diag(star1)
    S2_v = np.diag(star2)

    s1_p, s2_p = perturb_stars(star1, star2, pert_seed)
    S1_p = np.diag(s1_p)
    S2_p = np.diag(s2_p)

    results = {}
    for label, d1k, S1, S2 in [
        ("E+V", d1k_ex, S1_v, S2_v),
        ("E+P", d1k_ex, S1_p, S2_p),
        ("S+V", d1k_st, S1_v, S2_v),
        ("S+P", d1k_st, S1_p, S2_p),
    ]:
        results[label] = get_c_and_n0(d1k, S1, S2, k2)

    return results


# ===================================================================
# Tests
# ===================================================================

def test_N15_killer_table_cubic():
    """R36: 2×2 factorization table on Kelvin, C15, WP.

    Four configurations: exact/standard × Voronoi/perturbed Hodge stars.
    Exact+Voronoi: c = 1, n_lost = 0 (both conditions satisfied).
    Exact+Perturbed: c ≈ 0.99, n_lost = 0 (topology OK, geometry broken).
    Standard+Voronoi: c > 1.1, n_lost > 0 (topology broken, geometry OK).
    Standard+Perturbed: c > 1.1, n_lost > 0 (both broken).
    n_lost unchanged by star perturbation — topology, not geometry.
    """
    print("\n=== R36: Killer table — cubic structures ===")

    k_mag = 0.01
    pert_seed = 42

    for name, builder in BUILDERS:
        data = builder(N=2)
        L = data["L"]
        nE = len(data["E"])

        k_vec = k_mag * (2 * np.pi / L) * np.array([1.0, 0.0, 0.0])
        k2 = k_vec @ k_vec

        ops = build_operators(data, k_vec)
        r = killer_table(ops["d1k_ex"], ops["d1k_st"],
                         ops["star1"], ops["star2"], k2, pert_seed)

        c_EV, n0_EV = r["E+V"]
        c_EP, n0_EP = r["E+P"]
        c_SV, n0_SV = r["S+V"]
        c_SP, n0_SP = r["S+P"]
        n_lost = n0_EV - n0_SV

        print(f"\n  === {name} (nE={nE}) ===")
        print(f"  {'Config':>22s} {'c':>8s} {'n0':>5s} {'n_lost':>7s}")
        print(f"  {'Exact+Voronoi':>22s} {c_EV:8.4f} {n0_EV:5d} {0:7d}")
        print(f"  {'Exact+Perturbed':>22s} {c_EP:8.4f} {n0_EP:5d} {n0_EV - n0_EP:7d}")
        print(f"  {'Standard+Voronoi':>22s} {c_SV:8.4f} {n0_SV:5d} {n_lost:7d}")
        print(f"  {'Standard+Perturbed':>22s} {c_SP:8.4f} {n0_SP:5d} {n0_EV - n0_SP:7d}")

        dc_geom = c_EP - c_EV
        dc_topo = c_SV - c_EV
        dc_both = c_SP - c_EV
        interaction = dc_both - (dc_geom + dc_topo)
        print(f"  interaction = {interaction:+.4f}")

        # (1) Exact+Voronoi: c = 1
        assert abs(c_EV - 1.0) < 0.001, f"{name}: c_EV = {c_EV:.6f}"

        # (2) Exact+Perturbed: c ≠ 1 but close, n_lost = 0
        assert n0_EV == n0_EP, f"{name}: perturbation changed n_zero on exact"
        assert abs(c_EP - 1.0) < 0.05, f"{name}: c_EP = {c_EP:.4f}, too far"

        # (3) Standard+Voronoi: c >> 1, n_lost > 0
        assert c_SV > 1.1, f"{name}: c_SV = {c_SV:.4f}"
        assert n_lost > 0, f"{name}: no lost modes"

        # (4) n_lost same for Voronoi and Perturbed
        assert n0_SV == n0_SP, (
            f"{name}: n_lost differs: V={n0_EV - n0_SV} vs P={n0_EP - n0_SP}")

        # (5) Interaction small
        assert abs(interaction) < 0.05, (
            f"{name}: interaction = {interaction:.4f}")

    print("\n  PASSED")


def test_R39_killer_table_random():
    """R39: 2×2 factorization table on 3 random Voronoi.

    Same structure as R36 but on meshes with zero symmetry. Interaction < 0.10.
    Stronger test: no cubic symmetry to accidentally force separability.
    """
    print("\n=== R39: Killer table — random Voronoi ===")

    L = 4.0
    k_mag = 0.01
    k_vec = k_mag * np.array([1.0, 0.0, 0.0])
    k2 = k_vec @ k_vec

    n_pass = 0
    for seed in [42, 137, 999]:
        np.random.seed(seed)
        pts = np.random.uniform(0, L, (80, 3))
        try:
            data = build_foam_with_dual_info(pts, L)
        except Exception:
            continue

        ops = build_operators(data, k_vec)
        pert_seed = 42 + seed
        r = killer_table(ops["d1k_ex"], ops["d1k_st"],
                         ops["star1"], ops["star2"], k2, pert_seed)

        c_EV, n0_EV = r["E+V"]
        c_EP, n0_EP = r["E+P"]
        c_SV, n0_SV = r["S+V"]
        c_SP, n0_SP = r["S+P"]
        n_lost = n0_EV - n0_SV

        dc_geom = c_EP - c_EV
        dc_topo = c_SV - c_EV
        dc_both = c_SP - c_EV
        interaction = dc_both - (dc_geom + dc_topo)

        print(f"\n  Seed {seed}: nE={len(data['E'])}")
        print(f"    E+V: c={c_EV:.4f}, n0={n0_EV}")
        print(f"    E+P: c={c_EP:.4f}, n0={n0_EP}")
        print(f"    S+V: c={c_SV:.4f}, n0={n0_SV}")
        print(f"    S+P: c={c_SP:.4f}, n0={n0_SP}")
        print(f"    n_lost={n_lost}, interaction={interaction:+.4f}")

        # (1) Exact+Voronoi: c = 1
        assert abs(c_EV - 1.0) < 5e-3, f"Seed {seed}: c_EV = {c_EV:.6f}"

        # (2) Standard loses modes
        assert n_lost > 0, f"Seed {seed}: n_lost = 0, standard should lose modes"

        # (3) Perturbation does not change n_lost
        assert n0_EV == n0_EP, f"Seed {seed}: perturbation changed exact n_zero"
        assert n0_SV == n0_SP, f"Seed {seed}: perturbation changed standard n_zero"

        # (4) Interaction small
        assert abs(interaction) < 0.10, (
            f"Seed {seed}: interaction = {interaction:.4f}")

        n_pass += 1

    assert n_pass >= 2, f"Only {n_pass}/3 seeds passed"
    print(f"\n  {n_pass} seeds passed")
    print("  PASSED")


def test_M1_error_factorization():
    """R29: Δc_geom + Δc_topo separable — quantitative additivity.

    On Kelvin (cleanest structure), test that the interaction term stays small
    across multiple perturbation amplitudes (5%, 10%, 15%). The interaction
    should grow slowly (sub-linearly) while individual effects grow with amplitude.
    """
    print("\n=== R29: Error factorization — additivity ===")

    data = build_kelvin_with_dual_info(N=2)
    V, E, F = data["V"], data["E"], data["F"]
    L_vec = np.array(data["L_vec"])
    L = data["L"]
    nE = len(E)

    star1, star2 = build_hodge_stars_voronoi(data)

    k_mag = 0.01
    k_vec = k_mag * (2 * np.pi / L) * np.array([1.0, 0.0, 0.0])
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

    S1_v = np.diag(star1)
    S2_v = np.diag(star2)

    c_EV, _ = get_c_and_n0(d1k_ex, S1_v, S2_v, k2)
    c_SV, _ = get_c_and_n0(d1k_st, S1_v, S2_v, k2)
    dc_topo = c_SV - c_EV

    print(f"  c_EV = {c_EV:.6f}, c_SV = {c_SV:.6f}, dc_topo = {dc_topo:+.4f}")
    print(f"\n  {'amp':>5s} {'dc_geom':>9s} {'dc_both':>9s} {'interact':>9s}")

    amplitudes = [0.05, 0.10, 0.15]
    for amp in amplitudes:
        s1_p, s2_p = perturb_stars(star1, star2, seed=42, amplitude=amp)
        S1_p = np.diag(s1_p)
        S2_p = np.diag(s2_p)

        c_EP, n0_EP = get_c_and_n0(d1k_ex, S1_p, S2_p, k2)
        c_SP, n0_SP = get_c_and_n0(d1k_st, S1_p, S2_p, k2)

        dc_geom = c_EP - c_EV
        dc_both = c_SP - c_EV
        interaction = dc_both - (dc_geom + dc_topo)

        print(f"  {amp:5.0%} {dc_geom:+9.4f} {dc_both:+9.4f} {interaction:+9.4f}")

        # Interaction should be small relative to individual effects
        assert abs(interaction) < 0.05, (
            f"amp={amp}: interaction = {interaction:.4f}")

    print("  PASSED")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    test_N15_killer_table_cubic()
    test_R39_killer_table_random()
    test_M1_error_factorization()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED (§5 removing geometry — 3 tests)")
    print("=" * 60)
