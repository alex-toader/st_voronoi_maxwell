Dear Editor,

Please find enclosed our manuscript entitled "Geometry and topology of discrete Maxwell wave speed", which we submit for consideration for publication in the Journal of Computational Physics as a companion paper to JCOMP-D-26-00537.

The paper proves that on any periodic Voronoi tessellation equipped with an exact discrete cochain complex, the long-wavelength electromagnetic wave speed equals unity in natural units. The result is a geometric identity at finite resolution, not a convergence statement. The proof identifies two sufficient conditions: (A) the Voronoi metric tensors satisfy G = H = Vol I (proved via the divergence theorem), and (B) the Bloch-periodic cochain complex preserves exactness (provided by the companion paper JCOMP-D-26-00537). These combine through a Schur complement reduction and a discrete curl identity to yield the exact leading dispersion coefficient.

The necessity of both conditions is established: removing exactness (condition B) produces O(1) error in the wave speed with severe anisotropy and no convergence under refinement; removing the Voronoi metric (condition A) produces direction-dependent speed errors. The error contributions are approximately separable.

Numerical verification spans three ordered foams (Kelvin, Weaire-Phelan, C15) and five random Voronoi tessellations, supported by a structured test suite (7 files, 42 tests) publicly available at https://github.com/alex-toader/st_voronoi_maxwell.

Together with JCOMP-D-26-00537, which constructs the exact cochain complex, the present paper provides a complete discrete Maxwell theory on periodic Voronoi complexes: exact topology and exact long-wavelength wave speed.

This work is original, has not been published previously, and is not under consideration for publication elsewhere. The author declares no conflict of interest.

Thank you for your consideration.

Alex Toader
