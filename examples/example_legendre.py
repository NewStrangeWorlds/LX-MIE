"""
Compute Legendre series coefficients of the phase function for Cr particles.

Demonstrates both the single-wavelength and batch APIs.

Run from the lx-mie root directory:
    python examples/example_legendre.py
"""

import math
import lxmie


def read_refractive_index(path):
    """Read a refractive index file (wavelength in microns, n, k columns)."""
    wavelengths, n_real, n_imag = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            wavelengths.append(float(parts[0]))
            n_real.append(float(parts[1]))
            n_imag.append(float(parts[2]))
    return wavelengths, n_real, n_imag


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

refindex_file = "compilation/Cr.dat"
particle_radius = 10.0  # microns
nb_legendre_terms = 50

wavelengths, n_real, n_imag = read_refractive_index(refindex_file)
size_parameters = [2.0 * math.pi * particle_radius / wl for wl in wavelengths]

# ---------------------------------------------------------------------------
# Single-wavelength example
# ---------------------------------------------------------------------------

wl_idx = 0
result_single, moments_single = lxmie.mie_legendre(
    n_real[wl_idx], n_imag[wl_idx], size_parameters[wl_idx],
    nb_moments=nb_legendre_terms)

print(f"Single wavelength: {wavelengths[wl_idx]:.4e} um  "
      f"(x = {size_parameters[wl_idx]:.4f})")
print(f"  Q_ext = {result_single.q_ext:.6f}  "
      f"Q_sca = {result_single.q_sca:.6f}  "
      f"Q_abs = {result_single.q_abs:.6f}")
print(f"  First {nb_legendre_terms} Legendre coefficients:")
for i, c in enumerate(moments_single):
    print(f"    chi_{i:>3d} = {c:+.6e}")

# ---------------------------------------------------------------------------
# Batch example — all wavelengths at once
# ---------------------------------------------------------------------------

print(f"\nBatch: computing {len(wavelengths)} wavelengths in parallel...")

results_batch, all_moments = lxmie.mie_legendre(
    n_real, n_imag, size_parameters,
    nb_moments=nb_legendre_terms)

print(f"Done. Showing chi_1 (= asymmetry parameter) for each wavelength:\n")
print(f"{'wavelength (um)':>18}{'chi_1':>14}")
print("-" * 32)
for wl, moments in zip(wavelengths, all_moments):
    print(f"{wl:18.6e}{moments[1]:14.6f}")
