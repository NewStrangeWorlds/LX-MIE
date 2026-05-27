"""
Compute the phase function for Cr particles at a set of scattering angles.

The angular_grid passed to lxmie.mie_amplitudes() contains the COSINES of the
scattering angles (not angles in degrees or radians).

Demonstrates both the single-wavelength and batch APIs, and shows how to
convert the scattering amplitudes into a normalised phase function using
lxmie.phase_function().

Run from the lx-mie root directory:
    python examples/example_phase_function.py
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

wavelengths, n_real, n_imag = read_refractive_index(refindex_file)
size_parameters = [2.0 * math.pi * particle_radius / wl for wl in wavelengths]

# Angular grid: 37 angles from 0° to 180°, stored as cosines.
# cos(0°) = 1  →  cos(180°) = -1, so the grid runs from 1 to -1.
nb_angles = 37
angles_deg = [i * 180.0 / (nb_angles - 1) for i in range(nb_angles)]
cos_angles = [math.cos(math.radians(a)) for a in angles_deg]

# ---------------------------------------------------------------------------
# Single-wavelength example
# ---------------------------------------------------------------------------

wl_idx = 0
result_single, s1, s2 = lxmie.mie_amplitudes(
    n_real[wl_idx], n_imag[wl_idx], size_parameters[wl_idx],
    angles=cos_angles)

phase_function_single = lxmie.phase_function(
    s1, s2,
    size_parameter=size_parameters[wl_idx],
    q_sca=result_single.q_sca)

print(f"Single wavelength: {wavelengths[wl_idx]:.4e} um  "
      f"(x = {size_parameters[wl_idx]:.4f})")
print(f"  Q_ext = {result_single.q_ext:.6f}  "
      f"Q_sca = {result_single.q_sca:.6f}  "
      f"Q_abs = {result_single.q_abs:.6f}")
print(f"\n  {'angle (deg)':>12}{'cos(angle)':>14}{'phase function':>16}")
print("  " + "-" * 42)
for deg, cos_a, pf in zip(angles_deg, cos_angles, phase_function_single):
    print(f"  {deg:12.1f}{cos_a:14.6f}{pf:16.6e}")

# ---------------------------------------------------------------------------
# Batch example — all wavelengths at once
# ---------------------------------------------------------------------------

print(f"\nBatch: computing phase function for {len(wavelengths)} wavelengths in parallel...")

results_batch, all_s1, all_s2 = lxmie.mie_amplitudes(
    n_real, n_imag, size_parameters,
    angles=cos_angles)

all_phase_functions = [
    lxmie.phase_function(s1, s2, size_parameter=x, q_sca=r.q_sca)
    for s1, s2, x, r in zip(all_s1, all_s2, size_parameters, results_batch)
]

# Print the forward-scattering (0°) value for each wavelength
print(f"Done. Forward-scattering (0°) phase function per wavelength:\n")
print(f"{'wavelength (um)':>18}{'P(0°)':>16}")
print("-" * 34)
for wl, pf in zip(wavelengths, all_phase_functions):
    print(f"{wl:18.6e}{pf[0]:16.6e}")
