"""
Compute Mie scattering efficiencies and cross sections for Cr particles
across all wavelengths in a refractive index file.

Uses the batch API so all wavelengths are computed in parallel via OpenMP.

Run from the lx-mie root directory:
    python examples/example_efficiencies.py
"""

import math
import lxmie


def read_refractive_index(path):
    """Read a refractive index file (wavelength in microns, n, k columns).

    Lines starting with '#' are treated as comments/headers.
    Returns (wavelengths, n_real, n_imag) as lists of floats.
    """
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

# Size parameter x = 2 * pi * r / lambda for each wavelength
size_parameters = [2.0 * math.pi * particle_radius / wl for wl in wavelengths]

# Geometric cross section (cm^2), converting radius from microns to cm
r_cm = particle_radius * 1.0e-4
geom_cross_section = math.pi * r_cm ** 2

# ---------------------------------------------------------------------------
# Mie calculation — all wavelengths in one parallel batch call
# ---------------------------------------------------------------------------

results = lxmie.mie(n_real, n_imag, size_parameters)

# ---------------------------------------------------------------------------
# Print results
# ---------------------------------------------------------------------------

header = (
    f"{'wavelength (um)':>18}"
    f"{'size param':>14}"
    f"{'Q_ext':>12}"
    f"{'Q_sca':>12}"
    f"{'Q_abs':>12}"
    f"{'C_ext (cm^2)':>16}"
    f"{'asymmetry g':>14}"
)
print(header)
print("-" * len(header))

for wl, x, r in zip(wavelengths, size_parameters, results):
    print(
        f"{wl:18.6e}"
        f"{x:14.6e}"
        f"{r.q_ext:12.6f}"
        f"{r.q_sca:12.6f}"
        f"{r.q_abs:12.6f}"
        f"{r.q_ext * geom_cross_section:16.6e}"
        f"{r.asymmetry_parameter:14.6f}"
    )
