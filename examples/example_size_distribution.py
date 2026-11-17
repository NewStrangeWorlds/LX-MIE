"""
Compute optical properties of MgSiO3 particles integrated over a particle
size distribution.

Demonstrates the available size distributions, the size-integrated cross
sections for all wavelengths of a refractive index file, and the
size-averaged phase function and its Legendre series.

Run from the lx-mie root directory:
    python examples/example_size_distribution.py
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

refindex_file = "compilation/MgSiO3_amorph_glass.dat"

wavelengths, n_real, n_imag = read_refractive_index(refindex_file)

# ---------------------------------------------------------------------------
# Size distributions
#
# Radii have to be given in the same units as the wavelengths (here: microns).
# The size-integrated cross sections are then returned in microns^2.
# ---------------------------------------------------------------------------

# Log-normal distribution: median radius 1 micron, geometric standard deviation 1.5
distribution = lxmie.lognormal_distribution(1.0, 1.5)

# Other available distributions:
#
# Gamma distribution (Hansen & Travis 1974): effective radius, effective variance
#   distribution = lxmie.gamma_distribution(1.0, 0.1)
#
# Power law n(r) ~ r^-3.5 between 0.005 and 0.25 microns
#   distribution = lxmie.power_law_distribution(0.005, 0.25, 3.5)
#
# Tabulated distribution dN/dr on a grid of radii
#   distribution = lxmie.tabulated_distribution(radii, dn_dr)
#
# Discrete size bins: radii and the number of particles in each bin
#   distribution = lxmie.SizeDistribution(radii, numbers)

print(distribution)

# ---------------------------------------------------------------------------
# Size-integrated cross sections — all wavelengths in one parallel batch call
# ---------------------------------------------------------------------------

results = lxmie.mie_distribution(n_real, n_imag, wavelengths, distribution)

print(f"effective radius: {results[0].effective_radius:.4f} um")
print()

micron2_to_cm2 = 1.0e-8

header = (
    f"{'wavelength (um)':>18}"
    f"{'C_ext (cm^2)':>16}"
    f"{'C_sca (cm^2)':>16}"
    f"{'C_abs (cm^2)':>16}"
    f"{'Q_ext':>12}"
    f"{'asymmetry g':>14}"
)
print(header)
print("-" * len(header))

for wl, r in zip(wavelengths, results):
    print(
        f"{wl:18.6e}"
        f"{r.c_ext * micron2_to_cm2:16.6e}"
        f"{r.c_sca * micron2_to_cm2:16.6e}"
        f"{r.c_abs * micron2_to_cm2:16.6e}"
        f"{r.q_ext:12.6f}"
        f"{r.asymmetry_parameter:14.6f}"
    )

# ---------------------------------------------------------------------------
# Size-averaged phase function and Legendre series for a single wavelength
# ---------------------------------------------------------------------------

i = 0
cos_angles = [math.cos(math.radians(a)) for a in range(0, 181, 15)]

result, phase_function = lxmie.mie_distribution_phase_function(
    n_real[i], n_imag[i], wavelengths[i], distribution, cos_angles)

print()
print(f"Size-averaged phase function at {wavelengths[i]:.4e} um")
print(f"{'angle (deg)':>12}{'phase function':>18}")
for a, p in zip(range(0, 181, 15), phase_function):
    print(f"{a:12d}{p:18.6e}")

nb_moments = 10
result, moments = lxmie.mie_distribution_legendre(
    n_real[i], n_imag[i], wavelengths[i], distribution, nb_moments)

print()
print(f"First {nb_moments} Legendre coefficients (the second one equals g = {result.asymmetry_parameter:.6f})")
for l, m in enumerate(moments):
    print(f"{l:4d}{m:16.6e}")
