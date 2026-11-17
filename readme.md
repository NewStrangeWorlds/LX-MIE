# LX-MIE

**Author:** Daniel Kitzmann

[![PyPI version](https://img.shields.io/pypi/v/lxmie)](https://pypi.org/project/lxmie/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

LX-MIE is a Mie scattering code that computes:

- Extinction, scattering, and absorption efficiencies (or cross sections)
- The asymmetry parameter
- The scattering phase function at arbitrary angles
- The Legendre series representation of the phase function
- All of the above for single particle sizes or integrated over a particle size distribution

LX-MIE is optimised for arbitrarily large size parameters and has been tested up to x = 10⁹. Details of the numerical method are described in:

> Kitzmann & Heng (2018), *Optical properties of potential condensates in exoplanetary atmospheres*, MNRAS 475, 94–107. https://arxiv.org/abs/1710.04946

LX-MIE is available both as a standalone command-line tool and as a Python package.

---

## Requirements

- C++17-compliant compiler (GCC ≥ 7, Clang ≥ 5, MSVC ≥ 2017)
- CMake ≥ 3.15
- OpenMP (optional but recommended — enables parallel computation over wavelengths)

---

## Building the command-line tool

```bash
cmake -B build
cmake --build build
```

The executable `lx_mie` is placed in the root directory. With OpenMP available, the wavelength loop runs in parallel automatically; without it the code builds and runs correctly in single-threaded mode.

---

## Running the command-line tool

```bash
lx_mie <parameter_file>
```

Four example parameter files are provided in `input/`.

### Parameter file format

Lines beginning with `#` are treated as comments. Fields are read in the order shown below.

```
#file with refractive indices
compilation/Cr.dat

#particle radius (microns) or a size distribution, see below
10.0

#output filename
output/cross_sections.dat

#output format  (CS = cross sections, EF = efficiencies [default])
CS

#compute phase function / Legendre series  (omit or leave blank for neither)
#   PF  →  phase function at a set of angles
#   LS  →  Legendre series of the phase function
PF

#output filename for phase function / Legendre series
output/phase_function.dat

# For PF: list the cosines of the scattering angles, one per line.
# For LS: provide the number of Legendre terms on a single line.
0.9982
0.9659
...
-0.9982
```

### Particle size distributions

Instead of a single particle radius, the corresponding line of the parameter file can also describe a particle size distribution. All output quantities are then integrated over this distribution (see [Size-integrated quantities](#size-integrated-quantities)).

| Entry in the parameter file | Distribution |
|---|---|
| `10.0` | Single particle radius (microns) |
| `lognormal <r_m> <sigma_g> [N]` | Log-normal distribution, n(r) ∝ 1/r · exp(−ln²(r/r_m) / (2 ln²σ_g)), with the median radius r_m (microns) and the geometric standard deviation σ_g > 1 |
| `gamma <r_eff> <v_eff> [N]` | Gamma distribution of Hansen & Travis (1974), n(r) ∝ r^((1−3b)/b) · exp(−r/(ab)), with the effective radius a = r_eff (microns) and the effective variance b = v_eff (0 < b < 0.5) |
| `powerlaw <r_min> <r_max> <p> [N]` | Power law n(r) ∝ r^(−p) between r_min and r_max (microns) |
| `file <file name>` | Tabulated distribution read from a file |

The optional parameter `N` is the number of radii used to integrate over the distribution (default: 500). They are distributed equidistantly in ln(r). For the log-normal and gamma distributions the integration limits are chosen automatically.

A file with a tabulated distribution contains two columns: the particle radius in microns (strictly increasing) and the distribution n(r) = dN/dr in arbitrary units. Lines beginning with `#` are comments. The distribution is integrated with the trapezoidal rule.

Example (`input/parameters_size_distribution.input`):

```
#particle radius (microns) or a size distribution (lognormal, gamma, powerlaw, file)
lognormal 1.0 1.5
```

### Refractive index file format

Lines beginning with `#` are comments. Each data line contains three columns:

```
wavelength (microns)    n    k
```

where the complex refractive index is m = n + ik. The `compilation/` directory contains refractive index data for 33 condensate species relevant to planetary and exoplanetary atmospheres.

### Output file formats

**Efficiencies / cross sections** (`CS` or `EF`):

```
# wavelength (um)   size parameter   Q_ext   Q_sca   Q_abs   (cross section)   asymmetry g
```

For a size distribution, the cross sections are the mean cross sections per particle and the efficiencies are the effective efficiencies of the distribution. The size parameter is calculated with the effective radius.

**Phase function** (`PF`):

```
# cos(angle)   phase function
```

**Legendre series** (`LS`):

```
# term index   coefficient
```

### Size-integrated quantities

For a size distribution n(r), LX-MIE calculates the following quantities:

- Mean cross sections per particle: ⟨C⟩ = ∫ C(r) n(r) dr / ∫ n(r) dr, with C(r) = Q(r) π r²
- Effective efficiencies: Q_eff = ⟨C⟩ / ⟨π r²⟩
- Asymmetry parameter, phase function, and Legendre coefficients: averages weighted by the scattering cross section, e.g. g = ∫ g(r) C_sca(r) n(r) dr / ∫ C_sca(r) n(r) dr
- Effective radius: r_eff = ⟨r³⟩ / ⟨r²⟩

Multiplying the mean cross sections by the total particle number density yields the extinction, scattering, and absorption coefficients. The normalisation of n(r) is therefore irrelevant.

Mie efficiencies of weakly absorbing particles oscillate rapidly with the particle size. With the default number of 500 radii, the size-integrated quantities are typically accurate to about 10⁻⁴ – 10⁻³. For higher accuracy or very broad distributions, the number of radii can be increased.

---

## Python interface

### Installation

The easiest way to install the Python package of LX_Mie is directly from PyPI:

```bash
pip install lxmie
```

As an alternative, the Python package can also be installed directly from the repository with pip (requires a C++17 compiler and CMake):

```bash
pip install .
```



### API reference

All functions are in the `lxmie` module. The refractive index is always passed as **separate real and imaginary parts** (`n_real`, `n_imag`).

The sign of `n_imag` does not matter: LX-MIE internally uses the convention m = n - ik, so the imaginary part is always converted to a negative value. Both positive and negative values of `k` (e.g. taken directly from the files in `compilation/`) therefore give the same, physically correct result.

Every function accepts either scalar arguments (single particle) or equal-length lists (batch of particles). List inputs trigger an OpenMP-parallel loop internally, with the GIL released during computation.

---

#### `lxmie.mie`

```python
result = lxmie.mie(n_real, n_imag, size_parameter)
results = lxmie.mie(n_real_list, n_imag_list, size_parameter_list)
```

Returns a `MieResult` (or `list[MieResult]` for list inputs) with attributes:

| Attribute | Description |
|---|---|
| `q_ext` | Extinction efficiency |
| `q_sca` | Scattering efficiency |
| `q_abs` | Absorption efficiency |
| `asymmetry_parameter` | Asymmetry parameter g |

---

#### `lxmie.mie_legendre`

```python
result, moments = lxmie.mie_legendre(n_real, n_imag, size_parameter, nb_moments)
results, all_moments = lxmie.mie_legendre(n_real_list, n_imag_list, size_parameter_list, nb_moments)
```

Computes Mie efficiencies and the first `nb_moments` Legendre series coefficients of the phase function. The Legendre coefficients (`moments`) are returned as a `list[float]` (or `list[list[float]]` for batch inputs).

---

#### `lxmie.mie_amplitudes`

```python
result, s1, s2 = lxmie.mie_amplitudes(n_real, n_imag, size_parameter, angles)
results, all_s1, all_s2 = lxmie.mie_amplitudes(n_real_list, n_imag_list, size_parameter_list, angles)
```

Computes Mie efficiencies and the complex scattering amplitudes S₁ and S₂ at the given angles. `angles` is a list of the **cosines** of the scattering angles. `s1` and `s2` are returned as `list[complex]`.

---

#### `lxmie.phase_function`

```python
pf = lxmie.phase_function(s1, s2, size_parameter, q_sca)          # list[float]
pf = lxmie.phase_function(s1[i], s2[i], size_parameter, q_sca)    # float
```

Computes the normalised phase function from the scattering amplitudes. The vector overload accepts `list[complex]` for `s1` and `s2`; the scalar overload accepts individual `complex` values.

---

#### Size distributions

A size distribution is represented by a `SizeDistribution` object, which contains a list of radii (`radius`) and the number of particles for each radius (`weight`). The weights do not need to be normalised. A distribution can be created with one of the following functions:

```python
dist = lxmie.lognormal_distribution(median_radius, sigma_g, nb_points=500, r_min=0, r_max=0)
dist = lxmie.gamma_distribution(effective_radius, effective_variance, nb_points=500, r_min=0, r_max=0)
dist = lxmie.power_law_distribution(r_min, r_max, exponent, nb_points=500)
dist = lxmie.tabulated_distribution(radii, dn_dr)
dist = lxmie.SizeDistribution(radii, weights)
```

The first four correspond to the distributions described in [Particle size distributions](#particle-size-distributions). For the log-normal and gamma distributions, the integration limits are chosen automatically, unless `r_min` and/or `r_max` are set to non-zero values. `tabulated_distribution` takes the values of n(r) = dN/dr at a list of strictly increasing radii and integrates them with the trapezoidal rule. `SizeDistribution(radii, weights)` directly uses the given radii and weights, e.g. the size bins of a cloud model with the particle number (density) in each bin.

---

#### `lxmie.mie_distribution`

```python
result = lxmie.mie_distribution(n_real, n_imag, wavelength, dist)
results = lxmie.mie_distribution(n_real_list, n_imag_list, wavelength_list, dist)
```

Computes the optical properties integrated over the size distribution `dist` (see [Size-integrated quantities](#size-integrated-quantities)). In contrast to the single-particle functions, the **wavelength** is passed instead of the size parameter. It has to be given in the same units as the radii of the distribution.

Returns a `MieDistributionResult` (or `list[MieDistributionResult]` for list inputs) with attributes:

| Attribute | Description |
|---|---|
| `c_ext`, `c_sca`, `c_abs` | Mean extinction, scattering, and absorption cross sections per particle, in units of the squared length unit of the wavelength and radii (e.g. microns²) |
| `q_ext`, `q_sca`, `q_abs` | Effective efficiencies (mean cross sections divided by the mean geometric cross section) |
| `asymmetry_parameter` | Asymmetry parameter g of the distribution |
| `geometric_cross_section` | Mean geometric cross section ⟨π r²⟩ |
| `effective_radius` | Effective radius ⟨r³⟩ / ⟨r²⟩ |

For scalar arguments, the loop over the particle sizes runs in parallel; for list arguments, the loop over the wavelengths.

---

#### `lxmie.mie_distribution_legendre`

```python
result, moments = lxmie.mie_distribution_legendre(n_real, n_imag, wavelength, dist, nb_moments)
results, all_moments = lxmie.mie_distribution_legendre(n_real_list, n_imag_list, wavelength_list, dist, nb_moments)
```

Additionally returns the first `nb_moments` Legendre series coefficients of the size-averaged phase function, with the same normalisation as `lxmie.mie_legendre`.

---

#### `lxmie.mie_distribution_phase_function`

```python
result, pf = lxmie.mie_distribution_phase_function(n_real, n_imag, wavelength, dist, angles)
results, all_pf = lxmie.mie_distribution_phase_function(n_real_list, n_imag_list, wavelength_list, dist, angles)
```

Additionally returns the size-averaged phase function at the given angles. `angles` is a list of the **cosines** of the scattering angles. The phase function has the same normalisation as the one returned by `lxmie.phase_function`.

---

### Example scripts

Four example scripts in `examples/` demonstrate common workflows. Run them from the repository root directory.

**`examples/example_efficiencies.py`**
Reads `compilation/Cr.dat` and computes efficiencies and cross sections across all 273 wavelengths in a single parallel batch call.

```bash
python examples/example_efficiencies.py
```

**`examples/example_legendre.py`**
Computes Legendre series coefficients of the phase function, showing both single-wavelength and batch usage.

```bash
python examples/example_legendre.py
```

**`examples/example_phase_function.py`**
Computes the phase function at 37 angles from 0° to 180°. Demonstrates the conversion from degrees to cosines required by `mie_amplitudes`, and shows how to use `phase_function()` to normalise the scattering amplitudes.

```bash
python examples/example_phase_function.py
```

**`examples/example_size_distribution.py`**
Computes cross sections integrated over a log-normal size distribution for all wavelengths of `compilation/MgSiO3_amorph_glass.dat`, as well as the size-averaged phase function and its Legendre series.

```bash
python examples/example_size_distribution.py
```

### Quick-start example

```python
import math
import lxmie

# Read refractive index data (wavelength in microns, n, k)
wavelengths, n_real, n_imag = [], [], []
with open("compilation/Cr.dat") as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        wl, n, k = line.split()
        wavelengths.append(float(wl))
        n_real.append(float(n))
        n_imag.append(float(k))

particle_radius = 10.0  # microns
size_params = [2 * math.pi * particle_radius / wl for wl in wavelengths]

# Compute efficiencies for all wavelengths in parallel
results = lxmie.mie(n_real, n_imag, size_params)

for wl, r in zip(wavelengths, results):
    print(f"{wl:.4e} um   Q_ext={r.q_ext:.4f}   g={r.asymmetry_parameter:.4f}")

# Compute phase function at 37 angles for a single wavelength
cos_angles = [math.cos(math.radians(i * 5)) for i in range(37)]
result, s1, s2 = lxmie.mie_amplitudes(n_real[0], n_imag[0], size_params[0], cos_angles)
phase_function = lxmie.phase_function(s1, s2, size_params[0], result.q_sca)

# Integrate over a log-normal size distribution (median radius 1 micron, sigma_g = 1.5)
# Radii and wavelengths in microns, cross sections are returned in microns^2
distribution = lxmie.lognormal_distribution(1.0, 1.5)
results = lxmie.mie_distribution(n_real, n_imag, wavelengths, distribution)

for wl, r in zip(wavelengths, results):
    print(f"{wl:.4e} um   C_ext={r.c_ext * 1e-8:.4e} cm^2   g={r.asymmetry_parameter:.4f}")
```

---

## License

LX-MIE is free software distributed under the [GNU General Public License v3](licence.md). 
