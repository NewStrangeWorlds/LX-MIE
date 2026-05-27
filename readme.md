# LX-MIE

**Author:** Daniel Kitzmann

LX-MIE is a Mie scattering code that computes:

- Extinction, scattering, and absorption efficiencies (or cross sections)
- The asymmetry parameter
- The scattering phase function at arbitrary angles
- The Legendre series representation of the phase function

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

Three example parameter files are provided in `input/`.

### Parameter file format

Lines beginning with `#` are treated as comments. Fields are read in the order shown below.

```
#file with refractive indices
compilation/Cr.dat

#particle radius (microns)
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

**Phase function** (`PF`):

```
# cos(angle)   phase function
```

**Legendre series** (`LS`):

```
# term index   coefficient
```

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

### Example scripts

Three example scripts in `examples/` demonstrate common workflows. Run them from the repository root directory.

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
```

---

## License

LX-MIE is free software distributed under the [GNU General Public License v3](licence.md).
