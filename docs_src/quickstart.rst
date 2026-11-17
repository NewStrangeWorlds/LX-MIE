Quick start
===========

This page shows the most common calculations with the Python package. The
same calculations can be done with the :doc:`command-line tool <cli/index>`
and the :doc:`C++ library <cpp/index>`.

A single particle
-----------------

A Mie calculation requires the complex refractive index of the particle and
the size parameter :math:`x = 2 \pi r / \lambda`, where :math:`r` is the
particle radius and :math:`\lambda` the wavelength.

.. code-block:: python

   import math
   import lxmie

   wavelength = 0.55   # micron
   radius = 1.0        # micron
   size_parameter = 2.0 * math.pi * radius / wavelength

   result = lxmie.mie(1.5, 0.01, size_parameter)

   print(result.q_ext, result.q_sca, result.q_abs, result.asymmetry_parameter)

   # cross sections: efficiencies times the geometric cross section
   c_ext = result.q_ext * math.pi * radius**2   # micron^2

Many particles or wavelengths at once
-------------------------------------

Every function also accepts lists. The calculations are then run in parallel
with OpenMP. The following example calculates the efficiencies of chromium
particles for all wavelengths of a refractive index file from the
:doc:`compilation <background/refractive_indices>`:

.. code-block:: python

   import math
   import lxmie

   # read the refractive index data (wavelength in micron, n, k)
   wavelengths, n_real, n_imag = [], [], []

   with open("compilation/Cr.dat") as f:
       for line in f:
           if line.startswith("#") or not line.strip():
               continue
           wl, n, k = line.split()
           wavelengths.append(float(wl))
           n_real.append(float(n))
           n_imag.append(float(k))

   radius = 10.0  # micron
   size_parameters = [2.0 * math.pi * radius / wl for wl in wavelengths]

   results = lxmie.mie(n_real, n_imag, size_parameters)

   for wl, r in zip(wavelengths, results):
       print(f"{wl:.4e}  Q_ext = {r.q_ext:.4f}  g = {r.asymmetry_parameter:.4f}")

Phase function
--------------

The phase function is obtained from the scattering amplitudes. Angles are
always passed as the **cosines** of the scattering angles.

.. code-block:: python

   cos_angles = [math.cos(math.radians(theta)) for theta in range(0, 181, 5)]

   result, s1, s2 = lxmie.mie_amplitudes(1.5, 0.01, size_parameter, cos_angles)
   phase_function = lxmie.phase_function(s1, s2, size_parameter, result.q_sca)

Radiative transfer codes often require the phase function as a Legendre
series instead:

.. code-block:: python

   result, moments = lxmie.mie_legendre(1.5, 0.01, size_parameter, 32)

   # moments[0] = 1, moments[1] = asymmetry parameter

Particle size distributions
---------------------------

For a size distribution, the wavelength is passed instead of the size
parameter. It has to be given in the same units as the particle radii.

.. code-block:: python

   # log-normal distribution: median radius 1 micron, sigma_g = 1.5
   distribution = lxmie.lognormal_distribution(1.0, 1.5)

   results = lxmie.mie_distribution(n_real, n_imag, wavelengths, distribution)

   for wl, r in zip(wavelengths, results):
       # mean cross sections per particle in micron^2, converted to cm^2
       print(f"{wl:.4e}  C_ext = {r.c_ext * 1e-8:.4e} cm^2")

Where to go next
----------------

- :doc:`conventions` explains the definitions and normalisations used
  throughout LX-MIE.
- :doc:`python/api` is the complete reference of the Python package.
- :doc:`background/size_distributions` describes the available size
  distributions and the size-integrated quantities.
- :doc:`cli/index` describes the command-line tool.
