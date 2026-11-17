Output files
============

All output files are tab-separated text files with a header line that starts
with ``#``. The numbers are written in scientific notation.

Cross sections
--------------

With the output format ``CS``, the main output file contains one line per
wavelength with the following columns:

.. list-table::
   :header-rows: 1
   :widths: 10 40 50

   * - Column
     - Quantity
     - Unit
   * - 1
     - Wavelength
     - |micron|
   * - 2
     - Size parameter
     - –
   * - 3
     - Extinction cross section
     - cm\ :sup:`2`
   * - 4
     - Scattering cross section
     - cm\ :sup:`2`
   * - 5
     - Absorption cross section
     - cm\ :sup:`2`
   * - 6
     - Single scattering albedo
     - –
   * - 7
     - Asymmetry parameter
     - –

.. code-block:: text

   #wavelengths (mu)   size parameter     extinction cross section (cm^2)   scattering cross section (cm^2)   absorption cross section (cm^2)   single scattering albedo   asymmetry parameter
   4.2030000000e-02    1.4949286955e+03   6.3230336336e-06                  3.7216107144e-06                  2.6014229192e-06                  5.8857993332e-01           9.3407376153e-01
   4.3050000000e-02    1.4595087822e+03   6.3238877163e-06                  3.7434247712e-06                  2.5804629451e-06                  5.9194991106e-01           9.2986422962e-01

Efficiencies
------------

With any other output format (``EF``), the file contains the efficiencies:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Column
     - Quantity
   * - 1
     - Wavelength (|micron|)
   * - 2
     - Size parameter
   * - 3
     - Absorption efficiency :math:`Q_\mathrm{abs}`
   * - 4
     - Scattering efficiency :math:`Q_\mathrm{sca}`
   * - 5
     - Extinction efficiency :math:`Q_\mathrm{ext}`
   * - 6
     - Asymmetry parameter

.. note::

   The order of the columns differs between the two formats: cross sections
   are written in the order extinction – scattering – absorption, the
   efficiencies in the order absorption – scattering – extinction. The single
   scattering albedo is only included in the cross section format.

Size distributions
~~~~~~~~~~~~~~~~~~

For a particle size distribution, the columns have the following meaning:

- The cross sections are the mean cross sections per particle.
- The efficiencies are the effective efficiencies of the distribution, i.e.
  the mean cross sections divided by the mean geometric cross section.
- The size parameter is calculated with the effective radius of the
  distribution.
- The asymmetry parameter is the average weighted by the scattering cross
  sections.

The definitions are given in :ref:`size-integrated-quantities`.

Phase function
--------------

The phase function file (``PF``) contains one line for each angle. The first
column is the cosine of the scattering angle, followed by **one column per
wavelength** with the values of the phase function. The columns are in the
order of the wavelengths in the refractive index file.

.. code-block:: text

   #cos(angle)        phase function
   9.9820000000e-01   1.6571346029e+00   ...
   9.9070000000e-01   1.6434862967e+00   ...

The normalisation of the phase function is described in
:ref:`phase-function-normalisation`.

Legendre series
---------------

The Legendre series file (``LS``) contains one line for each coefficient. The
first column is the index :math:`l` of the coefficient, starting at zero,
followed by **one column per wavelength** with the coefficients
:math:`\chi_l`.

.. code-block:: text

   #Legendre series terms
   0   1.0000000000e+00   ...
   1   4.1650385436e-02   ...

The first coefficient is always 1, the second one is equal to the asymmetry
parameter (see :ref:`legendre-series`).

Reading the output with Python
------------------------------

.. code-block:: python

   import numpy as np

   # cross sections
   data = np.loadtxt("output/cross_sections.dat", ndmin=2)

   wavelength = data[:, 0]
   c_ext, c_sca, c_abs = data[:, 2], data[:, 3], data[:, 4]
   albedo, g = data[:, 5], data[:, 6]

   # phase function: one row per angle, one column per wavelength
   data = np.loadtxt("output/phase_function.dat", ndmin=2)

   cos_angles = data[:, 0]
   phase_function = data[:, 1:]     # phase_function[i_angle, i_wavelength]
