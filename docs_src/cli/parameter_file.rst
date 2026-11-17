Input files
===========

Parameter file
--------------

The parameter file is a plain text file that consists of a sequence of
blocks. Each block contains a comment line, starting with ``#``, followed by
a line with the value of the parameter. The blocks are separated by an empty
line.

.. code-block:: text
   :caption: Complete parameter file with all blocks

   #file with refractive indices
   compilation/Cr.dat

   #particle radius (microns)
   10.0

   #output filename
   output/cross_sections.dat

   #output format (cross section: CS, efficiencies: EF (default) )
   CS

   #compute phase function/Legendre series (No (default), phase function: PF, Legendre series: LS )
   PF

   #output filename for phase function/Legendre series
   output/phase_function.dat

   #cosine of phase function angles
   0.9982
   0.9659
   0.8660
   0.0
   -0.8660
   -0.9982

.. important::

   The parameter file is read **by position**, not by keywords. The blocks
   have to appear in exactly the order shown above, and each block has to
   consist of exactly one comment line, one line with the value, and one
   empty line. The text of the comment lines is ignored, but the lines
   themselves cannot be omitted. The safest way to create a new parameter
   file is to copy and modify one of the examples in the ``input/``
   directory.

   The only exception is the particle size: additional comment lines and
   empty lines are allowed directly in front of it.

The blocks of the parameter file are:

.. list-table::
   :header-rows: 1
   :widths: 6 28 66

   * - #
     - Parameter
     - Description
   * - 1
     - Refractive index file
     - Name of the file with the refractive indices, see
       :ref:`below <refractive-index-file>`. File names cannot contain
       spaces.
   * - 2
     - Particle size
     - Either a single particle radius in |micron| or the description of a
       size distribution, see :ref:`below <cli-size-distributions>`.
   * - 3
     - Output file
     - Name of the output file for the cross sections or efficiencies.
   * - 4
     - Output format
     - ``CS`` for cross sections. For any other entry (usually ``EF``),
       efficiencies are written.
   * - 5
     - Angular output (optional)
     - ``PF`` for the phase function at a set of angles, ``LS`` for the
       Legendre series of the phase function. For any other entry, neither of
       them is calculated.
   * - 6
     - Angular output file
     - Name of the output file for the phase function or the Legendre
       series. Only required for ``PF`` and ``LS``.
   * - 7
     - Angles or number of terms
     - For ``PF``: the cosines of the scattering angles, one per line. For
       ``LS``: the number of Legendre series coefficients. Only required for
       ``PF`` and ``LS``.

The blocks 5 – 7 are optional. If the phase function is not required, the
parameter file can end after the output format.

Phase function angles
~~~~~~~~~~~~~~~~~~~~~

For the phase function (``PF``), the last block lists the **cosines** of the
scattering angles, with one value per line. The values can be given in any
order. Every line after the comment line is interpreted as an angle. The list
must therefore not contain empty lines or comments.

.. warning::

   The last line of the parameter file has to be terminated by a line break.
   Otherwise, the last angle of the list is ignored.

Legendre series
~~~~~~~~~~~~~~~

For the Legendre series (``LS``), the last block contains the number of
coefficients:

.. code-block:: text

   #compute phase function/Legendre series (No (default), phase function: PF, Legendre series: LS )
   LS

   #output filename for phase function/Legendre series
   output/legendre_series.dat

   #Number of Legendre series terms
   30

The normalisation of the coefficients is described in
:ref:`legendre-series`.

.. _cli-size-distributions:

Particle size distributions
---------------------------

Instead of a single particle radius, the corresponding line of the parameter
file can describe a particle size distribution. All output quantities are
then integrated over this distribution.

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Entry in the parameter file
     - Distribution
   * - ``10.0``
     - Single particle radius (|micron|)
   * - ``lognormal <r_m> <sigma_g> [N]``
     - Log-normal distribution with the median radius :math:`r_m`
       (|micron|) and the geometric standard deviation
       :math:`\sigma_g > 1`
   * - ``gamma <r_eff> <v_eff> [N]``
     - Gamma distribution of Hansen & Travis (1974) with the effective
       radius :math:`r_\mathrm{eff}` (|micron|) and the effective variance
       :math:`0 < v_\mathrm{eff} < 0.5`
   * - ``powerlaw <r_min> <r_max> <p> [N]``
     - Power law :math:`n(r) \propto r^{-p}` between :math:`r_\mathrm{min}`
       and :math:`r_\mathrm{max}` (|micron|)
   * - ``file <file name>``
     - Tabulated distribution read from a file

The optional parameter ``N`` is the number of radii that are used to integrate
over the distribution (default: 500). The definitions of the distributions,
the integration grids, and the calculated size-integrated quantities are
described in :doc:`../background/size_distributions`.

.. code-block:: text
   :caption: Example: input/parameters_size_distribution.input

   #file with refractive indices
   compilation/MgSiO3_amorph_glass.dat

   #particle radius (microns) or a size distribution (lognormal, gamma, powerlaw, file)
   lognormal 1.0 1.5

   #output filename
   output/cross_sections.dat

   #output format (cross section: CS, efficiencies: EF (default) )
   CS

Further examples for the particle size entry:

.. code-block:: text

   lognormal 1.0 1.5 2000
   gamma 1.0 0.1
   powerlaw 0.005 0.25 3.5
   file input/my_distribution.dat

Tabulated distributions
~~~~~~~~~~~~~~~~~~~~~~~

A file with a tabulated distribution contains two columns:

1. the particle radius in |micron|, strictly increasing,
2. the distribution :math:`n(r) = \mathrm{d}N/\mathrm{d}r` in arbitrary units.

Empty lines and lines beginning with ``#`` are ignored. The file has to
contain at least two radii. The distribution is integrated with the
trapezoidal rule.

.. code-block:: text
   :caption: Example of a tabulated size distribution

   # radius (micron)   dN/dr
   0.10    0.00
   0.20    0.35
   0.40    1.00
   0.80    0.42
   1.60    0.05
   3.20    0.00

.. _refractive-index-file:

Refractive index file
---------------------

The refractive index file contains the wavelengths for which the calculations
are done and the corresponding complex refractive indices of the particle
material.

.. code-block:: text
   :caption: First lines of compilation/Cr.dat

   #Refractive index of Cr[s]; Lynch&Hunter (1991) in Palik: "Handbook of Optical Constants of Solids"; Rakic et al. (1998) Applied Optics Vol. 37, Issue 22
   #Data points: 273
   #Wavelength (micron) n  k
   4.203000e-02	7.800000e-01	2.100000e-01
   4.305000e-02	7.700000e-01	2.200000e-01
   ...

The file has the following format:

- The **first three lines** are header lines. They are always skipped,
  whether or not they begin with ``#``.
- Each of the following lines contains three columns, separated by spaces or
  tabs:

  1. the wavelength in |micron|,
  2. the real part :math:`n` of the refractive index,
  3. the imaginary part :math:`k` of the refractive index.

The sign of :math:`k` does not matter (see :doc:`../conventions`). Reading
stops at the first line that cannot be interpreted as three numbers. Comments
or empty lines within the data are therefore not allowed.

One Mie calculation is done for each wavelength of the file. LX-MIE does not
interpolate the refractive indices. To obtain results on a specific
wavelength grid, the refractive indices have to be interpolated to this grid
beforehand.

The ``compilation/`` directory of the repository contains refractive index
files for 32 condensate species, see :doc:`../background/refractive_indices`.
