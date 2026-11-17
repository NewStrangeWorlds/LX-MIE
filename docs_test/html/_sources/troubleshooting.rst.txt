Troubleshooting
===============

Python package
--------------

``AttributeError: module 'lxmie' has no attribute 'mie_distribution'``
   The installed version of the package does not include particle size
   distributions yet. Update it with ``pip install --upgrade lxmie``.

``TypeError: mie(): incompatible function arguments``
   The arguments do not match the scalar or the list version of the
   function. The most common reason is a mixture of scalars and lists, e.g.
   a scalar refractive index together with a list of size parameters. In a
   batch call, ``n_real``, ``n_imag``, and ``size_parameter`` (or
   ``wavelength``) all have to be lists (see :doc:`python/index`).

``RuntimeError: LX-MIE: complex continued fraction did not converge``
   The Mie calculation failed. This happens for a size parameter of zero or
   infinity. Check the size parameters and refractive indices that are
   passed to the function.

``ValueError: vector::_M_default_append`` (or a similar message about a vector length)
   The size parameter is negative or not a number (``nan``).

``ValueError: LX-MIE: ...``
   The parameters of a size distribution or the wavelength are invalid. The
   message describes the reason, e.g. a geometric standard deviation
   :math:`\sigma_g \leq 1` of a log-normal distribution.

The results of a batch call look wrong or change between runs
   Check that all lists have the same length. The single-particle functions
   do not verify this (see :doc:`python/index`).

The phase function looks wrong
   Angles have to be passed as the cosines of the scattering angles, not in
   degrees or radians.

The results for a size distribution are wrong by orders of magnitude
   The wavelength and the radii of the distribution have to be in the same
   units. Note that the functions for size distributions take the
   wavelength, not the size parameter.

``pip install .`` fails with a CMake error
   Building from source requires a C++17 compiler and CMake ≥ 3.15. If CMake
   cannot find ``pybind11``, update ``pip`` (``pip install --upgrade pip``)
   and try again.

Command-line tool
-----------------

``parameter file ... not found``, ``refractive index file ... not found``
   File names are interpreted relative to the directory from which ``lx_mie``
   is started. The example parameter files have to be run from the root
   directory of the repository.

``Could not create output file``
   The directory of the output file does not exist. LX-MIE does not create
   directories. For the example parameter files, run ``mkdir output`` first.

``Error reading the particle size from the parameter file``
   The particle size entry is neither a number nor a valid description of a
   size distribution. Check the keyword (``lognormal``, ``gamma``,
   ``powerlaw``, ``file``) and the number of its parameters (see
   :ref:`cli-size-distributions`).

The parameters printed by LX-MIE do not match the parameter file
   The parameter file is read by position. Every block has to consist of
   exactly one comment line, one line with the value, and one empty line (see
   :doc:`cli/parameter_file`). Compare the file with the examples in the
   ``input/`` directory.

The first wavelengths of the refractive index file are missing in the output
   The first three lines of a refractive index file are always treated as
   header lines. A file with fewer header lines loses its first data lines.

The output ends before the last wavelength of the refractive index file
   Reading of the refractive index file stops at the first line that does not
   consist of three numbers, e.g. a comment or an empty line within the
   data.

The last angle is missing in the phase function output
   The last line of the parameter file has to end with a line break.

The absorption cross section is slightly negative
   The absorption is calculated as the difference of extinction and
   scattering. For non-absorbing materials (:math:`k = 0`), the result
   is zero within the numerical round-off error, which can be negative.

Results
-------

The efficiencies oscillate strongly with wavelength or particle size
   This is a physical property of weakly absorbing spheres of a single size:
   the interference structure and the sharp resonances are part of the exact
   Mie solution. In nature, they are usually smoothed out by a distribution
   of particle sizes (see :doc:`background/size_distributions`).

The extinction efficiency approaches 2 for large particles, not 1
   This is the correct limit of Mie theory, known as the extinction paradox:
   in addition to the light that is geometrically blocked, the same amount of
   light is removed from the forward direction by diffraction.

The Legendre coefficients are zero beyond a certain index
   The Legendre series of a particle with the size parameter :math:`x` only
   contains about :math:`2x` non-zero terms. For small particles, the series
   is very short.

Differences to other Mie codes
   Check the conventions: the definition of the size parameter, the
   normalisation of the phase function and the Legendre coefficients
   (see :doc:`conventions`), and the sign of the imaginary part of the
   refractive index. The scattering amplitudes of LX-MIE are the complex
   conjugates of those in the convention of Bohren & Huffman (1983).

If a problem cannot be solved, please open an issue on
`GitHub <https://github.com/newstrangeworlds/lx-mie/issues>`_.
