Overview
========

The Python package ``lxmie`` is a thin interface to the C++ core of LX-MIE
(see :doc:`../installation`). All functions are located directly in the
``lxmie`` module:

.. code-block:: python

   import lxmie

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Function
     - Result
   * - :func:`lxmie.mie`
     - Efficiencies and asymmetry parameter
   * - :func:`lxmie.mie_legendre`
     - … plus the Legendre series of the phase function
   * - :func:`lxmie.mie_amplitudes`
     - … plus the scattering amplitudes at given angles
   * - :func:`lxmie.phase_function`
     - Phase function from the scattering amplitudes
   * - :func:`lxmie.lognormal_distribution`,
       :func:`lxmie.gamma_distribution`,
       :func:`lxmie.power_law_distribution`,
       :func:`lxmie.tabulated_distribution`,
       :class:`lxmie.SizeDistribution`
     - Particle size distributions
   * - :func:`lxmie.mie_distribution`
     - Cross sections, efficiencies, and asymmetry parameter integrated over
       a size distribution
   * - :func:`lxmie.mie_distribution_legendre`
     - … plus the Legendre series of the size-averaged phase function
   * - :func:`lxmie.mie_distribution_phase_function`
     - … plus the size-averaged phase function at given angles

General rules
-------------

A few rules apply to all functions of the package. The definitions of the
physical quantities are given in :doc:`../conventions`.

Refractive index
   The refractive index is always passed as two separate real numbers,
   ``n_real`` and ``n_imag``. The sign of ``n_imag`` does not matter.

Size parameter or wavelength
   The single-particle functions take the dimensionless **size parameter**
   :math:`x = 2 \pi r / \lambda`. The functions for size distributions take
   the **wavelength** instead, in the same units as the radii of the
   distribution.

Angles
   Angles are always given as the **cosines** of the scattering angles.

Scalar and batch calls
----------------------

Every calculation function exists in two versions that are selected
automatically from the type of the arguments:

- With **scalar** arguments, a single calculation is done and a single result
  is returned.
- With **list** arguments for the refractive index and the size parameter (or
  wavelength), one calculation is done for each element of the lists. The
  results are returned as lists in the same order.

.. code-block:: python

   # scalar call: a single MieResult
   result = lxmie.mie(1.5, 0.01, 5.0)

   # batch call: a list of MieResult objects
   results = lxmie.mie([1.5, 1.5, 1.6], [0.01, 0.02, 0.0], [5.0, 10.0, 20.0])

In a batch call, the lists ``n_real``, ``n_imag``, and ``size_parameter`` (or
``wavelength``) have to be of the **same length**. Scalars and lists cannot
be mixed. To calculate, for example, many size parameters for the same
refractive index, the refractive index has to be repeated:

.. code-block:: python

   x = [0.1, 1.0, 10.0, 100.0]
   results = lxmie.mie([1.5] * len(x), [0.01] * len(x), x)

Additional arguments, such as the list of angles, the number of Legendre
coefficients, or the size distribution, are shared by all elements of a batch.

.. warning::

   The functions for size distributions verify that the lists have the same
   length and raise a ``ValueError`` otherwise. The single-particle
   functions :func:`lxmie.mie`, :func:`lxmie.mie_legendre`, and
   :func:`lxmie.mie_amplitudes` currently do **not** perform this check. The
   number of calculations is given by the length of ``n_real``; shorter
   ``n_imag`` or ``size_parameter`` lists lead to undefined results.

Batch calls are run in parallel with OpenMP and release Python's global
interpreter lock (GIL) during the calculation. They are therefore much faster
than a Python loop over scalar calls. Details are given in
:doc:`../background/performance`.

NumPy arrays
------------

LX-MIE does not depend on NumPy, but one-dimensional NumPy arrays can be
passed wherever a list is expected. Results are always returned as Python
lists and can be converted to arrays if required:

.. code-block:: python

   import numpy as np
   import lxmie

   x = np.logspace(-1, 2, 1000)
   n = np.full_like(x, 1.5)
   k = np.full_like(x, 0.01)

   results = lxmie.mie(n, k, x)

   q_ext = np.array([r.q_ext for r in results])
   g = np.array([r.asymmetry_parameter for r in results])

Errors
------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Exception
     - Raised for
   * - ``ValueError``
     - Invalid parameters of a size distribution, an invalid wavelength, or
       lists of different lengths in the functions for size distributions.
   * - ``RuntimeError``
     - Failures within the Mie calculation, e.g. for a size parameter of
       zero. The size parameter has to be a positive, finite number.
   * - ``TypeError``
     - Arguments that do not match any of the two versions of a function,
       e.g. a mixture of scalars and lists.

If an error occurs in one element of a batch call, the exception is raised for
the entire call.
