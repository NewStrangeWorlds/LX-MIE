Performance and parallelisation
===============================

Computing time and memory
-------------------------

The number of terms in the Mie series is approximately equal to the size
parameter :math:`x` (see :doc:`theory`). Both the computing time and the
memory required for a single Mie calculation therefore increase linearly with
the size parameter.

.. list-table:: Typical values for a single calculation of the efficiencies
   on one core of a desktop computer
   :header-rows: 1
   :widths: 30 35 35

   * - Size parameter
     - Computing time
     - Memory
   * - :math:`10^2`
     - :math:`\sim 10^{-5}` s
     - negligible
   * - :math:`10^4`
     - :math:`\sim 5 \times 10^{-4}` s
     - 0.4 MB
   * - :math:`10^6`
     - :math:`\sim 0.05` s
     - 40 MB
   * - :math:`10^8`
     - :math:`\sim 5` s
     - 4 GB

The memory requirement is about 40 bytes per term of the Mie series. A
calculation with :math:`x = 10^9` therefore requires about 40 GB of memory.
In a parallel calculation, this memory is needed for each thread.

Additional costs arise for the angular quantities:

- The scattering amplitudes require a summation of the Mie series for each
  angle. The computing time is thus proportional to the number of angles
  times the size parameter.
- The complete Legendre series of the phase function contains about
  :math:`2x` terms. The computing time for the Legendre coefficients
  increases with both the number of requested coefficients and the size
  parameter. For large size parameters, a complete Legendre series is by far
  the most expensive calculation.

For a size distribution, one Mie calculation is done for each radius of the
distribution, i.e. 500 calculations per wavelength by default.

.. note::

   To avoid repeated memory allocations, LX-MIE reuses its internal work
   arrays in subsequent calls. The memory of the largest calculation
   therefore remains allocated for each thread, also after the calculation
   has finished.

Parallelisation
---------------

LX-MIE uses OpenMP to distribute independent Mie calculations over several
threads. Which loop is run in parallel depends on the interface:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Calculation
     - Parallel loop
   * - Command-line tool
     - Wavelengths
   * - Python: batch call for single particles, e.g. :func:`lxmie.mie` with
       lists
     - Elements of the lists
   * - Python: scalar call for single particles
     - None (a single Mie calculation is not parallelised)
   * - Python: :func:`lxmie.mie_distribution` with scalar arguments
     - Particle sizes of the distribution
   * - Python: :func:`lxmie.mie_distribution` with lists
     - Wavelengths
   * - C++: ``lxmie::Mie``
     - None
   * - C++: ``lxmie::MieDistribution``
     - Particle sizes, unless called from within a parallel region

The loops use a dynamic schedule, such that calculations with very different
size parameters are distributed evenly over the threads.

Number of threads
~~~~~~~~~~~~~~~~~

By default, OpenMP uses all available cores. The number of threads is
controlled by the environment variable ``OMP_NUM_THREADS``:

.. code-block:: bash

   export OMP_NUM_THREADS=4

   ./lx_mie input/parameters.input
   python my_script.py

For Python, the variable has to be set before the first parallel calculation,
either in the shell or at the beginning of the script:

.. code-block:: python

   import os
   os.environ["OMP_NUM_THREADS"] = "4"

   import lxmie

Recommendations for Python
--------------------------

**Use batch calls.** A batch call avoids the overhead of the Python
interpreter and runs in parallel. On a multi-core machine, calculating
10 000 size parameters with a single batch call can be an order of magnitude
faster than a Python loop over scalar calls:

.. code-block:: python

   # slow: Python loop
   results = [lxmie.mie(1.5, 0.01, x) for x in size_parameters]

   # fast: one parallel batch call
   n = len(size_parameters)
   results = lxmie.mie([1.5] * n, [0.01] * n, size_parameters)

**Loop over wavelengths for size distributions.** Pass all wavelengths in one
call of :func:`lxmie.mie_distribution`, rather than calling the function in a
Python loop.

**Create size distributions once.** A :class:`lxmie.SizeDistribution` object
does not depend on the wavelength or the refractive index and can be reused
for any number of calculations.

**Multi-threading.** The global interpreter lock (GIL) of Python is released
during all batch calculations and all calculations for size distributions.
Other Python threads can therefore continue to run during a long
calculation.

**Memory of large batches.** In a batch call of
:func:`lxmie.mie_amplitudes` or :func:`lxmie.mie_legendre`, the results for
all particles are stored and returned at the end. For a large number of
particles and angles (or Legendre coefficients), these results can require a
considerable amount of memory.
