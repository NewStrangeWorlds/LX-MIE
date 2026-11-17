Running LX-MIE
==============

The command-line tool ``lx_mie`` performs Mie calculations for all wavelengths
of a refractive index file and writes the results to text files. It is built
with CMake as described in :doc:`../installation`.

Usage
-----

.. code-block:: bash

   lx_mie <parameter_file>

The only argument is the name of a :doc:`parameter file <parameter_file>`
that specifies

- the file with the refractive indices,
- the particle radius or a particle size distribution,
- the output file and the output format (cross sections or efficiencies),
- optionally, the calculation of the phase function or its Legendre series.

LX-MIE prints the parameters it has read from the file and, after a
successful calculation, the message ``Model finished!``:

.. code-block:: text

   $ ./lx_mie input/parameters.input
   read from parameter file:
   refractive index: compilation/Cr.dat
   particle radius: 10
   output file: output/cross_sections.dat
   use cross section 1
   compute phase function 0
   Model finished!

.. important::

   - All file names in the parameter file are interpreted relative to the
     directory from which ``lx_mie`` is started, not relative to the location
     of the parameter file.
   - The directories of the output files have to exist. LX-MIE does not
     create them. The example parameter files write to the directory
     ``output/``, which has to be created first (``mkdir output``).

Example parameter files
-----------------------

Four example parameter files are provided in the ``input/`` directory. They
are meant to be run from the root directory of the repository.

.. list-table::
   :header-rows: 1
   :widths: 44 56

   * - File
     - Calculation
   * - ``parameters.input``
     - Cross sections of chromium particles with a radius of 10 |micron|
       for all wavelengths of ``compilation/Cr.dat``.
   * - ``parameters_phase_function.input``
     - Additionally calculates the phase function at a set of angles, for a
       single wavelength (``input/single_refractive_index.dat``).
   * - ``parameters_legendre_series.input``
     - Additionally calculates the first 30 coefficients of the Legendre
       series of the phase function.
   * - ``parameters_size_distribution.input``
     - Cross sections of MgSiO\ :sub:`3` particles integrated over a
       log-normal size distribution.

Input and output
----------------

The formats of all input files are described in :doc:`parameter_file`, those
of the output files in :doc:`output_files`. Wavelengths and particle radii are
always given in |micron|.

Parallelisation
---------------

If LX-MIE has been compiled with OpenMP, the loop over the wavelengths is run
in parallel. The number of threads can be set with the environment variable
``OMP_NUM_THREADS``:

.. code-block:: bash

   OMP_NUM_THREADS=4 ./lx_mie input/parameters.input

Exit status and error messages
------------------------------

``lx_mie`` returns the exit status 0 after a successful calculation and 1 if

- no parameter file was given,
- the parameter file, the refractive index file, or the file of a tabulated
  size distribution could not be read,
- the particle radius or the parameters of the size distribution are
  invalid,
- an error occurred during the Mie calculations.

.. warning::

   If an output file cannot be created, e.g. because its directory does not
   exist, LX-MIE prints the message ``Could not create output file``, but
   still ends with ``Model finished!`` and the exit status 0. In this case,
   the results are lost.
